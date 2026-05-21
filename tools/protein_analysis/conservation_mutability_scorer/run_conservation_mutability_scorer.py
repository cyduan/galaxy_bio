#!/usr/bin/env python3
"""Score residue conservation and mutability from a protein MSA."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS_AA = set("BJOUXZ")
GAP_CHARS = {"-", "."}

OUTPUT_COLUMNS = [
    "residue_index",
    "chain",
    "wt_aa",
    "alignment_column",
    "conservation_score",
    "mutability_score",
    "msa_entropy",
    "consensus_aa",
    "accepted_aas",
    "accepted_aas_with_frequency",
    "rate4site_score",
    "rate4site_score_normalized",
    "non_gap_count",
    "gap_fraction",
    "status",
    "warning",
]


@dataclass
class FastaRecord:
    identifier: str
    description: str
    sequence: str


@dataclass
class SiteStats:
    residue_index: int
    chain: str
    wt_aa: str
    alignment_column: int
    msa_entropy: float
    consensus_aa: str
    accepted_aas: list[str]
    accepted_aas_with_frequency: str
    non_gap_count: int
    gap_fraction: float
    status: str
    warning: str
    rate4site_score: float | None = None
    rate4site_score_normalized: float | None = None
    conservation_score: float | None = None
    mutability_score: float | None = None


class ToolError(RuntimeError):
    """User-facing error for Galaxy."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate conservation and mutability scores from an MSA, optionally using Rate4Site."
    )
    parser.add_argument("--msa-fasta", required=True, help="Input multiple sequence alignment FASTA.")
    parser.add_argument("--output-tsv", required=True, help="HotSpot workflow residue conservation TSV.")
    parser.add_argument("--rate4site-output", required=True, help="Preserved raw Rate4Site output.")
    parser.add_argument("--run-log", required=True, help="Run log output.")
    parser.add_argument(
        "--backend",
        choices=["rate4site", "auto", "entropy_only"],
        default="auto",
        help="Scoring backend. auto tries Rate4Site and falls back to MSA entropy.",
    )
    parser.add_argument(
        "--reference-mode",
        choices=["first_record", "sequence_id"],
        default="first_record",
        help="How to choose the target/reference sequence used for residue numbering.",
    )
    parser.add_argument("--reference-id", default="", help="Reference sequence ID when --reference-mode sequence_id.")
    parser.add_argument("--chain", default="A", help="Chain label written to the output TSV.")
    parser.add_argument(
        "--accepted-aa-min-frequency",
        type=float,
        default=0.05,
        help="Minimum non-gap column frequency for an amino acid to be reported as accepted.",
    )
    parser.add_argument(
        "--gap-warning-threshold",
        type=float,
        default=0.5,
        help="Warn when a column has at least this gap fraction.",
    )
    parser.add_argument("--rate4site-binary", default=os.environ.get("RATE4SITE_BINARY", "rate4site"))
    return parser.parse_args()


def read_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: str | None = None
    parts: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(record_from_parts(header, parts))
                header = line[1:].strip()
                parts = []
            else:
                parts.append("".join(line.split()).upper())
    if header is not None:
        records.append(record_from_parts(header, parts))
    if len(records) < 2:
        raise ToolError("MSA must contain at least two sequences for conservation scoring.")
    lengths = {len(record.sequence) for record in records}
    if len(lengths) != 1:
        raise ToolError("All MSA sequences must have the same aligned length.")
    return records


def record_from_parts(header: str, parts: list[str]) -> FastaRecord:
    identifier = header.split()[0] if header else "sequence"
    sequence = "".join(parts)
    if not sequence:
        raise ToolError(f"FASTA record {identifier!r} is empty.")
    return FastaRecord(identifier=identifier, description=header, sequence=sequence)


def choose_reference(records: list[FastaRecord], args: argparse.Namespace) -> FastaRecord:
    if args.reference_mode == "first_record":
        return records[0]
    wanted = args.reference_id.strip()
    if not wanted:
        raise ToolError("reference_id is required when reference-mode is sequence_id.")
    for record in records:
        if record.identifier == wanted or record.description == wanted:
            return record
    raise ToolError(f"Reference sequence {wanted!r} was not found in the MSA.")


def run_rate4site(args: argparse.Namespace, raw_output: Path, log_lines: list[str]) -> dict[int, float]:
    binary = shutil.which(args.rate4site_binary) or args.rate4site_binary
    command = [binary, "-s", str(Path(args.msa_fasta)), "-o", str(raw_output)]
    log_lines.append("$ " + " ".join(command))
    try:
        completed = subprocess.run(command, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise ToolError(f"Rate4Site executable not found: {args.rate4site_binary}") from exc

    log_lines.append(f"rate4site_exit_code={completed.returncode}")
    if completed.stdout:
        log_lines.append("Rate4Site STDOUT:\n" + completed.stdout)
    if completed.stderr:
        log_lines.append("Rate4Site STDERR:\n" + completed.stderr)

    if completed.returncode != 0:
        raise ToolError("Rate4Site failed; see run log for details.")
    if not raw_output.exists() or raw_output.stat().st_size == 0:
        raise ToolError("Rate4Site completed but did not produce a non-empty output file.")
    parsed = parse_rate4site_output(raw_output)
    if not parsed:
        raise ToolError("Rate4Site output could not be parsed into per-site scores.")
    return parsed


def parse_rate4site_output(path: Path) -> dict[int, float]:
    """Parse common Rate4Site/ConSurf-like per-site score tables.

    The parser is intentionally permissive: it looks for rows beginning with a
    residue position and uses the first floating-point value after the amino-acid
    token as the evolutionary-rate score.
    """
    scores: dict[int, float] = {}
    number_re = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?$")
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            tokens = line.split()
            if not tokens:
                continue
            try:
                pos = int(tokens[0])
            except ValueError:
                continue
            numeric_values: list[float] = []
            for token in tokens[1:]:
                cleaned = token.strip("[](),")
                if number_re.match(cleaned):
                    try:
                        numeric_values.append(float(cleaned))
                    except ValueError:
                        pass
            if numeric_values:
                scores[pos] = numeric_values[0]
    return scores


def shannon_entropy_normalized(counts: Counter[str]) -> float:
    total = sum(counts.values())
    if total <= 0:
        return 0.0
    entropy = 0.0
    for count in counts.values():
        p = count / total
        if p > 0:
            entropy -= p * math.log(p)
    return entropy / math.log(20)


def score_sites(records: list[FastaRecord], reference: FastaRecord, args: argparse.Namespace) -> list[SiteStats]:
    sites: list[SiteStats] = []
    residue_index = 0
    alignment_length = len(reference.sequence)
    for column_index in range(alignment_length):
        wt_aa = reference.sequence[column_index]
        if wt_aa in GAP_CHARS:
            continue
        residue_index += 1
        column_values = [record.sequence[column_index].upper() for record in records]
        residues = [
            aa
            for aa in column_values
            if aa not in GAP_CHARS and (aa in AMINO_ACIDS or aa in AMBIGUOUS_AA)
        ]
        counts = Counter(residues)
        standard_counts = Counter({aa: count for aa, count in counts.items() if aa in AMINO_ACIDS})
        non_gap_count = sum(counts.values())
        gap_fraction = 1.0 - (non_gap_count / len(records))
        entropy = shannon_entropy_normalized(standard_counts)
        consensus_aa = ""
        if standard_counts:
            consensus_aa = sorted(standard_counts.items(), key=lambda item: (-item[1], item[0]))[0][0]
        accepted = accepted_amino_acids(standard_counts, args.accepted_aa_min_frequency)
        accepted_with_frequency = ",".join(
            f"{aa}:{standard_counts[aa] / max(1, sum(standard_counts.values())):.3f}" for aa in accepted
        )
        warning_parts: list[str] = []
        status = "ok"
        if wt_aa not in AMINO_ACIDS:
            warning_parts.append("reference_residue_nonstandard")
            status = "warning"
        if gap_fraction >= args.gap_warning_threshold:
            warning_parts.append("high_gap_fraction")
            status = "warning"
        if non_gap_count == 0:
            warning_parts.append("empty_column")
            status = "warning"
        sites.append(
            SiteStats(
                residue_index=residue_index,
                chain=args.chain,
                wt_aa=wt_aa,
                alignment_column=column_index + 1,
                msa_entropy=entropy,
                consensus_aa=consensus_aa,
                accepted_aas=accepted,
                accepted_aas_with_frequency=accepted_with_frequency,
                non_gap_count=non_gap_count,
                gap_fraction=gap_fraction,
                status=status,
                warning=";".join(warning_parts),
            )
        )
    if not sites:
        raise ToolError("The selected reference sequence contains no ungapped residues.")
    return sites


def accepted_amino_acids(counts: Counter[str], min_frequency: float) -> list[str]:
    total = sum(counts.values())
    if total <= 0:
        return []
    accepted = [aa for aa, count in counts.items() if count / total >= min_frequency]
    return sorted(accepted, key=lambda aa: (-counts[aa], aa))


def assign_scores(sites: list[SiteStats], rate4site_scores: dict[int, float] | None, log_lines: list[str]) -> None:
    if rate4site_scores:
        values = [score for index, score in rate4site_scores.items() if 1 <= index <= len(sites)]
        if values:
            min_value = min(values)
            max_value = max(values)
            span = max(max_value - min_value, 1e-12)
            log_lines.append(
                "Rate4Site score normalization: lower evolutionary rate is treated as more conserved."
            )
            for site in sites:
                score = rate4site_scores.get(site.residue_index)
                site.rate4site_score = score
                if score is None:
                    conservation = max(0.0, min(1.0, 1.0 - site.msa_entropy))
                    site.warning = append_warning(site.warning, "missing_rate4site_score")
                else:
                    normalized_rate = (score - min_value) / span
                    site.rate4site_score_normalized = normalized_rate
                    conservation = 1.0 - normalized_rate
                site.conservation_score = max(0.0, min(1.0, conservation))
                site.mutability_score = mutability_from_conservation(site.conservation_score, site.gap_fraction)
            return

    for site in sites:
        site.conservation_score = max(0.0, min(1.0, 1.0 - site.msa_entropy))
        site.mutability_score = mutability_from_conservation(site.conservation_score, site.gap_fraction)


def mutability_from_conservation(conservation_score: float, gap_fraction: float) -> float:
    # High gap columns are less reliable candidates even when they look variable.
    return max(0.0, min(1.0, (1.0 - conservation_score) * (1.0 - gap_fraction)))


def append_warning(existing: str, warning: str) -> str:
    if not existing:
        return warning
    parts = existing.split(";")
    if warning not in parts:
        parts.append(warning)
    return ";".join(parts)


def write_output(sites: Iterable[SiteStats], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for site in sites:
            writer.writerow(
                {
                    "residue_index": site.residue_index,
                    "chain": site.chain,
                    "wt_aa": site.wt_aa,
                    "alignment_column": site.alignment_column,
                    "conservation_score": format_float(site.conservation_score),
                    "mutability_score": format_float(site.mutability_score),
                    "msa_entropy": format_float(site.msa_entropy),
                    "consensus_aa": site.consensus_aa,
                    "accepted_aas": ",".join(site.accepted_aas),
                    "accepted_aas_with_frequency": site.accepted_aas_with_frequency,
                    "rate4site_score": format_float(site.rate4site_score),
                    "rate4site_score_normalized": format_float(site.rate4site_score_normalized),
                    "non_gap_count": site.non_gap_count,
                    "gap_fraction": format_float(site.gap_fraction),
                    "status": site.status,
                    "warning": site.warning,
                }
            )


def format_float(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.6g}"


def main() -> int:
    args = parse_args()
    log_lines = [
        "Conservation / Mutability Scorer",
        f"backend={args.backend}",
        f"reference_mode={args.reference_mode}",
    ]
    try:
        records = read_fasta(Path(args.msa_fasta))
        reference = choose_reference(records, args)
        log_lines.append(f"sequence_count={len(records)}")
        log_lines.append(f"alignment_length={len(reference.sequence)}")
        log_lines.append(f"reference_id={reference.identifier}")
        sites = score_sites(records, reference, args)

        rate4site_scores: dict[int, float] | None = None
        raw_output = Path(args.rate4site_output)
        if args.backend in {"rate4site", "auto"}:
            try:
                rate4site_scores = run_rate4site(args, raw_output, log_lines)
                log_lines.append(f"rate4site_sites_parsed={len(rate4site_scores)}")
            except ToolError as exc:
                if args.backend == "rate4site":
                    raise
                log_lines.append(f"WARNING: {exc}")
                log_lines.append("Falling back to entropy-only conservation scoring.")
                raw_output.write_text(
                    "Rate4Site was not used successfully. See run_log.txt for details.\n",
                    encoding="utf-8",
                )
        else:
            raw_output.write_text(
                "Rate4Site was not requested; entropy-only backend was used.\n",
                encoding="utf-8",
            )

        assign_scores(sites, rate4site_scores, log_lines)
        write_output(sites, Path(args.output_tsv))
    except ToolError as exc:
        Path(args.run_log).write_text("\n".join(log_lines) + f"\nERROR: {exc}\n", encoding="utf-8")
        print(f"conservation_mutability_scorer: {exc}", file=sys.stderr)
        return 1

    Path(args.run_log).write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
