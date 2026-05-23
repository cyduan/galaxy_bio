#!/usr/bin/env python3
"""Design candidate hotspot mutations and compact smart-library tables."""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path


AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
GAP_CHARS = {"-", "."}

AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}

PROPERTY_GROUPS = {
    "hydrophobic": set("AVILM"),
    "aromatic": set("FWY"),
    "polar": set("STNQ"),
    "positive": set("KRH"),
    "negative": set("DE"),
    "special": set("CGP"),
}

CONSERVATIVE_SUBSTITUTIONS = {
    "A": "VGS",
    "V": "ILA",
    "I": "LVM",
    "L": "IVM",
    "M": "LI",
    "F": "YW",
    "Y": "FW",
    "W": "FY",
    "S": "TAN",
    "T": "SAV",
    "N": "QDS",
    "Q": "NEK",
    "D": "EN",
    "E": "DQ",
    "K": "RHQ",
    "R": "KH",
    "H": "KRY",
    "C": "SA",
    "G": "AS",
    "P": "A",
}

STANDARD_CODE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

IUPAC_BASES = {
    "A": set("A"),
    "C": set("C"),
    "G": set("G"),
    "T": set("T"),
    "R": set("AG"),
    "Y": set("CT"),
    "S": set("GC"),
    "W": set("AT"),
    "K": set("GT"),
    "M": set("AC"),
    "B": set("CGT"),
    "D": set("AGT"),
    "H": set("ACT"),
    "V": set("ACG"),
    "N": set("ACGT"),
}


CANDIDATE_COLUMNS = [
    "position",
    "residue_id",
    "chain",
    "wt_aa",
    "suggested_aa",
    "reason",
    "msa_frequency",
    "mutation_type",
    "source",
    "hotspot_rank",
    "hotspot_score",
    "hotspot_class",
    "consensus_aa",
    "accepted_aas",
    "notes",
]

SMART_LIBRARY_COLUMNS = [
    "library_id",
    "design_type",
    "positions",
    "residue_ids",
    "wt_aas",
    "designed_aas",
    "mutation_count",
    "theoretical_variant_count",
    "recommended_degenerate_codons",
    "coverage_fraction",
    "hotspot_classes",
    "rationale",
    "notes",
]

DEGENERATE_COLUMNS = [
    "position",
    "residue_id",
    "wt_aa",
    "designed_aas",
    "degenerate_codon",
    "covered_aas",
    "exact_match",
    "extra_aas",
    "missing_aas",
    "codon_count",
    "coverage_fraction",
    "notes",
]


@dataclass
class FastaRecord:
    identifier: str
    description: str
    sequence: str


@dataclass
class Hotspot:
    rank: int
    residue_id: str
    chain: str
    position: int
    wt_aa: str
    hotspot_score: float
    hotspot_class: str
    recommendation: str
    consensus_aa: str = ""
    accepted_aas: list[str] = field(default_factory=list)
    raw: dict[str, str] = field(default_factory=dict)


@dataclass
class MutationCandidate:
    hotspot: Hotspot
    suggested_aa: str
    reasons: set[str] = field(default_factory=set)
    sources: set[str] = field(default_factory=set)
    msa_frequency: float | None = None
    notes: list[str] = field(default_factory=list)


class ToolError(RuntimeError):
    """User-facing CLI error."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Recommend hotspot substitutions from hotspot ranking, MSA variation, consensus residues, "
        "conservative replacements, alanine scanning, and reduced amino-acid alphabets."
    )
    parser.add_argument("--hotspots-tsv", required=True, help="Tool 6 hotspot_residue_ranked.tsv.")
    parser.add_argument("--msa-fasta", required=True, help="Tool 3 MSA FASTA.")
    parser.add_argument("--protein-sequence", required=True, help="Ungapped target protein FASTA.")
    parser.add_argument("--candidate-mutations", required=True, help="Output candidate_mutations.tsv.")
    parser.add_argument("--smart-library", required=True, help="Output smart_library.tsv.")
    parser.add_argument("--degenerate-codons", required=True, help="Output degenerate_codons.tsv.")
    parser.add_argument("--run-log", required=True, help="Run log output.")
    parser.add_argument("--top-n-hotspots", type=int, default=20, help="Maximum hotspot rows to design.")
    parser.add_argument("--min-hotspot-score", type=float, default=0.0, help="Minimum hotspot_score to include.")
    parser.add_argument("--min-msa-frequency", type=float, default=0.05, help="Minimum MSA frequency for natural AAs.")
    parser.add_argument("--max-suggestions-per-site", type=int, default=8, help="Maximum candidate AAs per site.")
    parser.add_argument("--smart-library-max-sites", type=int, default=6, help="Maximum sites in pooled smart library row.")
    parser.add_argument("--reduced-alphabet", default="ADKSTVFY", help="Reduced amino acid alphabet for exploratory designs.")
    parser.add_argument("--max-reduced-suggestions-per-site", type=int, default=2)
    parser.add_argument("--msa-reference-id", default="", help="Optional MSA record ID matching the target sequence.")
    parser.add_argument("--include-natural", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--include-consensus", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--include-conservative", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--include-alanine", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--include-reduced-alphabet", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--include-wt-in-library", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--avoid-cysteine", action=argparse.BooleanOptionalAction, default=True)
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
    if not records:
        raise ToolError(f"No FASTA records were found in {path}.")
    return records


def record_from_parts(header: str, parts: list[str]) -> FastaRecord:
    identifier = header.split()[0] if header else "sequence"
    sequence = "".join(parts).upper()
    if not sequence:
        raise ToolError(f"FASTA record {identifier!r} is empty.")
    return FastaRecord(identifier=identifier, description=header, sequence=sequence)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ToolError(f"{path} is empty or does not contain a TSV header.")
        return [{key: (value or "") for key, value in row.items()} for row in reader]


def parse_float(value: str, default: float = 0.0) -> float:
    try:
        if value in {"", "None", "NA", "nan"}:
            return default
        return float(value)
    except ValueError:
        return default


def parse_int(value: str) -> int | None:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def normalize_aa(value: str) -> str:
    token = (value or "").strip().upper()
    if len(token) == 1 and token in AMINO_ACIDS:
        return token
    if len(token) == 3 and token in AA3_TO_1:
        return AA3_TO_1[token]
    return ""


def parse_accepted_aas(value: str) -> list[str]:
    accepted: list[str] = []
    for part in re.split(r"[,;|\s]+", value or ""):
        token = part.strip()
        if not token:
            continue
        aa = normalize_aa(token.split(":", 1)[0])
        if aa and aa not in accepted:
            accepted.append(aa)
    return accepted


def parse_hotspots(path: Path, args: argparse.Namespace) -> list[Hotspot]:
    rows = read_tsv(path)
    hotspots: list[Hotspot] = []
    for index, row in enumerate(rows, start=1):
        rank = parse_int(row.get("rank", "")) or index
        position = (
            parse_int(row.get("residue_number", ""))
            or parse_int(row.get("position", ""))
            or parse_position_from_residue_id(row.get("residue_id", ""))
        )
        if position is None:
            continue
        hotspot_score = parse_float(row.get("hotspot_score", "0"), 0.0)
        if hotspot_score < args.min_hotspot_score:
            continue
        wt_aa = normalize_aa(row.get("wt_aa", "")) or normalize_aa(row.get("residue_name", ""))
        if not wt_aa:
            continue
        hotspots.append(
            Hotspot(
                rank=rank,
                residue_id=row.get("residue_id", "") or f"{row.get('chain', 'A')}:{position}",
                chain=row.get("chain", "A") or "A",
                position=position,
                wt_aa=wt_aa,
                hotspot_score=hotspot_score,
                hotspot_class=row.get("hotspot_class", ""),
                recommendation=row.get("recommendation", ""),
                consensus_aa=normalize_aa(row.get("consensus_aa", "")),
                accepted_aas=parse_accepted_aas(
                    row.get("accepted_aas_with_frequency", "") or row.get("accepted_aas", "")
                ),
                raw=row,
            )
        )
    hotspots.sort(key=lambda item: (item.rank, -item.hotspot_score, item.position))
    if args.top_n_hotspots > 0:
        hotspots = hotspots[: args.top_n_hotspots]
    if not hotspots:
        raise ToolError("No usable hotspot rows were found after filtering.")
    return hotspots


def parse_position_from_residue_id(value: str) -> int | None:
    match = re.search(r"(-?\d+)", value or "")
    if not match:
        return None
    return parse_int(match.group(1))


def choose_msa_reference(records: list[FastaRecord], protein_sequence: str, reference_id: str, log: list[str]) -> FastaRecord:
    if reference_id:
        for record in records:
            if record.identifier == reference_id or record.description == reference_id:
                log.append(f"MSA reference selected by ID: {record.identifier}")
                return record
        raise ToolError(f"MSA reference ID {reference_id!r} was not found.")

    ungapped_target = strip_gaps(protein_sequence)
    for record in records:
        if strip_gaps(record.sequence) == ungapped_target:
            log.append(f"MSA reference selected by exact sequence match: {record.identifier}")
            return record
    log.append(
        "WARNING: no MSA record exactly matched protein_sequence.fasta; using the first MSA record "
        f"({records[0].identifier}) as reference."
    )
    return records[0]


def strip_gaps(sequence: str) -> str:
    return "".join(char for char in sequence.upper() if char not in GAP_CHARS)


def validate_msa(records: list[FastaRecord]) -> None:
    lengths = {len(record.sequence) for record in records}
    if len(lengths) != 1:
        raise ToolError("All MSA sequences must have the same aligned length.")
    if len(records) < 2:
        raise ToolError("MSA must contain at least two sequences.")


def build_position_to_column(reference: FastaRecord) -> dict[int, int]:
    mapping: dict[int, int] = {}
    residue_index = 0
    for column_index, aa in enumerate(reference.sequence, start=1):
        if aa in GAP_CHARS:
            continue
        residue_index += 1
        mapping[residue_index] = column_index
    return mapping


def msa_column_frequencies(
    records: list[FastaRecord], position_to_column: dict[int, int]
) -> dict[int, dict[str, float]]:
    frequencies: dict[int, dict[str, float]] = {}
    for position, column in position_to_column.items():
        counter: Counter[str] = Counter()
        for record in records:
            aa = record.sequence[column - 1].upper()
            if aa in AMINO_ACIDS:
                counter[aa] += 1
        total = sum(counter.values())
        frequencies[position] = {aa: count / total for aa, count in counter.items()} if total else {}
    return frequencies


def add_candidate(
    candidates: dict[tuple[int, str], MutationCandidate],
    hotspot: Hotspot,
    suggested_aa: str,
    reason: str,
    source: str,
    frequencies: dict[int, dict[str, float]],
    note: str = "",
) -> None:
    aa = normalize_aa(suggested_aa)
    if not aa or aa == hotspot.wt_aa:
        return
    key = (hotspot.position, aa)
    candidate = candidates.get(key)
    if candidate is None:
        candidate = MutationCandidate(hotspot=hotspot, suggested_aa=aa)
        candidates[key] = candidate
    candidate.reasons.add(reason)
    candidate.sources.add(source)
    if note:
        candidate.notes.append(note)
    freq = frequencies.get(hotspot.position, {}).get(aa)
    if freq is not None:
        candidate.msa_frequency = max(candidate.msa_frequency or 0.0, freq)


def same_property_group(left: str, right: str) -> bool:
    return any(left in group and right in group for group in PROPERTY_GROUPS.values())


def mutation_type(hotspot: Hotspot, suggested_aa: str, reasons: set[str]) -> str:
    if "alanine_scanning" in reasons:
        return "alanine_scan"
    if "physicochemical_conservative" in reasons or same_property_group(hotspot.wt_aa, suggested_aa):
        return "conservative"
    if "consensus" in reasons or "natural_observed" in reasons or "accepted_aa" in reasons:
        return "natural_or_consensus"
    return "exploratory"


def source_priority(candidate: MutationCandidate) -> tuple[int, int, float, str]:
    priorities = {
        "natural_observed": 0,
        "consensus": 1,
        "accepted_aa": 2,
        "conservative": 3,
        "alanine_scanning": 4,
        "reduced_alphabet": 5,
    }
    best = min((priorities.get(source, 9) for source in candidate.sources), default=9)
    return (
        candidate.hotspot.rank,
        best,
        -(candidate.msa_frequency or 0.0),
        candidate.suggested_aa,
    )


def design_candidates(
    hotspots: list[Hotspot],
    frequencies: dict[int, dict[str, float]],
    args: argparse.Namespace,
    log: list[str],
) -> list[MutationCandidate]:
    candidates: dict[tuple[int, str], MutationCandidate] = {}
    reduced_alphabet = [aa for aa in args.reduced_alphabet.upper() if aa in AMINO_ACIDS]
    if args.avoid_cysteine:
        reduced_alphabet = [aa for aa in reduced_alphabet if aa != "C"]

    for hotspot in hotspots:
        observed = frequencies.get(hotspot.position, {})
        if args.include_natural:
            for aa, freq in sorted(observed.items(), key=lambda item: (-item[1], item[0])):
                if aa == hotspot.wt_aa or freq < args.min_msa_frequency:
                    continue
                if args.avoid_cysteine and aa == "C":
                    continue
                add_candidate(
                    candidates,
                    hotspot,
                    aa,
                    "natural_observed",
                    "natural_observed",
                    frequencies,
                    f"MSA frequency {freq:.3f}",
                )
            for aa in hotspot.accepted_aas:
                if args.avoid_cysteine and aa == "C":
                    continue
                add_candidate(candidates, hotspot, aa, "accepted_aa", "accepted_aa", frequencies)

        if args.include_consensus and hotspot.consensus_aa:
            if not (args.avoid_cysteine and hotspot.consensus_aa == "C"):
                add_candidate(candidates, hotspot, hotspot.consensus_aa, "consensus", "consensus", frequencies)

        if args.include_conservative:
            for aa in CONSERVATIVE_SUBSTITUTIONS.get(hotspot.wt_aa, ""):
                if args.avoid_cysteine and aa == "C":
                    continue
                add_candidate(candidates, hotspot, aa, "physicochemical_conservative", "conservative", frequencies)

        if args.include_alanine:
            add_candidate(candidates, hotspot, "A", "alanine_scanning", "alanine_scanning", frequencies)

        if args.include_reduced_alphabet and reduced_alphabet:
            added = 0
            for aa in reduced_alphabet:
                if aa == hotspot.wt_aa:
                    continue
                add_candidate(candidates, hotspot, aa, "reduced_alphabet", "reduced_alphabet", frequencies)
                added += 1
                if added >= args.max_reduced_suggestions_per_site:
                    break

    grouped: dict[int, list[MutationCandidate]] = defaultdict(list)
    for candidate in candidates.values():
        grouped[candidate.hotspot.position].append(candidate)

    selected: list[MutationCandidate] = []
    for position, site_candidates in grouped.items():
        site_candidates.sort(key=source_priority)
        selected.extend(site_candidates[: args.max_suggestions_per_site])
        if len(site_candidates) > args.max_suggestions_per_site:
            log.append(
                f"Position {position}: kept {args.max_suggestions_per_site} of {len(site_candidates)} suggestions."
            )
    selected.sort(key=source_priority)
    return selected


def write_candidates(path: Path, candidates: list[MutationCandidate]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CANDIDATE_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for candidate in candidates:
            hotspot = candidate.hotspot
            writer.writerow(
                {
                    "position": hotspot.position,
                    "residue_id": hotspot.residue_id,
                    "chain": hotspot.chain,
                    "wt_aa": hotspot.wt_aa,
                    "suggested_aa": candidate.suggested_aa,
                    "reason": "/".join(sorted(candidate.reasons)),
                    "msa_frequency": format_float(candidate.msa_frequency),
                    "mutation_type": mutation_type(hotspot, candidate.suggested_aa, candidate.reasons),
                    "source": ",".join(sorted(candidate.sources)),
                    "hotspot_rank": hotspot.rank,
                    "hotspot_score": format_float(hotspot.hotspot_score),
                    "hotspot_class": hotspot.hotspot_class,
                    "consensus_aa": hotspot.consensus_aa,
                    "accepted_aas": ",".join(hotspot.accepted_aas),
                    "notes": "; ".join(dict.fromkeys(candidate.notes)),
                }
            )


def format_float(value: float | None) -> str:
    if value is None:
        return ""
    if math.isnan(value):
        return ""
    return f"{value:.6g}"


def make_library_rows(candidates: list[MutationCandidate], args: argparse.Namespace) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    by_position: dict[int, list[MutationCandidate]] = defaultdict(list)
    for candidate in candidates:
        by_position[candidate.hotspot.position].append(candidate)

    site_rows: list[dict[str, str]] = []
    degenerate_rows: list[dict[str, str]] = []

    for position in sorted(by_position, key=lambda pos: by_position[pos][0].hotspot.rank):
        site_candidates = sorted(by_position[position], key=source_priority)
        hotspot = site_candidates[0].hotspot
        designed_aas = sorted({candidate.suggested_aa for candidate in site_candidates})
        codon_targets = list(designed_aas)
        if args.include_wt_in_library and hotspot.wt_aa not in codon_targets:
            codon_targets.append(hotspot.wt_aa)
        codon_design = choose_degenerate_codon(set(codon_targets))
        degenerate_rows.append(degenerate_row(hotspot, sorted(codon_targets), codon_design))
        site_rows.append(
            {
                "library_id": f"site_{hotspot.chain}{hotspot.position}",
                "design_type": "single_site",
                "positions": str(hotspot.position),
                "residue_ids": hotspot.residue_id,
                "wt_aas": hotspot.wt_aa,
                "designed_aas": ",".join(codon_targets),
                "mutation_count": str(len(designed_aas)),
                "theoretical_variant_count": str(len(set(codon_design["covered_aas"]))),
                "recommended_degenerate_codons": codon_design["codon"],
                "coverage_fraction": format_float(codon_design["coverage_fraction"]),
                "hotspot_classes": hotspot.hotspot_class,
                "rationale": summarize_site_reasons(site_candidates),
                "notes": codon_design["notes"],
            }
        )

    pooled = make_pooled_library_row(site_rows, args.smart_library_max_sites)
    if pooled:
        site_rows.insert(0, pooled)
    return site_rows, degenerate_rows


def summarize_site_reasons(candidates: list[MutationCandidate]) -> str:
    reasons: list[str] = []
    for candidate in candidates:
        reasons.extend(sorted(candidate.reasons))
    return ",".join(dict.fromkeys(reasons))


def make_pooled_library_row(site_rows: list[dict[str, str]], max_sites: int) -> dict[str, str] | None:
    if not site_rows or max_sites <= 1:
        return None
    selected = site_rows[:max_sites]
    variant_count = 1
    for row in selected:
        variant_count *= max(1, parse_int(row.get("theoretical_variant_count", "1")) or 1)
    return {
        "library_id": f"pooled_top_{len(selected)}_sites",
        "design_type": "pooled_multi_site",
        "positions": ";".join(row["positions"] for row in selected),
        "residue_ids": ";".join(row["residue_ids"] for row in selected),
        "wt_aas": ";".join(row["wt_aas"] for row in selected),
        "designed_aas": ";".join(row["designed_aas"] for row in selected),
        "mutation_count": str(sum(parse_int(row["mutation_count"]) or 0 for row in selected)),
        "theoretical_variant_count": str(variant_count),
        "recommended_degenerate_codons": ";".join(row["recommended_degenerate_codons"] for row in selected),
        "coverage_fraction": "",
        "hotspot_classes": ";".join(row["hotspot_classes"] for row in selected),
        "rationale": "pooled design for top ranked hotspot sites",
        "notes": "Variant count is multiplicative across included sites; split into smaller libraries if too large.",
    }


def choose_degenerate_codon(target_aas: set[str]) -> dict[str, object]:
    if not target_aas:
        return {
            "codon": "",
            "covered_aas": [],
            "extra_aas": [],
            "missing_aas": [],
            "exact_match": "no",
            "codon_count": 0,
            "coverage_fraction": 0.0,
            "notes": "No target amino acids.",
        }

    best: dict[str, object] | None = None
    symbols = list(IUPAC_BASES)
    for codon in ("".join(parts) for parts in itertools.product(symbols, repeat=3)):
        concrete = expand_degenerate_codon(codon)
        aas = {STANDARD_CODE[triplet] for triplet in concrete}
        if "*" in aas:
            continue
        covered = aas & target_aas
        missing = target_aas - aas
        extra = aas - target_aas
        if not covered:
            continue
        score = (
            len(missing),
            len(extra),
            len(concrete),
            len(codon.replace("N", "")),
            codon,
        )
        candidate = {
            "codon": codon,
            "covered_aas": sorted(aas),
            "extra_aas": sorted(extra),
            "missing_aas": sorted(missing),
            "exact_match": "yes" if not extra and not missing else "no",
            "codon_count": len(concrete),
            "coverage_fraction": len(covered) / len(target_aas),
            "score": score,
        }
        if best is None or candidate["score"] < best["score"]:  # type: ignore[index]
            best = candidate

    if best is None:
        return {
            "codon": "",
            "covered_aas": [],
            "extra_aas": [],
            "missing_aas": sorted(target_aas),
            "exact_match": "no",
            "codon_count": 0,
            "coverage_fraction": 0.0,
            "notes": "No stop-free degenerate codon could cover the target set.",
        }

    notes = []
    if best["extra_aas"]:
        notes.append("degenerate codon introduces extra amino acids")
    if best["missing_aas"]:
        notes.append("degenerate codon does not cover all requested amino acids")
    return {
        **best,
        "notes": "; ".join(notes),
    }


def expand_degenerate_codon(codon: str) -> list[str]:
    bases = [IUPAC_BASES[base] for base in codon]
    return ["".join(parts) for parts in itertools.product(*bases)]


def degenerate_row(hotspot: Hotspot, designed_aas: list[str], design: dict[str, object]) -> dict[str, str]:
    return {
        "position": str(hotspot.position),
        "residue_id": hotspot.residue_id,
        "wt_aa": hotspot.wt_aa,
        "designed_aas": ",".join(designed_aas),
        "degenerate_codon": str(design["codon"]),
        "covered_aas": ",".join(design["covered_aas"]),  # type: ignore[arg-type]
        "exact_match": str(design["exact_match"]),
        "extra_aas": ",".join(design["extra_aas"]),  # type: ignore[arg-type]
        "missing_aas": ",".join(design["missing_aas"]),  # type: ignore[arg-type]
        "codon_count": str(design["codon_count"]),
        "coverage_fraction": format_float(float(design["coverage_fraction"])),
        "notes": str(design["notes"]),
    }


def write_rows(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    args = parse_args()
    log: list[str] = []
    try:
        hotspots = parse_hotspots(Path(args.hotspots_tsv), args)
        msa_records = read_fasta(Path(args.msa_fasta))
        validate_msa(msa_records)
        protein_records = read_fasta(Path(args.protein_sequence))
        protein_sequence = strip_gaps(protein_records[0].sequence)
        reference = choose_msa_reference(msa_records, protein_sequence, args.msa_reference_id.strip(), log)
        position_to_column = build_position_to_column(reference)
        frequencies = msa_column_frequencies(msa_records, position_to_column)

        for hotspot in hotspots:
            if hotspot.position <= len(protein_sequence) and protein_sequence[hotspot.position - 1] != hotspot.wt_aa:
                log.append(
                    f"WARNING: hotspot {hotspot.residue_id} has wt_aa={hotspot.wt_aa}, "
                    f"but protein_sequence position {hotspot.position} is {protein_sequence[hotspot.position - 1]}."
                )
            if hotspot.position not in frequencies:
                log.append(f"WARNING: no MSA frequency could be mapped for position {hotspot.position}.")

        candidates = design_candidates(hotspots, frequencies, args, log)
        library_rows, degenerate_rows = make_library_rows(candidates, args)

        write_candidates(Path(args.candidate_mutations), candidates)
        write_rows(Path(args.smart_library), SMART_LIBRARY_COLUMNS, library_rows)
        write_rows(Path(args.degenerate_codons), DEGENERATE_COLUMNS, degenerate_rows)

        log.extend(
            [
                f"hotspots_used={len(hotspots)}",
                f"candidate_mutations={len(candidates)}",
                f"smart_library_rows={len(library_rows)}",
                f"degenerate_codon_rows={len(degenerate_rows)}",
                f"msa_sequences={len(msa_records)}",
                f"msa_reference={reference.identifier}",
            ]
        )
        Path(args.run_log).write_text("\n".join(log) + "\n", encoding="utf-8")
        return 0
    except ToolError as exc:
        Path(args.run_log).write_text("\n".join(log + [f"ERROR: {exc}"]) + "\n", encoding="utf-8")
        print(f"mutation_designer: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
