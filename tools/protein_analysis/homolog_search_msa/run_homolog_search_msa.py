#!/usr/bin/env python3
"""Find homologous protein sequences, reduce redundancy, and build an MSA."""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


SEARCH_COLUMNS = [
    "query_id",
    "subject_id",
    "percent_identity",
    "alignment_length",
    "mismatches",
    "gap_openings",
    "query_start",
    "query_end",
    "subject_start",
    "subject_end",
    "evalue",
    "bitscore",
    "query_length",
    "subject_length",
    "query_coverage",
    "subject_coverage",
    "source",
]

SUMMARY_COLUMNS = [
    "sequence_id",
    "source",
    "best_query_id",
    "best_evalue",
    "best_bitscore",
    "best_percent_identity",
    "best_query_coverage",
    "best_subject_coverage",
    "best_alignment_length",
    "subject_length",
    "passed_search_filter",
    "in_filtered_homologs",
    "included_in_msa",
    "warning",
]


@dataclass
class FastaRecord:
    identifier: str
    description: str
    sequence: str


@dataclass
class SearchHit:
    query_id: str
    subject_id: str
    percent_identity: float | None
    alignment_length: int | None
    mismatches: int | None
    gap_openings: int | None
    query_start: int | None
    query_end: int | None
    subject_start: int | None
    subject_end: int | None
    evalue: float | None
    bitscore: float | None
    query_length: int | None
    subject_length: int | None
    query_coverage: float | None
    subject_coverage: float | None
    source: str


class ToolError(RuntimeError):
    """User-facing runtime error."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Find homologs with BLAST+/MMseqs2, remove redundancy, and build an MSA."
    )
    parser.add_argument("--query-fasta", required=True, help="Input query protein FASTA.")
    parser.add_argument("--subject-fasta", required=True, help="Subject/database FASTA containing candidate homologs.")
    parser.add_argument("--homologs-fasta", required=True, help="Output FASTA of homologs passing search filters.")
    parser.add_argument(
        "--filtered-homologs-fasta", required=True, help="Output FASTA after redundancy filtering."
    )
    parser.add_argument("--msa-fasta", required=True, help="Output multiple-sequence alignment FASTA.")
    parser.add_argument("--summary-tsv", required=True, help="Output homolog summary TSV.")
    parser.add_argument("--raw-search-output", required=True, help="Output standardized raw search TSV.")
    parser.add_argument("--cluster-report", required=True, help="Output clustering report.")
    parser.add_argument("--msa-log", required=True, help="Output MSA tool log.")
    parser.add_argument("--run-log", required=True, help="Output run log.")

    parser.add_argument(
        "--search-backend",
        choices=["blastp", "mmseqs", "all_subjects"],
        default="blastp",
        help="Homolog search backend. all_subjects treats the subject FASTA as a curated homolog set.",
    )
    parser.add_argument(
        "--filter-backend",
        choices=["cd-hit", "mmseqs", "internal"],
        default="cd-hit",
        help="Redundancy filtering backend.",
    )
    parser.add_argument(
        "--msa-backend",
        choices=["mafft", "muscle", "copy"],
        default="mafft",
        help="MSA backend. copy writes unaligned FASTA and is intended only for smoke tests or curated prealigned data.",
    )

    parser.add_argument("--evalue", type=float, default=1e-5, help="Maximum search E-value.")
    parser.add_argument("--max-target-seqs", type=int, default=500, help="Maximum targets to request from search backend.")
    parser.add_argument("--max-homologs", type=int, default=500, help="Maximum homolog records to keep after filtering hits.")
    parser.add_argument("--min-identity", type=float, default=20.0, help="Minimum percent identity for search hits.")
    parser.add_argument("--min-query-coverage", type=float, default=40.0, help="Minimum query coverage percentage.")
    parser.add_argument("--min-subject-coverage", type=float, default=0.0, help="Minimum subject coverage percentage.")
    parser.add_argument("--min-alignment-length", type=int, default=20, help="Minimum alignment length.")
    parser.add_argument("--identity-threshold", type=float, default=0.9, help="CD-HIT/MMseqs sequence identity threshold.")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads for external tools.")
    parser.add_argument("--include-query-in-msa", action="store_true", help="Include query sequences in the MSA input.")
    parser.add_argument("--exclude-self-hits", action="store_true", help="Drop hits where query id equals subject id.")

    parser.add_argument("--blastp-binary", default=os.environ.get("BLASTP_BINARY", "blastp"))
    parser.add_argument("--makeblastdb-binary", default=os.environ.get("MAKEBLASTDB_BINARY", "makeblastdb"))
    parser.add_argument("--mmseqs-binary", default=os.environ.get("MMSEQS_BINARY", "mmseqs"))
    parser.add_argument("--cdhit-binary", default=os.environ.get("CDHIT_BINARY", "cd-hit"))
    parser.add_argument("--mafft-binary", default=os.environ.get("MAFFT_BINARY", "mafft"))
    parser.add_argument("--muscle-binary", default=os.environ.get("MUSCLE_BINARY", "muscle"))
    return parser.parse_args()


def read_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    current_header: str | None = None
    sequence_parts: list[str] = []

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    records.append(record_from_parts(current_header, sequence_parts))
                current_header = line[1:].strip()
                sequence_parts = []
            else:
                sequence_parts.append("".join(line.split()))

    if current_header is not None:
        records.append(record_from_parts(current_header, sequence_parts))

    if not records:
        raise ToolError(f"No FASTA records found in {path}")
    return records


def record_from_parts(header: str, sequence_parts: list[str]) -> FastaRecord:
    identifier = header.split()[0] if header else "sequence"
    sequence = "".join(sequence_parts).upper()
    if not sequence:
        raise ToolError(f"FASTA record {identifier!r} has an empty sequence.")
    return FastaRecord(identifier=identifier, description=header, sequence=sequence)


def write_fasta(records: Iterable[FastaRecord], path: Path, width: int = 80) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            header = record.description or record.identifier
            if not header.startswith(record.identifier):
                header = record.identifier
            handle.write(f">{header}\n")
            for index in range(0, len(record.sequence), width):
                handle.write(record.sequence[index : index + width] + "\n")


def run_command(command: list[str], log_lines: list[str], stdout_path: Path | None = None) -> subprocess.CompletedProcess[str]:
    log_lines.append("$ " + " ".join(command))
    try:
        if stdout_path is None:
            completed = subprocess.run(command, capture_output=True, text=True)
        else:
            with stdout_path.open("w", encoding="utf-8", newline="\n") as stdout_handle:
                completed = subprocess.run(command, stdout=stdout_handle, stderr=subprocess.PIPE, text=True)
    except FileNotFoundError as exc:
        raise ToolError(f"Executable not found: {command[0]}. Configure the binary path or install the tool.") from exc
    log_lines.append(f"exit_code={completed.returncode}")
    if completed.stdout:
        log_lines.append("STDOUT:\n" + completed.stdout)
    if completed.stderr:
        log_lines.append("STDERR:\n" + completed.stderr)
    return completed


def run_blastp(args: argparse.Namespace, raw_path: Path, tmp_dir: Path, log_lines: list[str]) -> list[SearchHit]:
    db_prefix = tmp_dir / "blast_subject_db"
    makeblastdb_cmd = [
        args.makeblastdb_binary,
        "-in",
        str(Path(args.subject_fasta)),
        "-dbtype",
        "prot",
        "-out",
        str(db_prefix),
    ]
    completed = run_command(makeblastdb_cmd, log_lines)
    if completed.returncode != 0:
        raise ToolError("makeblastdb failed; see run log for details.")

    blast_tmp = tmp_dir / "blast_raw.tsv"
    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs"
    blast_cmd = [
        args.blastp_binary,
        "-query",
        str(Path(args.query_fasta)),
        "-db",
        str(db_prefix),
        "-out",
        str(blast_tmp),
        "-outfmt",
        outfmt,
        "-evalue",
        str(args.evalue),
        "-max_target_seqs",
        str(args.max_target_seqs),
        "-num_threads",
        str(args.threads),
    ]
    completed = run_command(blast_cmd, log_lines)
    if completed.returncode != 0:
        raise ToolError("blastp failed; see run log for details.")

    hits: list[SearchHit] = []
    with blast_tmp.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 15:
                continue
            alignment_length = parse_int(parts[3])
            subject_length = parse_int(parts[13])
            subject_coverage = None
            if alignment_length is not None and subject_length:
                subject_coverage = 100.0 * alignment_length / subject_length
            hits.append(
                SearchHit(
                    query_id=parts[0],
                    subject_id=parts[1],
                    percent_identity=parse_float(parts[2]),
                    alignment_length=alignment_length,
                    mismatches=parse_int(parts[4]),
                    gap_openings=parse_int(parts[5]),
                    query_start=parse_int(parts[6]),
                    query_end=parse_int(parts[7]),
                    subject_start=parse_int(parts[8]),
                    subject_end=parse_int(parts[9]),
                    evalue=parse_float(parts[10]),
                    bitscore=parse_float(parts[11]),
                    query_length=parse_int(parts[12]),
                    subject_length=subject_length,
                    query_coverage=parse_float(parts[14]),
                    subject_coverage=subject_coverage,
                    source="blastp",
                )
            )
    write_search_hits(hits, raw_path)
    return hits


def run_mmseqs_search(args: argparse.Namespace, raw_path: Path, tmp_dir: Path, log_lines: list[str]) -> list[SearchHit]:
    mmseqs_raw = tmp_dir / "mmseqs_search.tsv"
    mmseqs_tmp = tmp_dir / "mmseqs_tmp"
    format_output = ",".join(
        [
            "query",
            "target",
            "pident",
            "alnlen",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "tstart",
            "tend",
            "evalue",
            "bits",
            "qlen",
            "tlen",
            "qcov",
            "tcov",
        ]
    )
    command = [
        args.mmseqs_binary,
        "easy-search",
        str(Path(args.query_fasta)),
        str(Path(args.subject_fasta)),
        str(mmseqs_raw),
        str(mmseqs_tmp),
        "-e",
        str(args.evalue),
        "--threads",
        str(args.threads),
        "--max-seqs",
        str(args.max_target_seqs),
        "--format-output",
        format_output,
    ]
    completed = run_command(command, log_lines)
    if completed.returncode != 0:
        raise ToolError("mmseqs easy-search failed; see run log for details.")

    hits: list[SearchHit] = []
    with mmseqs_raw.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 16:
                continue
            pident = parse_float(parts[2])
            query_coverage = normalize_fraction_or_percent(parse_float(parts[14]))
            subject_coverage = normalize_fraction_or_percent(parse_float(parts[15]))
            hits.append(
                SearchHit(
                    query_id=parts[0],
                    subject_id=parts[1],
                    percent_identity=normalize_fraction_or_percent(pident),
                    alignment_length=parse_int(parts[3]),
                    mismatches=parse_int(parts[4]),
                    gap_openings=parse_int(parts[5]),
                    query_start=parse_int(parts[6]),
                    query_end=parse_int(parts[7]),
                    subject_start=parse_int(parts[8]),
                    subject_end=parse_int(parts[9]),
                    evalue=parse_float(parts[10]),
                    bitscore=parse_float(parts[11]),
                    query_length=parse_int(parts[12]),
                    subject_length=parse_int(parts[13]),
                    query_coverage=query_coverage,
                    subject_coverage=subject_coverage,
                    source="mmseqs",
                )
            )
    write_search_hits(hits, raw_path)
    return hits


def use_all_subjects(args: argparse.Namespace, raw_path: Path, subject_records: list[FastaRecord]) -> list[SearchHit]:
    hits = [
        SearchHit(
            query_id="curated_subject_set",
            subject_id=record.identifier,
            percent_identity=None,
            alignment_length=len(record.sequence),
            mismatches=None,
            gap_openings=None,
            query_start=None,
            query_end=None,
            subject_start=1,
            subject_end=len(record.sequence),
            evalue=None,
            bitscore=None,
            query_length=None,
            subject_length=len(record.sequence),
            query_coverage=None,
            subject_coverage=100.0,
            source="all_subjects",
        )
        for record in subject_records
    ]
    write_search_hits(hits, raw_path)
    return hits


def parse_float(value: str | None) -> float | None:
    if value is None or value in {"", "-", "nan", "NaN"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_int(value: str | None) -> int | None:
    parsed = parse_float(value)
    if parsed is None:
        return None
    return int(parsed)


def normalize_fraction_or_percent(value: float | None) -> float | None:
    if value is None:
        return None
    if 0 <= value <= 1:
        return value * 100.0
    return value


def write_search_hits(hits: list[SearchHit], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SEARCH_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for hit in hits:
            writer.writerow(hit_to_row(hit))


def hit_to_row(hit: SearchHit) -> dict[str, str | int | float | None]:
    return {
        "query_id": hit.query_id,
        "subject_id": hit.subject_id,
        "percent_identity": format_optional_float(hit.percent_identity),
        "alignment_length": hit.alignment_length,
        "mismatches": hit.mismatches,
        "gap_openings": hit.gap_openings,
        "query_start": hit.query_start,
        "query_end": hit.query_end,
        "subject_start": hit.subject_start,
        "subject_end": hit.subject_end,
        "evalue": format_optional_float(hit.evalue),
        "bitscore": format_optional_float(hit.bitscore),
        "query_length": hit.query_length,
        "subject_length": hit.subject_length,
        "query_coverage": format_optional_float(hit.query_coverage),
        "subject_coverage": format_optional_float(hit.subject_coverage),
        "source": hit.source,
    }


def format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.6g}"


def passes_filters(hit: SearchHit, args: argparse.Namespace) -> bool:
    if args.exclude_self_hits and hit.query_id == hit.subject_id:
        return False
    if hit.percent_identity is not None and hit.percent_identity < args.min_identity:
        return False
    if hit.query_coverage is not None and hit.query_coverage < args.min_query_coverage:
        return False
    if hit.subject_coverage is not None and hit.subject_coverage < args.min_subject_coverage:
        return False
    if hit.alignment_length is not None and hit.alignment_length < args.min_alignment_length:
        return False
    if hit.evalue is not None and hit.evalue > args.evalue:
        return False
    return True


def sort_hits(hits: list[SearchHit]) -> list[SearchHit]:
    return sorted(
        hits,
        key=lambda hit: (
            -(hit.bitscore if hit.bitscore is not None else -1.0),
            hit.evalue if hit.evalue is not None else float("inf"),
            -(hit.percent_identity if hit.percent_identity is not None else -1.0),
            hit.subject_id,
        ),
    )


def choose_best_hits(hits: list[SearchHit], args: argparse.Namespace) -> dict[str, SearchHit]:
    best: dict[str, SearchHit] = {}
    for hit in sort_hits([h for h in hits if passes_filters(h, args)]):
        if hit.subject_id not in best:
            best[hit.subject_id] = hit
        if len(best) >= args.max_homologs:
            break
    return best


def run_clustering(
    args: argparse.Namespace,
    homologs_path: Path,
    filtered_path: Path,
    cluster_report: Path,
    tmp_dir: Path,
    log_lines: list[str],
) -> None:
    homolog_records = read_fasta_or_empty(homologs_path)
    if not homolog_records:
        filtered_path.write_text("", encoding="utf-8")
        cluster_report.write_text("No homolog sequences passed search filters.\n", encoding="utf-8")
        return

    if args.filter_backend == "internal":
        internal_deduplicate(homolog_records, filtered_path, cluster_report)
        return

    if args.filter_backend == "cd-hit":
        word_size = cdhit_word_size(args.identity_threshold)
        command = [
            args.cdhit_binary,
            "-i",
            str(homologs_path),
            "-o",
            str(filtered_path),
            "-c",
            str(args.identity_threshold),
            "-n",
            str(word_size),
            "-T",
            str(args.threads),
            "-M",
            "0",
        ]
        completed = run_command(command, log_lines)
        if completed.returncode != 0:
            raise ToolError("CD-HIT failed; see run log for details.")
        clstr = Path(str(filtered_path) + ".clstr")
        if clstr.exists():
            shutil.copyfile(clstr, cluster_report)
        else:
            cluster_report.write_text("CD-HIT completed, but no .clstr file was produced.\n", encoding="utf-8")
        return

    if args.filter_backend == "mmseqs":
        prefix = tmp_dir / "mmseqs_cluster"
        mmseqs_tmp = tmp_dir / "mmseqs_cluster_tmp"
        command = [
            args.mmseqs_binary,
            "easy-cluster",
            str(homologs_path),
            str(prefix),
            str(mmseqs_tmp),
            "--min-seq-id",
            str(args.identity_threshold),
            "--threads",
            str(args.threads),
        ]
        completed = run_command(command, log_lines)
        if completed.returncode != 0:
            raise ToolError("mmseqs easy-cluster failed; see run log for details.")
        rep_fasta = Path(str(prefix) + "_rep_seq.fasta")
        cluster_tsv = Path(str(prefix) + "_cluster.tsv")
        if not rep_fasta.exists():
            raise ToolError("mmseqs easy-cluster completed but representative FASTA was not found.")
        shutil.copyfile(rep_fasta, filtered_path)
        if cluster_tsv.exists():
            shutil.copyfile(cluster_tsv, cluster_report)
        else:
            cluster_report.write_text("MMseqs clustering completed; cluster TSV was not produced.\n", encoding="utf-8")
        return

    raise ToolError(f"Unsupported filter backend: {args.filter_backend}")


def read_fasta_or_empty(path: Path) -> list[FastaRecord]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    return read_fasta(path)


def internal_deduplicate(records: list[FastaRecord], filtered_path: Path, cluster_report: Path) -> None:
    seen: dict[str, FastaRecord] = {}
    members: dict[str, list[str]] = {}
    for record in records:
        if record.sequence not in seen:
            seen[record.sequence] = record
            members[record.identifier] = [record.identifier]
        else:
            representative = seen[record.sequence].identifier
            members.setdefault(representative, [representative]).append(record.identifier)
    write_fasta(seen.values(), filtered_path)
    with cluster_report.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("representative_id\tmember_ids\n")
        for representative, member_ids in members.items():
            handle.write(f"{representative}\t{','.join(member_ids)}\n")


def cdhit_word_size(identity: float) -> int:
    if identity >= 0.7:
        return 5
    if identity >= 0.6:
        return 4
    if identity >= 0.5:
        return 3
    return 2


def build_msa_input(query_records: list[FastaRecord], filtered_records: list[FastaRecord], args: argparse.Namespace) -> list[FastaRecord]:
    if not args.include_query_in_msa:
        return filtered_records

    existing_ids = {record.identifier for record in filtered_records}
    msa_records: list[FastaRecord] = []
    for query in query_records:
        if query.identifier in existing_ids:
            msa_records.append(
                FastaRecord(
                    identifier=f"{query.identifier}_query",
                    description=f"{query.identifier}_query original_query",
                    sequence=query.sequence,
                )
            )
        else:
            msa_records.append(query)
    msa_records.extend(filtered_records)
    return msa_records


def run_msa(
    args: argparse.Namespace,
    msa_input_records: list[FastaRecord],
    msa_path: Path,
    msa_log: Path,
    tmp_dir: Path,
    log_lines: list[str],
) -> None:
    if not msa_input_records:
        msa_path.write_text("", encoding="utf-8")
        msa_log.write_text("No records available for MSA.\n", encoding="utf-8")
        return

    msa_input = tmp_dir / "msa_input.fasta"
    write_fasta(msa_input_records, msa_input)

    if args.msa_backend == "copy" or len(msa_input_records) == 1:
        shutil.copyfile(msa_input, msa_path)
        msa_log.write_text("MSA backend copy was used, or only one sequence was available.\n", encoding="utf-8")
        return

    if args.msa_backend == "mafft":
        command = [args.mafft_binary, "--auto", "--thread", str(args.threads), str(msa_input)]
        completed = run_command(command, log_lines, stdout_path=msa_path)
        msa_log.write_text(completed.stderr or "MAFFT completed.\n", encoding="utf-8")
        if completed.returncode != 0:
            raise ToolError("MAFFT failed; see MSA log and run log for details.")
        return

    if args.msa_backend == "muscle":
        v5_command = [
            args.muscle_binary,
            "-align",
            str(msa_input),
            "-output",
            str(msa_path),
            "-threads",
            str(args.threads),
        ]
        completed = run_command(v5_command, log_lines)
        if completed.returncode == 0 and msa_path.exists():
            msa_log.write_text(completed.stderr or completed.stdout or "MUSCLE completed.\n", encoding="utf-8")
            return

        v3_command = [args.muscle_binary, "-in", str(msa_input), "-out", str(msa_path)]
        completed = run_command(v3_command, log_lines)
        msa_log.write_text(completed.stderr or completed.stdout or "MUSCLE completed.\n", encoding="utf-8")
        if completed.returncode != 0:
            raise ToolError("MUSCLE failed with both v5 and v3 command syntaxes; see logs for details.")
        return

    raise ToolError(f"Unsupported MSA backend: {args.msa_backend}")


def write_summary(
    path: Path,
    subject_records: list[FastaRecord],
    best_hits: dict[str, SearchHit],
    filtered_ids: set[str],
    msa_ids: set[str],
    selected_ids: set[str],
) -> None:
    subject_by_id = {record.identifier: record for record in subject_records}
    # Avoid writing a row for every sequence in a large subject database; the
    # summary is intentionally limited to homologs selected for downstream use.
    all_ids = sorted(best_hits)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for sequence_id in all_ids:
            hit = best_hits.get(sequence_id)
            record = subject_by_id.get(sequence_id)
            warning = ""
            if hit is None:
                warning = "not_selected_or_failed_filters"
            if record is None:
                warning = "selected_hit_missing_from_subject_fasta"
            writer.writerow(
                {
                    "sequence_id": sequence_id,
                    "source": hit.source if hit else "",
                    "best_query_id": hit.query_id if hit else "",
                    "best_evalue": format_optional_float(hit.evalue if hit else None),
                    "best_bitscore": format_optional_float(hit.bitscore if hit else None),
                    "best_percent_identity": format_optional_float(hit.percent_identity if hit else None),
                    "best_query_coverage": format_optional_float(hit.query_coverage if hit else None),
                    "best_subject_coverage": format_optional_float(hit.subject_coverage if hit else None),
                    "best_alignment_length": hit.alignment_length if hit else "",
                    "subject_length": len(record.sequence) if record else (hit.subject_length if hit else ""),
                    "passed_search_filter": "yes" if sequence_id in selected_ids else "no",
                    "in_filtered_homologs": "yes" if sequence_id in filtered_ids else "no",
                    "included_in_msa": "yes" if sequence_id in msa_ids else "no",
                    "warning": warning,
                }
            )


def main() -> int:
    args = parse_args()
    log_lines: list[str] = [
        "Homolog Search and MSA",
        f"search_backend={args.search_backend}",
        f"filter_backend={args.filter_backend}",
        f"msa_backend={args.msa_backend}",
    ]

    try:
        query_records = read_fasta(Path(args.query_fasta))
        subject_records = read_fasta(Path(args.subject_fasta))
        subject_by_id: dict[str, FastaRecord] = {}
        for record in subject_records:
            if record.identifier not in subject_by_id:
                subject_by_id[record.identifier] = record
            else:
                log_lines.append(f"WARNING: duplicate subject id ignored after first occurrence: {record.identifier}")

        with tempfile.TemporaryDirectory(prefix="homolog_search_msa_") as tmp_name:
            tmp_dir = Path(tmp_name)
            raw_search_path = Path(args.raw_search_output)
            if args.search_backend == "blastp":
                hits = run_blastp(args, raw_search_path, tmp_dir, log_lines)
            elif args.search_backend == "mmseqs":
                hits = run_mmseqs_search(args, raw_search_path, tmp_dir, log_lines)
            else:
                hits = use_all_subjects(args, raw_search_path, subject_records)

            best_hits = choose_best_hits(hits, args)
            selected_ids = set(best_hits)
            homolog_records: list[FastaRecord] = []
            for sequence_id in sorted(
                selected_ids,
                key=lambda seq_id: (
                    -(best_hits[seq_id].bitscore if best_hits[seq_id].bitscore is not None else -1.0),
                    best_hits[seq_id].evalue if best_hits[seq_id].evalue is not None else float("inf"),
                    seq_id,
                ),
            ):
                record = subject_by_id.get(sequence_id)
                if record is None:
                    log_lines.append(f"WARNING: selected subject id not found in subject FASTA: {sequence_id}")
                    continue
                homolog_records.append(record)

            homologs_path = Path(args.homologs_fasta)
            write_fasta(homolog_records, homologs_path)
            log_lines.append(f"homologs_written={len(homolog_records)}")

            filtered_path = Path(args.filtered_homologs_fasta)
            run_clustering(args, homologs_path, filtered_path, Path(args.cluster_report), tmp_dir, log_lines)
            filtered_records = read_fasta_or_empty(filtered_path)
            filtered_ids = {record.identifier for record in filtered_records}
            log_lines.append(f"filtered_homologs_written={len(filtered_records)}")

            msa_records = build_msa_input(query_records, filtered_records, args)
            run_msa(args, msa_records, Path(args.msa_fasta), Path(args.msa_log), tmp_dir, log_lines)
            msa_ids = {record.identifier for record in msa_records}
            log_lines.append(f"msa_records={len(msa_records)}")

            write_summary(
                Path(args.summary_tsv),
                subject_records,
                best_hits,
                filtered_ids,
                msa_ids,
                selected_ids,
            )
    except ToolError as exc:
        Path(args.run_log).write_text("\n".join(log_lines) + f"\nERROR: {exc}\n", encoding="utf-8")
        print(f"homolog_search_msa: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:  # pragma: no cover - defensive guard for Galaxy users.
        Path(args.run_log).write_text("\n".join(log_lines) + f"\nUNEXPECTED ERROR: {exc}\n", encoding="utf-8")
        raise

    Path(args.run_log).write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
