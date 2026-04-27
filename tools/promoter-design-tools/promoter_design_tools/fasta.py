from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


VALID_DNA = set("ACGTN")


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    sequence: str
    description: str = ""


@dataclass(frozen=True)
class NormalizationResult:
    record: FastaRecord
    original_length: int
    normalized_length: int
    status: str
    warning: str


def read_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    for record in SeqIO.parse(str(path), "fasta"):
        sequence = str(record.seq).strip().upper().replace("U", "T")
        if not sequence:
            raise ValueError(f"FASTA record '{record.id}' is empty.")
        records.append(FastaRecord(record.id, sequence, record.description))
    if not records:
        raise ValueError(f"No FASTA records were found in {path}.")
    return records


def write_fasta(records: list[FastaRecord], path: Path, line_width: int = 80) -> None:
    seq_records = [
        SeqRecord(Seq(record.sequence), id=record.identifier, description="")
        for record in records
    ]
    with path.open("w", encoding="utf-8") as handle:
        SeqIO.write(seq_records, handle, "fasta")


def sanitize_dna(sequence: str, padding_base: str = "N") -> tuple[str, str]:
    padding_base = padding_base.upper()
    if padding_base not in VALID_DNA:
        raise ValueError("padding_base must be one of A, C, G, T, or N.")
    warnings: list[str] = []
    cleaned: list[str] = []
    for base in sequence.upper().replace("U", "T"):
        if base in VALID_DNA:
            cleaned.append(base)
        elif base.isspace():
            continue
        else:
            cleaned.append(padding_base)
            warnings.append(f"invalid_base_{base}_replaced")
    return "".join(cleaned), ";".join(sorted(set(warnings)))


def normalize_record(
    record: FastaRecord,
    target_length: int = 81,
    tss_offset: int = 60,
    padding_base: str = "N",
) -> NormalizationResult:
    if target_length <= 0:
        raise ValueError("target_length must be positive.")
    if tss_offset < 0:
        raise ValueError("tss_offset must be non-negative.")

    cleaned, warning = sanitize_dna(record.sequence, padding_base)
    original_length = len(cleaned)
    status_parts: list[str] = []
    if warning:
        status_parts.append("warning")

    if original_length > target_length:
        normalized = cleaned[:target_length]
        status_parts.append("truncated")
    elif original_length < target_length:
        normalized = cleaned + padding_base.upper() * (target_length - original_length)
        status_parts.append("padded")
    else:
        normalized = cleaned
        status_parts.append("ok")

    if tss_offset >= target_length:
        warning = append_warning(warning, "tss_offset_outside_normalized_window")
    if tss_offset >= original_length:
        warning = append_warning(warning, "tss_offset_outside_original_sequence")

    status = "ok" if status_parts == ["ok"] and not warning else ";".join(status_parts)
    return NormalizationResult(
        record=FastaRecord(record.identifier, normalized, record.description),
        original_length=original_length,
        normalized_length=len(normalized),
        status=status,
        warning=warning,
    )


def append_warning(existing: str, new_warning: str) -> str:
    if not existing:
        return new_warning
    warnings = set(existing.split(";"))
    warnings.add(new_warning)
    return ";".join(sorted(warnings))


def normalize_records(
    records: list[FastaRecord],
    target_length: int = 81,
    tss_offset: int = 60,
    padding_base: str = "N",
) -> list[NormalizationResult]:
    return [
        normalize_record(record, target_length, tss_offset, padding_base)
        for record in records
    ]

