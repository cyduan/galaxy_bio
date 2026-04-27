from __future__ import annotations

import csv
import math
import shutil
import subprocess
import uuid
from dataclasses import dataclass
from pathlib import Path

from Bio.Seq import Seq

from promoter_design_tools.fasta import FastaRecord, read_fasta


BASES = ["A", "C", "G", "T"]


@dataclass(frozen=True)
class MotifModel:
    name: str
    consensus: str
    pwm: list[dict[str, float]] | None = None

    @property
    def length(self) -> int:
        return len(self.pwm) if self.pwm is not None else len(self.consensus)


@dataclass(frozen=True)
class MotifHit:
    sequence_id: str
    motif_name: str
    start: int
    end: int
    strand: str
    matched_sequence: str
    raw_score: float
    normalized_score: float
    p_value: str
    q_value: str
    source: str


def model_from_consensus(name: str, consensus: str) -> MotifModel:
    consensus = consensus.upper().replace("U", "T")
    if not consensus or any(base not in BASES for base in consensus):
        raise ValueError(f"Consensus for {name} must contain only A/C/G/T.")
    return MotifModel(name=name, consensus=consensus)


def read_first_meme_motif(path: Path, fallback_name: str) -> MotifModel:
    motif_name = fallback_name
    matrix: list[dict[str, float]] = []
    reading = False
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                if reading and matrix:
                    break
                continue
            if line.startswith("MOTIF"):
                parts = line.split()
                if len(parts) >= 2:
                    motif_name = parts[1]
            elif line.startswith("letter-probability matrix"):
                reading = True
            elif reading:
                parts = line.split()
                if len(parts) < 4:
                    if matrix:
                        break
                    continue
                try:
                    values = [float(value) for value in parts[:4]]
                except ValueError:
                    if matrix:
                        break
                    continue
                matrix.append(dict(zip(BASES, values, strict=True)))
    if not matrix:
        raise ValueError(f"No MEME letter-probability matrix was found in {path}.")
    consensus = "".join(max(row, key=row.get) for row in matrix)
    return MotifModel(name=motif_name or fallback_name, consensus=consensus, pwm=matrix)


def score_window(sequence: str, model: MotifModel) -> tuple[float, float]:
    sequence = sequence.upper().replace("U", "T")
    if len(sequence) != model.length:
        raise ValueError("Window length does not match motif length.")
    if model.pwm is None:
        matches = sum(1 for observed, expected in zip(sequence, model.consensus, strict=True) if observed == expected)
        raw_score = float(matches)
        return raw_score, raw_score / model.length

    raw_score = 0.0
    max_score = 0.0
    for base, row in zip(sequence, model.pwm, strict=True):
        raw_score += row.get(base, 0.0)
        max_score += max(row.values())
    normalized = raw_score / max_score if max_score else 0.0
    return raw_score, normalized


def scan_record(record: FastaRecord, models: list[MotifModel], scan_reverse_complement: bool = False) -> list[MotifHit]:
    hits: list[MotifHit] = []
    sequence = record.sequence.upper().replace("U", "T")
    strands = [("+", sequence)]
    if scan_reverse_complement:
        strands.append(("-", str(Seq(sequence).reverse_complement())))

    for strand, scan_sequence in strands:
        for model in models:
            if len(scan_sequence) < model.length:
                continue
            for zero_based in range(0, len(scan_sequence) - model.length + 1):
                matched = scan_sequence[zero_based : zero_based + model.length]
                raw_score, normalized_score = score_window(matched, model)
                if strand == "+":
                    start = zero_based + 1
                    end = zero_based + model.length
                    matched_sequence = matched
                else:
                    start = len(sequence) - (zero_based + model.length) + 1
                    end = len(sequence) - zero_based
                    matched_sequence = sequence[start - 1 : end]
                hits.append(
                    MotifHit(
                        sequence_id=record.identifier,
                        motif_name=model.name,
                        start=start,
                        end=end,
                        strand=strand,
                        matched_sequence=matched_sequence,
                        raw_score=round(raw_score, 6),
                        normalized_score=round(normalized_score, 6),
                        p_value="NA",
                        q_value="NA",
                        source="simple",
                    )
                )
    return hits


def scan_simple(records: list[FastaRecord], models: list[MotifModel], scan_reverse_complement: bool = False) -> list[MotifHit]:
    hits: list[MotifHit] = []
    for record in records:
        hits.extend(scan_record(record, models, scan_reverse_complement))
    return hits


def write_hits_tsv(hits: list[MotifHit], path: Path) -> None:
    fieldnames = [
        "sequence_id",
        "motif_name",
        "start",
        "end",
        "strand",
        "matched_sequence",
        "raw_score",
        "normalized_score",
        "p_value",
        "q_value",
        "source",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for hit in hits:
            writer.writerow(hit.__dict__)


def run_fimo(fasta: Path, meme_file: Path, fimo_binary: str = "fimo") -> list[MotifHit]:
    if not meme_file.exists():
        raise FileNotFoundError(f"MEME motif file does not exist: {meme_file}")
    tmp_path = Path.cwd() / f"promoter_fimo_{uuid.uuid4().hex}"
    tmp_path.mkdir(parents=True, exist_ok=False)
    try:
        command = [fimo_binary, "--oc", str(tmp_path), str(meme_file), str(fasta)]
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
        if completed.returncode != 0:
            raise RuntimeError(
                "FIMO failed with exit code "
                f"{completed.returncode}\nSTDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
            )
        fimo_tsv = tmp_path / "fimo.tsv"
        if not fimo_tsv.exists():
            raise RuntimeError("FIMO completed but fimo.tsv was not created.")
        return parse_fimo_tsv(fimo_tsv)
    finally:
        shutil.rmtree(tmp_path, ignore_errors=True)


def parse_fimo_tsv(path: Path) -> list[MotifHit]:
    hits: list[MotifHit] = []
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t")
        for row in reader:
            if not row or row.get("motif_id") is None:
                continue
            raw_score = safe_float(row.get("score", "0"))
            hits.append(
                MotifHit(
                    sequence_id=row.get("sequence_name", ""),
                    motif_name=row.get("motif_id", ""),
                    start=int(float(row.get("start", "0"))),
                    end=int(float(row.get("stop", row.get("end", "0")))),
                    strand=row.get("strand", "+"),
                    matched_sequence=row.get("matched_sequence", ""),
                    raw_score=raw_score,
                    normalized_score=raw_score,
                    p_value=row.get("p-value", "NA"),
                    q_value=row.get("q-value", "NA"),
                    source="fimo",
                )
            )
    return hits


def safe_float(value: str | None) -> float:
    try:
        return float(value or 0)
    except ValueError:
        return math.nan
