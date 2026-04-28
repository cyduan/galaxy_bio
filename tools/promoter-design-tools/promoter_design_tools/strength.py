from __future__ import annotations

import csv
import math
import os
import shutil
import shlex
import subprocess
import uuid
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from promoter_design_tools.fasta import FastaRecord, write_fasta
from promoter_design_tools.motifs import model_from_consensus, scan_record
from promoter_design_tools.pairing import spacer_score


@dataclass(frozen=True)
class StrengthPrediction:
    sequence_id: str
    predicted_strength: float
    log10_strength: float
    backend: str
    model_version: str
    warning: str


def gc_content(sequence: str) -> float:
    sequence = sequence.upper()
    return (sequence.count("G") + sequence.count("C")) / len(sequence)


def at_fraction(sequence: str) -> float:
    sequence = sequence.upper()
    return (sequence.count("A") + sequence.count("T")) / len(sequence) if sequence else 0.0


def predict_heuristic(records: list[FastaRecord], tss_offset: int = 60) -> list[StrengthPrediction]:
    minus35 = model_from_consensus("minus35", "TTGACA")
    minus10 = model_from_consensus("minus10", "TATAAT")

    predictions: list[StrengthPrediction] = []
    for record in records:
        sequence = record.sequence.upper()
        pair = best_pair_for_record(record, minus35, minus10)
        conservation = pair["conservation_score"] if pair else 0.0
        gc_score = max(0.0, 1.0 - abs(gc_content(sequence) - 0.45) / 0.45)
        if pair:
            minus10_start = int(pair["minus10_start"])
            window = sequence[max(0, minus10_start - 4) : min(len(sequence), minus10_start + 8)]
        else:
            center = max(0, min(len(sequence), tss_offset - 10))
            window = sequence[max(0, center - 6) : min(len(sequence), center + 6)]
        at_score = at_fraction(window)
        score_0_1 = 0.55 * conservation + 0.20 * gc_score + 0.20 * at_score + 0.05
        predicted = max(0.001, round(100.0 * min(1.0, score_0_1), 6))
        predictions.append(
            StrengthPrediction(
                sequence_id=record.identifier,
                predicted_strength=predicted,
                log10_strength=round(math.log10(predicted), 6),
                backend="heuristic",
                model_version="heuristic_sigma70_v0.1",
                warning="relative_score_not_absolute_transcription_rate",
            )
        )
    return predictions


def best_pair_for_record(record: FastaRecord, minus35, minus10) -> dict | None:
    hits = scan_record(record, [minus35, minus10], False)
    left_hits = [hit for hit in hits if hit.motif_name == "minus35"]
    right_hits = [hit for hit in hits if hit.motif_name == "minus10"]
    best: dict | None = None
    for left in left_hits:
        for right in right_hits:
            if left.start >= right.start:
                continue
            spacer = right.start - left.end - 1
            s_score = spacer_score(spacer)
            conservation = 0.4 * left.normalized_score + 0.4 * right.normalized_score + 0.2 * s_score
            candidate = {
                "minus10_start": right.start,
                "conservation_score": conservation,
            }
            if best is None or conservation > best["conservation_score"]:
                best = candidate
    return best


def predict_promoter_calculator(
    records: list[FastaRecord],
    organism: str,
    threads: int,
    command_template: str | None = None,
) -> list[StrengthPrediction]:
    command_template = command_template or os.environ.get("PROMOTER_CALCULATOR_COMMAND", "")
    if not command_template.strip():
        raise RuntimeError(
            "backend=promoter_calculator requires PROMOTER_CALCULATOR_COMMAND or --promoter-calculator-command. "
            "Use placeholders {input_fasta}, {output_tsv}, {organism}, and {threads}."
        )
    tmp_path = Path.cwd() / f"promoter_calculator_{uuid.uuid4().hex}"
    tmp_path.mkdir(parents=True, exist_ok=False)
    try:
        predictions: list[StrengthPrediction] = []
        for index, record in enumerate(records, start=1):
            # The Barrick Lab CLI reads FASTA files as one concatenated sequence,
            # so run it once per record to preserve Galaxy sequence IDs.
            safe_id = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in record.identifier) or f"record_{index}"
            input_fasta = tmp_path / f"{safe_id}.fasta"
            output_tsv = tmp_path / f"{safe_id}.csv"
            write_fasta([record], input_fasta)
            rendered = command_template.format(
                input_fasta=input_fasta,
                output_tsv=output_tsv,
                organism=organism,
                threads=threads,
            )
            completed = subprocess.run(shlex.split(rendered, posix=(os.name != "nt")), capture_output=True, text=True)
            if completed.returncode != 0:
                raise RuntimeError(
                    f"Promoter Calculator command failed for {record.identifier} with exit code "
                    f"{completed.returncode}\nSTDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
                )
            if not output_tsv.exists():
                if completed.stdout.strip():
                    output_tsv.write_text(completed.stdout, encoding="utf-8")
                else:
                    predictions.append(
                        StrengthPrediction(
                            sequence_id=record.identifier,
                            predicted_strength=0.0,
                            log10_strength=round(math.log10(1e-9), 6),
                            backend="promoter_calculator",
                            model_version="barricklab_promotercalculator",
                            warning="no_promoter_calculator_output",
                        )
                    )
                    continue
            predictions.extend(parse_promoter_calculator_output(output_tsv, fallback_sequence_id=record.identifier))
        return predictions
    finally:
        shutil.rmtree(tmp_path, ignore_errors=True)


def parse_promoter_calculator_output(path: Path, fallback_sequence_id: str | None = None) -> list[StrengthPrediction]:
    df = pd.read_csv(path, sep=None, engine="python")
    df.columns = [str(column).strip() for column in df.columns]
    if "sequence_id" not in df.columns:
        for candidate in ["id", "name", "promoter", "variant"]:
            if candidate in df.columns:
                df = df.rename(columns={candidate: "sequence_id"})
                break
    if "predicted_strength" not in df.columns:
        lower_to_original = {str(column).lower(): column for column in df.columns}
        for candidate in [
            "predicted_strength",
            "strength",
            "expression",
            "tx_rate",
            "tx rate",
            "txrate",
            "transcription_rate",
            "transcription initiation rate",
            "transcription_initiation_rate",
        ]:
            original_column = lower_to_original.get(candidate)
            if original_column is not None:
                df = df.rename(columns={original_column: "predicted_strength"})
                break
    if "predicted_strength" in df.columns and fallback_sequence_id and "sequence_id" not in df.columns:
        df.insert(0, "sequence_id", fallback_sequence_id)
    if "sequence_id" not in df.columns or "predicted_strength" not in df.columns:
        raise RuntimeError(
            "Promoter Calculator output must contain sequence_id and predicted_strength columns, "
            "or a supported strength column such as Tx_rate."
        )
    if fallback_sequence_id:
        df["sequence_id"] = fallback_sequence_id
        df["predicted_strength"] = pd.to_numeric(df["predicted_strength"], errors="coerce")
        df = df.dropna(subset=["predicted_strength"])
        if df.empty:
            return [
                StrengthPrediction(
                    sequence_id=fallback_sequence_id,
                    predicted_strength=0.0,
                    log10_strength=round(math.log10(1e-9), 6),
                    backend="promoter_calculator",
                    model_version="barricklab_promotercalculator",
                    warning="no_numeric_tx_rate_found",
                )
            ]
        # One Galaxy input sequence should yield one row. Use the strongest
        # predicted promoter hit reported by the external calculator.
        df = df.sort_values("predicted_strength", ascending=False).head(1)
    predictions: list[StrengthPrediction] = []
    for _, row in df.iterrows():
        predicted = float(row["predicted_strength"])
        predictions.append(
            StrengthPrediction(
                sequence_id=str(row["sequence_id"]),
                predicted_strength=round(predicted, 6),
                log10_strength=round(math.log10(max(predicted, 1e-9)), 6),
                backend="promoter_calculator",
                model_version=str(row.get("model_version", "barricklab_promotercalculator")),
                warning=str(row.get("warning", "")),
            )
        )
    return predictions


def write_strength(predictions: list[StrengthPrediction], path: Path) -> None:
    fieldnames = list(StrengthPrediction.__dataclass_fields__.keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for prediction in predictions:
            writer.writerow(prediction.__dict__)
