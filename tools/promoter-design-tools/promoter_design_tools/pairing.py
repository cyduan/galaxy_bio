from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from promoter_design_tools.fasta import FastaRecord


@dataclass(frozen=True)
class PairScore:
    sequence_id: str
    minus35_seq: str
    minus35_start: str
    minus35_end: str
    minus35_score: str
    minus10_seq: str
    minus10_start: str
    minus10_end: str
    minus10_score: str
    spacer_len: str
    spacer_score: str
    conservation_score: str
    status: str


def spacer_score(spacer: int, min_spacer: int = 16, max_spacer: int = 18, optimal_spacer: int = 17) -> float:
    if spacer == optimal_spacer:
        return 1.0
    if min_spacer <= spacer <= max_spacer:
        distance = abs(spacer - optimal_spacer)
        scale = max(optimal_spacer - min_spacer, max_spacer - optimal_spacer, 1)
        return max(0.5, 1.0 - 0.5 * distance / scale)
    distance = min(abs(spacer - min_spacer), abs(spacer - max_spacer))
    return max(0.0, 0.25 - 0.05 * distance)


def score_pairs(
    records: list[FastaRecord],
    hits_tsv: Path,
    minus35_motif_name: str = "minus35",
    minus10_motif_name: str = "minus10",
    min_spacer: int = 16,
    max_spacer: int = 18,
    optimal_spacer: int = 17,
    weight_minus35: float = 0.4,
    weight_minus10: float = 0.4,
    weight_spacer: float = 0.2,
    strict: bool = False,
) -> list[PairScore]:
    hits = pd.read_csv(hits_tsv, sep="\t")
    if hits.empty:
        hits = pd.DataFrame(
            columns=[
                "sequence_id",
                "motif_name",
                "start",
                "end",
                "matched_sequence",
                "normalized_score",
            ]
        )
    results: list[PairScore] = []
    for record in records:
        seq_hits = hits[hits["sequence_id"] == record.identifier].copy()
        minus35_hits = seq_hits[seq_hits["motif_name"].astype(str) == minus35_motif_name]
        minus10_hits = seq_hits[seq_hits["motif_name"].astype(str) == minus10_motif_name]
        best: tuple[float, dict] | None = None
        for _, left in minus35_hits.iterrows():
            for _, right in minus10_hits.iterrows():
                if int(left["start"]) >= int(right["start"]):
                    continue
                spacer = int(right["start"]) - int(left["end"]) - 1
                if strict and not (min_spacer <= spacer <= max_spacer):
                    continue
                s_score = spacer_score(spacer, min_spacer, max_spacer, optimal_spacer)
                minus35_score = float(left["normalized_score"])
                minus10_score = float(right["normalized_score"])
                conservation = (
                    weight_minus35 * minus35_score
                    + weight_minus10 * minus10_score
                    + weight_spacer * s_score
                )
                row = {
                    "minus35_seq": left["matched_sequence"],
                    "minus35_start": int(left["start"]),
                    "minus35_end": int(left["end"]),
                    "minus35_score": minus35_score,
                    "minus10_seq": right["matched_sequence"],
                    "minus10_start": int(right["start"]),
                    "minus10_end": int(right["end"]),
                    "minus10_score": minus10_score,
                    "spacer_len": spacer,
                    "spacer_score": s_score,
                    "conservation_score": conservation,
                    "status": "ok" if min_spacer <= spacer <= max_spacer else "suboptimal_spacer",
                }
                if best is None or conservation > best[0]:
                    best = (conservation, row)
        if best is None:
            results.append(empty_pair_score(record.identifier, "no_valid_pair"))
        else:
            row = best[1]
            results.append(
                PairScore(
                    sequence_id=record.identifier,
                    minus35_seq=str(row["minus35_seq"]),
                    minus35_start=str(row["minus35_start"]),
                    minus35_end=str(row["minus35_end"]),
                    minus35_score=f"{row['minus35_score']:.6f}",
                    minus10_seq=str(row["minus10_seq"]),
                    minus10_start=str(row["minus10_start"]),
                    minus10_end=str(row["minus10_end"]),
                    minus10_score=f"{row['minus10_score']:.6f}",
                    spacer_len=str(row["spacer_len"]),
                    spacer_score=f"{row['spacer_score']:.6f}",
                    conservation_score=f"{row['conservation_score']:.6f}",
                    status=str(row["status"]),
                )
            )
    return results


def empty_pair_score(sequence_id: str, status: str) -> PairScore:
    return PairScore(
        sequence_id=sequence_id,
        minus35_seq="",
        minus35_start="",
        minus35_end="",
        minus35_score="0.000000",
        minus10_seq="",
        minus10_start="",
        minus10_end="",
        minus10_score="0.000000",
        spacer_len="",
        spacer_score="0.000000",
        conservation_score="0.000000",
        status=status,
    )


def write_pair_scores(results: list[PairScore], path: Path) -> None:
    fieldnames = list(PairScore.__dataclass_fields__.keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for result in results:
            writer.writerow(result.__dict__)

