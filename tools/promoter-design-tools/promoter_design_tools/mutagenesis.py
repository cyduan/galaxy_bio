from __future__ import annotations

import csv
import random
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from promoter_design_tools.fasta import FastaRecord


DNA_BASES = ["A", "C", "G", "T"]


@dataclass(frozen=True)
class MutationEvent:
    mutant_id: str
    parent_id: str
    mutation_index: int
    mutation_type: str
    position_1based: int
    ref: str
    alt: str
    region: str
    seed: int


def load_conservation(path: Path | None) -> pd.DataFrame | None:
    if path is None:
        return None
    df = pd.read_csv(path, sep="\t")
    return df if not df.empty else None


def allowed_positions(
    sequence: str,
    sequence_id: str,
    region_mode: str,
    conservation: pd.DataFrame | None = None,
    custom_start: int | None = None,
    custom_end: int | None = None,
    allow_n: bool = False,
) -> tuple[list[int], str]:
    length = len(sequence)
    region = region_mode
    intervals: list[tuple[int, int]] = []

    if region_mode == "full":
        intervals = [(1, length)]
    elif region_mode == "custom":
        if custom_start is None or custom_end is None:
            raise ValueError("custom_start and custom_end are required when region_mode=custom.")
        intervals = [(custom_start, custom_end)]
    else:
        row = conservation_row(sequence_id, conservation)
        if row is None:
            intervals = []
        elif region_mode == "minus35":
            intervals = [(int(row["minus35_start"]), int(row["minus35_end"]))]
        elif region_mode == "minus10":
            intervals = [(int(row["minus10_start"]), int(row["minus10_end"]))]
        elif region_mode == "spacer":
            intervals = [(int(row["minus35_end"]) + 1, int(row["minus10_start"]) - 1)]
        elif region_mode == "non_core":
            core = {
                pos
                for start, end in [
                    (int(row["minus35_start"]), int(row["minus35_end"])),
                    (int(row["minus10_start"]), int(row["minus10_end"])),
                    (int(row["minus35_end"]) + 1, int(row["minus10_start"]) - 1),
                ]
                for pos in range(max(1, start), min(length, end) + 1)
            }
            positions = [pos for pos in range(1, length + 1) if pos not in core]
            return filter_positions(sequence, positions, allow_n), region
        else:
            raise ValueError(f"Unsupported region_mode: {region_mode}")

    positions: list[int] = []
    for start, end in intervals:
        positions.extend(range(max(1, start), min(length, end) + 1))
    return filter_positions(sequence, sorted(set(positions)), allow_n), region


def conservation_row(sequence_id: str, conservation: pd.DataFrame | None) -> pd.Series | None:
    if conservation is None:
        return None
    rows = conservation[conservation["sequence_id"] == sequence_id]
    if rows.empty:
        return None
    row = rows.iloc[0]
    if str(row.get("status", "")) == "no_valid_pair":
        return None
    return row


def filter_positions(sequence: str, positions: list[int], allow_n: bool) -> list[int]:
    if allow_n:
        return positions
    return [pos for pos in positions if sequence[pos - 1] in DNA_BASES]


def choose_mutation_count(mode: str, mutation_count: int, mutation_rate: float, allowed_count: int) -> int:
    if mode == "fixed":
        return max(0, mutation_count)
    if mode == "rate":
        return max(1, round(allowed_count * mutation_rate))
    raise ValueError(f"Unsupported mutation_count_mode: {mode}")


def generate_mutants(
    records: list[FastaRecord],
    conservation_path: Path | None,
    mutants_per_sequence: int = 10,
    mutation_mode: str = "snv",
    mutation_count_mode: str = "fixed",
    mutation_count: int = 1,
    mutation_rate: float = 0.01,
    region_mode: str = "full",
    custom_start: int | None = None,
    custom_end: int | None = None,
    keep_length: bool = True,
    allow_n: bool = False,
    seed: int = 1,
) -> tuple[list[FastaRecord], list[MutationEvent]]:
    if keep_length and mutation_mode in {"insertion", "deletion"}:
        raise ValueError("keep_length=true does not allow insertion or deletion mutation modes.")
    if keep_length and mutation_mode == "mixed":
        mutation_mode = "snv"
    conservation = load_conservation(conservation_path)
    rng = random.Random(seed)
    mutant_records: list[FastaRecord] = []
    events: list[MutationEvent] = []

    for record in records:
        positions, region = allowed_positions(
            record.sequence,
            record.identifier,
            region_mode,
            conservation,
            custom_start,
            custom_end,
            allow_n,
        )
        if not positions:
            raise ValueError(f"No mutable positions are available for {record.identifier} with region_mode={region_mode}.")
        count = choose_mutation_count(mutation_count_mode, mutation_count, mutation_rate, len(positions))
        for mutant_index in range(1, mutants_per_sequence + 1):
            mutable = list(record.sequence)
            local_events: list[MutationEvent] = []
            current_positions = positions.copy()
            for event_index in range(1, count + 1):
                mutation_type = choose_mutation_type(rng, mutation_mode)
                position = rng.choice(current_positions)
                zero_based = position - 1
                ref = mutable[zero_based]
                if mutation_type == "snv":
                    alt = rng.choice([base for base in DNA_BASES if base != ref])
                    mutable[zero_based] = alt
                elif mutation_type == "insertion":
                    alt = rng.choice(DNA_BASES)
                    mutable.insert(zero_based, alt)
                elif mutation_type == "deletion":
                    alt = "-"
                    del mutable[zero_based]
                    current_positions = [pos for pos in current_positions if pos != position]
                    if not current_positions:
                        break
                else:
                    raise ValueError(f"Unsupported mutation type: {mutation_type}")
                mutant_id = f"{record.identifier}_mut{mutant_index:06d}"
                local_events.append(
                    MutationEvent(
                        mutant_id=mutant_id,
                        parent_id=record.identifier,
                        mutation_index=event_index,
                        mutation_type=mutation_type,
                        position_1based=position,
                        ref=ref,
                        alt=alt,
                        region=region,
                        seed=seed,
                    )
                )
            mutant_id = f"{record.identifier}_mut{mutant_index:06d}"
            mutant_records.append(FastaRecord(mutant_id, "".join(mutable)))
            events.extend(local_events)
    return mutant_records, events


def choose_mutation_type(rng: random.Random, mode: str) -> str:
    if mode == "mixed":
        return rng.choice(["snv", "insertion", "deletion"])
    if mode in {"snv", "insertion", "deletion"}:
        return mode
    raise ValueError(f"Unsupported mutation_mode: {mode}")


def write_mutations(events: list[MutationEvent], path: Path) -> None:
    fieldnames = list(MutationEvent.__dataclass_fields__.keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for event in events:
            writer.writerow(event.__dict__)

