#!/usr/bin/env python

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute ProtParam-style protein properties for FASTA records.")
    parser.add_argument("--input-fasta", required=True)
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--ignore-invalid", action="store_true")
    return parser.parse_args()


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    identifier: str | None = None
    sequence_lines: list[str] = []

    def flush() -> None:
        nonlocal identifier, sequence_lines
        if identifier is None:
            return
        sequence = "".join(sequence_lines).replace(" ", "").upper()
        if not sequence:
            raise ValueError(f"FASTA record '{identifier}' is empty.")
        records.append((identifier, sequence))

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                identifier = line[1:].strip() or f"record_{len(records) + 1}"
                sequence_lines = []
            else:
                if identifier is None:
                    raise ValueError("Input is not a valid FASTA file: sequence appears before a header.")
                sequence_lines.append(line)
    flush()
    if not records:
        raise ValueError("No FASTA records were found.")
    return records


def validate_sequence(identifier: str, sequence: str, ignore_invalid: bool) -> str:
    invalid = sorted(set(sequence) - VALID_AA)
    if invalid and not ignore_invalid:
        raise ValueError(
            f"Record '{identifier}' contains unsupported amino-acid letters: {', '.join(invalid)}. "
            "Use canonical amino acids or enable 'Ignore invalid residues'."
        )
    if invalid:
        return "".join(aa for aa in sequence if aa in VALID_AA)
    return sequence


def import_protein_analysis():
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
    except Exception as exc:  # pragma: no cover - exercised on servers missing Biopython
        raise RuntimeError(
            "ProtParam requires Biopython in the Galaxy job Python environment. "
            "Install it with: python -m pip install biopython"
        ) from exc
    return ProteinAnalysis


def analyse_record(identifier: str, sequence: str, protein_analysis_cls) -> dict:
    analysis = protein_analysis_cls(sequence)
    helix, turn, sheet = analysis.secondary_structure_fraction()
    extinction_reduced, extinction_oxidized = analysis.molar_extinction_coefficient()
    counts = analysis.count_amino_acids()
    percentages = analysis.amino_acids_percent
    row = {
        "name": identifier,
        "length": len(sequence),
        "molecular_weight": round(analysis.molecular_weight(), 4),
        "aromaticity": round(analysis.aromaticity(), 6),
        "instability_index": round(analysis.instability_index(), 6),
        "instability_class": "unstable" if analysis.instability_index() > 40 else "stable",
        "isoelectric_point": round(analysis.isoelectric_point(), 6),
        "gravy": round(analysis.gravy(), 6),
        "charge_at_pH7": round(analysis.charge_at_pH(7.0), 6),
        "helix_fraction": round(helix, 6),
        "turn_fraction": round(turn, 6),
        "sheet_fraction": round(sheet, 6),
        "extinction_reduced": extinction_reduced,
        "extinction_oxidized": extinction_oxidized,
    }
    for aa in sorted(VALID_AA):
        row[f"count_{aa}"] = counts.get(aa, 0)
        row[f"percent_{aa}"] = round(percentages.get(aa, 0.0), 6)
    return row


def main() -> int:
    args = parse_args()
    protein_analysis_cls = import_protein_analysis()
    records = read_fasta(Path(args.input_fasta))
    rows: list[dict] = []
    for identifier, sequence in records:
        cleaned = validate_sequence(identifier, sequence, args.ignore_invalid)
        rows.append(analyse_record(identifier, cleaned, protein_analysis_cls))

    fieldnames = list(rows[0].keys())
    with Path(args.output_csv).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "tool": "ProtParam",
        "record_count": len(rows),
        "fields": fieldnames,
        "biopython_required": True,
    }
    Path(args.summary_json).write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
