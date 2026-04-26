#!/usr/bin/env python

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shlex
import subprocess
import sys
import zipfile
from pathlib import Path


VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
CHARGED = set("DEKRH")
HYDROPHOBIC = set("AILMFWYV")
POLAR = set("STNQCY")
POSITIVE = set("KRH")
NEGATIVE = set("DE")
AROMATIC = set("FWY")
HYDROPATHY = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Protein-Sol or compute Protein-Sol-style sequence features.")
    parser.add_argument("--input-fasta", required=True)
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--external-output", required=True)
    parser.add_argument("--archive", required=True)
    parser.add_argument("--proteinsol-command", default=os.environ.get("PROTEINSOL_COMMAND", ""))
    parser.add_argument("--mode", choices=["features", "external"], default="features")
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


def clean_sequence(sequence: str) -> str:
    cleaned = "".join(aa for aa in sequence if aa in VALID_AA)
    if not cleaned:
        raise ValueError("Sequence contains no canonical amino-acid residues.")
    return cleaned


def fraction(sequence: str, alphabet: set[str]) -> float:
    return sum(1 for aa in sequence if aa in alphabet) / len(sequence)


def sequence_entropy(sequence: str) -> float:
    entropy = 0.0
    for aa in VALID_AA:
        count = sequence.count(aa)
        if count:
            p = count / len(sequence)
            entropy -= p * math.log2(p)
    return entropy


def feature_score(sequence: str) -> tuple[float, str]:
    charge = fraction(sequence, CHARGED)
    hydrophobic = fraction(sequence, HYDROPHOBIC)
    aromatic = fraction(sequence, AROMATIC)
    gravy = sum(HYDROPATHY[aa] for aa in sequence) / len(sequence)
    score = 0.45 + 0.75 * charge - 0.35 * hydrophobic - 0.08 * max(gravy, 0) - 0.25 * aromatic
    score = max(0.0, min(1.0, score))
    label = "above_population_average_proxy" if score >= 0.45 else "below_population_average_proxy"
    return score, label


def features_for_record(identifier: str, raw_sequence: str) -> dict:
    sequence = clean_sequence(raw_sequence)
    score, label = feature_score(sequence)
    return {
        "name": identifier,
        "length": len(sequence),
        "protein_sol_proxy_score": round(score, 6),
        "proxy_class": label,
        "charged_fraction": round(fraction(sequence, CHARGED), 6),
        "positive_fraction": round(fraction(sequence, POSITIVE), 6),
        "negative_fraction": round(fraction(sequence, NEGATIVE), 6),
        "hydrophobic_fraction": round(fraction(sequence, HYDROPHOBIC), 6),
        "polar_fraction": round(fraction(sequence, POLAR), 6),
        "aromatic_fraction": round(fraction(sequence, AROMATIC), 6),
        "gravy": round(sum(HYDROPATHY[aa] for aa in sequence) / len(sequence), 6),
        "sequence_entropy": round(sequence_entropy(sequence), 6),
        "method_note": "offline feature summary; configure official Protein-Sol for QuerySol scores",
    }


def write_features(records: list[tuple[str, str]], output_csv: Path) -> list[dict]:
    rows = [features_for_record(identifier, sequence) for identifier, sequence in records]
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def run_external(command_template: str, input_fasta: Path, output_dir: Path) -> subprocess.CompletedProcess[str]:
    if not command_template.strip():
        raise RuntimeError(
            "Protein-Sol external mode requires PROTEINSOL_COMMAND or --proteinsol-command. "
            "Use placeholders {input_fasta} and {output_dir} in the command template."
        )
    rendered = command_template.format(input_fasta=input_fasta, output_dir=output_dir)
    return subprocess.run(shlex.split(rendered, posix=(os.name != "nt")), capture_output=True, text=True)


def archive_outputs(archive: Path, paths: list[Path]) -> None:
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as zip_handle:
        for source in paths:
            if source.is_file():
                zip_handle.write(source, source.name)
            elif source.is_dir():
                for path in sorted(source.rglob("*")):
                    if path.is_file():
                        zip_handle.write(path, source.name + "/" + str(path.relative_to(source)))


def main() -> int:
    args = parse_args()
    input_fasta = Path(args.input_fasta)
    output_csv = Path(args.output_csv)
    summary_json = Path(args.summary_json)
    external_output = Path(args.external_output)
    archive = Path(args.archive)
    external_dir = Path.cwd() / "protein_sol_external_output"
    external_dir.mkdir(parents=True, exist_ok=True)

    records = read_fasta(input_fasta)
    rows = write_features(records, output_csv)
    external_status: dict | None = None

    if args.mode == "external":
        completed = run_external(args.proteinsol_command, input_fasta, external_dir)
        external_output.write_text(
            "\n".join(
                [
                    "Protein-Sol external command output",
                    f"returncode: {completed.returncode}",
                    "",
                    "[stdout]",
                    completed.stdout or "",
                    "",
                    "[stderr]",
                    completed.stderr or "",
                ]
            ),
            encoding="utf-8",
        )
        if completed.returncode != 0:
            sys.stderr.write(external_output.read_text(encoding="utf-8"))
            return 1
        external_status = {"returncode": completed.returncode, "command": args.proteinsol_command}
    else:
        external_output.write_text(
            "Protein-Sol external command was not run. The CSV contains offline sequence features only.\n",
            encoding="utf-8",
        )

    archive_outputs(archive, [output_csv, external_output, external_dir])
    summary = {
        "tool": "Protein-Sol",
        "mode": args.mode,
        "record_count": len(rows),
        "external_status": external_status,
        "note": (
            "The built-in features output is not the official QuerySol score. "
            "Use external mode with the official Protein-Sol package for official scores."
        ),
    }
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
