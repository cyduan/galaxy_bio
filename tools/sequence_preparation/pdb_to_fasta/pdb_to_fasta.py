#!/usr/bin/env python3
"""Extract protein chain sequences from PDB files and write FASTA."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path


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
    # Common modified or ambiguous amino acids.
    "MSE": "M",
    "SEC": "U",
    "PYL": "O",
    "ASX": "B",
    "GLX": "Z",
    "XLE": "J",
    "UNK": "X",
    "SEP": "S",
    "TPO": "T",
    "PTR": "Y",
    "HYP": "P",
    "CSO": "C",
    "CSD": "C",
    "CME": "C",
    "CYX": "C",
    "HID": "H",
    "HIE": "H",
    "HIP": "H",
    "HSD": "H",
    "HSE": "H",
    "HSP": "H",
    "MLY": "K",
    "MLZ": "K",
    "LLP": "K",
    "PCA": "E",
}

STANDARD_AA3 = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}


@dataclass
class Residue:
    chain_id: str
    residue_name: str
    residue_number: str
    insertion_code: str
    record_type: str
    one_letter: str


@dataclass
class ChainSummary:
    chain_id: str
    sequence_id: str
    model: str
    residues_total: int = 0
    sequence_length: int = 0
    hetatm_residues_used: int = 0
    nonstandard_residues_used: int = 0
    unknown_residues_as_x: int = 0
    skipped_hetatm_records: int = 0
    skipped_duplicate_atom_records: int = 0
    warning: list[str] = field(default_factory=list)


class PdbToFastaError(RuntimeError):
    """User-facing error."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert PDB protein chains to FASTA.")
    parser.add_argument("--input-pdb", required=True, help="Input PDB file.")
    parser.add_argument("--output-fasta", required=True, help="Output FASTA file.")
    parser.add_argument("--report-tsv", required=True, help="Output extraction report TSV.")
    parser.add_argument(
        "--chain-mode",
        choices=["all", "selected"],
        default="all",
        help="Extract all chains or only one selected chain.",
    )
    parser.add_argument("--chain-id", default="", help="Chain ID used when --chain-mode selected.")
    parser.add_argument(
        "--record-mode",
        choices=["separate_chains", "concatenate_chains"],
        default="separate_chains",
        help="Write one FASTA record per chain or concatenate extracted chains into one sequence.",
    )
    parser.add_argument(
        "--model",
        default="1",
        help="PDB MODEL number to extract. Use 'all' to read all models; default is first model (1).",
    )
    parser.add_argument(
        "--include-hetatm",
        action="store_true",
        help="Use HETATM records for known modified amino acids such as MSE.",
    )
    parser.add_argument(
        "--unknown-to-x",
        action="store_true",
        help="Represent unknown ATOM residues as X instead of skipping them.",
    )
    parser.add_argument(
        "--sequence-id-prefix",
        default="",
        help="Optional FASTA ID prefix. If omitted, the input file stem is used.",
    )
    return parser.parse_args()


def normalize_chain_id(chain_id: str) -> str:
    chain_id = chain_id.strip()
    return chain_id if chain_id else "_"


def selected_model(line: str, current_model: str | None, requested_model: str) -> tuple[str | None, bool]:
    record = line[0:6].strip()
    if record == "MODEL":
        model_value = line[10:14].strip() or line.split(maxsplit=1)[1].strip() if len(line.split()) > 1 else "1"
        return model_value, requested_model == "all" or model_value == requested_model
    if record == "ENDMDL":
        return None, requested_model == "all"
    if current_model is None:
        # PDB files without MODEL records are treated as model 1.
        return current_model, requested_model in {"1", "all"}
    return current_model, requested_model == "all" or current_model == requested_model


def parse_pdb_atom_line(line: str) -> tuple[str, str, str, str, str] | None:
    record_type = line[0:6].strip()
    if record_type not in {"ATOM", "HETATM"}:
        return None

    if len(line) >= 27:
        residue_name = line[17:20].strip().upper()
        chain_id = normalize_chain_id(line[21:22])
        residue_number = line[22:26].strip()
        insertion_code = line[26:27].strip()
        if residue_name and residue_number:
            return record_type, residue_name, chain_id, residue_number, insertion_code

    # Fallback for non-strict PDB-like files with whitespace separated columns.
    parts = line.split()
    if len(parts) >= 6:
        residue_name = parts[3].upper()
        chain_id = normalize_chain_id(parts[4])
        residue_number = parts[5]
        return record_type, residue_name, chain_id, residue_number, ""
    return None


def read_pdb_sequences(args: argparse.Namespace) -> tuple[OrderedDict[str, list[Residue]], dict[str, ChainSummary]]:
    input_path = Path(args.input_pdb)
    requested_chain = normalize_chain_id(args.chain_id) if args.chain_mode == "selected" else None
    sequence_prefix = sanitize_id(args.sequence_id_prefix or input_path.stem)
    requested_model = str(args.model).strip() or "1"
    if requested_model.lower() == "all":
        requested_model = "all"

    chains: OrderedDict[str, list[Residue]] = OrderedDict()
    summaries: dict[str, ChainSummary] = {}
    seen_residues: set[tuple[str, str, str, str, str]] = set()
    current_model: str | None = None
    any_atom_records = False

    with input_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            current_model, model_is_active = selected_model(line, current_model, requested_model)
            parsed = parse_pdb_atom_line(line)
            if parsed is None:
                continue
            any_atom_records = True
            if not model_is_active:
                continue

            record_type, residue_name, chain_id, residue_number, insertion_code = parsed
            if requested_chain is not None and chain_id != requested_chain:
                continue

            sequence_id = f"{sequence_prefix}_chain_{chain_id}"
            summary = summaries.setdefault(
                chain_id,
                ChainSummary(
                    chain_id=chain_id,
                    sequence_id=sequence_id,
                    model=requested_model,
                ),
            )

            key_model = current_model or "1"
            residue_key = (key_model, chain_id, residue_number, insertion_code, residue_name)
            if residue_key in seen_residues:
                summary.skipped_duplicate_atom_records += 1
                continue
            seen_residues.add(residue_key)

            one_letter = AA3_TO_1.get(residue_name)
            if record_type == "HETATM":
                if not args.include_hetatm or one_letter is None:
                    summary.skipped_hetatm_records += 1
                    continue
                summary.hetatm_residues_used += 1

            if one_letter is None:
                if record_type == "ATOM" and args.unknown_to_x:
                    one_letter = "X"
                    summary.unknown_residues_as_x += 1
                else:
                    summary.warning.append(f"skipped_unknown_residue:{residue_name}{residue_number}{insertion_code}")
                    continue

            if residue_name not in STANDARD_AA3:
                summary.nonstandard_residues_used += 1

            residue = Residue(
                chain_id=chain_id,
                residue_name=residue_name,
                residue_number=residue_number,
                insertion_code=insertion_code,
                record_type=record_type,
                one_letter=one_letter,
            )
            chains.setdefault(chain_id, []).append(residue)
            summary.residues_total += 1
            summary.sequence_length += 1

    if not any_atom_records:
        raise PdbToFastaError("No ATOM/HETATM records were found. Please upload a PDB coordinate file.")
    if requested_chain is not None and requested_chain not in summaries:
        raise PdbToFastaError(f"Selected chain {requested_chain!r} was not found in the requested model.")
    if not chains:
        raise PdbToFastaError("No protein residues could be extracted from the PDB file.")

    for summary in summaries.values():
        if summary.sequence_length == 0:
            summary.warning.append("no_sequence_extracted")
        if summary.unknown_residues_as_x:
            summary.warning.append("unknown_residues_written_as_X")
        if summary.nonstandard_residues_used:
            summary.warning.append("nonstandard_residues_mapped")
    return chains, summaries


def sanitize_id(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.:-]+", "_", value.strip())
    return value.strip("_") or "pdb_sequence"


def write_outputs(
    chains: OrderedDict[str, list[Residue]],
    summaries: dict[str, ChainSummary],
    output_fasta: Path,
    report_tsv: Path,
    record_mode: str,
    sequence_prefix: str,
) -> None:
    with output_fasta.open("w", encoding="utf-8", newline="\n") as fasta_handle:
        if record_mode == "concatenate_chains":
            ordered_chains = [chain_id for chain_id, residues in chains.items() if residues]
            concatenated = "".join("".join(residue.one_letter for residue in chains[chain_id]) for chain_id in ordered_chains)
            chain_label = ",".join(ordered_chains)
            fasta_handle.write(f">{sanitize_id(sequence_prefix)}_chains_{sanitize_id(chain_label)} chains={chain_label}\n")
            write_wrapped_sequence(fasta_handle, concatenated)
        else:
            for chain_id, residues in chains.items():
                if not residues:
                    continue
                sequence = "".join(residue.one_letter for residue in residues)
                summary = summaries[chain_id]
                fasta_handle.write(f">{summary.sequence_id} chain={chain_id} length={len(sequence)}\n")
                write_wrapped_sequence(fasta_handle, sequence)

    with report_tsv.open("w", encoding="utf-8", newline="") as report_handle:
        fieldnames = [
            "chain_id",
            "sequence_id",
            "model",
            "sequence_length",
            "hetatm_residues_used",
            "nonstandard_residues_used",
            "unknown_residues_as_x",
            "skipped_hetatm_records",
            "skipped_duplicate_atom_records",
            "warning",
        ]
        writer = csv.DictWriter(report_handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for chain_id in chains:
            summary = summaries[chain_id]
            writer.writerow(
                {
                    "chain_id": summary.chain_id,
                    "sequence_id": summary.sequence_id,
                    "model": summary.model,
                    "sequence_length": summary.sequence_length,
                    "hetatm_residues_used": summary.hetatm_residues_used,
                    "nonstandard_residues_used": summary.nonstandard_residues_used,
                    "unknown_residues_as_x": summary.unknown_residues_as_x,
                    "skipped_hetatm_records": summary.skipped_hetatm_records,
                    "skipped_duplicate_atom_records": summary.skipped_duplicate_atom_records,
                    "warning": ";".join(dict.fromkeys(summary.warning)),
                }
            )


def write_wrapped_sequence(handle, sequence: str, width: int = 80) -> None:
    for index in range(0, len(sequence), width):
        handle.write(sequence[index : index + width] + "\n")


def main() -> int:
    args = parse_args()
    try:
        chains, summaries = read_pdb_sequences(args)
        sequence_prefix = sanitize_id(args.sequence_id_prefix or Path(args.input_pdb).stem)
        write_outputs(
            chains,
            summaries,
            Path(args.output_fasta),
            Path(args.report_tsv),
            args.record_mode,
            sequence_prefix,
        )
    except PdbToFastaError as exc:
        print(f"pdb_to_fasta: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
