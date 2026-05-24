#!/usr/bin/env python3
"""Validate residue numbering consistency across PDB, FASTA, MSA, and TSV outputs."""

from __future__ import annotations

import argparse
import csv
import html
import sys
from dataclasses import dataclass
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
    "MSE": "M",
    "SEC": "U",
    "PYL": "O",
}

GAP_CHARS = {"-", "."}

VALIDATION_COLUMNS = [
    "sequence_position",
    "pdb_chain",
    "pdb_residue_number",
    "pdb_insertion_code",
    "pdb_residue_name",
    "pdb_aa",
    "fasta_aa",
    "msa_column",
    "msa_aa",
    "structure_feature_aa",
    "conservation_aa",
    "pdb_matches_fasta",
    "msa_matches_fasta",
    "structure_features_present",
    "conservation_present",
    "numbering_matches_sequence_position",
    "status",
    "warning",
]

SUMMARY_COLUMNS = ["metric", "value", "severity", "note"]


@dataclass
class PdbResidue:
    sequence_position: int
    chain: str
    residue_number: str
    insertion_code: str
    residue_name: str
    aa: str


@dataclass
class FastaRecord:
    identifier: str
    sequence: str


class ToolError(RuntimeError):
    """User-facing error."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Check that PDB residue numbering, protein FASTA positions, MSA columns, "
            "structure-feature TSV rows, and conservation TSV rows refer to the same residues."
        )
    )
    parser.add_argument("--input-pdb", required=True, help="Input PDB structure.")
    parser.add_argument("--protein-fasta", required=True, help="Target protein FASTA.")
    parser.add_argument("--msa-fasta", default="", help="Optional MSA FASTA.")
    parser.add_argument("--structure-features", default="", help="Optional residue_structure_features.tsv.")
    parser.add_argument("--conservation-tsv", default="", help="Optional residue_conservation.tsv.")
    parser.add_argument("--chain-id", default="", help="PDB chain to validate. Empty/auto uses the first protein chain.")
    parser.add_argument("--msa-reference-id", default="", help="Optional MSA record ID corresponding to the target protein.")
    parser.add_argument("--validation-tsv", required=True, help="Output row-level validation TSV.")
    parser.add_argument("--summary-tsv", required=True, help="Output summary TSV.")
    parser.add_argument("--report-html", required=True, help="Output HTML report.")
    parser.add_argument("--run-log", required=True, help="Output log file.")
    return parser.parse_args()


def normalize_chain(chain: str) -> str:
    return (chain or "").strip() or "_"


def parse_pdb_atom_line(line: str) -> tuple[str, str, str, str] | None:
    if not line.startswith("ATOM"):
        return None
    try:
        residue_name = line[17:20].strip().upper()
        chain = normalize_chain(line[21:22])
        residue_number = line[22:26].strip()
        insertion_code = line[26:27].strip()
        if residue_name and residue_number:
            return chain, residue_number, insertion_code, residue_name
    except Exception:
        pass

    parts = line.split()
    if len(parts) >= 6:
        residue_name = parts[3].upper()
        chain = normalize_chain(parts[4])
        residue_number = parts[5]
        return chain, residue_number, "", residue_name
    return None


def read_pdb_residues(path: Path, chain_id: str) -> tuple[str, list[PdbResidue]]:
    seen: set[tuple[str, str, str]] = set()
    by_chain: dict[str, list[tuple[str, str, str, str]]] = {}
    chain_order: list[str] = []

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parsed = parse_pdb_atom_line(line)
            if parsed is None:
                continue
            chain, residue_number, insertion_code, residue_name = parsed
            if residue_name not in AA3_TO_1:
                continue
            key = (chain, residue_number, insertion_code)
            if key in seen:
                continue
            seen.add(key)
            if chain not in by_chain:
                by_chain[chain] = []
                chain_order.append(chain)
            by_chain[chain].append((chain, residue_number, insertion_code, residue_name))

    if not by_chain:
        raise ToolError("No protein ATOM residues were parsed from the PDB.")

    requested = (chain_id or "").strip()
    if requested.lower() in {"", "auto", "first"}:
        selected_chain = chain_order[0]
    else:
        selected_chain = normalize_chain(requested)
        if selected_chain not in by_chain:
            available = ", ".join(chain_order)
            raise ToolError(f"Chain '{selected_chain}' was not found. Available chains: {available}")

    residues = [
        PdbResidue(
            sequence_position=i + 1,
            chain=chain,
            residue_number=residue_number,
            insertion_code=insertion_code,
            residue_name=residue_name,
            aa=AA3_TO_1.get(residue_name, "X"),
        )
        for i, (chain, residue_number, insertion_code, residue_name) in enumerate(by_chain[selected_chain])
    ]
    return selected_chain, residues


def read_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    current_id = ""
    parts: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    records.append(FastaRecord(current_id, "".join(parts).upper()))
                current_id = line[1:].split()[0] or f"record_{len(records) + 1}"
                parts = []
            else:
                parts.append(line.replace(" ", "").upper())
    if current_id:
        records.append(FastaRecord(current_id, "".join(parts).upper()))
    if not records:
        raise ToolError(f"No FASTA records found in {path}")
    return records


def ungap(sequence: str) -> str:
    return "".join(aa for aa in sequence.upper() if aa not in GAP_CHARS)


def choose_msa_reference(records: list[FastaRecord], target_sequence: str, requested_id: str) -> tuple[FastaRecord, str]:
    if requested_id:
        for record in records:
            if record.identifier == requested_id:
                return record, "requested_id"
        raise ToolError(f"MSA reference ID '{requested_id}' was not found.")

    for record in records:
        if ungap(record.sequence) == target_sequence:
            return record, "exact_sequence_match"
    return records[0], "first_record_fallback"


def map_msa_columns(reference: FastaRecord) -> dict[int, tuple[int, str]]:
    mapping: dict[int, tuple[int, str]] = {}
    sequence_position = 0
    for column_index, aa in enumerate(reference.sequence.upper(), start=1):
        if aa in GAP_CHARS:
            continue
        sequence_position += 1
        mapping[sequence_position] = (column_index, aa)
    return mapping


def read_tsv(path: str) -> list[dict[str, str]]:
    if not path:
        return []
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return []
    with p.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            return []
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def index_structure_features(rows: list[dict[str, str]]) -> dict[tuple[str, str, str], dict[str, str]]:
    indexed: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in rows:
        chain = normalize_chain(row.get("chain_id") or row.get("chain") or "")
        number = row.get("residue_number") or row.get("residue_index") or ""
        insertion_code = row.get("insertion_code") or ""
        if number:
            indexed[(chain, number, insertion_code)] = row
    return indexed


def index_conservation(rows: list[dict[str, str]], chain: str) -> dict[int, dict[str, str]]:
    indexed: dict[int, dict[str, str]] = {}
    for row in rows:
        row_chain = normalize_chain(row.get("chain") or chain)
        if row_chain != chain:
            continue
        try:
            position = int(row.get("residue_index") or row.get("sequence_position") or "")
        except ValueError:
            continue
        indexed[position] = row
    return indexed


def is_int_string(value: str) -> bool:
    try:
        int(value)
        return True
    except ValueError:
        return False


def numbering_matches_sequence_position(residue: PdbResidue) -> bool:
    return residue.insertion_code == "" and is_int_string(residue.residue_number) and int(residue.residue_number) == residue.sequence_position


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def build_validation_rows(
    chain: str,
    pdb_residues: list[PdbResidue],
    fasta_sequence: str,
    msa_mapping: dict[int, tuple[int, str]],
    structure_features: dict[tuple[str, str, str], dict[str, str]],
    conservation: dict[int, dict[str, str]],
) -> list[dict[str, str]]:
    max_len = max(len(pdb_residues), len(fasta_sequence), max(conservation.keys(), default=0))
    rows: list[dict[str, str]] = []
    for sequence_position in range(1, max_len + 1):
        pdb = pdb_residues[sequence_position - 1] if sequence_position <= len(pdb_residues) else None
        fasta_aa = fasta_sequence[sequence_position - 1] if sequence_position <= len(fasta_sequence) else ""
        msa_column, msa_aa = msa_mapping.get(sequence_position, ("", ""))
        feature_row = structure_features.get((chain, pdb.residue_number, pdb.insertion_code)) if pdb else None
        conservation_row = conservation.get(sequence_position)
        feature_aa = ""
        if feature_row:
            feature_aa = (
                feature_row.get("amino_acid")
                or feature_row.get("wt_aa")
                or feature_row.get("aa")
                or feature_row.get("residue_aa")
                or ""
            ).upper()
        conservation_aa = (conservation_row.get("wt_aa") or conservation_row.get("amino_acid") or "").upper() if conservation_row else ""

        warnings: list[str] = []
        errors: list[str] = []

        if pdb is None:
            errors.append("missing_pdb_residue_for_sequence_position")
        elif fasta_aa and pdb.aa != fasta_aa:
            errors.append("pdb_fasta_amino_acid_mismatch")

        if not fasta_aa:
            warnings.append("missing_fasta_position_for_pdb_or_conservation_row")

        if msa_mapping and not msa_aa:
            warnings.append("missing_msa_column")
        elif msa_aa and fasta_aa and msa_aa != fasta_aa:
            errors.append("msa_fasta_amino_acid_mismatch")

        if structure_features and feature_row is None:
            warnings.append("missing_structure_feature_row")
        elif feature_aa and fasta_aa and feature_aa != fasta_aa:
            errors.append("structure_feature_amino_acid_mismatch")

        if conservation and conservation_row is None:
            warnings.append("missing_conservation_row")
        elif conservation_aa and fasta_aa and conservation_aa != fasta_aa:
            errors.append("conservation_amino_acid_mismatch")

        if pdb and not numbering_matches_sequence_position(pdb):
            warnings.append("pdb_residue_number_not_equal_sequence_position")

        status = "error" if errors else "warning" if warnings else "ok"
        rows.append(
            {
                "sequence_position": str(sequence_position),
                "pdb_chain": pdb.chain if pdb else chain,
                "pdb_residue_number": pdb.residue_number if pdb else "",
                "pdb_insertion_code": pdb.insertion_code if pdb else "",
                "pdb_residue_name": pdb.residue_name if pdb else "",
                "pdb_aa": pdb.aa if pdb else "",
                "fasta_aa": fasta_aa,
                "msa_column": str(msa_column),
                "msa_aa": msa_aa,
                "structure_feature_aa": feature_aa,
                "conservation_aa": conservation_aa,
                "pdb_matches_fasta": yes_no(bool(pdb and fasta_aa and pdb.aa == fasta_aa)),
                "msa_matches_fasta": yes_no(bool(msa_aa and fasta_aa and msa_aa == fasta_aa)) if msa_mapping else "not_checked",
                "structure_features_present": yes_no(feature_row is not None) if structure_features else "not_checked",
                "conservation_present": yes_no(conservation_row is not None) if conservation else "not_checked",
                "numbering_matches_sequence_position": yes_no(bool(pdb and numbering_matches_sequence_position(pdb))),
                "status": status,
                "warning": ";".join(errors + warnings),
            }
        )
    return rows


def summarize(
    rows: list[dict[str, str]],
    chain: str,
    pdb_residues: list[PdbResidue],
    fasta_record: FastaRecord,
    msa_records: list[FastaRecord],
    msa_reference: str,
    msa_reference_mode: str,
) -> list[dict[str, str]]:
    error_count = sum(1 for row in rows if row["status"] == "error")
    warning_count = sum(1 for row in rows if row["status"] == "warning")
    numbering_warning_count = sum(
        1 for row in rows if "pdb_residue_number_not_equal_sequence_position" in row["warning"].split(";")
    )
    aa_mismatch_count = sum(1 for row in rows if "amino_acid_mismatch" in row["warning"])
    summary = [
        {"metric": "selected_chain", "value": chain, "severity": "info", "note": "PDB chain used for validation."},
        {"metric": "pdb_residue_count", "value": str(len(pdb_residues)), "severity": "info", "note": "Protein residues parsed from the selected PDB chain."},
        {"metric": "protein_fasta_id", "value": fasta_record.identifier, "severity": "info", "note": "Target sequence record."},
        {"metric": "protein_fasta_length", "value": str(len(fasta_record.sequence)), "severity": "info", "note": "Ungapped target sequence length."},
        {"metric": "msa_record_count", "value": str(len(msa_records)), "severity": "info", "note": "Zero means no MSA file was provided."},
        {"metric": "msa_reference_id", "value": msa_reference, "severity": "info", "note": f"Reference selected by {msa_reference_mode}."},
        {"metric": "row_error_count", "value": str(error_count), "severity": "error" if error_count else "ok", "note": "Amino-acid or length mismatches that can invalidate downstream residue mapping."},
        {"metric": "row_warning_count", "value": str(warning_count), "severity": "warning" if warning_count else "ok", "note": "Missing optional rows or numbering differences that may require attention."},
        {"metric": "amino_acid_mismatch_count", "value": str(aa_mismatch_count), "severity": "error" if aa_mismatch_count else "ok", "note": "PDB/FASTA/MSA/TSV residue identity mismatches."},
        {"metric": "pdb_numbering_not_sequence_position_count", "value": str(numbering_warning_count), "severity": "warning" if numbering_warning_count else "ok", "note": "Non-1-based or insertion-coded PDB numbering. Downstream merges by residue number may need a mapping table."},
    ]

    merge_risk = "high" if error_count else "moderate" if numbering_warning_count else "low"
    risk_note = {
        "low": "PDB numbering and sequence positions look directly mergeable.",
        "moderate": "Residue identities match, but PDB residue numbers differ from sequence positions. Use the validation table when merging Tool 2 and Tool 4 outputs.",
        "high": "Residue identities or lengths disagree. Fix chain selection, FASTA sequence, or structure before trusting hotspot ranking.",
    }[merge_risk]
    summary.append({"metric": "hotspot_merge_risk", "value": merge_risk, "severity": "error" if merge_risk == "high" else "warning" if merge_risk == "moderate" else "ok", "note": risk_note})
    return summary


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def render_table(rows: list[dict[str, str]], columns: list[str], limit: int | None = None) -> str:
    display_rows = rows if limit is None else rows[:limit]
    header = "".join(f"<th>{html.escape(col)}</th>" for col in columns)
    body = []
    for row in display_rows:
        status = row.get("status") or row.get("severity") or ""
        css = "error" if status == "error" else "warning" if status == "warning" else "ok" if status == "ok" else ""
        body.append(
            f"<tr class='{css}'>" + "".join(f"<td>{html.escape(str(row.get(col, '')))}</td>" for col in columns) + "</tr>"
        )
    if not body:
        body.append(f"<tr><td colspan='{len(columns)}'>No rows.</td></tr>")
    return f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def render_report(path: Path, summary_rows: list[dict[str, str]], validation_rows: list[dict[str, str]]) -> None:
    risk = next((row["value"] for row in summary_rows if row["metric"] == "hotspot_merge_risk"), "unknown")
    html_text = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>Residue Numbering Validator Report</title>
  <style>
    body {{ font-family: system-ui, sans-serif; margin: 2rem; color: #1f2937; background: #f8fafc; }}
    table {{ border-collapse: collapse; width: 100%; background: white; margin-bottom: 1.5rem; }}
    th, td {{ border: 1px solid #d1d5db; padding: 6px 8px; font-size: 13px; text-align: left; }}
    th {{ background: #111827; color: white; }}
    .card {{ display: inline-block; padding: 14px 18px; border-radius: 12px; background: white; border: 1px solid #d1d5db; margin: 0 12px 12px 0; }}
    .ok td {{ background: #ecfdf5; }}
    .warning td {{ background: #fffbeb; }}
    .error td {{ background: #fef2f2; }}
    code {{ background: #e5e7eb; padding: 2px 4px; border-radius: 4px; }}
  </style>
</head>
<body>
  <h1>Residue Numbering Validator Report</h1>
  <div class="card"><strong>HotSpot merge risk</strong><br>{html.escape(risk)}</div>
  <p>This report checks whether PDB residues, FASTA sequence positions, MSA reference columns, Tool 2 structure features, and Tool 4 conservation rows point to the same biological residues.</p>
  <h2>Summary</h2>
  {render_table(summary_rows, SUMMARY_COLUMNS)}
  <h2>Residue Mapping Preview</h2>
  {render_table(validation_rows, VALIDATION_COLUMNS, limit=200)}
  <p>If PDB numbering differs from sequence positions but amino acids match, the workflow can still be biologically valid; use this table as the mapping bridge when reviewing hotspot rows.</p>
</body>
</html>"""
    path.write_text(html_text, encoding="utf-8")


def main() -> int:
    args = parse_args()
    try:
        selected_chain, pdb_residues = read_pdb_residues(Path(args.input_pdb), args.chain_id)
        fasta_record = read_fasta(Path(args.protein_fasta))[0]
        fasta_sequence = ungap(fasta_record.sequence)

        msa_records: list[FastaRecord] = []
        msa_mapping: dict[int, tuple[int, str]] = {}
        msa_reference_id = ""
        msa_reference_mode = "not_provided"
        if args.msa_fasta:
            msa_records = read_fasta(Path(args.msa_fasta))
            reference, msa_reference_mode = choose_msa_reference(msa_records, fasta_sequence, args.msa_reference_id)
            msa_reference_id = reference.identifier
            msa_mapping = map_msa_columns(reference)

        structure_feature_rows = read_tsv(args.structure_features)
        conservation_rows = read_tsv(args.conservation_tsv)
        structure_features = index_structure_features(structure_feature_rows)
        conservation = index_conservation(conservation_rows, selected_chain)

        validation_rows = build_validation_rows(
            selected_chain,
            pdb_residues,
            fasta_sequence,
            msa_mapping,
            structure_features,
            conservation,
        )
        summary_rows = summarize(
            validation_rows,
            selected_chain,
            pdb_residues,
            fasta_record,
            msa_records,
            msa_reference_id,
            msa_reference_mode,
        )

        write_tsv(Path(args.validation_tsv), VALIDATION_COLUMNS, validation_rows)
        write_tsv(Path(args.summary_tsv), SUMMARY_COLUMNS, summary_rows)
        render_report(Path(args.report_html), summary_rows, validation_rows)
        Path(args.run_log).write_text(
            "\n".join(f"{row['metric']}={row['value']}" for row in summary_rows) + "\n",
            encoding="utf-8",
        )
        return 0
    except Exception as exc:
        Path(args.run_log).write_text(f"ERROR: {exc}\n", encoding="utf-8")
        print(f"residue_numbering_validator: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
