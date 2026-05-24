#!/usr/bin/env python3
"""Annotate functional residues from PDB ligands/metals and optional TSV evidence."""

from __future__ import annotations

import argparse
import csv
import html
import math
import sys
from collections import defaultdict
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
    "MSE": "M",
    "SEC": "U",
    "PYL": "O",
}

WATER_NAMES = {"HOH", "WAT", "DOD", "H2O"}
METAL_NAMES = {
    "ZN",
    "MG",
    "MN",
    "FE",
    "FE2",
    "CU",
    "CO",
    "NI",
    "CA",
    "NA",
    "K",
    "CD",
    "HG",
    "MO",
    "W",
}

ANNOTATION_COLUMNS = [
    "residue_id",
    "chain",
    "residue_number",
    "insertion_code",
    "residue_name",
    "wt_aa",
    "functional_role",
    "protect_from_mutation",
    "evidence_sources",
    "evidence_count",
    "nearest_ligand_id",
    "nearest_ligand_distance",
    "ligand_ids",
    "ligand_contact_count",
    "metal_contact",
    "annotations",
    "notes",
]

ACTIVE_SITE_COLUMNS = [
    "site_id",
    "chain",
    "residue_number",
    "residue_name",
    "role",
    "x",
    "y",
    "z",
    "evidence_source",
    "note",
]

CONTACT_COLUMNS = [
    "ligand_id",
    "ligand_name",
    "ligand_chain",
    "ligand_residue_number",
    "ligand_atom",
    "residue_id",
    "chain",
    "residue_number",
    "residue_name",
    "atom_name",
    "distance",
    "contact_type",
]


@dataclass
class Atom:
    record: str
    atom_name: str
    residue_name: str
    chain: str
    residue_number: str
    insertion_code: str
    x: float
    y: float
    z: float
    element: str = ""


@dataclass
class Residue:
    chain: str
    number: str
    insertion_code: str
    name: str
    atoms: list[Atom] = field(default_factory=list)
    roles: set[str] = field(default_factory=set)
    sources: set[str] = field(default_factory=set)
    ligand_ids: set[str] = field(default_factory=set)
    annotations: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    nearest_ligand_id: str = ""
    nearest_ligand_distance: float | None = None
    ligand_contact_count: int = 0
    metal_contact: bool = False
    protect: bool = False

    @property
    def residue_id(self) -> str:
        return f"{self.chain}:{self.name}{self.number}{self.insertion_code}"

    @property
    def wt_aa(self) -> str:
        return AA3_TO_1.get(self.name, "")


class ToolError(RuntimeError):
    """User-facing error."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Annotate functional residues from structure ligands/metals and TSV evidence.")
    parser.add_argument("--input-pdb", required=True)
    parser.add_argument("--functional-annotation-tsv", default="")
    parser.add_argument("--residue-annotation-tsv", required=True)
    parser.add_argument("--active-site-tsv", required=True)
    parser.add_argument("--ligand-contacts-tsv", required=True)
    parser.add_argument("--report-html", required=True)
    parser.add_argument("--run-log", required=True)
    parser.add_argument("--ligand-distance-cutoff", type=float, default=4.5)
    parser.add_argument("--metal-distance-cutoff", type=float, default=3.0)
    parser.add_argument("--ignore-ligands", default="HOH,WAT,DOD,H2O")
    parser.add_argument(
        "--protect-roles",
        default="catalytic,active_site,metal_binding,essential,cofactor_binding",
        help="Comma-separated role substrings that should be marked as protected.",
    )
    parser.add_argument("--include-ligand-contacts-in-active-site", action="store_true")
    return parser.parse_args()


def normalize_chain(chain: str) -> str:
    return chain.strip() or "_"


def parse_pdb_atom_line(line: str) -> Atom | None:
    record = line[0:6].strip()
    if record not in {"ATOM", "HETATM"}:
        return None
    try:
        if len(line) >= 54:
            return Atom(
                record=record,
                atom_name=line[12:16].strip(),
                residue_name=line[17:20].strip().upper(),
                chain=normalize_chain(line[21:22]),
                residue_number=line[22:26].strip(),
                insertion_code=line[26:27].strip(),
                x=float(line[30:38]),
                y=float(line[38:46]),
                z=float(line[46:54]),
                element=line[76:78].strip().upper() if len(line) >= 78 else "",
            )
    except ValueError:
        pass

    parts = line.split()
    if len(parts) >= 9:
        try:
            return Atom(
                record=record,
                atom_name=parts[2],
                residue_name=parts[3].upper(),
                chain=normalize_chain(parts[4]),
                residue_number=parts[5],
                insertion_code="",
                x=float(parts[6]),
                y=float(parts[7]),
                z=float(parts[8]),
                element=parts[-1].upper() if parts else "",
            )
        except ValueError:
            return None
    return None


def residue_key(atom: Atom) -> tuple[str, str, str]:
    return atom.chain, atom.residue_number, atom.insertion_code


def ligand_id(atom: Atom) -> str:
    return f"{atom.residue_name}:{atom.chain}:{atom.residue_number}{atom.insertion_code}"


def is_protein_atom(atom: Atom) -> bool:
    return atom.record == "ATOM" and atom.residue_name in AA3_TO_1


def is_metal(atom: Atom) -> bool:
    element = atom.element or atom.residue_name
    return element.upper() in METAL_NAMES or atom.residue_name.upper() in METAL_NAMES


def read_structure(path: Path, ignore_ligands: set[str]) -> tuple[dict[tuple[str, str, str], Residue], list[Atom]]:
    residues: dict[tuple[str, str, str], Residue] = {}
    ligands: list[Atom] = []
    seen_residue_atoms: set[tuple[str, str, str, str]] = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            atom = parse_pdb_atom_line(line)
            if atom is None:
                continue
            if is_protein_atom(atom):
                key = residue_key(atom)
                residue = residues.setdefault(
                    key,
                    Residue(
                        chain=atom.chain,
                        number=atom.residue_number,
                        insertion_code=atom.insertion_code,
                        name=atom.residue_name,
                    ),
                )
                atom_key = (*key, atom.atom_name)
                if atom_key not in seen_residue_atoms:
                    residue.atoms.append(atom)
                    seen_residue_atoms.add(atom_key)
            elif atom.record == "HETATM" and atom.residue_name not in ignore_ligands and atom.residue_name not in AA3_TO_1:
                ligands.append(atom)
    if not residues:
        raise ToolError("No protein ATOM residues were parsed from the PDB.")
    return residues, ligands


def distance(left: Atom, right: Atom) -> float:
    return math.sqrt((left.x - right.x) ** 2 + (left.y - right.y) ** 2 + (left.z - right.z) ** 2)


def residue_centroid(residue: Residue) -> tuple[float, float, float]:
    if not residue.atoms:
        return 0.0, 0.0, 0.0
    n = len(residue.atoms)
    return (
        sum(atom.x for atom in residue.atoms) / n,
        sum(atom.y for atom in residue.atoms) / n,
        sum(atom.z for atom in residue.atoms) / n,
    )


def annotate_contacts(
    residues: dict[tuple[str, str, str], Residue],
    ligands: list[Atom],
    ligand_cutoff: float,
    metal_cutoff: float,
) -> list[dict[str, str]]:
    contacts: list[dict[str, str]] = []
    ligand_atoms_by_id: dict[str, list[Atom]] = defaultdict(list)
    for ligand in ligands:
        ligand_atoms_by_id[ligand_id(ligand)].append(ligand)

    for residue in residues.values():
        nearest_distance: float | None = None
        nearest_id = ""
        for ligand_key, ligand_atoms in ligand_atoms_by_id.items():
            ligand_name = ligand_atoms[0].residue_name
            cutoff = metal_cutoff if any(is_metal(atom) for atom in ligand_atoms) else ligand_cutoff
            best_pair: tuple[Atom, Atom, float] | None = None
            for res_atom in residue.atoms:
                for lig_atom in ligand_atoms:
                    dist = distance(res_atom, lig_atom)
                    if nearest_distance is None or dist < nearest_distance:
                        nearest_distance = dist
                        nearest_id = ligand_key
                    if dist <= cutoff and (best_pair is None or dist < best_pair[2]):
                        best_pair = (res_atom, lig_atom, dist)
            if best_pair is None:
                continue
            res_atom, lig_atom, dist = best_pair
            contact_type = "metal_binding" if any(is_metal(atom) for atom in ligand_atoms) else "ligand_contact"
            residue.roles.add(contact_type)
            residue.sources.add("pdb_hetatm_distance")
            residue.ligand_ids.add(ligand_key)
            residue.ligand_contact_count += 1
            if contact_type == "metal_binding":
                residue.metal_contact = True
            contacts.append(
                {
                    "ligand_id": ligand_key,
                    "ligand_name": ligand_name,
                    "ligand_chain": lig_atom.chain,
                    "ligand_residue_number": lig_atom.residue_number,
                    "ligand_atom": lig_atom.atom_name,
                    "residue_id": residue.residue_id,
                    "chain": residue.chain,
                    "residue_number": residue.number,
                    "residue_name": residue.name,
                    "atom_name": res_atom.atom_name,
                    "distance": f"{dist:.3f}",
                    "contact_type": contact_type,
                }
            )
        residue.nearest_ligand_id = nearest_id
        residue.nearest_ligand_distance = nearest_distance
    return contacts


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path or not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            return []
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def apply_external_annotations(residues: dict[tuple[str, str, str], Residue], path: str) -> int:
    if not path:
        return 0
    rows = read_tsv(Path(path))
    applied = 0
    for row in rows:
        chain = row.get("chain") or row.get("chain_id") or "A"
        number = row.get("residue_number") or row.get("residue_index") or row.get("resnum") or ""
        insertion = row.get("insertion_code", "")
        if not number:
            continue
        key = (normalize_chain(chain), number, insertion)
        residue = residues.get(key)
        if residue is None:
            name = (row.get("residue_name") or row.get("resname") or "").upper()
            residue = residues.setdefault(key, Residue(chain=key[0], number=number, insertion_code=insertion, name=name))
        role = row.get("role") or row.get("type") or row.get("feature") or "functional_annotation"
        source = row.get("source") or row.get("database") or "uploaded_annotation"
        evidence = row.get("evidence") or row.get("note") or row.get("description") or ""
        residue.roles.add(role)
        residue.sources.add(source)
        if evidence:
            residue.annotations.append(evidence)
        if row.get("protect_from_mutation", "").lower() in {"yes", "true", "1"}:
            residue.protect = True
        applied += 1
    return applied


def should_protect(residue: Residue, protect_terms: list[str]) -> bool:
    if residue.protect:
        return True
    role_text = ",".join(residue.roles).lower()
    return any(term and term.lower() in role_text for term in protect_terms)


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_annotation_rows(residues: dict[tuple[str, str, str], Residue], protect_terms: list[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for residue in sorted(residues.values(), key=lambda item: (item.chain, safe_int(item.number), item.insertion_code)):
        if not residue.roles:
            continue
        residue.protect = should_protect(residue, protect_terms)
        rows.append(
            {
                "residue_id": residue.residue_id,
                "chain": residue.chain,
                "residue_number": residue.number,
                "insertion_code": residue.insertion_code,
                "residue_name": residue.name,
                "wt_aa": residue.wt_aa,
                "functional_role": ",".join(sorted(residue.roles)),
                "protect_from_mutation": "yes" if residue.protect else "no",
                "evidence_sources": ",".join(sorted(residue.sources)),
                "evidence_count": str(len(residue.sources) + len(residue.ligand_ids) + len(residue.annotations)),
                "nearest_ligand_id": residue.nearest_ligand_id,
                "nearest_ligand_distance": "" if residue.nearest_ligand_distance is None else f"{residue.nearest_ligand_distance:.3f}",
                "ligand_ids": ",".join(sorted(residue.ligand_ids)),
                "ligand_contact_count": str(residue.ligand_contact_count),
                "metal_contact": "yes" if residue.metal_contact else "no",
                "annotations": "; ".join(dict.fromkeys(residue.annotations)),
                "notes": "; ".join(dict.fromkeys(residue.notes)),
            }
        )
    return rows


def build_active_site_rows(rows: list[dict[str, str]], include_ligand_contacts: bool) -> list[dict[str, str]]:
    active_rows: list[dict[str, str]] = []
    for row in rows:
        roles = row["functional_role"].split(",")
        if row["protect_from_mutation"] != "yes" and not include_ligand_contacts:
            continue
        if not include_ligand_contacts and not any("catal" in role or "active" in role or "metal" in role for role in roles):
            continue
        active_rows.append(
            {
                "site_id": row["residue_id"],
                "chain": row["chain"],
                "residue_number": row["residue_number"],
                "residue_name": row["residue_name"],
                "role": row["functional_role"],
                "x": "",
                "y": "",
                "z": "",
                "evidence_source": row["evidence_sources"],
                "note": row["annotations"] or row["ligand_ids"],
            }
        )
    return active_rows


def safe_int(value: str) -> int:
    try:
        return int(value)
    except ValueError:
        return 0


def render_report(path: Path, annotation_rows: list[dict[str, str]], contact_rows: list[dict[str, str]], summary: dict[str, str]) -> None:
    html_text = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Functional Residue Annotator Report</title>
<style>
body {{ font-family: system-ui, sans-serif; margin: 2rem; background:#f6f8fb; color:#172033; }}
table {{ border-collapse: collapse; width: 100%; background: white; }}
th, td {{ border: 1px solid #d8dee9; padding: 6px 8px; text-align: left; font-size: 13px; }}
th {{ background: #172033; color: white; }}
.card {{ display:inline-block; background:white; border:1px solid #d8dee9; border-radius:12px; padding:14px; margin: 0 10px 12px 0; }}
</style></head><body>
<h1>Functional Residue Annotator Report</h1>
{''.join(f'<div class="card"><strong>{html.escape(k)}</strong><br>{html.escape(str(v))}</div>' for k, v in summary.items())}
<h2>Functional Residues</h2>
{table_html(annotation_rows, ANNOTATION_COLUMNS)}
<h2>Ligand / Metal Contacts</h2>
{table_html(contact_rows[:100], CONTACT_COLUMNS)}
</body></html>"""
    path.write_text(html_text, encoding="utf-8")


def table_html(rows: list[dict[str, str]], columns: list[str]) -> str:
    header = "".join(f"<th>{html.escape(col)}</th>" for col in columns)
    body = []
    for row in rows:
        body.append("<tr>" + "".join(f"<td>{html.escape(row.get(col, ''))}</td>" for col in columns) + "</tr>")
    if not body:
        body.append(f"<tr><td colspan='{len(columns)}'>No rows.</td></tr>")
    return f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def main() -> int:
    args = parse_args()
    log: list[str] = []
    try:
        ignore = {item.strip().upper() for item in args.ignore_ligands.split(",") if item.strip()}
        protect_terms = [item.strip() for item in args.protect_roles.split(",") if item.strip()]
        residues, ligands = read_structure(Path(args.input_pdb), ignore)
        contact_rows = annotate_contacts(residues, ligands, args.ligand_distance_cutoff, args.metal_distance_cutoff)
        external_count = apply_external_annotations(residues, args.functional_annotation_tsv)
        annotation_rows = build_annotation_rows(residues, protect_terms)
        active_rows = build_active_site_rows(annotation_rows, args.include_ligand_contacts_in_active_site)
        write_tsv(Path(args.residue_annotation_tsv), ANNOTATION_COLUMNS, annotation_rows)
        write_tsv(Path(args.active_site_tsv), ACTIVE_SITE_COLUMNS, active_rows)
        write_tsv(Path(args.ligand_contacts_tsv), CONTACT_COLUMNS, contact_rows)
        summary = {
            "protein_residues": str(len(residues)),
            "ligand_atoms": str(len(ligands)),
            "functional_residues": str(len(annotation_rows)),
            "active_site_rows": str(len(active_rows)),
            "ligand_contacts": str(len(contact_rows)),
            "external_annotations": str(external_count),
        }
        render_report(Path(args.report_html), annotation_rows, contact_rows, summary)
        Path(args.run_log).write_text("\n".join([f"{k}={v}" for k, v in summary.items()]) + "\n", encoding="utf-8")
        return 0
    except Exception as exc:
        Path(args.run_log).write_text(f"ERROR: {exc}\n", encoding="utf-8")
        print(f"functional_residue_annotator: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
