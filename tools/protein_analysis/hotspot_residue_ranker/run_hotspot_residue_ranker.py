#!/usr/bin/env python3
"""Rank hotspot residues by integrating structure, conservation, and pocket data."""

from __future__ import annotations

import argparse
import csv
import html
import math
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


OUTPUT_COLUMNS = [
    "rank",
    "residue_id",
    "chain",
    "residue_number",
    "insertion_code",
    "wt_aa",
    "residue_name",
    "hotspot_score",
    "hotspot_class",
    "recommendation",
    "pocket_or_tunnel_score",
    "mutability_score",
    "surface_accessibility_score",
    "flexibility_score",
    "active_site_distance_score",
    "catalytic_core_penalty",
    "highly_conserved_penalty",
    "buried_core_penalty",
    "conservation_score",
    "msa_entropy",
    "consensus_aa",
    "accepted_aas",
    "relative_asa",
    "freesasa_relative_asa",
    "exposure_class",
    "secondary_structure_class",
    "avg_b_factor",
    "in_pocket",
    "pocket_ids",
    "in_tunnel",
    "tunnel_ids",
    "nearest_active_site_id",
    "nearest_active_site_distance",
    "quality_flags",
    "notes",
]

AA1_TO_3 = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
    "U": "SEC",
    "O": "PYL",
    "X": "UNK",
}

AA3_TO_1 = {value: key for key, value in AA1_TO_3.items()}


@dataclass
class ResidueRecord:
    residue_id: str
    chain: str
    residue_number: str
    insertion_code: str = ""
    wt_aa: str = ""
    residue_name: str = ""
    structure: dict[str, str] | None = None
    conservation: dict[str, str] | None = None
    pocket_annotation: dict[str, str] | None = None
    nearest_active_site_id: str = ""
    nearest_active_site_distance: float | None = None
    active_site_role: str = ""


class ToolError(RuntimeError):
    """User-facing error."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Rank hotspot residues from integrated HotSpot workflow outputs.")
    parser.add_argument("--structure-features", required=True, help="Tool 2 residue_structure_features.tsv.")
    parser.add_argument("--conservation-tsv", required=True, help="Tool 4 residue_conservation.tsv.")
    parser.add_argument("--pockets-tsv", required=True, help="Tool 5 pockets.tsv.")
    parser.add_argument("--tunnels-tsv", required=True, help="Tool 5 tunnels.tsv.")
    parser.add_argument("--output-tsv", required=True, help="Ranked hotspot table.")
    parser.add_argument("--report-html", required=True, help="HTML summary report.")
    parser.add_argument("--viewer-html", required=True, help="HTML viewer scaffold.")
    parser.add_argument("--run-log", required=True, help="Run log.")
    parser.add_argument("--active-site-tsv", default="", help="Optional active/catalytic site TSV.")
    parser.add_argument("--structure-pdb", default="", help="Optional PDB path for Mol* viewer embedding.")

    parser.add_argument("--pocket-or-tunnel-weight", type=float, default=0.25)
    parser.add_argument("--mutability-weight", type=float, default=0.25)
    parser.add_argument("--surface-accessibility-weight", type=float, default=0.15)
    parser.add_argument("--flexibility-weight", type=float, default=0.1)
    parser.add_argument("--active-site-distance-weight", type=float, default=0.1)
    parser.add_argument("--catalytic-core-penalty", type=float, default=0.25)
    parser.add_argument("--highly-conserved-penalty", type=float, default=0.2)
    parser.add_argument("--buried-core-penalty", type=float, default=0.15)

    parser.add_argument("--high-conservation-threshold", type=float, default=0.85)
    parser.add_argument("--buried-threshold", type=float, default=0.09)
    parser.add_argument("--flexibility-b-factor-threshold", type=float, default=70.0)
    parser.add_argument("--active-site-near-distance", type=float, default=6.0)
    parser.add_argument("--active-site-caution-distance", type=float, default=10.0)
    parser.add_argument("--top-n-report", type=int, default=50)
    return parser.parse_args()


def read_tsv(path: Path, required: bool = True) -> list[dict[str, str]]:
    if not path or not path.exists() or path.stat().st_size == 0:
        if required:
            raise ToolError(f"Missing or empty TSV file: {path}")
        return []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return []
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def norm_float(value: str | float | None, default: float | None = None) -> float | None:
    if value is None or value == "":
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def clamp(value: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, value))


def make_residue_id(chain: str, residue_number: str, residue_name: str = "", insertion_code: str = "") -> str:
    if residue_name and len(residue_name) == 3:
        return f"{chain}:{residue_name}{residue_number}{insertion_code}"
    return f"{chain}:{residue_number}{insertion_code}"


def parse_residue_id(residue_id: str) -> tuple[str, str, str, str]:
    residue_id = residue_id.strip()
    match = None
    for pattern in [
        r"^([^:]+):([A-Za-z0-9]{3})(-?\d+)([A-Za-z]?)$",
        r"^([^:]+):(-?\d+)([A-Za-z]?)$",
    ]:
        match = __import__("re").match(pattern, residue_id)
        if match:
            break
    if match is None:
        return "", "", "", ""
    if len(match.groups()) == 4:
        return match.group(1), match.group(3), match.group(4), match.group(2).upper()
    return match.group(1), match.group(2), match.group(3), ""


def canonical_key(chain: str, number: str, insertion_code: str = "") -> str:
    return f"{chain}:{number}{insertion_code}"


def build_records(args: argparse.Namespace, log_lines: list[str]) -> dict[str, ResidueRecord]:
    records: dict[str, ResidueRecord] = {}

    for row in read_tsv(Path(args.structure_features)):
        chain = row.get("chain") or row.get("chain_id") or "_"
        number = row.get("residue_number") or row.get("resnum") or ""
        insertion = row.get("insertion_code", "")
        aa = row.get("amino_acid") or row.get("wt_aa") or ""
        residue_name = AA1_TO_3.get(aa.upper(), aa.upper()[:3] if aa else "")
        key = canonical_key(chain, number, insertion)
        record = records.setdefault(
            key,
            ResidueRecord(
                residue_id=row.get("residue_id") or make_residue_id(chain, number, "", insertion),
                chain=chain,
                residue_number=number,
                insertion_code=insertion,
                wt_aa=aa,
                residue_name=residue_name,
            ),
        )
        record.structure = row
        if aa and not record.wt_aa:
            record.wt_aa = aa
        if residue_name and not record.residue_name:
            record.residue_name = residue_name

    for row in read_tsv(Path(args.conservation_tsv)):
        chain = row.get("chain") or "_"
        number = row.get("residue_index") or row.get("residue_number") or ""
        insertion = row.get("insertion_code", "")
        aa = row.get("wt_aa", "")
        key = canonical_key(chain, number, insertion)
        record = records.setdefault(
            key,
            ResidueRecord(
                residue_id=make_residue_id(chain, number, "", insertion),
                chain=chain,
                residue_number=number,
                insertion_code=insertion,
                wt_aa=aa,
                residue_name=AA1_TO_3.get(aa.upper(), ""),
            ),
        )
        record.conservation = row
        if aa and not record.wt_aa:
            record.wt_aa = aa
            record.residue_name = AA1_TO_3.get(aa.upper(), "")

    pocket_annotation = build_pocket_annotation(Path(args.pockets_tsv), Path(args.tunnels_tsv))
    for residue_id, annotation in pocket_annotation.items():
        chain, number, insertion, residue_name = parse_residue_id(residue_id)
        if not chain:
            continue
        key = canonical_key(chain, number, insertion)
        record = records.setdefault(
            key,
            ResidueRecord(
                residue_id=residue_id,
                chain=chain,
                residue_number=number,
                insertion_code=insertion,
                residue_name=residue_name,
                wt_aa=AA3_TO_1.get(residue_name, ""),
            ),
        )
        record.pocket_annotation = annotation
        if residue_name and not record.residue_name:
            record.residue_name = residue_name
            record.wt_aa = AA3_TO_1.get(residue_name, record.wt_aa)

    active_sites = read_active_sites(Path(args.active_site_tsv)) if args.active_site_tsv else []
    assign_active_site_distances(records, active_sites)
    log_lines.append(f"residue_records={len(records)}")
    log_lines.append(f"active_site_entries={len(active_sites)}")
    return records


def build_pocket_annotation(pockets_path: Path, tunnels_path: Path) -> dict[str, dict[str, str]]:
    annotations: dict[str, dict[str, str]] = defaultdict(lambda: {"pocket_ids": "", "tunnel_ids": ""})
    for row in read_tsv(pockets_path, required=False):
        pocket_id = row.get("pocket_id", "")
        for residue_id in split_ids(row.get("residue_ids", "")):
            ann = annotations[residue_id]
            ann["pocket_ids"] = append_csv_value(ann.get("pocket_ids", ""), pocket_id)
    for row in read_tsv(tunnels_path, required=False):
        if row.get("status") == "not_computed":
            continue
        tunnel_id = row.get("tunnel_id", "")
        for residue_id in split_ids(row.get("residue_ids", "")):
            ann = annotations[residue_id]
            ann["tunnel_ids"] = append_csv_value(ann.get("tunnel_ids", ""), tunnel_id)
    return dict(annotations)


def split_ids(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def append_csv_value(existing: str, value: str) -> str:
    if not value:
        return existing
    values = split_ids(existing)
    if value not in values:
        values.append(value)
    return ",".join(values)


def read_active_sites(path: Path) -> list[dict[str, str]]:
    rows = read_tsv(path, required=False)
    normalized: list[dict[str, str]] = []
    for row in rows:
        chain = row.get("chain") or row.get("chain_id") or "_"
        number = row.get("residue_number") or row.get("residue_index") or row.get("resnum") or ""
        residue_id = row.get("residue_id") or make_residue_id(chain, number, row.get("residue_name", ""), row.get("insertion_code", ""))
        x = row.get("x") or row.get("coord_x") or row.get("center_x") or ""
        y = row.get("y") or row.get("coord_y") or row.get("center_y") or ""
        z = row.get("z") or row.get("coord_z") or row.get("center_z") or ""
        normalized.append(
            {
                "site_id": row.get("site_id") or row.get("id") or residue_id,
                "residue_id": residue_id,
                "role": row.get("role") or row.get("type") or "active_site",
                "x": x,
                "y": y,
                "z": z,
            }
        )
    return normalized


def assign_active_site_distances(records: dict[str, ResidueRecord], active_sites: list[dict[str, str]]) -> None:
    active_by_residue = {site["residue_id"]: site for site in active_sites if site.get("residue_id")}
    active_coords = [
        (site, (norm_float(site.get("x")), norm_float(site.get("y")), norm_float(site.get("z"))))
        for site in active_sites
        if site.get("x") and site.get("y") and site.get("z")
    ]
    active_coords = [(site, coords) for site, coords in active_coords if None not in coords]

    for record in records.values():
        by_exact = active_by_residue.get(record.residue_id)
        if not by_exact and record.residue_name:
            by_exact = active_by_residue.get(make_residue_id(record.chain, record.residue_number, record.residue_name, record.insertion_code))
        if by_exact:
            record.nearest_active_site_id = by_exact.get("site_id", "")
            record.nearest_active_site_distance = 0.0
            record.active_site_role = by_exact.get("role", "")
            continue
        # Coordinate distance is only possible if future tools provide residue centroids.
        # We keep the schema here; without residue coordinates this remains blank.
        if active_coords:
            continue


def score_record(record: ResidueRecord, args: argparse.Namespace) -> dict[str, str | float]:
    structure = record.structure or {}
    conservation = record.conservation or {}
    pocket = record.pocket_annotation or {}

    pocket_ids = pocket.get("pocket_ids", "")
    tunnel_ids = pocket.get("tunnel_ids", "")
    in_pocket = bool(split_ids(pocket_ids))
    in_tunnel = bool(split_ids(tunnel_ids))
    pocket_or_tunnel_score = 1.0 if in_tunnel else (0.8 if in_pocket else 0.0)

    mutability_score = norm_float(conservation.get("mutability_score"), 0.0) or 0.0
    conservation_score = norm_float(conservation.get("conservation_score"), 0.0) or 0.0
    entropy = norm_float(conservation.get("msa_entropy"), 0.0) or 0.0

    rel_asa = first_float(
        structure.get("freesasa_relative_asa"),
        structure.get("relative_asa"),
        default=0.0,
    )
    surface_accessibility_score = clamp(rel_asa / 0.36) if rel_asa is not None else 0.0

    avg_b = norm_float(structure.get("avg_b_factor"), None)
    flexibility_score = clamp(avg_b / args.flexibility_b_factor_threshold) if avg_b is not None else 0.0

    distance_score = 0.0
    catalytic_penalty = 0.0
    distance = record.nearest_active_site_distance
    if distance is not None:
        if distance <= args.active_site_near_distance:
            distance_score = 1.0
        elif distance <= args.active_site_caution_distance:
            distance_score = 0.5
        if distance == 0.0 or (record.active_site_role and "catal" in record.active_site_role.lower()):
            catalytic_penalty = 1.0

    highly_conserved_penalty = 1.0 if conservation_score >= args.high_conservation_threshold else 0.0
    buried_core_penalty = 0.0
    exposure = structure.get("freesasa_exposure_class") or structure.get("exposure_class") or ""
    if exposure == "buried" or (rel_asa is not None and rel_asa < args.buried_threshold):
        buried_core_penalty = 1.0

    score = (
        args.pocket_or_tunnel_weight * pocket_or_tunnel_score
        + args.mutability_weight * mutability_score
        + args.surface_accessibility_weight * surface_accessibility_score
        + args.flexibility_weight * flexibility_score
        + args.active_site_distance_weight * distance_score
        - args.catalytic_core_penalty * catalytic_penalty
        - args.highly_conserved_penalty * highly_conserved_penalty
        - args.buried_core_penalty * buried_core_penalty
    )
    hotspot_score = clamp(score)
    hotspot_class, recommendation, notes = classify_hotspot(
        in_pocket,
        in_tunnel,
        mutability_score,
        conservation_score,
        buried_core_penalty,
        catalytic_penalty,
        structure,
        conservation,
    )
    return {
        "residue_id": record.residue_id,
        "chain": record.chain,
        "residue_number": record.residue_number,
        "insertion_code": record.insertion_code,
        "wt_aa": record.wt_aa,
        "residue_name": record.residue_name,
        "hotspot_score": hotspot_score,
        "hotspot_class": hotspot_class,
        "recommendation": recommendation,
        "pocket_or_tunnel_score": pocket_or_tunnel_score,
        "mutability_score": mutability_score,
        "surface_accessibility_score": surface_accessibility_score,
        "flexibility_score": flexibility_score,
        "active_site_distance_score": distance_score,
        "catalytic_core_penalty": catalytic_penalty,
        "highly_conserved_penalty": highly_conserved_penalty,
        "buried_core_penalty": buried_core_penalty,
        "conservation_score": conservation_score,
        "msa_entropy": entropy,
        "consensus_aa": conservation.get("consensus_aa", ""),
        "accepted_aas": conservation.get("accepted_aas", ""),
        "relative_asa": structure.get("relative_asa", ""),
        "freesasa_relative_asa": structure.get("freesasa_relative_asa", ""),
        "exposure_class": exposure,
        "secondary_structure_class": structure.get("secondary_structure_class", ""),
        "avg_b_factor": structure.get("avg_b_factor", ""),
        "in_pocket": "yes" if in_pocket else "no",
        "pocket_ids": pocket_ids,
        "in_tunnel": "yes" if in_tunnel else "no",
        "tunnel_ids": tunnel_ids,
        "nearest_active_site_id": record.nearest_active_site_id,
        "nearest_active_site_distance": "" if record.nearest_active_site_distance is None else record.nearest_active_site_distance,
        "quality_flags": structure.get("quality_flags", ""),
        "notes": ";".join(notes),
    }


def first_float(*values: str | None, default: float | None = None) -> float | None:
    for value in values:
        parsed = norm_float(value, None)
        if parsed is not None:
            return parsed
    return default


def classify_hotspot(
    in_pocket: bool,
    in_tunnel: bool,
    mutability_score: float,
    conservation_score: float,
    buried_core_penalty: float,
    catalytic_penalty: float,
    structure: dict[str, str],
    conservation: dict[str, str],
) -> tuple[str, str, list[str]]:
    notes: list[str] = []
    if catalytic_penalty > 0 or conservation_score >= 0.9:
        notes.append("close_to_catalytic_or_highly_conserved")
        return "Risky hotspot", "caution", notes
    if buried_core_penalty > 0 and conservation_score >= 0.65:
        notes.append("buried_and_conserved")
        return "Risky hotspot", "avoid_or_validate", notes
    if (in_pocket or in_tunnel) and conservation_score < 0.85:
        notes.append("pocket_or_tunnel_nearby")
        return "Functional hotspot", "prioritize_for_function_screen", notes
    consensus = conservation.get("consensus_aa", "")
    wt = conservation.get("wt_aa", "")
    if consensus and wt and consensus != wt and conservation_score >= 0.5:
        notes.append("differs_from_family_consensus")
        return "Consensus hotspot", "consider_back_to_consensus", notes
    exposure = structure.get("freesasa_exposure_class") or structure.get("exposure_class") or ""
    if mutability_score >= 0.3 and exposure in {"exposed", "intermediate"}:
        notes.append("mutable_and_accessible")
        return "Stability hotspot", "consider_stability_library", notes
    return "Background", "lower_priority", notes


def rank_records(records: dict[str, ResidueRecord], args: argparse.Namespace) -> list[dict[str, str | float]]:
    scored = [score_record(record, args) for record in records.values()]
    scored.sort(
        key=lambda row: (
            -float(row["hotspot_score"]),
            str(row["hotspot_class"]) == "Risky hotspot",
            str(row["chain"]),
            int(row["residue_number"]) if str(row["residue_number"]).lstrip("-").isdigit() else 999999,
        )
    )
    for index, row in enumerate(scored, start=1):
        row["rank"] = index
    return scored


def write_ranked(rows: list[dict[str, str | float]], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: format_value(row.get(column, "")) for column in OUTPUT_COLUMNS})


def format_value(value: str | float | int) -> str | int:
    if isinstance(value, float):
        return f"{value:.6g}"
    return value


def write_report(rows: list[dict[str, str | float]], path: Path, top_n: int) -> None:
    counts = defaultdict(int)
    for row in rows:
        counts[str(row["hotspot_class"])] += 1
    top_rows = rows[:top_n]
    table_rows = "\n".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(str(format_value(row.get(column, ''))))}</td>"
            for column in [
                "rank",
                "residue_id",
                "hotspot_score",
                "hotspot_class",
                "recommendation",
                "conservation_score",
                "mutability_score",
                "in_pocket",
                "in_tunnel",
                "consensus_aa",
                "accepted_aas",
            ]
        )
        + "</tr>"
        for row in top_rows
    )
    summary_items = "\n".join(f"<li>{html.escape(key)}: {value}</li>" for key, value in sorted(counts.items()))
    path.write_text(
        f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>Hotspot Residue Ranker Report</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2rem; color: #172033; }}
    table {{ border-collapse: collapse; width: 100%; font-size: 0.9rem; }}
    th, td {{ border: 1px solid #d7dde8; padding: 0.35rem 0.5rem; text-align: left; }}
    th {{ background: #edf2f7; }}
    .card {{ border: 1px solid #d7dde8; padding: 1rem; border-radius: 8px; margin-bottom: 1rem; }}
  </style>
</head>
<body>
  <h1>Hotspot Residue Ranker Report</h1>
  <div class="card">
    <h2>Summary</h2>
    <p>Total ranked residues: {len(rows)}</p>
    <ul>{summary_items}</ul>
  </div>
  <h2>Top {min(top_n, len(rows))} Candidates</h2>
  <table>
    <thead>
      <tr><th>Rank</th><th>Residue</th><th>Score</th><th>Class</th><th>Recommendation</th><th>Conservation</th><th>Mutability</th><th>Pocket</th><th>Tunnel</th><th>Consensus</th><th>Accepted AAs</th></tr>
    </thead>
    <tbody>{table_rows}</tbody>
  </table>
</body>
</html>
""",
        encoding="utf-8",
    )


def write_viewer(rows: list[dict[str, str | float]], path: Path, structure_pdb: str) -> None:
    top_residues = [str(row["residue_id"]) for row in rows[:50]]
    residue_items = "\n".join(f"<li>{html.escape(residue)}</li>" for residue in top_residues)
    structure_note = (
        f"Structure dataset supplied: {html.escape(Path(structure_pdb).name)}"
        if structure_pdb
        else "No structure PDB was supplied to embed directly. Open the ranked table together with the original PDB in Galaxy visualization."
    )
    path.write_text(
        f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>Hotspot Viewer</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2rem; color: #172033; }}
    .viewer {{ border: 1px dashed #9aa7b7; border-radius: 8px; padding: 1rem; background: #f8fafc; }}
    code {{ background: #edf2f7; padding: 0.1rem 0.25rem; }}
  </style>
</head>
<body>
  <h1>Hotspot Viewer</h1>
  <div class="viewer">
    <p>{structure_note}</p>
    <p>This lightweight HTML lists top hotspot residues. For interactive structure coloring, open the PDB in Galaxy's structure viewer and use these residue IDs as selections.</p>
  </div>
  <h2>Top hotspot residues</h2>
  <ol>{residue_items}</ol>
</body>
</html>
""",
        encoding="utf-8",
    )


def main() -> int:
    args = parse_args()
    log_lines = ["Hotspot Residue Ranker"]
    try:
        records = build_records(args, log_lines)
        rows = rank_records(records, args)
        write_ranked(rows, Path(args.output_tsv))
        write_report(rows, Path(args.report_html), args.top_n_report)
        write_viewer(rows, Path(args.viewer_html), args.structure_pdb)
    except ToolError as exc:
        Path(args.run_log).write_text("\n".join(log_lines) + f"\nERROR: {exc}\n", encoding="utf-8")
        print(f"hotspot_residue_ranker: {exc}", file=sys.stderr)
        return 1
    Path(args.run_log).write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
