#!/usr/bin/env python3
"""Run fpocket and summarize pocket/cavity residues for HotSpot workflows."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from collections import OrderedDict, defaultdict
from dataclasses import dataclass, field
from pathlib import Path


POCKET_COLUMNS = [
    "pocket_id",
    "source",
    "score",
    "druggability_score",
    "number_of_alpha_spheres",
    "total_sasa",
    "polar_sasa",
    "apolar_sasa",
    "volume",
    "mean_local_hydrophobic_density",
    "mean_alpha_sphere_radius",
    "residue_count",
    "residue_ids",
    "pocket_file",
    "extra_attributes_json",
    "status",
    "warning",
]

TUNNEL_COLUMNS = [
    "tunnel_id",
    "source",
    "linked_pocket_id",
    "start_x",
    "start_y",
    "start_z",
    "bottleneck_radius",
    "length",
    "curvature",
    "residue_ids",
    "status",
    "warning",
]

ANNOTATION_COLUMNS = [
    "chain",
    "residue_number",
    "insertion_code",
    "residue_name",
    "residue_id",
    "in_pocket",
    "pocket_ids",
    "pocket_count",
    "in_tunnel",
    "tunnel_ids",
    "source",
]


@dataclass
class PocketSummary:
    pocket_id: str
    attributes: dict[str, str] = field(default_factory=dict)
    residue_ids: list[str] = field(default_factory=list)
    pocket_file: str = ""
    status: str = "ok"
    warning: str = ""


@dataclass(frozen=True)
class ResidueKey:
    chain: str
    number: str
    insertion_code: str
    residue_name: str

    @property
    def residue_id(self) -> str:
        insertion = self.insertion_code if self.insertion_code else ""
        return f"{self.chain}:{self.residue_name}{self.number}{insertion}"


class ToolError(RuntimeError):
    """User-facing runtime error."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Find protein pockets/cavities with fpocket.")
    parser.add_argument("--input-pdb", required=True, help="Input PDB structure.")
    parser.add_argument("--pockets-tsv", required=True, help="Output pocket summary TSV.")
    parser.add_argument("--tunnels-tsv", required=True, help="Output tunnel summary TSV placeholder.")
    parser.add_argument("--residue-annotation-tsv", required=True, help="Output residue-pocket annotation TSV.")
    parser.add_argument("--pocket-structure-pdb", required=True, help="Output PDB with fpocket pocket annotations.")
    parser.add_argument("--fpocket-info", required=True, help="Preserved fpocket info text output.")
    parser.add_argument("--raw-archive", required=True, help="ZIP archive of complete fpocket output directory.")
    parser.add_argument("--run-log", required=True, help="Run log.")
    parser.add_argument("--fpocket-binary", default=os.environ.get("FPOCKET_BINARY", "fpocket"))
    parser.add_argument("--min-alpha-spheres", type=int, default=0, help="Optional minimum alpha spheres filter.")
    parser.add_argument("--top-n-pockets", type=int, default=0, help="Keep top N pockets in summary; 0 keeps all.")
    return parser.parse_args()


def run_command(command: list[str], cwd: Path, log_lines: list[str]) -> subprocess.CompletedProcess[str]:
    log_lines.append("$ " + " ".join(command))
    try:
        completed = subprocess.run(command, cwd=str(cwd), capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise ToolError(f"Executable not found: {command[0]}. Install fpocket or configure FPOCKET_BINARY.") from exc
    log_lines.append(f"exit_code={completed.returncode}")
    if completed.stdout:
        log_lines.append("STDOUT:\n" + completed.stdout)
    if completed.stderr:
        log_lines.append("STDERR:\n" + completed.stderr)
    return completed


def copy_input_pdb(input_path: Path, work_dir: Path) -> Path:
    destination = work_dir / "input.pdb"
    with input_path.open("r", encoding="utf-8", errors="replace") as src, destination.open(
        "w", encoding="utf-8", newline="\n"
    ) as dst:
        for raw_line in src:
            dst.write(raw_line.rstrip("\n") + "\n")
    return destination


def run_fpocket(args: argparse.Namespace, work_dir: Path, log_lines: list[str]) -> Path:
    input_pdb = copy_input_pdb(Path(args.input_pdb), work_dir)
    completed = run_command([args.fpocket_binary, "-f", str(input_pdb.name)], work_dir, log_lines)
    if completed.returncode != 0:
        raise ToolError("fpocket failed; see run log for details.")
    output_dir = work_dir / "input_out"
    if not output_dir.exists():
        candidates = list(work_dir.glob("*_out"))
        if not candidates:
            raise ToolError("fpocket completed but no output directory was found.")
        output_dir = candidates[0]
    return output_dir


def parse_info_file(info_path: Path) -> OrderedDict[str, PocketSummary]:
    pockets: OrderedDict[str, PocketSummary] = OrderedDict()
    current: PocketSummary | None = None
    key_value_re = re.compile(r"^\s*([^:]+?)\s*:\s*(.*?)\s*$")
    pocket_re = re.compile(r"^\s*Pocket\s+(\d+)\s*:?\s*$", re.IGNORECASE)

    if not info_path.exists():
        return pockets

    with info_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            pocket_match = pocket_re.match(line)
            if pocket_match:
                pocket_id = pocket_match.group(1)
                current = pockets.setdefault(pocket_id, PocketSummary(pocket_id=pocket_id))
                continue
            if current is None:
                continue
            match = key_value_re.match(line)
            if match:
                key = normalize_key(match.group(1))
                value = match.group(2).strip()
                current.attributes[key] = value
    return pockets


def normalize_key(key: str) -> str:
    key = key.strip().lower()
    key = key.replace("(", "").replace(")", "")
    key = re.sub(r"[^a-z0-9]+", "_", key)
    return key.strip("_")


def find_fpocket_info(output_dir: Path) -> Path | None:
    candidates = sorted(output_dir.glob("*_info.txt"))
    if candidates:
        return candidates[0]
    candidates = sorted(output_dir.glob("*info*.txt"))
    return candidates[0] if candidates else None


def find_pocket_files(output_dir: Path) -> list[Path]:
    pocket_dir = output_dir / "pockets"
    candidates: list[Path] = []
    if pocket_dir.exists():
        candidates.extend(sorted(pocket_dir.glob("pocket*_atm.pdb")))
        candidates.extend(path for path in sorted(pocket_dir.glob("pocket*.pdb")) if path not in candidates)
    return candidates


def pocket_id_from_path(path: Path) -> str:
    match = re.search(r"pocket(\d+)", path.name, re.IGNORECASE)
    return match.group(1) if match else path.stem


def parse_pocket_residues(path: Path) -> list[ResidueKey]:
    residues: OrderedDict[tuple[str, str, str, str], ResidueKey] = OrderedDict()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            parsed = parse_pdb_residue(line)
            if parsed is None:
                continue
            residues[(parsed.chain, parsed.number, parsed.insertion_code, parsed.residue_name)] = parsed
    return list(residues.values())


def parse_pdb_residue(line: str) -> ResidueKey | None:
    if len(line) >= 27:
        residue_name = line[17:20].strip().upper()
        chain = line[21:22].strip() or "_"
        number = line[22:26].strip()
        insertion_code = line[26:27].strip()
        if residue_name and number:
            return ResidueKey(chain=chain, number=number, insertion_code=insertion_code, residue_name=residue_name)
    parts = line.split()
    if len(parts) >= 6:
        return ResidueKey(chain=parts[4] or "_", number=parts[5], insertion_code="", residue_name=parts[3].upper())
    return None


def enrich_pockets_from_files(pockets: OrderedDict[str, PocketSummary], output_dir: Path) -> OrderedDict[str, PocketSummary]:
    for pocket_file in find_pocket_files(output_dir):
        pocket_id = pocket_id_from_path(pocket_file)
        pocket = pockets.setdefault(pocket_id, PocketSummary(pocket_id=pocket_id))
        pocket.pocket_file = str(pocket_file.relative_to(output_dir))
        residues = parse_pocket_residues(pocket_file)
        pocket.residue_ids = [residue.residue_id for residue in residues]
    return pockets


def numeric_attr(attributes: dict[str, str], *keys: str) -> str:
    for key in keys:
        value = attributes.get(key)
        if value not in {None, ""}:
            return str(value)
    return ""


def pocket_sort_key(pocket: PocketSummary) -> tuple[float, int]:
    score = parse_float(
        numeric_attr(
            pocket.attributes,
            "score",
            "pocket_score",
            "druggability_score",
            "druggability_score_score",
        )
    )
    try:
        pocket_num = int(pocket.pocket_id)
    except ValueError:
        pocket_num = 999999
    return (-(score if score is not None else -1e9), pocket_num)


def filter_pockets(pockets: OrderedDict[str, PocketSummary], min_alpha_spheres: int, top_n: int) -> list[PocketSummary]:
    pocket_list = list(pockets.values())
    if min_alpha_spheres > 0:
        pocket_list = [
            pocket
            for pocket in pocket_list
            if (parse_float(numeric_attr(pocket.attributes, "number_of_alpha_spheres", "alpha_spheres")) or 0)
            >= min_alpha_spheres
        ]
    pocket_list = sorted(pocket_list, key=pocket_sort_key)
    if top_n > 0:
        pocket_list = pocket_list[:top_n]
    return pocket_list


def parse_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def write_pockets_tsv(pockets: list[PocketSummary], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=POCKET_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for pocket in pockets:
            attrs = pocket.attributes
            known_keys = {
                "score",
                "pocket_score",
                "druggability_score",
                "number_of_alpha_spheres",
                "alpha_spheres",
                "total_sasa",
                "polar_sasa",
                "apolar_sasa",
                "volume",
                "mean_local_hydrophobic_density",
                "mean_alpha_sphere_radius",
            }
            extras = {key: value for key, value in attrs.items() if key not in known_keys}
            writer.writerow(
                {
                    "pocket_id": pocket.pocket_id,
                    "source": "fpocket",
                    "score": numeric_attr(attrs, "score", "pocket_score"),
                    "druggability_score": numeric_attr(attrs, "druggability_score"),
                    "number_of_alpha_spheres": numeric_attr(attrs, "number_of_alpha_spheres", "alpha_spheres"),
                    "total_sasa": numeric_attr(attrs, "total_sasa"),
                    "polar_sasa": numeric_attr(attrs, "polar_sasa"),
                    "apolar_sasa": numeric_attr(attrs, "apolar_sasa"),
                    "volume": numeric_attr(attrs, "volume"),
                    "mean_local_hydrophobic_density": numeric_attr(attrs, "mean_local_hydrophobic_density"),
                    "mean_alpha_sphere_radius": numeric_attr(attrs, "mean_alpha_sphere_radius"),
                    "residue_count": len(pocket.residue_ids),
                    "residue_ids": ",".join(pocket.residue_ids),
                    "pocket_file": pocket.pocket_file,
                    "extra_attributes_json": json.dumps(extras, ensure_ascii=False, sort_keys=True),
                    "status": pocket.status,
                    "warning": pocket.warning,
                }
            )


def write_tunnels_placeholder(path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TUNNEL_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(
            {
                "tunnel_id": "",
                "source": "fpocket",
                "linked_pocket_id": "",
                "start_x": "",
                "start_y": "",
                "start_z": "",
                "bottleneck_radius": "",
                "length": "",
                "curvature": "",
                "residue_ids": "",
                "status": "not_computed",
                "warning": "fpocket detects pockets/cavities, not tunnels; CAVER integration will populate this table later.",
            }
        )


def write_residue_annotation(pockets: list[PocketSummary], path: Path) -> None:
    residue_to_pockets: dict[str, list[str]] = defaultdict(list)
    residue_name_by_id: dict[str, ResidueKey] = {}
    for pocket in pockets:
        for residue_id in pocket.residue_ids:
            residue_to_pockets[residue_id].append(pocket.pocket_id)
            parsed = residue_id_to_key(residue_id)
            if parsed:
                residue_name_by_id[residue_id] = parsed

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ANNOTATION_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for residue_id in sorted(residue_to_pockets, key=residue_sort_key):
            key = residue_name_by_id.get(residue_id)
            pocket_ids = sorted(residue_to_pockets[residue_id], key=lambda value: int(value) if value.isdigit() else 999999)
            writer.writerow(
                {
                    "chain": key.chain if key else "",
                    "residue_number": key.number if key else "",
                    "insertion_code": key.insertion_code if key else "",
                    "residue_name": key.residue_name if key else "",
                    "residue_id": residue_id,
                    "in_pocket": "yes",
                    "pocket_ids": ",".join(pocket_ids),
                    "pocket_count": len(pocket_ids),
                    "in_tunnel": "no",
                    "tunnel_ids": "",
                    "source": "fpocket",
                }
            )


def residue_id_to_key(residue_id: str) -> ResidueKey | None:
    match = re.match(r"([^:]+):([A-Za-z0-9]{3})(-?\d+)([A-Za-z]?)$", residue_id)
    if not match:
        return None
    return ResidueKey(
        chain=match.group(1),
        residue_name=match.group(2),
        number=match.group(3),
        insertion_code=match.group(4),
    )


def residue_sort_key(residue_id: str) -> tuple[str, int, str, str]:
    key = residue_id_to_key(residue_id)
    if key is None:
        return ("", 999999, "", residue_id)
    try:
        number = int(key.number)
    except ValueError:
        number = 999999
    return (key.chain, number, key.insertion_code, key.residue_name)


def copy_pocket_structure(output_dir: Path, destination: Path) -> None:
    candidates = sorted(output_dir.glob("*_out.pdb")) + sorted(output_dir.glob("*.pdb"))
    if candidates:
        shutil.copyfile(candidates[0], destination)
        return

    pocket_files = find_pocket_files(output_dir)
    if pocket_files:
        with destination.open("w", encoding="utf-8", newline="\n") as out:
            for pocket_file in pocket_files:
                out.write(f"MODEL     {pocket_id_from_path(pocket_file)}\n")
                with pocket_file.open("r", encoding="utf-8", errors="replace") as handle:
                    for line in handle:
                        if line.startswith(("ATOM", "HETATM", "TER")):
                            out.write(line.rstrip("\n") + "\n")
                out.write("ENDMDL\n")
            out.write("END\n")
        return
    destination.write_text("HEADER    NO FPOCKET POCKET STRUCTURE AVAILABLE\nEND\n", encoding="utf-8")


def zip_directory(source_dir: Path, destination_zip: Path) -> None:
    with zipfile.ZipFile(destination_zip, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(source_dir.rglob("*")):
            if path.is_file():
                archive.write(path, arcname=str(path.relative_to(source_dir.parent)))


def main() -> int:
    args = parse_args()
    log_lines = [
        "Pocket and Tunnel Finder",
        "enabled_backend=fpocket",
        f"fpocket_binary={args.fpocket_binary}",
    ]
    try:
        with tempfile.TemporaryDirectory(prefix="fpocket_run_") as tmp_name:
            work_dir = Path(tmp_name)
            output_dir = run_fpocket(args, work_dir, log_lines)
            info_path = find_fpocket_info(output_dir)
            pockets: OrderedDict[str, PocketSummary] = OrderedDict()
            if info_path:
                shutil.copyfile(info_path, args.fpocket_info)
                pockets = parse_info_file(info_path)
            else:
                Path(args.fpocket_info).write_text("No fpocket info file was found.\n", encoding="utf-8")
                log_lines.append("WARNING: no fpocket info file was found.")

            pockets = enrich_pockets_from_files(pockets, output_dir)
            filtered_pockets = filter_pockets(pockets, args.min_alpha_spheres, args.top_n_pockets)
            log_lines.append(f"pockets_detected={len(pockets)}")
            log_lines.append(f"pockets_reported={len(filtered_pockets)}")

            write_pockets_tsv(filtered_pockets, Path(args.pockets_tsv))
            write_tunnels_placeholder(Path(args.tunnels_tsv))
            write_residue_annotation(filtered_pockets, Path(args.residue_annotation_tsv))
            copy_pocket_structure(output_dir, Path(args.pocket_structure_pdb))
            zip_directory(output_dir, Path(args.raw_archive))
    except ToolError as exc:
        Path(args.run_log).write_text("\n".join(log_lines) + f"\nERROR: {exc}\n", encoding="utf-8")
        print(f"pocket_tunnel_finder: {exc}", file=sys.stderr)
        return 1

    Path(args.run_log).write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
