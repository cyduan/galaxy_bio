#!/usr/bin/env python3
"""Run fpocket/P2Rank/CAVER and summarize pocket/tunnel residues."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shlex
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
    "center_x",
    "center_y",
    "center_z",
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
    "tunnel_file",
    "extra_attributes_json",
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
    "tunnel_count",
    "source",
]


@dataclass
class PocketSummary:
    pocket_id: str
    source: str
    attributes: dict[str, str] = field(default_factory=dict)
    residue_ids: list[str] = field(default_factory=list)
    pocket_file: str = ""
    status: str = "ok"
    warning: str = ""


@dataclass
class TunnelSummary:
    tunnel_id: str
    source: str = "caver"
    linked_pocket_id: str = ""
    start_x: str = ""
    start_y: str = ""
    start_z: str = ""
    bottleneck_radius: str = ""
    length: str = ""
    curvature: str = ""
    residue_ids: list[str] = field(default_factory=list)
    tunnel_file: str = ""
    attributes: dict[str, str] = field(default_factory=dict)
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
    parser = argparse.ArgumentParser(description="Find protein pockets/cavities/tunnels.")
    parser.add_argument("--input-pdb", required=True, help="Input PDB structure.")
    parser.add_argument("--pockets-tsv", required=True, help="Output pocket summary TSV.")
    parser.add_argument("--tunnels-tsv", required=True, help="Output tunnel summary TSV.")
    parser.add_argument("--residue-annotation-tsv", required=True, help="Output residue pocket/tunnel annotation TSV.")
    parser.add_argument("--pocket-structure-pdb", required=True, help="Output PDB with pocket/tunnel annotations.")
    parser.add_argument("--fpocket-info", required=True, help="Preserved fpocket info text output.")
    parser.add_argument("--p2rank-predictions", required=True, help="Preserved P2Rank predictions CSV.")
    parser.add_argument("--p2rank-residues", required=True, help="Preserved P2Rank residues CSV.")
    parser.add_argument("--caver-config", required=True, help="CAVER config used for the run.")
    parser.add_argument("--raw-archive", required=True, help="ZIP archive of complete raw output directories.")
    parser.add_argument("--run-log", required=True, help="Run log.")

    parser.add_argument("--run-fpocket", action="store_true", help="Run fpocket.")
    parser.add_argument("--run-p2rank", action="store_true", help="Run P2Rank.")
    parser.add_argument("--run-caver", action="store_true", help="Run CAVER.")

    parser.add_argument("--fpocket-binary", default=os.environ.get("FPOCKET_BINARY", "fpocket"))
    parser.add_argument("--p2rank-binary", default=os.environ.get("P2RANK_BINARY", "prank"))
    parser.add_argument("--caver-home", default=os.environ.get("CAVER_HOME", "/data/tools/caver_3.0/caver"))
    parser.add_argument("--caver-jar", default=os.environ.get("CAVER_JAR", ""))
    parser.add_argument("--java-binary", default=os.environ.get("JAVA_BINARY", "java"))
    parser.add_argument("--caver-java-mem", default=os.environ.get("CAVER_JAVA_MEM", "4000m"))
    parser.add_argument("--caver-command-template", default=os.environ.get("CAVER_COMMAND", ""))

    parser.add_argument("--threads", type=int, default=1, help="Threads for P2Rank when supported.")
    parser.add_argument("--p2rank-config", default="", help="Optional P2Rank config/profile, e.g. alphafold.")
    parser.add_argument("--min-alpha-spheres", type=int, default=0, help="Optional fpocket minimum alpha spheres filter.")
    parser.add_argument("--top-n-pockets", type=int, default=0, help="Keep top N pockets in summary; 0 keeps all.")
    parser.add_argument("--caver-start-mode", choices=["manual", "auto_centroid"], default="auto_centroid")
    parser.add_argument("--caver-start-x", default="")
    parser.add_argument("--caver-start-y", default="")
    parser.add_argument("--caver-start-z", default="")
    parser.add_argument("--caver-probe-radius", default="0.9")
    parser.add_argument("--caver-shell-radius", default="3")
    parser.add_argument("--caver-shell-depth", default="4")
    return parser.parse_args()


def run_command(command: list[str], cwd: Path, log_lines: list[str], label: str) -> subprocess.CompletedProcess[str]:
    log_lines.append(f"[{label}] $ " + " ".join(command))
    try:
        completed = subprocess.run(command, cwd=str(cwd), capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise ToolError(f"{label}: executable not found: {command[0]}") from exc
    log_lines.append(f"[{label}] exit_code={completed.returncode}")
    if completed.stdout:
        log_lines.append(f"[{label}] STDOUT:\n{completed.stdout}")
    if completed.stderr:
        log_lines.append(f"[{label}] STDERR:\n{completed.stderr}")
    return completed


def copy_input_pdb(input_path: Path, work_dir: Path, name: str = "input.pdb") -> Path:
    destination = work_dir / name
    with input_path.open("r", encoding="utf-8", errors="replace") as src, destination.open(
        "w", encoding="utf-8", newline="\n"
    ) as dst:
        for raw_line in src:
            dst.write(raw_line.rstrip("\n") + "\n")
    return destination


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


def normalize_key(key: str) -> str:
    key = key.strip().lower()
    key = key.replace("(", "").replace(")", "")
    key = re.sub(r"[^a-z0-9]+", "_", key)
    return key.strip("_")


def parse_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def numeric_attr(attributes: dict[str, str], *keys: str) -> str:
    for key in keys:
        value = attributes.get(key)
        if value not in {None, ""}:
            return str(value)
    return ""


def prefixed_id(source: str, local_id: str) -> str:
    return f"{source}:{local_id}"


# ----------------------------- fpocket -----------------------------


def run_fpocket(args: argparse.Namespace, root_work_dir: Path, log_lines: list[str]) -> tuple[list[PocketSummary], Path | None]:
    work_dir = root_work_dir / "fpocket"
    work_dir.mkdir()
    input_pdb = copy_input_pdb(Path(args.input_pdb), work_dir)
    completed = run_command([args.fpocket_binary, "-f", str(input_pdb.name)], work_dir, log_lines, "fpocket")
    if completed.returncode != 0:
        raise ToolError("fpocket failed; see run log for details.")
    output_dir = work_dir / "input_out"
    if not output_dir.exists():
        candidates = list(work_dir.glob("*_out"))
        if not candidates:
            raise ToolError("fpocket completed but no output directory was found.")
        output_dir = candidates[0]
    info_path = find_fpocket_info(output_dir)
    pockets: OrderedDict[str, PocketSummary] = OrderedDict()
    if info_path:
        shutil.copyfile(info_path, args.fpocket_info)
        pockets = parse_fpocket_info_file(info_path)
    else:
        Path(args.fpocket_info).write_text("No fpocket info file was found.\n", encoding="utf-8")
        log_lines.append("WARNING: no fpocket info file was found.")
    enrich_fpocket_pockets_from_files(pockets, output_dir)
    filtered = filter_pockets(pockets, args.min_alpha_spheres, args.top_n_pockets)
    log_lines.append(f"fpocket_pockets_reported={len(filtered)}")
    return filtered, output_dir


def parse_fpocket_info_file(info_path: Path) -> OrderedDict[str, PocketSummary]:
    pockets: OrderedDict[str, PocketSummary] = OrderedDict()
    current: PocketSummary | None = None
    key_value_re = re.compile(r"^\s*([^:]+?)\s*:\s*(.*?)\s*$")
    pocket_re = re.compile(r"^\s*Pocket\s+(\d+)\s*:?\s*$", re.IGNORECASE)
    with info_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            pocket_match = pocket_re.match(line)
            if pocket_match:
                local_id = pocket_match.group(1)
                pocket_id = prefixed_id("fpocket", local_id)
                current = pockets.setdefault(pocket_id, PocketSummary(pocket_id=pocket_id, source="fpocket"))
                continue
            if current is None:
                continue
            match = key_value_re.match(line)
            if match:
                current.attributes[normalize_key(match.group(1))] = match.group(2).strip()
    return pockets


def find_fpocket_info(output_dir: Path) -> Path | None:
    candidates = sorted(output_dir.glob("*_info.txt"))
    if candidates:
        return candidates[0]
    candidates = sorted(output_dir.glob("*info*.txt"))
    return candidates[0] if candidates else None


def find_fpocket_pocket_files(output_dir: Path) -> list[Path]:
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


def enrich_fpocket_pockets_from_files(pockets: OrderedDict[str, PocketSummary], output_dir: Path) -> None:
    for pocket_file in find_fpocket_pocket_files(output_dir):
        local_id = pocket_id_from_path(pocket_file)
        pocket_id = prefixed_id("fpocket", local_id)
        pocket = pockets.setdefault(pocket_id, PocketSummary(pocket_id=pocket_id, source="fpocket"))
        pocket.pocket_file = str(pocket_file.relative_to(output_dir))
        pocket.residue_ids = [residue.residue_id for residue in parse_pocket_residues(pocket_file)]


def pocket_sort_key(pocket: PocketSummary) -> tuple[float, str]:
    score = parse_float(
        numeric_attr(pocket.attributes, "score", "pocket_score", "druggability_score", "probability")
    )
    return (-(score if score is not None else -1e9), pocket.pocket_id)


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


# ----------------------------- P2Rank -----------------------------


def run_p2rank(args: argparse.Namespace, root_work_dir: Path, log_lines: list[str]) -> tuple[list[PocketSummary], Path | None]:
    work_dir = root_work_dir / "p2rank"
    out_dir = work_dir / "p2rank_output"
    work_dir.mkdir()
    out_dir.mkdir()
    input_pdb = copy_input_pdb(Path(args.input_pdb), work_dir)
    command = [args.p2rank_binary, "predict", "-f", str(input_pdb), "-o", str(out_dir), "-threads", str(args.threads)]
    if args.p2rank_config:
        command.extend(["-c", args.p2rank_config])
    completed = run_command(command, work_dir, log_lines, "p2rank")
    if completed.returncode != 0:
        raise ToolError("P2Rank failed; see run log for details.")
    predictions = find_first(out_dir, ["*_predictions.csv", "*_prediction.csv", "predictions.csv"])
    residues = find_first(out_dir, ["*_residues.csv", "residues.csv"])
    if predictions:
        shutil.copyfile(predictions, args.p2rank_predictions)
    else:
        Path(args.p2rank_predictions).write_text("No P2Rank predictions CSV was found.\n", encoding="utf-8")
    if residues:
        shutil.copyfile(residues, args.p2rank_residues)
    else:
        Path(args.p2rank_residues).write_text("No P2Rank residues CSV was found.\n", encoding="utf-8")
    pockets = parse_p2rank_predictions(predictions) if predictions else []
    residue_map = parse_p2rank_residues(residues) if residues else {}
    for pocket in pockets:
        local_id = pocket.pocket_id.split(":", 1)[-1]
        if local_id in residue_map:
            pocket.residue_ids = residue_map[local_id]
    log_lines.append(f"p2rank_pockets_reported={len(pockets)}")
    return pockets, out_dir


def find_first(directory: Path, patterns: list[str]) -> Path | None:
    for pattern in patterns:
        candidates = sorted(directory.rglob(pattern))
        if candidates:
            return candidates[0]
    return None


def read_csv_dicts(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delimiter = "," if sample.count(",") >= sample.count("\t") else "\t"
        reader = csv.DictReader(handle, delimiter=delimiter, skipinitialspace=True)
        return [{normalize_key(k): (v or "").strip() for k, v in row.items()} for row in reader]


def parse_p2rank_predictions(path: Path) -> list[PocketSummary]:
    pockets: list[PocketSummary] = []
    for row in read_csv_dicts(path):
        local_id = first_value(row, "rank", "pocket", "pocket_id", "name") or str(len(pockets) + 1)
        pocket = PocketSummary(pocket_id=prefixed_id("p2rank", local_id), source="p2rank")
        pocket.attributes = dict(row)
        residues = first_value(row, "residue_ids", "residues", "residue_labels", "sas_points")
        pocket.residue_ids = parse_residue_list(residues)
        pocket.warning = "" if pocket.residue_ids else "no_residue_ids_parsed_from_predictions"
        pockets.append(pocket)
    return sorted(pockets, key=pocket_sort_key)


def parse_p2rank_residues(path: Path) -> dict[str, list[str]]:
    pocket_residues: dict[str, list[str]] = defaultdict(list)
    for row in read_csv_dicts(path):
        pocket_value = first_value(row, "pocket", "pocket_id", "prediction", "predicted_pocket", "rank")
        residue_id = residue_from_p2rank_row(row)
        if not pocket_value or not residue_id or pocket_value in {"0", "-", "none", "None"}:
            continue
        for local_id in re.split(r"[;,| ]+", pocket_value):
            local_id = local_id.strip()
            if local_id and local_id not in {"0", "-", "none", "None"}:
                pocket_residues[local_id].append(residue_id)
    return {pocket: sorted(set(residues), key=residue_sort_key) for pocket, residues in pocket_residues.items()}


def first_value(row: dict[str, str], *keys: str) -> str:
    for key in keys:
        value = row.get(key)
        if value not in {None, ""}:
            return str(value)
    return ""


def residue_from_p2rank_row(row: dict[str, str]) -> str:
    raw = first_value(row, "residue_id", "residue", "residue_label", "label")
    parsed = parse_residue_list(raw)
    if parsed:
        return parsed[0]
    chain = first_value(row, "chain", "chain_id")
    number = first_value(row, "residue_number", "res_seq", "residue_index", "resid", "residue")
    name = first_value(row, "residue_name", "res_name", "aa", "amino_acid") or "UNK"
    if number:
        return ResidueKey(chain=chain or "_", number=number, insertion_code="", residue_name=name[:3].upper()).residue_id
    return ""


def parse_residue_list(raw: str) -> list[str]:
    if not raw:
        return []
    candidates = re.split(r"[,;| ]+", raw.strip().strip("[]"))
    residues: list[str] = []
    for candidate in candidates:
        candidate = candidate.strip()
        if not candidate:
            continue
        parsed = parse_residue_token(candidate)
        if parsed:
            residues.append(parsed)
    return sorted(set(residues), key=residue_sort_key)


def parse_residue_token(token: str) -> str:
    token = token.strip()
    if re.match(r"^[^:]+:[A-Za-z0-9]{3}-?\d+[A-Za-z]?$", token):
        return token
    match = re.match(r"^([A-Za-z0-9_])[_:/.-]?([A-Za-z]{3})?[_:/.-]?(-?\d+)([A-Za-z]?)$", token)
    if match:
        chain = match.group(1) or "_"
        residue_name = (match.group(2) or "UNK").upper()
        number = match.group(3)
        insertion = match.group(4) or ""
        return ResidueKey(chain=chain, number=number, insertion_code=insertion, residue_name=residue_name).residue_id
    match = re.match(r"^([A-Za-z]{3})(-?\d+)([A-Za-z]?)$", token)
    if match:
        return ResidueKey(chain="_", residue_name=match.group(1).upper(), number=match.group(2), insertion_code=match.group(3)).residue_id
    return ""


# ----------------------------- CAVER -----------------------------


def run_caver(args: argparse.Namespace, root_work_dir: Path, log_lines: list[str]) -> tuple[list[TunnelSummary], Path | None]:
    work_dir = root_work_dir / "caver"
    pdb_dir = work_dir / "pdb"
    out_dir = work_dir / "caver_output"
    work_dir.mkdir()
    pdb_dir.mkdir()
    out_dir.mkdir()
    input_pdb = copy_input_pdb(Path(args.input_pdb), pdb_dir)
    start_x, start_y, start_z = caver_start_coordinates(args, input_pdb)
    config_path = work_dir / "config.txt"
    write_caver_config(args, config_path, start_x, start_y, start_z)
    shutil.copyfile(config_path, args.caver_config)

    if args.caver_command_template:
        command = shlex.split(
            args.caver_command_template.format(
                caver_home=args.caver_home,
                caver_jar=args.caver_jar or str(Path(args.caver_home) / "caver.jar"),
                java=args.java_binary,
                java_mem=args.caver_java_mem,
                pdb_dir=str(pdb_dir),
                config=str(config_path),
                out_dir=str(out_dir),
            )
        )
    else:
        caver_jar = args.caver_jar or str(Path(args.caver_home) / "caver.jar")
        command = [
            args.java_binary,
            f"-Xmx{args.caver_java_mem}",
            "-cp",
            str(Path(args.caver_home) / "lib"),
            "-jar",
            caver_jar,
            "-home",
            args.caver_home,
            "-pdb",
            str(pdb_dir),
            "-conf",
            str(config_path),
            "-out",
            str(out_dir),
        ]
    completed = run_command(command, work_dir, log_lines, "caver")
    if completed.returncode != 0:
        raise ToolError("CAVER failed; see run log for details.")
    tunnels = parse_caver_outputs(out_dir, start_x, start_y, start_z)
    if not tunnels:
        tunnels = [
            TunnelSummary(
                tunnel_id="caver:none",
                start_x=start_x,
                start_y=start_y,
                start_z=start_z,
                status="no_tunnels_parsed",
                warning="CAVER completed but no tunnel summary CSV/PDB residues were parsed.",
            )
        ]
    log_lines.append(f"caver_tunnels_reported={len(tunnels)}")
    return tunnels, out_dir


def caver_start_coordinates(args: argparse.Namespace, input_pdb: Path) -> tuple[str, str, str]:
    if args.caver_start_mode == "manual":
        if not (args.caver_start_x and args.caver_start_y and args.caver_start_z):
            raise ToolError("CAVER manual start mode requires x, y, and z coordinates.")
        return args.caver_start_x, args.caver_start_y, args.caver_start_z
    coords: list[tuple[float, float, float]] = []
    with input_pdb.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                continue
    if not coords:
        raise ToolError("Could not calculate automatic CAVER start point; no atom coordinates parsed.")
    return tuple(f"{sum(values) / len(values):.3f}" for values in zip(*coords))  # type: ignore[return-value]


def write_caver_config(args: argparse.Namespace, path: Path, x: str, y: str, z: str) -> None:
    path.write_text(
        "\n".join(
            [
                "# Auto-generated by Galaxy Pocket and Tunnel Finder",
                "load_tunnels no",
                "load_cluster_tree no",
                "stop_after never",
                "time_sparsity 1",
                "first_frame 1",
                "last_frame 1",
                f"starting_point_coordinates {x} {y} {z}",
                f"probe_radius {args.caver_probe_radius}",
                f"shell_radius {args.caver_shell_radius}",
                f"shell_depth {args.caver_shell_depth}",
                "generate_tunnel_characteristics yes",
                "generate_tunnel_profiles yes",
                "generate_histograms no",
                "",
            ]
        ),
        encoding="utf-8",
    )


def parse_caver_outputs(out_dir: Path, x: str, y: str, z: str) -> list[TunnelSummary]:
    tunnels: list[TunnelSummary] = []
    characteristics = find_first(out_dir, ["tunnel_characteristics.csv", "*tunnel*characteristics*.csv"])
    if characteristics:
        tunnels.extend(parse_caver_characteristics(characteristics, x, y, z))
    add_caver_residue_files(out_dir, tunnels)
    return tunnels


def parse_caver_characteristics(path: Path, x: str, y: str, z: str) -> list[TunnelSummary]:
    tunnels: list[TunnelSummary] = []
    rows = read_csv_dicts(path)
    for index, row in enumerate(rows, start=1):
        local_id = first_value(row, "tunnel", "tunnel_id", "cluster", "cluster_id", "id") or str(index)
        tunnel = TunnelSummary(tunnel_id=prefixed_id("caver", local_id), start_x=x, start_y=y, start_z=z)
        tunnel.bottleneck_radius = first_value(row, "bottleneck_radius", "bottleneck", "min_radius", "radius")
        tunnel.length = first_value(row, "length", "tunnel_length")
        tunnel.curvature = first_value(row, "curvature")
        tunnel.attributes = row
        tunnels.append(tunnel)
    return tunnels


def add_caver_residue_files(out_dir: Path, tunnels: list[TunnelSummary]) -> None:
    tunnel_map = {tunnel.tunnel_id.split(":", 1)[-1]: tunnel for tunnel in tunnels}
    for path in sorted(out_dir.rglob("*.pdb")):
        if "tunnel" not in path.name.lower() and "cl_" not in path.name.lower():
            continue
        residues = [residue.residue_id for residue in parse_pocket_residues(path)]
        if not residues:
            continue
        match = re.search(r"(?:tunnel|cl_?)(\d+)", path.name, re.IGNORECASE)
        local_id = match.group(1).lstrip("0") if match else str(len(tunnels) + 1)
        local_id = local_id or "0"
        tunnel = tunnel_map.get(local_id)
        if tunnel is None:
            tunnel = TunnelSummary(tunnel_id=prefixed_id("caver", local_id))
            tunnels.append(tunnel)
            tunnel_map[local_id] = tunnel
        tunnel.residue_ids = sorted(set(tunnel.residue_ids + residues), key=residue_sort_key)
        tunnel.tunnel_file = str(path.relative_to(out_dir))


# ----------------------------- output -----------------------------


def write_pockets_tsv(pockets: list[PocketSummary], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=POCKET_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for pocket in sorted(pockets, key=pocket_sort_key):
            attrs = pocket.attributes
            known_keys = {
                "score",
                "pocket_score",
                "druggability_score",
                "probability",
                "number_of_alpha_spheres",
                "alpha_spheres",
                "total_sasa",
                "polar_sasa",
                "apolar_sasa",
                "volume",
                "mean_local_hydrophobic_density",
                "mean_alpha_sphere_radius",
                "center_x",
                "center_y",
                "center_z",
                "center",
            }
            extras = {key: value for key, value in attrs.items() if key not in known_keys}
            writer.writerow(
                {
                    "pocket_id": pocket.pocket_id,
                    "source": pocket.source,
                    "score": numeric_attr(attrs, "score", "pocket_score"),
                    "druggability_score": numeric_attr(attrs, "druggability_score", "probability"),
                    "number_of_alpha_spheres": numeric_attr(attrs, "number_of_alpha_spheres", "alpha_spheres"),
                    "total_sasa": numeric_attr(attrs, "total_sasa"),
                    "polar_sasa": numeric_attr(attrs, "polar_sasa"),
                    "apolar_sasa": numeric_attr(attrs, "apolar_sasa"),
                    "volume": numeric_attr(attrs, "volume"),
                    "mean_local_hydrophobic_density": numeric_attr(attrs, "mean_local_hydrophobic_density"),
                    "mean_alpha_sphere_radius": numeric_attr(attrs, "mean_alpha_sphere_radius"),
                    "center_x": numeric_attr(attrs, "center_x", "x"),
                    "center_y": numeric_attr(attrs, "center_y", "y"),
                    "center_z": numeric_attr(attrs, "center_z", "z"),
                    "residue_count": len(pocket.residue_ids),
                    "residue_ids": ",".join(pocket.residue_ids),
                    "pocket_file": pocket.pocket_file,
                    "extra_attributes_json": json.dumps(extras, ensure_ascii=False, sort_keys=True),
                    "status": pocket.status,
                    "warning": pocket.warning,
                }
            )


def write_tunnels_tsv(tunnels: list[TunnelSummary], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TUNNEL_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        if not tunnels:
            writer.writerow(
                {
                    "tunnel_id": "",
                    "source": "caver",
                    "linked_pocket_id": "",
                    "start_x": "",
                    "start_y": "",
                    "start_z": "",
                    "bottleneck_radius": "",
                    "length": "",
                    "curvature": "",
                    "residue_ids": "",
                    "tunnel_file": "",
                    "extra_attributes_json": "{}",
                    "status": "not_computed",
                    "warning": "CAVER was not enabled.",
                }
            )
            return
        for tunnel in tunnels:
            writer.writerow(
                {
                    "tunnel_id": tunnel.tunnel_id,
                    "source": tunnel.source,
                    "linked_pocket_id": tunnel.linked_pocket_id,
                    "start_x": tunnel.start_x,
                    "start_y": tunnel.start_y,
                    "start_z": tunnel.start_z,
                    "bottleneck_radius": tunnel.bottleneck_radius,
                    "length": tunnel.length,
                    "curvature": tunnel.curvature,
                    "residue_ids": ",".join(tunnel.residue_ids),
                    "tunnel_file": tunnel.tunnel_file,
                    "extra_attributes_json": json.dumps(tunnel.attributes, ensure_ascii=False, sort_keys=True),
                    "status": tunnel.status,
                    "warning": tunnel.warning,
                }
            )


def write_residue_annotation(pockets: list[PocketSummary], tunnels: list[TunnelSummary], path: Path) -> None:
    residue_to_pockets: dict[str, list[str]] = defaultdict(list)
    residue_to_tunnels: dict[str, list[str]] = defaultdict(list)
    residue_name_by_id: dict[str, ResidueKey] = {}
    for pocket in pockets:
        for residue_id in pocket.residue_ids:
            residue_to_pockets[residue_id].append(pocket.pocket_id)
            parsed = residue_id_to_key(residue_id)
            if parsed:
                residue_name_by_id[residue_id] = parsed
    for tunnel in tunnels:
        for residue_id in tunnel.residue_ids:
            residue_to_tunnels[residue_id].append(tunnel.tunnel_id)
            parsed = residue_id_to_key(residue_id)
            if parsed:
                residue_name_by_id[residue_id] = parsed
    residue_ids = sorted(set(residue_to_pockets) | set(residue_to_tunnels), key=residue_sort_key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ANNOTATION_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for residue_id in residue_ids:
            key = residue_name_by_id.get(residue_id)
            pocket_ids = sorted(set(residue_to_pockets.get(residue_id, [])))
            tunnel_ids = sorted(set(residue_to_tunnels.get(residue_id, [])))
            sources = sorted(
                set([pocket_id.split(":", 1)[0] for pocket_id in pocket_ids] + [tunnel_id.split(":", 1)[0] for tunnel_id in tunnel_ids])
            )
            writer.writerow(
                {
                    "chain": key.chain if key else "",
                    "residue_number": key.number if key else "",
                    "insertion_code": key.insertion_code if key else "",
                    "residue_name": key.residue_name if key else "",
                    "residue_id": residue_id,
                    "in_pocket": "yes" if pocket_ids else "no",
                    "pocket_ids": ",".join(pocket_ids),
                    "pocket_count": len(pocket_ids),
                    "in_tunnel": "yes" if tunnel_ids else "no",
                    "tunnel_ids": ",".join(tunnel_ids),
                    "tunnel_count": len(tunnel_ids),
                    "source": ",".join(sources),
                }
            )


def residue_id_to_key(residue_id: str) -> ResidueKey | None:
    match = re.match(r"([^:]+):([A-Za-z0-9]{3})(-?\d+)([A-Za-z]?)$", residue_id)
    if not match:
        return None
    return ResidueKey(match.group(1), match.group(3), match.group(4), match.group(2).upper())


def residue_sort_key(residue_id: str) -> tuple[str, int, str, str]:
    key = residue_id_to_key(residue_id)
    if key is None:
        return ("", 999999, "", residue_id)
    try:
        number = int(key.number)
    except ValueError:
        number = 999999
    return (key.chain, number, key.insertion_code, key.residue_name)


def write_combined_structure(raw_dirs: list[Path], destination: Path) -> None:
    with destination.open("w", encoding="utf-8", newline="\n") as out:
        model_index = 1
        wrote_any = False
        for raw_dir in raw_dirs:
            for pdb_file in sorted(raw_dir.rglob("*.pdb")):
                if pdb_file.name.lower() == "input.pdb":
                    continue
                out.write(f"MODEL     {model_index}\n")
                with pdb_file.open("r", encoding="utf-8", errors="replace") as handle:
                    for line in handle:
                        if line.startswith(("ATOM", "HETATM", "TER")):
                            out.write(line.rstrip("\n") + "\n")
                out.write("ENDMDL\n")
                wrote_any = True
                model_index += 1
        if not wrote_any:
            out.write("HEADER    NO POCKET OR TUNNEL STRUCTURE AVAILABLE\n")
        out.write("END\n")


def zip_directories(source_dirs: list[Path], destination_zip: Path) -> None:
    with zipfile.ZipFile(destination_zip, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for source_dir in source_dirs:
            if not source_dir or not source_dir.exists():
                continue
            for path in sorted(source_dir.rglob("*")):
                if path.is_file():
                    archive.write(path, arcname=str(path.relative_to(source_dir.parent)))


def main() -> int:
    args = parse_args()
    log_lines = [
        "Pocket and Tunnel Finder",
        f"run_fpocket={args.run_fpocket}",
        f"run_p2rank={args.run_p2rank}",
        f"run_caver={args.run_caver}",
    ]
    if not (args.run_fpocket or args.run_p2rank or args.run_caver):
        args.run_fpocket = True
        log_lines.append("No backend enabled by CLI; fpocket enabled by default.")

    pockets: list[PocketSummary] = []
    tunnels: list[TunnelSummary] = []
    raw_dirs: list[Path] = []
    try:
        with tempfile.TemporaryDirectory(prefix="pocket_tunnel_run_") as tmp_name:
            root_work_dir = Path(tmp_name)
            if args.run_fpocket:
                fpocket_pockets, fpocket_dir = run_fpocket(args, root_work_dir, log_lines)
                pockets.extend(fpocket_pockets)
                if fpocket_dir:
                    raw_dirs.append(fpocket_dir)
            else:
                Path(args.fpocket_info).write_text("fpocket was not enabled.\n", encoding="utf-8")

            if args.run_p2rank:
                p2rank_pockets, p2rank_dir = run_p2rank(args, root_work_dir, log_lines)
                pockets.extend(p2rank_pockets)
                if p2rank_dir:
                    raw_dirs.append(p2rank_dir)
            else:
                Path(args.p2rank_predictions).write_text("P2Rank was not enabled.\n", encoding="utf-8")
                Path(args.p2rank_residues).write_text("P2Rank was not enabled.\n", encoding="utf-8")

            if args.run_caver:
                caver_tunnels, caver_dir = run_caver(args, root_work_dir, log_lines)
                tunnels.extend(caver_tunnels)
                if caver_dir:
                    raw_dirs.append(caver_dir)
            else:
                Path(args.caver_config).write_text("CAVER was not enabled.\n", encoding="utf-8")

            write_pockets_tsv(pockets, Path(args.pockets_tsv))
            write_tunnels_tsv(tunnels, Path(args.tunnels_tsv))
            write_residue_annotation(pockets, tunnels, Path(args.residue_annotation_tsv))
            write_combined_structure(raw_dirs, Path(args.pocket_structure_pdb))
            zip_directories(raw_dirs, Path(args.raw_archive))
            log_lines.append(f"total_pockets={len(pockets)}")
            log_lines.append(f"total_tunnels={len(tunnels)}")
    except ToolError as exc:
        Path(args.run_log).write_text("\n".join(log_lines) + f"\nERROR: {exc}\n", encoding="utf-8")
        print(f"pocket_tunnel_finder: {exc}", file=sys.stderr)
        return 1

    Path(args.run_log).write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
