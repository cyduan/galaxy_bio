#!/usr/bin/env python

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path
from typing import Any


INPUT_PLACEHOLDER = "__INPUT_STRUCTURE__"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Foundry RFDiffusion3 from Galaxy.")
    parser.add_argument("--mode", required=True)
    parser.add_argument("--config")
    parser.add_argument("--structure")
    parser.add_argument("--run-name", default="galaxy_rfd3_design")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--prepared-input", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--archive", required=True)
    parser.add_argument("--structures-dir", required=True)
    parser.add_argument("--rfd3-command", default=os.environ.get("RFD3_BINARY", "rfd3"))
    parser.add_argument("--ckpt-path", default=os.environ.get("RFD3_CKPT_PATH"))
    parser.add_argument("--n-batches", type=int)
    parser.add_argument("--diffusion-batch-size", type=int)
    parser.add_argument("--num-timesteps", type=int)
    parser.add_argument("--step-scale", type=float)
    parser.add_argument("--gamma-0", type=float)
    parser.add_argument("--partial-t", type=int)
    parser.add_argument("--low-memory-mode", action="store_true")
    parser.add_argument("--dump-trajectories", action="store_true")
    parser.add_argument("--additional-overrides", default="")

    parser.add_argument("--length")
    parser.add_argument("--contig")
    parser.add_argument("--ligand")
    parser.add_argument("--unindex")
    parser.add_argument("--ori-token")
    parser.add_argument("--infer-ori-strategy")
    parser.add_argument("--is-non-loopy", choices=["true", "false"])
    parser.add_argument("--redesign-motif-sidechains", choices=["true", "false"])
    parser.add_argument("--select-fixed-atoms-json")
    parser.add_argument("--select-hotspots-json")
    parser.add_argument("--select-buried-json")
    parser.add_argument("--select-exposed-json")
    parser.add_argument("--select-hbond-donors-json")
    parser.add_argument("--select-hbond-acceptors-json")
    parser.add_argument("--select-unfixed-sequence")
    return parser.parse_args()


def safe_name(text: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("._")
    return sanitized or "galaxy_rfd3_design"


def load_json_object(text: str | None, label: str) -> Any:
    if text is None or not text.strip():
        return None
    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{label} must be valid JSON: {exc}") from exc


def bool_value(value: str | None) -> bool | None:
    if value is None:
        return None
    return value == "true"


def copy_structure(path: str | None, work_dir: Path) -> str | None:
    if not path:
        return None
    source = Path(path)
    if not source.exists():
        raise FileNotFoundError(f"Input structure does not exist: {source}")
    suffix = "".join(source.suffixes) or ".pdb"
    destination = work_dir / f"input_structure{suffix}"
    shutil.copyfile(source, destination)
    return str(destination)


def replace_placeholder(value: Any, replacement: str | None) -> Any:
    if replacement is None:
        return value
    if isinstance(value, str):
        return value.replace(INPUT_PLACEHOLDER, replacement)
    if isinstance(value, list):
        return [replace_placeholder(item, replacement) for item in value]
    if isinstance(value, dict):
        return {key: replace_placeholder(item, replacement) for key, item in value.items()}
    return value


def write_raw_config(config_path: Path, output_path: Path, input_structure: str | None) -> None:
    raw_text = config_path.read_text(encoding="utf-8")
    try:
        data = json.loads(raw_text)
    except json.JSONDecodeError:
        if input_structure:
            raw_text = raw_text.replace(INPUT_PLACEHOLDER, input_structure)
        output_path.write_text(raw_text, encoding="utf-8")
        return
    data = replace_placeholder(data, input_structure)
    output_path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def add_if_present(target: dict[str, Any], key: str, value: Any) -> None:
    if value is not None and value != "":
        target[key] = value


def build_spec(args: argparse.Namespace, input_structure: str | None) -> dict[str, Any]:
    spec: dict[str, Any] = {}
    if args.mode != "de_novo":
        if not input_structure:
            raise ValueError(f"Mode '{args.mode}' requires an input structure.")
        spec["input"] = input_structure

    add_if_present(spec, "length", args.length)
    add_if_present(spec, "contig", args.contig)
    add_if_present(spec, "ligand", args.ligand)
    add_if_present(spec, "unindex", args.unindex)
    add_if_present(spec, "ori_token", args.ori_token)
    add_if_present(spec, "infer_ori_strategy", args.infer_ori_strategy)
    add_if_present(spec, "partial_T", args.partial_t)

    for key, value in (
        ("select_fixed_atoms", load_json_object(args.select_fixed_atoms_json, "Fixed atom selection")),
        ("select_hotspots", load_json_object(args.select_hotspots_json, "Hotspot selection")),
        ("select_buried", load_json_object(args.select_buried_json, "Buried selection")),
        ("select_exposed", load_json_object(args.select_exposed_json, "Exposed selection")),
        ("select_hbond_donors", load_json_object(args.select_hbond_donors_json, "H-bond donor selection")),
        ("select_hbond_acceptors", load_json_object(args.select_hbond_acceptors_json, "H-bond acceptor selection")),
        ("select_unfixed_sequence", load_json_object(args.select_unfixed_sequence, "Unfixed sequence selection")),
    ):
        add_if_present(spec, key, value)

    is_non_loopy = bool_value(args.is_non_loopy)
    redesign_motif_sidechains = bool_value(args.redesign_motif_sidechains)
    add_if_present(spec, "is_non_loopy", is_non_loopy)
    add_if_present(spec, "redesign_motif_sidechains", redesign_motif_sidechains)

    return {safe_name(args.run_name): spec}


def prepare_input(args: argparse.Namespace, work_dir: Path, prepared_input: Path) -> None:
    input_structure = copy_structure(args.structure, work_dir)
    if args.mode == "raw_config":
        if not args.config:
            raise ValueError("Raw configuration mode requires a JSON/YAML input specification.")
        write_raw_config(Path(args.config), prepared_input, input_structure)
        return
    spec = build_spec(args, input_structure)
    prepared_input.write_text(json.dumps(spec, indent=2, sort_keys=True), encoding="utf-8")


def split_overrides(text: str) -> list[str]:
    if not text.strip():
        return []
    return shlex.split(text, posix=(os.name != "nt"))


def build_command(args: argparse.Namespace, prepared_input: Path, out_dir: Path) -> list[str]:
    if args.rfd3_command == "mock-rfd3":
        return ["mock-rfd3"]
    command = shlex.split(args.rfd3_command, posix=(os.name != "nt"))
    command.extend(["design", f"out_dir={out_dir}", f"inputs={prepared_input}"])
    add_scalar_overrides(command, args)
    command.extend(split_overrides(args.additional_overrides))
    return command


def add_scalar_overrides(command: list[str], args: argparse.Namespace) -> None:
    mapping = {
        "ckpt_path": args.ckpt_path,
        "n_batches": args.n_batches,
        "diffusion_batch_size": args.diffusion_batch_size,
        "num_timesteps": args.num_timesteps,
        "step_scale": args.step_scale,
        "gamma_0": args.gamma_0,
        "low_memory_mode": args.low_memory_mode if args.low_memory_mode else None,
        "dump_trajectories": args.dump_trajectories if args.dump_trajectories else None,
    }
    for key, value in mapping.items():
        if value is not None and value != "":
            command.append(f"{key}={value}")


def run_mock(out_dir: Path, prepared_input: Path) -> subprocess.CompletedProcess[str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    cif_text = """data_mock_rfd3
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA GLY A 1 0.000 0.000 0.000
#
"""
    with gzip.open(out_dir / "mock_design_0.cif.gz", "wt", encoding="utf-8") as handle:
        handle.write(cif_text)
    (out_dir / "mock_scores.json").write_text(
        json.dumps({"mock": True, "input": str(prepared_input)}, indent=2), encoding="utf-8"
    )
    return subprocess.CompletedProcess(["mock-rfd3"], 0, "mock RFD3 completed\n", "")


def collect_outputs(out_dir: Path, structures_dir: Path, archive: Path) -> list[dict[str, str]]:
    structures_dir.mkdir(parents=True, exist_ok=True)
    collected: list[dict[str, str]] = []
    for path in sorted(out_dir.rglob("*")):
        if not path.is_file():
            continue
        if path.name.endswith(".cif.gz"):
            target = structures_dir / path.name[:-3]
            with gzip.open(path, "rt", encoding="utf-8") as source, target.open("w", encoding="utf-8") as dest:
                shutil.copyfileobj(source, dest)
            collected.append({"source": str(path), "galaxy_structure": str(target), "type": "cif.gz"})
        elif path.suffix.lower() in {".cif", ".pdb"}:
            target = structures_dir / path.name
            shutil.copyfile(path, target)
            collected.append({"source": str(path), "galaxy_structure": str(target), "type": path.suffix.lower()[1:]})

    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as zip_handle:
        for path in sorted(out_dir.rglob("*")):
            if path.is_file():
                zip_handle.write(path, path.relative_to(out_dir.parent))
    return collected


def main() -> int:
    args = parse_args()
    work_dir = Path.cwd() / ".rfd3_galaxy_work"
    out_dir = Path(args.out_dir)
    prepared_input = Path(args.prepared_input)
    summary_json = Path(args.summary_json)
    archive = Path(args.archive)
    structures_dir = Path(args.structures_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    prepare_input(args, work_dir, prepared_input)
    command = build_command(args, prepared_input, out_dir)

    if command == ["mock-rfd3"]:
        completed = run_mock(out_dir, prepared_input)
    else:
        completed = subprocess.run(command, capture_output=True, text=True)

    sys.stdout.write(completed.stdout or "")
    sys.stderr.write(completed.stderr or "")
    if completed.returncode != 0:
        raise RuntimeError(f"RFDiffusion3 exited with status {completed.returncode}: {' '.join(command)}")

    structures = collect_outputs(out_dir, structures_dir, archive)
    summary = {
        "tool": "RFDiffusion3",
        "mode": args.mode,
        "command": command,
        "prepared_input": str(prepared_input),
        "output_archive": str(archive),
        "structures": structures,
    }
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
