#!/usr/bin/env python

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path


STRUCTURE_PLACEHOLDER = "__STRUCTURE__"
OUTPUT_PLACEHOLDER = "__OUTPUT_DIR__"
CONFIG_PLACEHOLDER = "__CONFIG__"
WEIGHTS_PLACEHOLDER = "__WEIGHTS__"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the Foundry ProteinMPNN CLI from Galaxy.")
    parser.add_argument("--structure", required=True)
    parser.add_argument("--config")
    parser.add_argument("--weights")
    parser.add_argument("--arguments", required=True)
    parser.add_argument("--mpnn-command", default=os.environ.get("MPNN_BINARY", "mpnn"))
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--archive", required=True)
    parser.add_argument("--sequences-fasta", required=True)
    return parser.parse_args()


def copy_input(path: str, destination: Path) -> Path:
    source = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    suffix = "".join(source.suffixes) or ".dat"
    target = destination.with_suffix(suffix)
    shutil.copyfile(source, target)
    return target


def render_arguments(template: str, structure: Path, out_dir: Path, config: Path | None, weights: str | None) -> list[str]:
    replacements = {
        STRUCTURE_PLACEHOLDER: str(structure),
        OUTPUT_PLACEHOLDER: str(out_dir),
        CONFIG_PLACEHOLDER: str(config or ""),
        WEIGHTS_PLACEHOLDER: weights or "",
    }
    rendered = template
    for placeholder, value in replacements.items():
        rendered = rendered.replace(placeholder, value)
    return shlex.split(rendered, posix=(os.name != "nt"))


def run_mock(out_dir: Path) -> subprocess.CompletedProcess[str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "mock_sequences.fasta").write_text(">mock_design_1\nMKTAYIAK\n", encoding="utf-8")
    (out_dir / "mock_scores.json").write_text(json.dumps({"mock": True}, indent=2), encoding="utf-8")
    return subprocess.CompletedProcess(["mock-mpnn"], 0, "mock MPNN completed\n", "")


def collect_fasta(out_dir: Path, sequences_fasta: Path) -> list[str]:
    fasta_files = sorted(path for path in out_dir.rglob("*") if path.suffix.lower() in {".fa", ".fasta"})
    with sequences_fasta.open("w", encoding="utf-8") as output:
        for path in fasta_files:
            text = path.read_text(encoding="utf-8").strip()
            if text:
                output.write(text)
                output.write("\n")
    return [str(path) for path in fasta_files]


def make_archive(out_dir: Path, archive: Path) -> None:
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as zip_handle:
        for path in sorted(out_dir.rglob("*")):
            if path.is_file():
                zip_handle.write(path, path.relative_to(out_dir.parent))


def main() -> int:
    args = parse_args()
    work_dir = Path.cwd() / ".mpnn_galaxy_work"
    out_dir = Path(args.out_dir)
    summary_json = Path(args.summary_json)
    archive = Path(args.archive)
    sequences_fasta = Path(args.sequences_fasta)
    work_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    structure = copy_input(args.structure, work_dir / "input_structure")
    config = copy_input(args.config, work_dir / "input_config") if args.config else None

    if args.mpnn_command == "mock-mpnn":
        completed = run_mock(out_dir)
        command = ["mock-mpnn"]
    else:
        command = shlex.split(args.mpnn_command, posix=(os.name != "nt"))
        command.extend(render_arguments(args.arguments, structure, out_dir, config, args.weights))
        completed = subprocess.run(command, capture_output=True, text=True)

    sys.stdout.write(completed.stdout or "")
    sys.stderr.write(completed.stderr or "")
    if completed.returncode != 0:
        raise RuntimeError(f"Foundry ProteinMPNN exited with status {completed.returncode}: {' '.join(command)}")

    fasta_files = collect_fasta(out_dir, sequences_fasta)
    make_archive(out_dir, archive)
    summary = {
        "tool": "Foundry ProteinMPNN",
        "command": command,
        "structure": str(structure),
        "config": str(config) if config else None,
        "weights": args.weights,
        "fasta_files": fasta_files,
        "output_archive": str(archive),
    }
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
