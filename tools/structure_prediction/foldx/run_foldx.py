#!/usr/bin/env python

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
import zipfile
from pathlib import Path


STABILITY_COLUMNS = [
    "Pdb",
    "Total Energy",
    "Backbone Hbond",
    "Sidechain Hbond",
    "Van der Waals",
    "Electrostatics",
    "Solvation Polar",
    "Solvation Hydrophobic",
    "Van der Waals clashes",
    "Entropy Sidechain",
    "Entropy Mainchain",
    "Sloop Entropy",
    "Mloop Entropy",
    "Cis Bond",
    "Torsional Clash",
    "Backbone Clash",
    "Helix Dipole",
    "Water Bridge",
    "Disulfide",
    "Electrostatic Kon",
    "Partial Covalent Bonds",
    "Energy Ionisation",
    "Entropy Complex",
    "Number of Residues",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run FoldX stability and mutation-energy analyses.")
    parser.add_argument("--mode", choices=["stability", "buildmodel"], required=True)
    parser.add_argument("--pdb", required=True)
    parser.add_argument("--foldx-command", default=os.environ.get("FOLDX_BINARY", "foldx"))
    parser.add_argument("--repair", action="store_true")
    parser.add_argument("--pH", default="7.0")
    parser.add_argument("--temperature", default="298")
    parser.add_argument("--ionic-strength", default="0.05")
    parser.add_argument("--mutations")
    parser.add_argument("--mutation-file")
    parser.add_argument("--number-of-runs", type=int, default=5)
    parser.add_argument("--stability-output")
    parser.add_argument("--ddg-output")
    parser.add_argument("--repaired-pdb")
    parser.add_argument("--mutant-structures-dir")
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--archive", required=True)
    parser.add_argument("--run-log", required=True)
    return parser.parse_args()


def safe_name(text: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("._")
    return sanitized or "foldx_input"


def resolve_foldx(command: str) -> list[str]:
    if command == "mock-foldx":
        return ["mock-foldx"]
    parts = shlex.split(command, posix=(os.name != "nt"))
    if not parts:
        raise ValueError("FoldX command is empty.")
    return parts


def copy_input_pdb(source: Path, work_dir: Path) -> Path:
    # Galaxy stores datasets as .dat files on disk even when the datatype is PDB.
    # FoldX expects the copied working file to look like a PDB file.
    target = work_dir / f"{safe_name(source.stem)}.pdb"
    shutil.copyfile(source, target)
    return target


def run_command(command: list[str], cwd: Path, log_lines: list[str]) -> subprocess.CompletedProcess[str]:
    log_lines.append("$ " + " ".join(command))
    if command[0] == "mock-foldx":
        completed = run_mock(command, cwd)
        log_lines.append(f"[returncode] {completed.returncode}")
        if completed.stdout:
            log_lines.append("[stdout]\n" + completed.stdout.rstrip())
        if completed.stderr:
            log_lines.append("[stderr]\n" + completed.stderr.rstrip())
        return completed
    completed = subprocess.run(command, cwd=cwd, capture_output=True, text=True)
    log_lines.append(f"[returncode] {completed.returncode}")
    if completed.stdout:
        log_lines.append("[stdout]\n" + completed.stdout.rstrip())
    if completed.stderr:
        log_lines.append("[stderr]\n" + completed.stderr.rstrip())
    return completed


def foldx_failure_message(step: str, completed: subprocess.CompletedProcess[str]) -> str:
    message = [f"FoldX {step} failed with status {completed.returncode}."]
    if completed.stdout:
        message.append("FoldX stdout:\n" + completed.stdout.strip())
    if completed.stderr:
        message.append("FoldX stderr:\n" + completed.stderr.strip())
    return "\n\n".join(message)


def run_mock(command: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    command_text = " ".join(command)
    pdb_name = next((item.split("=", 1)[1] for item in command if item.startswith("--pdb=")), "input.pdb")
    pdb_path = cwd / pdb_name
    pdb_stem = Path(pdb_name).stem
    if "--command=RepairPDB" in command_text:
        repaired = cwd / f"{pdb_stem}_Repair.pdb"
        shutil.copyfile(pdb_path, repaired)
        return subprocess.CompletedProcess(command, 0, f"RepairPDB wrote {repaired.name}\n", "")
    if "--command=Stability" in command_text:
        (cwd / "stability_ST.fxout").write_text(
            "\n".join(
                [
                    "# FoldX mock Stability output",
                    "Pdb\tTotal Energy",
                    f"{pdb_name}\t-12.345",
                    "",
                ]
            ),
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(command, 0, "Stability completed\n", "")
    if "--command=BuildModel" in command_text:
        (cwd / "Dif_buildmodel.fxout").write_text(
            "\n".join(
                [
                    "# FoldX mock BuildModel output",
                    "Pdb\tMutations\tddG",
                    f"{pdb_name}\tGA1V\t1.234",
                    "",
                ]
            ),
            encoding="utf-8",
        )
        (cwd / "Average_buildmodel.fxout").write_text("Pdb\tAverage\nmock\t1.234\n", encoding="utf-8")
        (cwd / f"{pdb_stem}_1.pdb").write_text(
            "\n".join(
                [
                    "HEADER    MOCK FOLDX MUTANT",
                    "ATOM      1  N   VAL A   1       0.000   0.000   0.000  1.00 20.00           N",
                    "TER",
                    "END",
                    "",
                ]
            ),
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(command, 0, "BuildModel completed\n", "")
    return subprocess.CompletedProcess(command, 1, "", "Unknown mock FoldX command\n")


def find_repaired_pdb(work_dir: Path, input_pdb: Path) -> Path | None:
    expected = work_dir / f"{input_pdb.stem}_Repair.pdb"
    if expected.exists():
        return expected
    repaired = sorted(work_dir.glob("*_Repair.pdb"))
    return repaired[0] if repaired else None


def first_matching_file(work_dir: Path, patterns: list[str]) -> Path | None:
    matches: list[Path] = []
    for pattern in patterns:
        matches.extend(sorted(work_dir.glob(pattern)))
    return matches[0] if matches else None


def non_comment_rows(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            rows.append(re.split(r"\s+", line))
    return rows


def convert_stability_to_csv(source: Path, destination: str) -> None:
    rows = non_comment_rows(source)
    if not rows:
        raise RuntimeError(f"FoldX stability output is empty: {source}")
    data_rows = [row for row in rows if is_stability_data_row(row)]
    if not data_rows:
        data_rows = rows
    with Path(destination).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(STABILITY_COLUMNS)
        for row in data_rows:
            if len(row) == len(STABILITY_COLUMNS):
                writer.writerow(row)
            elif len(row) == len(STABILITY_COLUMNS) - 1:
                writer.writerow([source.name, *row])
            else:
                padded = row[: len(STABILITY_COLUMNS)]
                padded.extend([""] * (len(STABILITY_COLUMNS) - len(padded)))
                writer.writerow(padded)


def is_stability_data_row(row: list[str]) -> bool:
    return len(row) >= 2 and all(looks_numeric(value) for value in row[1:])


def convert_foldx_table_to_csv(source: Path, destination: str) -> None:
    rows = non_comment_rows(source)
    with Path(destination).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        if not rows:
            writer.writerow(["FoldX output"])
            return
        first_row = rows[0]
        has_header = any(not looks_numeric(value) for value in first_row[1:])
        if has_header:
            writer.writerow(first_row)
            writer.writerows(rows[1:])
        else:
            max_columns = max(len(row) for row in rows)
            writer.writerow([f"column_{index}" for index in range(1, max_columns + 1)])
            for row in rows:
                padded = row[:max_columns]
                padded.extend([""] * (max_columns - len(padded)))
                writer.writerow(padded)


def looks_numeric(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


def collect_mutant_structures(work_dir: Path, structures_dir: Path, excluded_stems: set[str]) -> list[str]:
    structures_dir.mkdir(parents=True, exist_ok=True)
    copied: list[str] = []
    for path in sorted(work_dir.glob("*.pdb")):
        if path.name.endswith("_Repair.pdb"):
            continue
        if path.stem in excluded_stems:
            continue
        target = structures_dir / path.name
        shutil.copyfile(path, target)
        copied.append(str(target))
    return copied


def read_mutation_lines(args: argparse.Namespace) -> list[str]:
    lines: list[str] = []
    if args.mutation_file:
        lines.extend(Path(args.mutation_file).read_text(encoding="utf-8").splitlines())
    if args.mutations:
        lines.extend(args.mutations.splitlines())
    cleaned = [line.strip() for line in lines if line.strip() and not line.strip().startswith("#")]
    if args.mode == "buildmodel" and not cleaned:
        raise ValueError("FoldX BuildModel requires at least one mutation line.")
    return [line if line.endswith(";") else f"{line};" for line in cleaned]


def write_mutation_file(args: argparse.Namespace, work_dir: Path) -> Path | None:
    mutation_lines = read_mutation_lines(args)
    if not mutation_lines:
        return None
    mutation_file = work_dir / "individual_list.txt"
    mutation_file.write_text("\n".join(mutation_lines) + "\n", encoding="utf-8")
    return mutation_file


def build_base_options(args: argparse.Namespace) -> list[str]:
    return [
        f"--pH={args.pH}",
        f"--temperature={args.temperature}",
        f"--ionStrength={args.ionic_strength}",
    ]


def run_repair(args: argparse.Namespace, foldx: list[str], input_pdb: Path, work_dir: Path, log_lines: list[str]) -> Path:
    command = [
        *foldx,
        "--command=RepairPDB",
        f"--pdb={input_pdb.name}",
        *build_base_options(args),
    ]
    completed = run_command(command, work_dir, log_lines)
    if completed.returncode != 0:
        raise RuntimeError(foldx_failure_message("RepairPDB", completed))
    repaired = find_repaired_pdb(work_dir, input_pdb)
    if repaired is None:
        raise RuntimeError("FoldX RepairPDB completed but no *_Repair.pdb file was found.")
    return repaired


def run_stability(args: argparse.Namespace, foldx: list[str], pdb: Path, work_dir: Path, log_lines: list[str]) -> Path:
    command = [
        *foldx,
        "--command=Stability",
        f"--pdb={pdb.name}",
        "--output-file=stability",
        *build_base_options(args),
    ]
    completed = run_command(command, work_dir, log_lines)
    if completed.returncode != 0:
        raise RuntimeError(foldx_failure_message("Stability", completed))
    result = first_matching_file(work_dir, ["stability_ST.fxout", "*_ST.fxout", "*.fxout"])
    if result is None:
        raise RuntimeError("FoldX Stability completed but no stability .fxout file was found.")
    return result


def run_buildmodel(
    args: argparse.Namespace, foldx: list[str], pdb: Path, mutation_file: Path, work_dir: Path, log_lines: list[str]
) -> Path:
    command = [
        *foldx,
        "--command=BuildModel",
        f"--pdb={pdb.name}",
        f"--mutant-file={mutation_file.name}",
        f"--numberOfRuns={args.number_of_runs}",
        "--output-file=buildmodel",
        *build_base_options(args),
    ]
    completed = run_command(command, work_dir, log_lines)
    if completed.returncode != 0:
        raise RuntimeError(foldx_failure_message("BuildModel", completed))
    result = first_matching_file(work_dir, ["Dif_*.fxout", "*Dif*.fxout", "Average_*.fxout", "*.fxout"])
    if result is None:
        raise RuntimeError("FoldX BuildModel completed but no .fxout file was found.")
    return result


def make_archive(work_dir: Path, archive: Path) -> None:
    archive.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as zip_handle:
        for path in sorted(work_dir.rglob("*")):
            if path.is_file():
                zip_handle.write(path, path.relative_to(work_dir))


def main() -> int:
    args = parse_args()
    input_pdb = Path(args.pdb)
    work_dir = Path.cwd() / "foldx_work"
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True)
    log_lines = ["FoldX Galaxy run log", f"Mode: {args.mode}", f"Input: {input_pdb}", ""]

    foldx = resolve_foldx(args.foldx_command)
    copied_pdb = copy_input_pdb(input_pdb, work_dir)
    active_pdb = copied_pdb
    repaired_pdb_path: str | None = None

    try:
        mutation_file = write_mutation_file(args, work_dir)
        if args.repair:
            active_pdb = run_repair(args, foldx, copied_pdb, work_dir, log_lines)
            repaired_pdb_path = str(active_pdb)
            if args.repaired_pdb:
                shutil.copyfile(active_pdb, args.repaired_pdb)

        if args.mode == "stability":
            result = run_stability(args, foldx, active_pdb, work_dir, log_lines)
            if not args.stability_output:
                raise ValueError("--stability-output is required in stability mode.")
            convert_stability_to_csv(result, args.stability_output)
            primary_output = args.stability_output
            mutant_structures: list[str] = []
        else:
            if mutation_file is None:
                raise ValueError("--mutation-file or --mutations is required in buildmodel mode.")
            result = run_buildmodel(args, foldx, active_pdb, mutation_file, work_dir, log_lines)
            if not args.ddg_output:
                raise ValueError("--ddg-output is required in buildmodel mode.")
            convert_foldx_table_to_csv(result, args.ddg_output)
            primary_output = args.ddg_output
            mutant_structures = (
                collect_mutant_structures(work_dir, Path(args.mutant_structures_dir), {copied_pdb.stem, active_pdb.stem})
                if args.mutant_structures_dir
                else []
            )

        log_lines.append("")
        log_lines.append("Status: success")
        log_lines.append(f"Primary output: {primary_output}")
        Path(args.run_log).write_text("\n".join(log_lines) + "\n", encoding="utf-8")
        make_archive(work_dir, Path(args.archive))
        summary = {
            "tool": "FoldX",
            "mode": args.mode,
            "status": "success",
            "foldx_command": args.foldx_command,
            "input_pdb": str(input_pdb),
            "active_pdb": str(active_pdb),
            "repaired_pdb": repaired_pdb_path,
            "primary_output": primary_output,
            "mutant_structures": mutant_structures,
            "archive": args.archive,
            "pH": args.pH,
            "temperature": args.temperature,
            "ionic_strength": args.ionic_strength,
            "number_of_runs": args.number_of_runs if args.mode == "buildmodel" else None,
            "mutations": read_mutation_lines(args) if args.mode == "buildmodel" else [],
        }
        Path(args.summary_json).write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    except Exception as exc:
        log_lines.append("")
        log_lines.append("Status: failed")
        log_lines.append(str(exc))
        Path(args.run_log).write_text("\n".join(log_lines) + "\n", encoding="utf-8")
        sys.stderr.write("\n".join(log_lines) + "\n")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
