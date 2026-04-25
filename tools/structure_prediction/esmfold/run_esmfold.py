#!/usr/bin/env python

from __future__ import annotations

import argparse
import importlib.util
import json
import os
from os.path import expanduser
import re
import shlex
import shutil
import subprocess
import sys
import uuid
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the esm-fold CLI on a single FASTA record.")
    parser.add_argument("--input-fasta", required=True)
    parser.add_argument("--output-pdb", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--esmfold-command", default=os.environ.get("ESMFOLD_BINARY", "esm-fold"))
    parser.add_argument("--num-recycles", type=int, default=4)
    parser.add_argument("--max-tokens-per-batch", type=int)
    parser.add_argument("--chunk-size", type=int)
    parser.add_argument("--output-name")
    parser.add_argument("--cpu-only", action="store_true")
    parser.add_argument("--cpu-offload", action="store_true")
    args = parser.parse_args()
    if args.cpu_only and args.cpu_offload:
        parser.error("--cpu-only and --cpu-offload are mutually exclusive")
    return args


def read_single_fasta(path: Path) -> tuple[str, str]:
    identifier = None
    sequence_lines: list[str] = []
    record_count = 0
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                record_count += 1
                if record_count > 1:
                    raise ValueError(
                        "ESMFold wrapper currently accepts exactly one FASTA record. "
                        "Use a collection and map over it for batches."
                    )
                identifier = line[1:].strip() or "query"
            else:
                if identifier is None:
                    raise ValueError("Input is not a valid FASTA file: sequence data found before a header line.")
                sequence_lines.append(line)
    if identifier is None:
        raise ValueError("No FASTA records were found in the input dataset.")
    sequence = "".join(sequence_lines).strip()
    if not sequence:
        raise ValueError("The FASTA record is empty.")
    return identifier, sequence


def sanitize_identifier(identifier: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9_.-]+", "_", identifier).strip("._")
    return sanitized or "query"


def split_command(command: str) -> list[str]:
    return shlex.split(command, posix=(os.name != "nt"))


def iter_esmfold_binary_candidates() -> list[Path]:
    candidates: list[Path] = []
    seen: set[str] = set()

    def add_candidate(path: Path) -> None:
        resolved = str(path)
        if resolved not in seen:
            seen.add(resolved)
            candidates.append(path)

    current_python = Path(sys.executable).resolve()
    add_candidate(current_python.with_name("esm-fold"))

    env_name = os.environ.get("ESMFOLD_ENV_NAME", "esmfold_official")

    current_parts = list(current_python.parts)
    if "envs" in current_parts:
        envs_index = current_parts.index("envs")
        envs_root = Path(*current_parts[: envs_index + 1])
        add_candidate(envs_root / env_name / "bin" / "esm-fold")

    conda_exe = os.environ.get("CONDA_EXE")
    if conda_exe:
        conda_path = Path(conda_exe).resolve()
        conda_root = conda_path.parent.parent
        add_candidate(conda_root / "envs" / env_name / "bin" / "esm-fold")

    add_candidate(Path(expanduser(f"~/miniconda3/envs/{env_name}/bin/esm-fold")))
    add_candidate(Path(expanduser(f"~/anaconda3/envs/{env_name}/bin/esm-fold")))
    add_candidate(Path(expanduser(f"~/.conda/envs/{env_name}/bin/esm-fold")))
    add_candidate(Path(f"/data/conda_envs/{env_name}/bin/esm-fold"))
    return candidates


def resolve_command(command: str) -> list[str]:
    if command == "mock-esmfold":
        mock_cli = Path(__file__).resolve().parent / "test-data" / "mock_esmfold_cli.py"
        return [sys.executable, str(mock_cli)]
    if command == "esm-fold":
        for candidate in iter_esmfold_binary_candidates():
            if candidate.exists():
                return [str(candidate)]
        try:
            has_fold_module = importlib.util.find_spec("esm.scripts.fold") is not None
        except ModuleNotFoundError:
            has_fold_module = False
        if has_fold_module:
            return [sys.executable, "-m", "esm.scripts.fold"]
        esmfold_python = os.environ.get("ESMFOLD_PYTHON")
        if esmfold_python:
            helper_cli = Path(__file__).resolve().parent / "esmfold_api_cli.py"
            return [esmfold_python, str(helper_cli)]
    return split_command(command)


def find_single_pdb(directory: Path) -> Path:
    pdb_files = sorted(directory.rglob("*.pdb"))
    if not pdb_files:
        raise RuntimeError("ESMFold did not produce any PDB files.")
    if len(pdb_files) > 1:
        names = ", ".join(path.name for path in pdb_files)
        raise RuntimeError(f"ESMFold produced multiple PDB files unexpectedly: {names}")
    return pdb_files[0]


def parse_pdb_metrics(path: Path) -> dict:
    atom_count = 0
    residue_ids = set()
    chain_ids = set()
    b_factors: list[float] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")):
                atom_count += 1
                chain_id = line[21].strip() or "_"
                residue_number = line[22:26].strip()
                insertion_code = line[26].strip()
                chain_ids.add(chain_id)
                residue_ids.add((chain_id, residue_number, insertion_code))
                b_text = line[60:66].strip()
                if b_text:
                    try:
                        b_factors.append(float(b_text))
                    except ValueError:
                        pass
    mean_b_factor = round(sum(b_factors) / len(b_factors), 3) if b_factors else None
    return {
        "atom_count": atom_count,
        "residue_count": len(residue_ids),
        "chain_ids": sorted(chain_ids),
        "mean_b_factor": mean_b_factor,
    }


def dependency_error_hint(stderr: str) -> str | None:
    env_name = os.environ.get("ESMFOLD_ENV_NAME", "esmfold_official")
    esmfold_python = os.environ.get("ESMFOLD_PYTHON")
    if not esmfold_python:
        esmfold_python = str(Path(f"/data/conda_envs/{env_name}/bin/python"))

    if "No module named 'modelcif'" in stderr:
        return (
            "The ESMFold runtime is missing the Python package 'modelcif'. "
            f"Install it into the ESMFold environment, for example: {esmfold_python} -m pip install modelcif"
        )
    if "No module named 'torch._six'" in stderr:
        return (
            "The ESMFold runtime has an incompatible DeepSpeed/PyTorch combination. "
            "For the official ESMFold CLI, prefer rebuilding the environment to match Meta's archived pins "
            "(environment.yml / README), including the pinned OpenFold commit "
            "4b41059694619831a7db195b7e0988fc4ff3a307."
        )
    if "linear_kv_points.linear" in stderr and "linear_q_points.linear" in stderr and "are missing" in stderr:
        return (
            "The ESMFold weights do not match the installed OpenFold implementation. "
            "This usually happens when a newer OpenFold tree (for example an editable OpenFold 2.x checkout) "
            "is installed instead of Meta's pinned ESMFold dependency. Rebuild the ESMFold environment with "
            f"the official pinned OpenFold commit: {esmfold_python} -m pip install "
            "'openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307' "
            "and remove any conflicting editable openfold installation first."
        )
    return None


def main() -> int:
    args = parse_args()
    input_path = Path(args.input_fasta)
    output_pdb = Path(args.output_pdb)
    summary_json = Path(args.summary_json)
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    summary_json.parent.mkdir(parents=True, exist_ok=True)

    identifier, sequence = read_single_fasta(input_path)
    output_name = sanitize_identifier(args.output_name or identifier)

    temp_dir = output_pdb.parent / f".esmfold_tmp_{uuid.uuid4().hex}"
    temp_dir.mkdir(parents=True, exist_ok=False)
    try:
        temp_fasta = temp_dir / f"{output_name}.fasta"
        output_dir = temp_dir / "esmfold_output"
        output_dir.mkdir()
        temp_fasta.write_text(f">{output_name}\n{sequence}\n", encoding="utf-8")

        command = resolve_command(args.esmfold_command)
        command.extend(["-i", str(temp_fasta), "-o", str(output_dir), "--num-recycles", str(args.num_recycles)])
        if args.max_tokens_per_batch is not None:
            command.extend(["--max-tokens-per-batch", str(args.max_tokens_per_batch)])
        if args.chunk_size is not None:
            command.extend(["--chunk-size", str(args.chunk_size)])
        if args.cpu_only:
            command.append("--cpu-only")
        elif args.cpu_offload:
            command.append("--cpu-offload")

        try:
            completed = subprocess.run(command, capture_output=True, text=True)
        except FileNotFoundError as exc:
            attempted_command = command[0] if command else args.esmfold_command
            candidate_hint = ""
            if args.esmfold_command == "esm-fold":
                candidates = [str(path) for path in iter_esmfold_binary_candidates()]
                if candidates:
                    candidate_hint = " Attempted auto-discovery in: " + ", ".join(candidates)
            raise RuntimeError(
                "Could not find the ESMFold executable "
                f"'{attempted_command}'. Install the official 'esm-fold' CLI in Galaxy's job environment, "
                "or set ESMFOLD_BINARY / --esmfold-command to the absolute executable path."
                f"{candidate_hint}"
            ) from exc
        if completed.returncode != 0:
            sys.stderr.write(completed.stderr)
            sys.stdout.write(completed.stdout)
            hint = dependency_error_hint(completed.stderr)
            if hint:
                raise RuntimeError(hint)
            raise RuntimeError(f"esm-fold exited with status {completed.returncode}")

        pdb_path = find_single_pdb(output_dir)
        shutil.copyfile(pdb_path, output_pdb)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

    metrics = parse_pdb_metrics(output_pdb)
    summary = {
        "tool": "ESMFold",
        "input_identifier": identifier,
        "output_name": output_name,
        "sequence_length": len(sequence.replace(":", "")),
        "multimer_chain_count": sequence.count(":") + 1,
        "num_recycles": args.num_recycles,
        "execution_mode": "cpu_only" if args.cpu_only else "cpu_offload" if args.cpu_offload else "gpu",
        "esmfold_command": args.esmfold_command,
        **metrics,
    }
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
