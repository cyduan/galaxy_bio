#!/usr/bin/env python

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
import uuid
from os.path import expanduser
from pathlib import Path
from typing import NoReturn


TABULAR_OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"


def stop_err(message: str) -> NoReturn:
    raise SystemExit(f"ERROR: {message}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run a local BLAST search against a FASTA subject dataset.")
    parser.add_argument("--program", required=True, choices=["blastn", "blastp"])
    parser.add_argument("--query", required=True)
    parser.add_argument("--subject", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--output-format", required=True, choices=["tabular", "pairwise"])
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--max-target-seqs", type=int, default=10)
    parser.add_argument("--num-threads", type=int, default=1)
    parser.add_argument("--blast-task", default="")
    parser.add_argument("--blast-binary", default="")
    parser.add_argument("--makeblastdb-binary", default="")
    return parser.parse_args()


def split_command(command: str) -> list[str]:
    return shlex.split(command, posix=(os.name != "nt"))


def iter_binary_candidates(binary_name: str) -> list[Path]:
    candidates: list[Path] = []
    seen: set[str] = set()

    def add_candidate(path: Path) -> None:
        resolved = str(path)
        if resolved not in seen:
            seen.add(resolved)
            candidates.append(path)

    current_python = Path(sys.executable).resolve()
    add_candidate(current_python.with_name(binary_name))

    env_name = os.environ.get("BLAST_ENV_NAME", "blastenv")
    current_parts = list(current_python.parts)
    if "envs" in current_parts:
        envs_index = current_parts.index("envs")
        envs_root = Path(*current_parts[: envs_index + 1])
        add_candidate(envs_root / env_name / "bin" / binary_name)

    conda_exe = os.environ.get("CONDA_EXE")
    if conda_exe:
        conda_path = Path(conda_exe).resolve()
        conda_root = conda_path.parent.parent
        add_candidate(conda_root / "envs" / env_name / "bin" / binary_name)

    add_candidate(Path(expanduser(f"~/miniconda3/envs/{env_name}/bin/{binary_name}")))
    add_candidate(Path(expanduser(f"~/anaconda3/envs/{env_name}/bin/{binary_name}")))
    add_candidate(Path(expanduser(f"~/.conda/envs/{env_name}/bin/{binary_name}")))
    return candidates


def resolve_binary(command: str, env_var: str, default_binary: str) -> list[str]:
    if command == "mock-makeblastdb":
        mock_cli = Path(__file__).resolve().parent / "test-data" / "mock_makeblastdb_cli.py"
        return [sys.executable, str(mock_cli)]
    if command == "mock-blastn":
        mock_cli = Path(__file__).resolve().parent / "test-data" / "mock_blast_cli.py"
        return [sys.executable, str(mock_cli), "--program", "blastn"]
    if command == "mock-blastp":
        mock_cli = Path(__file__).resolve().parent / "test-data" / "mock_blast_cli.py"
        return [sys.executable, str(mock_cli), "--program", "blastp"]

    configured = command.strip() if command else os.environ.get(env_var, default_binary)
    if configured == default_binary:
        for candidate in iter_binary_candidates(default_binary):
            if candidate.exists():
                return [str(candidate)]
    return split_command(configured)


def run_checked(command: list[str], missing_hint: str) -> subprocess.CompletedProcess[str]:
    try:
        completed = subprocess.run(command, capture_output=True, text=True)
    except FileNotFoundError as exc:
        stop_err(missing_hint)
    if completed.returncode != 0:
        sys.stderr.write(completed.stderr)
        sys.stdout.write(completed.stdout)
        stop_err(f"Command failed with exit code {completed.returncode}: {' '.join(command)}")
    return completed


def subject_dbtype(program: str) -> str:
    return "nucl" if program == "blastn" else "prot"


def main() -> int:
    args = parse_args()
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    blast_binary = resolve_binary(args.blast_binary, f"{args.program.upper()}_BINARY", args.program)
    makeblastdb_binary = resolve_binary(args.makeblastdb_binary, "MAKEBLASTDB_BINARY", "makeblastdb")

    query_path = Path(args.query)
    subject_path = Path(args.subject)

    temp_dir = output_path.parent / f".blast_tmp_{uuid.uuid4().hex}"
    temp_dir.mkdir(parents=True, exist_ok=False)
    try:
        db_prefix = temp_dir / "subject_db"
        makeblastdb_cmd = [
            *makeblastdb_binary,
            "-in",
            str(subject_path),
            "-dbtype",
            subject_dbtype(args.program),
            "-out",
            str(db_prefix),
        ]
        run_checked(
            makeblastdb_cmd,
            "Could not find 'makeblastdb'. Install NCBI BLAST+ in the BLAST environment or set MAKEBLASTDB_BINARY.",
        )

        blast_cmd = [
            *blast_binary,
            "-query",
            str(query_path),
            "-db",
            str(db_prefix),
            "-out",
            str(output_path),
            "-evalue",
            str(args.evalue),
            "-max_target_seqs",
            str(args.max_target_seqs),
            "-num_threads",
            str(args.num_threads),
        ]
        if args.program == "blastn" and args.blast_task.strip():
            blast_cmd.extend(["-task", args.blast_task.strip()])
        if args.output_format == "tabular":
            blast_cmd.extend(["-outfmt", TABULAR_OUTFMT])
        run_checked(
            blast_cmd,
            f"Could not find '{args.program}'. Install NCBI BLAST+ in the BLAST environment or set {args.program.upper()}_BINARY.",
        )
        if not output_path.exists():
            stop_err("BLAST completed without creating an output file.")
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
