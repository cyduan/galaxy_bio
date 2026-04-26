#!/usr/bin/env python

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import uuid
import zipfile
from pathlib import Path

from run_esmfold import dependency_error_hint, parse_pdb_metrics, resolve_command, sanitize_identifier


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the esm-fold CLI on a multi-record FASTA file.")
    parser.add_argument("--input-fasta", required=True)
    parser.add_argument("--structures-dir", required=True)
    parser.add_argument("--summaries-dir", required=True)
    parser.add_argument("--combined-summary-json", required=True)
    parser.add_argument("--run-log", required=True)
    parser.add_argument("--output-archive")
    parser.add_argument("--esmfold-command", default=os.environ.get("ESMFOLD_BINARY", "esm-fold"))
    parser.add_argument("--num-recycles", type=int, default=4)
    parser.add_argument("--max-tokens-per-batch", type=int)
    parser.add_argument("--chunk-size", type=int)
    parser.add_argument("--cpu-only", action="store_true")
    parser.add_argument("--cpu-offload", action="store_true")
    args = parser.parse_args()
    if args.cpu_only and args.cpu_offload:
        parser.error("--cpu-only and --cpu-offload are mutually exclusive")
    return args


def read_fasta_records(path: Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    current_identifier: str | None = None
    sequence_lines: list[str] = []

    def flush_record() -> None:
        if current_identifier is None:
            return
        sequence = "".join(sequence_lines).strip()
        if not sequence:
            raise ValueError(f"FASTA record '{current_identifier}' is empty.")
        records.append({"identifier": current_identifier, "sequence": sequence})

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush_record()
                current_identifier = line[1:].strip() or f"record_{len(records) + 1}"
                sequence_lines = []
            else:
                if current_identifier is None:
                    raise ValueError("Input is not a valid FASTA file: sequence data found before a header line.")
                sequence_lines.append(line)
    flush_record()

    if not records:
        raise ValueError("No FASTA records were found in the input dataset.")
    return records


def assign_output_names(records: list[dict[str, str]]) -> None:
    seen: dict[str, int] = {}
    for index, record in enumerate(records, start=1):
        base = sanitize_identifier(record["identifier"]) or f"record_{index}"
        count = seen.get(base, 0) + 1
        seen[base] = count
        record["output_name"] = base if count == 1 else f"{base}_{count}"


def write_batch_fasta(records: list[dict[str, str]], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{record['output_name']}\n{record['sequence']}\n")


def build_command(args: argparse.Namespace, input_fasta: Path, output_dir: Path) -> list[str]:
    command = resolve_command(args.esmfold_command)
    command.extend(["-i", str(input_fasta), "-o", str(output_dir), "--num-recycles", str(args.num_recycles)])
    if args.max_tokens_per_batch is not None:
        command.extend(["--max-tokens-per-batch", str(args.max_tokens_per_batch)])
    if args.chunk_size is not None:
        command.extend(["--chunk-size", str(args.chunk_size)])
    if args.cpu_only:
        command.append("--cpu-only")
    elif args.cpu_offload:
        command.append("--cpu-offload")
    return command


def summarize_record(
    record: dict[str, str],
    raw_record_dir: Path,
    structures_dir: Path,
    summaries_dir: Path,
    args: argparse.Namespace,
) -> dict:
    output_name = record["output_name"]
    pdb_files = sorted(raw_record_dir.rglob("*.pdb"))
    if not pdb_files:
        raise RuntimeError("esm-fold completed but did not produce a PDB file.")
    pdb_path = pdb_files[0]
    structure_path = structures_dir / f"{output_name}.pdb"
    shutil.copyfile(pdb_path, structure_path)
    metrics = parse_pdb_metrics(structure_path)
    summary = {
        "tool": "ESMFold Batch",
        "status": "success",
        "input_identifier": record["identifier"],
        "output_name": output_name,
        "sequence_length": len(record["sequence"].replace(":", "")),
        "multimer_chain_count": record["sequence"].count(":") + 1,
        "num_recycles": args.num_recycles,
        "execution_mode": "cpu_only" if args.cpu_only else "cpu_offload" if args.cpu_offload else "gpu",
        "structure": str(structure_path),
        "raw_output_archive_path": f"esmfold_raw_output/{output_name}",
        **metrics,
    }
    (summaries_dir / f"{output_name}.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    return summary


def write_failed_summary(
    record: dict[str, str],
    summaries_dir: Path,
    args: argparse.Namespace,
    command: list[str],
    returncode: int | None,
    stderr: str,
    stdout: str,
) -> dict:
    output_name = record["output_name"]
    summary = {
        "tool": "ESMFold Batch",
        "status": "failed",
        "input_identifier": record["identifier"],
        "output_name": output_name,
        "sequence_length": len(record["sequence"].replace(":", "")),
        "multimer_chain_count": record["sequence"].count(":") + 1,
        "num_recycles": args.num_recycles,
        "execution_mode": "cpu_only" if args.cpu_only else "cpu_offload" if args.cpu_offload else "gpu",
        "command": command,
        "returncode": returncode,
        "stdout": stdout[-4000:],
        "stderr": stderr[-4000:],
        "hint": dependency_error_hint(stderr),
    }
    (summaries_dir / f"{output_name}.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    return summary


def write_archive(archive_path: Path, paths: list[tuple[Path, str]]) -> None:
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED) as zip_handle:
        for source, prefix in paths:
            if source.is_file():
                zip_handle.write(source, Path(prefix) / source.name)
            elif source.is_dir():
                for path in sorted(source.rglob("*")):
                    if path.is_file():
                        zip_handle.write(path, Path(prefix) / path.relative_to(source))


def main() -> int:
    args = parse_args()
    input_path = Path(args.input_fasta)
    structures_dir = Path(args.structures_dir)
    summaries_dir = Path(args.summaries_dir)
    combined_summary_json = Path(args.combined_summary_json)
    run_log = Path(args.run_log)
    output_archive = Path(args.output_archive) if args.output_archive else None

    records = read_fasta_records(input_path)
    assign_output_names(records)
    structures_dir.mkdir(parents=True, exist_ok=True)
    summaries_dir.mkdir(parents=True, exist_ok=True)

    temp_dir = Path.cwd() / f".esmfold_batch_tmp_{uuid.uuid4().hex}"
    temp_dir.mkdir(parents=True, exist_ok=False)
    try:
        raw_output_dir = temp_dir / "esmfold_raw_output"
        raw_output_dir.mkdir()
        batch_input_copy = temp_dir / "batch_input.fasta"
        write_batch_fasta(records, batch_input_copy)
        summaries: list[dict] = []
        log_lines = [
            "ESMFold Batch run log",
            f"Input FASTA: {input_path}",
            f"Record count: {len(records)}",
            "",
        ]
        for index, record in enumerate(records, start=1):
            output_name = record["output_name"]
            structure_path = structures_dir / f"{output_name}.pdb"
            summary_path = summaries_dir / f"{output_name}.json"
            if structure_path.exists() and summary_path.exists():
                summary = json.loads(summary_path.read_text(encoding="utf-8"))
                summaries.append(summary)
                log_lines.append(f"[{index}/{len(records)}] SKIP {output_name}: existing outputs found")
                continue

            prepared_fasta = temp_dir / f"{output_name}.fasta"
            write_batch_fasta([record], prepared_fasta)
            raw_record_dir = raw_output_dir / output_name
            raw_record_dir.mkdir(parents=True, exist_ok=True)
            command = build_command(args, prepared_fasta, raw_record_dir)
            log_lines.append(f"[{index}/{len(records)}] START {output_name}: {' '.join(command)}")
            try:
                completed = subprocess.run(command, capture_output=True, text=True)
            except FileNotFoundError as exc:
                raise RuntimeError(
                    "Could not find the ESMFold executable. Set ESMFOLD_BINARY or --esmfold-command "
                    "to the absolute esm-fold path."
                ) from exc

            sys.stdout.write(completed.stdout or "")
            sys.stderr.write(completed.stderr or "")
            if completed.returncode == 0:
                try:
                    summary = summarize_record(record, raw_record_dir, structures_dir, summaries_dir, args)
                    summaries.append(summary)
                    log_lines.append(f"[{index}/{len(records)}] SUCCESS {output_name}: {summary['structure']}")
                except Exception as exc:
                    summary = write_failed_summary(
                        record,
                        summaries_dir,
                        args,
                        command,
                        completed.returncode,
                        str(exc),
                        completed.stdout or "",
                    )
                    summaries.append(summary)
                    log_lines.append(f"[{index}/{len(records)}] FAILED {output_name}: {exc}")
            else:
                summary = write_failed_summary(
                    record,
                    summaries_dir,
                    args,
                    command,
                    completed.returncode,
                    completed.stderr or "",
                    completed.stdout or "",
                )
                summaries.append(summary)
                log_lines.append(f"[{index}/{len(records)}] FAILED {output_name}: exit code {completed.returncode}")

        combined = {
            "tool": "ESMFold Batch",
            "record_count": len(records),
            "success_count": sum(1 for summary in summaries if summary.get("status") == "success"),
            "failed_count": sum(1 for summary in summaries if summary.get("status") == "failed"),
            "model_count": sum(1 for summary in summaries if summary.get("status") == "success"),
            "num_recycles": args.num_recycles,
            "esmfold_command": args.esmfold_command,
            "records": summaries,
        }
        combined_summary_json.write_text(json.dumps(combined, indent=2, sort_keys=True), encoding="utf-8")
        log_lines.extend(
            [
                "",
                f"Success count: {combined['success_count']}",
                f"Failed count: {combined['failed_count']}",
                f"Combined summary: {combined_summary_json}",
            ]
        )
        run_log.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
        if output_archive:
            write_archive(
                output_archive,
                [
                    (batch_input_copy, "prepared_input"),
                    (raw_output_dir, "esmfold_raw_output"),
                    (structures_dir, "structures"),
                    (summaries_dir, "summaries"),
                    (combined_summary_json, "summaries"),
                    (run_log, "logs"),
                ],
            )
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
