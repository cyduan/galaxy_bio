#!/usr/bin/env python

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Galaxy compatibility runner for ESMFold via the Python API.")
    parser.add_argument("-i", "--fasta", required=True, help="Path to the input FASTA file.")
    parser.add_argument("-o", "--pdb", required=True, help="Directory where PDB outputs will be written.")
    parser.add_argument("--num-recycles", type=int, default=4)
    parser.add_argument("--max-tokens-per-batch", type=int)
    parser.add_argument("--chunk-size", type=int)
    parser.add_argument("--cpu-only", action="store_true")
    parser.add_argument("--cpu-offload", action="store_true")
    args = parser.parse_args()
    if args.cpu_only and args.cpu_offload:
        parser.error("--cpu-only and --cpu-offload are mutually exclusive")
    return args


def sanitize_identifier(identifier: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9_.-]+", "_", identifier).strip("._")
    return sanitized or "query"


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    identifier: str | None = None
    sequence_lines: list[str] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records.append((identifier, "".join(sequence_lines)))
                identifier = line[1:].strip() or f"sequence_{len(records) + 1}"
                sequence_lines = []
            else:
                if identifier is None:
                    raise ValueError("Input is not a valid FASTA file: sequence data found before a header line.")
                sequence_lines.append(line)
    if identifier is not None:
        records.append((identifier, "".join(sequence_lines)))
    if not records:
        raise ValueError("No FASTA records were found in the input dataset.")
    return records


def main() -> int:
    args = parse_args()
    if args.cpu_offload:
        raise RuntimeError(
            "The Galaxy API fallback runner supports GPU and CPU-only execution, "
            "but not the official CPU offload mode. Use GPU/CPU-only or install the official esm-fold CLI."
        )

    try:
        import esm
        import torch
    except ImportError as exc:
        raise RuntimeError(
            "Could not import the official fair-esm runtime in the configured ESMFOLD_PYTHON environment."
        ) from exc

    records = read_fasta(Path(args.fasta))
    output_dir = Path(args.pdb)
    output_dir.mkdir(parents=True, exist_ok=True)

    model = esm.pretrained.esmfold_v1()
    model = model.eval()
    if args.chunk_size is not None:
        model.set_chunk_size(args.chunk_size)

    use_cuda = torch.cuda.is_available() and not args.cpu_only
    if use_cuda:
        model = model.cuda()
    else:
        sys.stderr.write("Running ESMFold on CPU.\n")
        model = model.cpu()

    with torch.no_grad():
        for identifier, sequence in records:
            output_path = output_dir / f"{sanitize_identifier(identifier)}.pdb"
            try:
                pdb_text = model.infer_pdb(sequence, num_recycles=args.num_recycles)
            except TypeError:
                pdb_text = model.infer_pdb(sequence)
            output_path.write_text(pdb_text, encoding="utf-8")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
