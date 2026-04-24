#!/usr/bin/env python

from __future__ import annotations

import argparse
from pathlib import Path


def read_first_identifier(path: Path) -> str:
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith(">"):
                return line[1:].strip() or "query"
    return "query"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output-dir", required=True)
    parser.add_argument("--num-recycles")
    parser.add_argument("--max-tokens-per-batch")
    parser.add_argument("--chunk-size")
    parser.add_argument("--cpu-only", action="store_true")
    parser.add_argument("--cpu-offload", action="store_true")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    identifier = read_first_identifier(input_path)
    output_path = output_dir / f"{identifier}.pdb"
    output_path.write_text(
        "\n".join(
            [
                "HEADER    MOCK ESMFOLD OUTPUT",
                "TITLE     GENERATED FOR LOCAL WRAPPER VERIFICATION",
                "ATOM      1  N   MET A   1      11.104  13.207   9.317  1.00 91.25           N",
                "ATOM      2  CA  MET A   1      12.560  13.425   9.188  1.00 90.75           C",
                "ATOM      3  C   MET A   1      13.197  12.140   8.621  1.00 89.50           C",
                "ATOM      4  N   GLY B   2      14.104  11.907   7.917  1.00 87.25           N",
                "ATOM      5  CA  GLY B   2      14.560  10.625   7.388  1.00 86.75           C",
                "TER",
                "END",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
