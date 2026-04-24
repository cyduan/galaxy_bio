#!/usr/bin/env python

from __future__ import annotations

import argparse
from pathlib import Path


def read_first_identifier(path: Path) -> str:
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                return line[1:].strip() or "query"
    return "query"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--program", required=True, choices=["blastn", "blastp"])
    parser.add_argument("-query", dest="query", required=True)
    parser.add_argument("-db", dest="db", required=True)
    parser.add_argument("-out", dest="out", required=True)
    parser.add_argument("-evalue", dest="evalue", required=False)
    parser.add_argument("-max_target_seqs", dest="max_target_seqs", required=False)
    parser.add_argument("-num_threads", dest="num_threads", required=False)
    parser.add_argument("-task", dest="task", required=False)
    parser.add_argument("-outfmt", dest="outfmt", required=False)
    args = parser.parse_args()

    query_name = read_first_identifier(Path(args.query))
    subject_fasta = Path(f"{args.db}.fasta")
    subject_name = read_first_identifier(subject_fasta)
    out_path = Path(args.out)
    if args.outfmt:
        out_path.write_text(
            f"{query_name}\t{subject_name}\t100.000\t8\t0\t0\t1\t8\t1\t8\t1.0e-50\t42.0\n",
            encoding="utf-8",
        )
    else:
        out_path.write_text(
            "\n".join(
                [
                    f"{args.program.upper()} mock report",
                    f"Query= {query_name}",
                    f"Subject= {subject_name}",
                    "Score = 42.0 bits, Expect = 1.0e-50",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
