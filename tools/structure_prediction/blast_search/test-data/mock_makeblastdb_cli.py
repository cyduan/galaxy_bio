#!/usr/bin/env python

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("-in", dest="input_path", required=True)
    parser.add_argument("-dbtype", required=True)
    parser.add_argument("-out", dest="out_prefix", required=True)
    args = parser.parse_args()
    marker = Path(f"{args.out_prefix}.mockdb")
    marker.write_text(f"input={args.input_path}\ndbtype={args.dbtype}\n", encoding="utf-8")
    shutil.copyfile(args.input_path, f"{args.out_prefix}.fasta")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
