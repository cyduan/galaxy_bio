#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from promoter_design_tools.fasta import read_fasta
from promoter_design_tools.strength import predict_heuristic, predict_promoter_calculator, write_strength


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Predict relative promoter strength.")
    parser.add_argument("--input-fasta", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("--backend", choices=["heuristic", "promoter_calculator"], default="heuristic")
    parser.add_argument("--organism", default="Escherichia_coli")
    parser.add_argument("--tss-offset", type=int, default=60)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--promoter-calculator-command", default="")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    log_lines: list[str] = []
    try:
        records = read_fasta(Path(args.input_fasta))
        if args.backend == "heuristic":
            predictions = predict_heuristic(records, args.tss_offset)
        else:
            predictions = predict_promoter_calculator(
                records,
                organism=args.organism,
                threads=args.threads,
                command_template=args.promoter_calculator_command or None,
            )
        write_strength(predictions, Path(args.output_tsv))
        log_lines.append(f"Wrote {len(predictions)} strength prediction(s) with backend={args.backend}.")
        return 0
    except Exception as exc:
        log_lines.append(f"ERROR: {exc}")
        print(f"promoter_strength_predictor: {exc}", file=sys.stderr)
        return 1
    finally:
        Path(args.log).write_text("\n".join(log_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())

