#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from promoter_design_tools.fasta import read_fasta
from promoter_design_tools.pairing import score_pairs, write_pair_scores


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Pair -35 and -10 motif hits and score sigma70 promoter conservation.")
    parser.add_argument("--motif-hits", required=True)
    parser.add_argument("--promoters-fasta", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("--minus35-motif-name", default="minus35")
    parser.add_argument("--minus10-motif-name", default="minus10")
    parser.add_argument("--min-spacer", type=int, default=16)
    parser.add_argument("--max-spacer", type=int, default=18)
    parser.add_argument("--optimal-spacer", type=int, default=17)
    parser.add_argument("--weight-minus35", type=float, default=0.4)
    parser.add_argument("--weight-minus10", type=float, default=0.4)
    parser.add_argument("--weight-spacer", type=float, default=0.2)
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    log_lines: list[str] = []
    try:
        records = read_fasta(Path(args.promoters_fasta))
        results = score_pairs(
            records,
            Path(args.motif_hits),
            args.minus35_motif_name,
            args.minus10_motif_name,
            args.min_spacer,
            args.max_spacer,
            args.optimal_spacer,
            args.weight_minus35,
            args.weight_minus10,
            args.weight_spacer,
            args.strict,
        )
        write_pair_scores(results, Path(args.output_tsv))
        log_lines.append(f"Scored motif pairs for {len(results)} sequence(s).")
        return 0
    except Exception as exc:
        log_lines.append(f"ERROR: {exc}")
        print(f"promoter_motif_pair_score: {exc}", file=sys.stderr)
        return 1
    finally:
        Path(args.log).write_text("\n".join(log_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())

