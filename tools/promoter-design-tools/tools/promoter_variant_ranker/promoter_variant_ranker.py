#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from promoter_design_tools.ranking import rank_variants, write_plots, write_summary_html


def optional_float(value: str) -> float | None:
    if value in {"", "None", "none", "NA"}:
        return None
    return float(value)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merge, rank, and visualize promoter variant analysis results.")
    parser.add_argument("--wt-conservation", required=True)
    parser.add_argument("--mutant-conservation", required=True)
    parser.add_argument("--wt-strength", required=True)
    parser.add_argument("--mutant-strength", required=True)
    parser.add_argument("--mutations", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--summary-html", required=True)
    parser.add_argument("--strength-distribution-png", required=True)
    parser.add_argument("--conservation-vs-strength-png", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("--rank-by", choices=["predicted_strength", "fold_change", "conservation_score", "custom"], default="predicted_strength")
    parser.add_argument("--min-conservation-score", type=optional_float, default=None)
    parser.add_argument("--max-conservation-score", type=optional_float, default=None)
    parser.add_argument("--min-fold-change", type=optional_float, default=None)
    parser.add_argument("--max-fold-change", type=optional_float, default=None)
    parser.add_argument("--keep-core-motif", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    log_lines: list[str] = []
    try:
        ranked = rank_variants(
            wt_conservation=Path(args.wt_conservation),
            mutant_conservation=Path(args.mutant_conservation),
            wt_strength=Path(args.wt_strength),
            mutant_strength=Path(args.mutant_strength),
            mutations_tsv=Path(args.mutations),
            rank_by=args.rank_by,
            min_conservation_score=args.min_conservation_score,
            max_conservation_score=args.max_conservation_score,
            min_fold_change=args.min_fold_change,
            max_fold_change=args.max_fold_change,
            keep_core_motif=args.keep_core_motif,
        )
        ranked.to_csv(args.output_tsv, sep="\t", index=False)
        write_summary_html(ranked, Path(args.summary_html))
        write_plots(ranked, Path(args.strength_distribution_png), Path(args.conservation_vs_strength_png))
        log_lines.append(f"Ranked {len(ranked)} promoter variant(s).")
        return 0
    except Exception as exc:
        log_lines.append(f"ERROR: {exc}")
        print(f"promoter_variant_ranker: {exc}", file=sys.stderr)
        return 1
    finally:
        Path(args.log).write_text("\n".join(log_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())

