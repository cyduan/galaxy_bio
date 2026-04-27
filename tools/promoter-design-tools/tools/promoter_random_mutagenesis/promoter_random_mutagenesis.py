#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from promoter_design_tools.fasta import read_fasta, write_fasta
from promoter_design_tools.mutagenesis import generate_mutants, write_mutations


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate random promoter mutant libraries.")
    parser.add_argument("--input-fasta", required=True)
    parser.add_argument("--conservation-tsv", default="")
    parser.add_argument("--output-fasta", required=True)
    parser.add_argument("--mutations-tsv", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("--mutants-per-sequence", type=int, default=10)
    parser.add_argument("--mutation-mode", choices=["snv", "insertion", "deletion", "mixed"], default="snv")
    parser.add_argument("--mutation-count-mode", choices=["fixed", "rate"], default="fixed")
    parser.add_argument("--mutation-count", type=int, default=1)
    parser.add_argument("--mutation-rate", type=float, default=0.01)
    parser.add_argument("--region-mode", choices=["full", "minus35", "minus10", "spacer", "non_core", "custom"], default="full")
    parser.add_argument("--custom-start", type=int)
    parser.add_argument("--custom-end", type=int)
    parser.add_argument("--keep-length", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--avoid-stop-codons", action="store_true", help="Reserved for CDS workflows; ignored for promoters.")
    parser.add_argument("--allow-n", action="store_true")
    parser.add_argument("--seed", type=int, default=1)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    log_lines: list[str] = []
    try:
        records = read_fasta(Path(args.input_fasta))
        conservation_path = Path(args.conservation_tsv) if args.conservation_tsv else None
        mutants, events = generate_mutants(
            records=records,
            conservation_path=conservation_path,
            mutants_per_sequence=args.mutants_per_sequence,
            mutation_mode=args.mutation_mode,
            mutation_count_mode=args.mutation_count_mode,
            mutation_count=args.mutation_count,
            mutation_rate=args.mutation_rate,
            region_mode=args.region_mode,
            custom_start=args.custom_start,
            custom_end=args.custom_end,
            keep_length=args.keep_length,
            allow_n=args.allow_n,
            seed=args.seed,
        )
        write_fasta(mutants, Path(args.output_fasta))
        write_mutations(events, Path(args.mutations_tsv))
        log_lines.append(f"Generated {len(mutants)} mutant sequence(s) and {len(events)} mutation event(s).")
        if args.avoid_stop_codons:
            log_lines.append("avoid_stop_codons was provided but ignored because promoters are non-coding DNA.")
        return 0
    except Exception as exc:
        log_lines.append(f"ERROR: {exc}")
        print(f"promoter_random_mutagenesis: {exc}", file=sys.stderr)
        return 1
    finally:
        Path(args.log).write_text("\n".join(log_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())

