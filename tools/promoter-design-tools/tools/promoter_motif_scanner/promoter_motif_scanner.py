#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from promoter_design_tools.fasta import read_fasta
from promoter_design_tools.motifs import model_from_consensus, read_first_meme_motif, run_fimo, scan_simple, write_hits_tsv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Scan promoter FASTA records for -35 and -10 motif candidates.")
    parser.add_argument("--input-fasta", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("--mode", choices=["simple", "fimo"], default="simple")
    parser.add_argument("--minus35-consensus", default="TTGACA")
    parser.add_argument("--minus10-consensus", default="TATAAT")
    parser.add_argument("--minus35-meme", default="")
    parser.add_argument("--minus10-meme", default="")
    parser.add_argument("--meme-file", default="", help="Combined MEME motif file for FIMO mode.")
    parser.add_argument("--fimo-binary", default="fimo")
    parser.add_argument("--scan-reverse-complement", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    log_lines: list[str] = []
    try:
        if args.mode == "fimo":
            if not args.meme_file:
                raise ValueError("--meme-file is required when --mode=fimo.")
            hits = run_fimo(Path(args.input_fasta), Path(args.meme_file), args.fimo_binary)
        else:
            records = read_fasta(Path(args.input_fasta))
            minus35 = (
                read_first_meme_motif(Path(args.minus35_meme), "minus35")
                if args.minus35_meme
                else model_from_consensus("minus35", args.minus35_consensus)
            )
            minus10 = (
                read_first_meme_motif(Path(args.minus10_meme), "minus10")
                if args.minus10_meme
                else model_from_consensus("minus10", args.minus10_consensus)
            )
            minus35 = type(minus35)("minus35", minus35.consensus, minus35.pwm)
            minus10 = type(minus10)("minus10", minus10.consensus, minus10.pwm)
            hits = scan_simple(records, [minus35, minus10], args.scan_reverse_complement)
        write_hits_tsv(hits, Path(args.output_tsv))
        log_lines.append(f"Wrote {len(hits)} motif hit(s).")
        return 0
    except Exception as exc:
        log_lines.append(f"ERROR: {exc}")
        print(f"promoter_motif_scanner: {exc}", file=sys.stderr)
        return 1
    finally:
        Path(args.log).write_text("\n".join(log_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())

