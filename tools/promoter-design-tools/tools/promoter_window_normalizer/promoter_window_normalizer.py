#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from promoter_design_tools.fasta import normalize_records, read_fasta, write_fasta


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Normalize promoter windows to a fixed TSS-aligned length.")
    parser.add_argument("--input-fasta", required=True, help="Input promoter FASTA.")
    parser.add_argument("--output-fasta", required=True, help="Normalized promoter FASTA.")
    parser.add_argument("--report-tsv", required=True, help="Normalization report TSV.")
    parser.add_argument("--log", required=True, help="Run log path.")
    parser.add_argument("--target-length", type=int, default=81, help="Output sequence length. Default: 81.")
    parser.add_argument("--tss-offset", type=int, default=60, help="0-based +1 TSS offset. Default: 60.")
    parser.add_argument("--padding-base", default="N", help="Padding/replacement base. Default: N.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    log_lines: list[str] = []
    try:
        records = read_fasta(Path(args.input_fasta))
        results = normalize_records(records, args.target_length, args.tss_offset, args.padding_base)
        write_fasta([result.record for result in results], Path(args.output_fasta))
        with Path(args.report_tsv).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["sequence_id", "original_length", "normalized_length", "status", "warning"],
                delimiter="\t",
            )
            writer.writeheader()
            for result in results:
                writer.writerow(
                    {
                        "sequence_id": result.record.identifier,
                        "original_length": result.original_length,
                        "normalized_length": result.normalized_length,
                        "status": result.status,
                        "warning": result.warning,
                    }
                )
        log_lines.append(f"Normalized {len(results)} promoter sequence(s).")
        return 0
    except Exception as exc:
        log_lines.append(f"ERROR: {exc}")
        print(f"promoter_window_normalizer: {exc}", file=sys.stderr)
        return 1
    finally:
        Path(args.log).write_text("\n".join(log_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())

