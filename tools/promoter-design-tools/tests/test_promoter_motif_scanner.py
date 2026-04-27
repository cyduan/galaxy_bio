from __future__ import annotations

import pandas as pd

from conftest import ROOT, example_fasta, run_cli


def test_promoter_motif_scanner_simple(tmp_path):
    hits = tmp_path / "hits.tsv"
    log = tmp_path / "run.log"
    run_cli(
        [
            "tools/promoter_motif_scanner/promoter_motif_scanner.py",
            "--input-fasta",
            str(example_fasta()),
            "--output-tsv",
            str(hits),
            "--log",
            str(log),
        ],
        cwd=ROOT,
    )
    df = pd.read_csv(hits, sep="\t")
    assert {"minus35", "minus10"}.issubset(set(df["motif_name"]))
    assert ((df["sequence_id"] == "promoter_A") & (df["matched_sequence"] == "TTGACA")).any()

