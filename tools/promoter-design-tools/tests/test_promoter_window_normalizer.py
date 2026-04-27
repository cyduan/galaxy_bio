from __future__ import annotations

import pandas as pd

from conftest import ROOT, example_fasta, run_cli


def test_promoter_window_normalizer(tmp_path):
    output_fasta = tmp_path / "normalized.fasta"
    report = tmp_path / "report.tsv"
    log = tmp_path / "run.log"
    run_cli(
        [
            "tools/promoter_window_normalizer/promoter_window_normalizer.py",
            "--input-fasta",
            str(example_fasta()),
            "--output-fasta",
            str(output_fasta),
            "--report-tsv",
            str(report),
            "--log",
            str(log),
        ],
        cwd=ROOT,
    )
    assert "TTGACA" in output_fasta.read_text()
    df = pd.read_csv(report, sep="\t")
    assert set(df["normalized_length"]) == {81}

