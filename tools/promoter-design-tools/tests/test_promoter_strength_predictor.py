from __future__ import annotations

import pandas as pd

from conftest import ROOT, example_fasta, run_cli


def test_promoter_strength_predictor_heuristic(tmp_path):
    output = tmp_path / "strength.tsv"
    log = tmp_path / "run.log"
    run_cli(
        [
            "tools/promoter_strength_predictor/promoter_strength_predictor.py",
            "--input-fasta",
            str(example_fasta()),
            "--output-tsv",
            str(output),
            "--log",
            str(log),
        ],
        cwd=ROOT,
    )
    df = pd.read_csv(output, sep="\t")
    assert "predicted_strength" in df.columns
    assert df["predicted_strength"].min() > 0

