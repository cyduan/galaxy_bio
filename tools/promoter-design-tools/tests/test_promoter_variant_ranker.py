from __future__ import annotations

import pandas as pd

from conftest import ROOT, run_cli


def test_promoter_variant_ranker(tmp_path):
    ranked = tmp_path / "ranked.tsv"
    summary = tmp_path / "summary.html"
    strength_png = tmp_path / "strength.png"
    scatter_png = tmp_path / "scatter.png"
    log = tmp_path / "run.log"
    run_cli(
        [
            "tools/promoter_variant_ranker/promoter_variant_ranker.py",
            "--wt-conservation",
            "tools/promoter_variant_ranker/test-data/wt_conservation.tsv",
            "--mutant-conservation",
            "tools/promoter_variant_ranker/test-data/mutant_conservation.tsv",
            "--wt-strength",
            "tools/promoter_variant_ranker/test-data/wt_strength.tsv",
            "--mutant-strength",
            "tools/promoter_variant_ranker/test-data/mutant_strength.tsv",
            "--mutations",
            "tools/promoter_variant_ranker/test-data/mutations.tsv",
            "--output-tsv",
            str(ranked),
            "--summary-html",
            str(summary),
            "--strength-distribution-png",
            str(strength_png),
            "--conservation-vs-strength-png",
            str(scatter_png),
            "--log",
            str(log),
        ],
        cwd=ROOT,
    )
    df = pd.read_csv(ranked, sep="\t")
    assert df.iloc[0]["mutant_id"] == "promoter_A_mut000001"
    assert summary.exists()
    assert strength_png.exists()
    assert scatter_png.exists()

