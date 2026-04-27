from __future__ import annotations

import pandas as pd

from conftest import ROOT, run_cli


def test_promoter_motif_pair_score(tmp_path):
    output = tmp_path / "conservation.tsv"
    log = tmp_path / "run.log"
    run_cli(
        [
            "tools/promoter_motif_pair_score/promoter_motif_pair_score.py",
            "--motif-hits",
            "tools/promoter_motif_pair_score/test-data/motif_hits.tsv",
            "--promoters-fasta",
            "tools/promoter_motif_pair_score/test-data/example_promoters.fasta",
            "--output-tsv",
            str(output),
            "--log",
            str(log),
        ],
        cwd=ROOT,
    )
    df = pd.read_csv(output, sep="\t")
    row = df[df["sequence_id"] == "promoter_A"].iloc[0]
    assert row["spacer_len"] == 17
    assert row["conservation_score"] == 1.0

