from __future__ import annotations

import pandas as pd
from Bio import SeqIO

from conftest import ROOT, run_cli


def test_promoter_random_mutagenesis_reproducible(tmp_path):
    mutants = tmp_path / "mutants.fasta"
    mutations = tmp_path / "mutations.tsv"
    log = tmp_path / "run.log"
    args = [
        "tools/promoter_random_mutagenesis/promoter_random_mutagenesis.py",
        "--input-fasta",
        "tools/promoter_random_mutagenesis/test-data/example_promoters.fasta",
        "--conservation-tsv",
        "tools/promoter_random_mutagenesis/test-data/promoter_conservation.tsv",
        "--output-fasta",
        str(mutants),
        "--mutations-tsv",
        str(mutations),
        "--log",
        str(log),
        "--mutants-per-sequence",
        "3",
        "--seed",
        "42",
    ]
    run_cli(args, cwd=ROOT)
    records = list(SeqIO.parse(mutants, "fasta"))
    df = pd.read_csv(mutations, sep="\t")
    assert len(records) == 3
    assert len(df) == 3
    assert df.iloc[0]["mutant_id"] == "promoter_A_mut000001"

