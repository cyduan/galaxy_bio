from __future__ import annotations

import sys

import pandas as pd

from conftest import ROOT, example_fasta, run_cli

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from promoter_design_tools.strength import parse_promoter_calculator_output


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


def test_promoter_calculator_tx_rate_output_is_normalized(tmp_path):
    output = tmp_path / "promoter_calculator.csv"
    output.write_text(
        "TSS_position,Tx_rate,UP_position,hex35,spacer,hex10,disc\n"
        "61,12.5,10,TTGACA,17,TATAAT,AAAA\n"
        "63,7.0,11,TTGACA,18,TATAAT,AAAA\n",
        encoding="utf-8",
    )
    predictions = parse_promoter_calculator_output(output, fallback_sequence_id="promoter_A")
    assert len(predictions) == 1
    assert predictions[0].sequence_id == "promoter_A"
    assert predictions[0].predicted_strength == 12.5
    assert predictions[0].backend == "promoter_calculator"
