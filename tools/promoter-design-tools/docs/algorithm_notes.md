# Algorithm Notes

## Normalization

Sequences are upper-cased, `U` is converted to `T`, and non-ACGTN characters are
replaced by the selected padding base. Long sequences are truncated to the
target length; short sequences are padded on the right. The report records
lengths, status, and warnings.

## Motif Scanning

Simple consensus mode scores every window by match fraction. PWM mode reads the
first MEME letter-probability matrix and scores each window by observed
probability divided by the maximum possible PWM score.

## Pair Scoring

The spacer is defined as:

```text
minus10_start - minus35_end - 1
```

The conservation score is:

```text
0.4 * minus35_score + 0.4 * minus10_score + 0.2 * spacer_score
```

Weights and spacer preferences are user configurable. Outside-range spacers are
strongly penalized but retained unless `--strict` is used.

## Strength Prediction

The heuristic backend combines motif conservation, spacer quality, GC content,
and AT richness around the -10 region. It produces a relative score for ranking
within a design batch and should not be treated as an experimentally calibrated
transcription rate.

## Random Mutagenesis

All random choices are controlled by `--seed`. SNVs never mutate to the same
base. Positions containing `N` are skipped by default.

