# Usage

## Minimal FASTA Input

```fasta
>promoter_A
ACGTACGTACGTACGTACGTTGACAGCGTATGCGTATGCGTATATAATCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGT
```

The default window is interpreted as -60 to +20 relative to TSS, with +1 at
offset 60.

## Tool Examples

Normalize windows:

```bash
python tools/promoter_window_normalizer/promoter_window_normalizer.py \
  --input-fasta input.fasta \
  --output-fasta normalized_promoters.fasta \
  --report-tsv normalization_report.tsv \
  --log normalizer.log
```

Scan motifs:

```bash
python tools/promoter_motif_scanner/promoter_motif_scanner.py \
  --input-fasta normalized_promoters.fasta \
  --output-tsv motif_hits.tsv \
  --log scanner.log \
  --minus35-consensus TTGACA \
  --minus10-consensus TATAAT
```

Score motif pairs:

```bash
python tools/promoter_motif_pair_score/promoter_motif_pair_score.py \
  --motif-hits motif_hits.tsv \
  --promoters-fasta normalized_promoters.fasta \
  --output-tsv promoter_conservation.tsv \
  --log pair.log
```

Generate a library:

```bash
python tools/promoter_random_mutagenesis/promoter_random_mutagenesis.py \
  --input-fasta normalized_promoters.fasta \
  --conservation-tsv promoter_conservation.tsv \
  --output-fasta mutants.fasta \
  --mutations-tsv mutations.tsv \
  --log mutagenesis.log \
  --mutants-per-sequence 20 \
  --region-mode non_core \
  --seed 42
```

Predict relative strength:

```bash
python tools/promoter_strength_predictor/promoter_strength_predictor.py \
  --input-fasta mutants.fasta \
  --output-tsv promoter_strength.tsv \
  --log strength.log
```

Rank variants:

```bash
python tools/promoter_variant_ranker/promoter_variant_ranker.py \
  --wt-conservation wt_conservation.tsv \
  --mutant-conservation mutant_conservation.tsv \
  --wt-strength wt_strength.tsv \
  --mutant-strength mutant_strength.tsv \
  --mutations mutations.tsv \
  --output-tsv ranked_promoter_variants.tsv \
  --summary-html summary.html \
  --strength-distribution-png strength_distribution.png \
  --conservation-vs-strength-png conservation_vs_strength.png \
  --log ranker.log
```

## Output Examples

`promoter_conservation.tsv` includes:

```text
sequence_id minus35_seq minus35_start minus35_end minus35_score minus10_seq minus10_start minus10_end minus10_score spacer_len spacer_score conservation_score status
```

`ranked_promoter_variants.tsv` includes:

```text
mutant_id parent_id mutation_summary mutated_regions parent_conservation_score mutant_conservation_score delta_conservation_score parent_strength mutant_strength fold_change_vs_parent rank recommendation
```

