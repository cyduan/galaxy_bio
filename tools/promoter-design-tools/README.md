# Promoter Design Tools

Galaxy tools for sigma70-style promoter design workflows, focused on E. coli
and common engineered bacteria. The package is intentionally split into small
tools so a Galaxy workflow can expose intermediate files instead of hiding the
entire process in one black-box script.

## Tools

- `promoter_window_normalizer`: normalize promoter FASTA windows around TSS.
- `promoter_motif_scanner`: scan -35 and -10 candidates with consensus/PWM or FIMO.
- `promoter_motif_pair_score`: pair -35/-10 hits and calculate conservation.
- `promoter_random_mutagenesis`: generate reproducible promoter libraries.
- `promoter_strength_predictor`: estimate relative promoter strength.
- `promoter_variant_ranker`: merge, rank, and visualize variant results.

Default assumptions:

- -35 consensus: `TTGACA`
- -10 consensus: `TATAAT`
- Recommended spacer: 16-18 bp
- Optimal spacer: 17 bp
- Default promoter window: -60 to +20, 81 bp total

## Quick Start

```bash
cd tools/promoter-design-tools
conda env create -f environment.yml
conda activate promoter-design-tools
python -m pytest
```

Run a minimal command-line workflow:

```bash
mkdir -p demo
cp tools/promoter_window_normalizer/test-data/example_promoters.fasta demo/input.fasta

python tools/promoter_window_normalizer/promoter_window_normalizer.py \
  --input-fasta demo/input.fasta \
  --output-fasta demo/normalized_promoters.fasta \
  --report-tsv demo/normalization_report.tsv \
  --log demo/normalizer.log

python tools/promoter_motif_scanner/promoter_motif_scanner.py \
  --input-fasta demo/normalized_promoters.fasta \
  --output-tsv demo/motif_hits.tsv \
  --log demo/scanner.log

python tools/promoter_motif_pair_score/promoter_motif_pair_score.py \
  --motif-hits demo/motif_hits.tsv \
  --promoters-fasta demo/normalized_promoters.fasta \
  --output-tsv demo/wt_conservation.tsv \
  --log demo/pair.log

python tools/promoter_random_mutagenesis/promoter_random_mutagenesis.py \
  --input-fasta demo/normalized_promoters.fasta \
  --conservation-tsv demo/wt_conservation.tsv \
  --output-fasta demo/mutants.fasta \
  --mutations-tsv demo/mutations.tsv \
  --log demo/mutagenesis.log \
  --mutants-per-sequence 5 \
  --seed 42

python tools/promoter_motif_scanner/promoter_motif_scanner.py \
  --input-fasta demo/mutants.fasta \
  --output-tsv demo/mutant_hits.tsv \
  --log demo/mutant_scanner.log

python tools/promoter_motif_pair_score/promoter_motif_pair_score.py \
  --motif-hits demo/mutant_hits.tsv \
  --promoters-fasta demo/mutants.fasta \
  --output-tsv demo/mutant_conservation.tsv \
  --log demo/mutant_pair.log

python tools/promoter_strength_predictor/promoter_strength_predictor.py \
  --input-fasta demo/normalized_promoters.fasta \
  --output-tsv demo/wt_strength.tsv \
  --log demo/wt_strength.log

python tools/promoter_strength_predictor/promoter_strength_predictor.py \
  --input-fasta demo/mutants.fasta \
  --output-tsv demo/mutant_strength.tsv \
  --log demo/mutant_strength.log

python tools/promoter_variant_ranker/promoter_variant_ranker.py \
  --wt-conservation demo/wt_conservation.tsv \
  --mutant-conservation demo/mutant_conservation.tsv \
  --wt-strength demo/wt_strength.tsv \
  --mutant-strength demo/mutant_strength.tsv \
  --mutations demo/mutations.tsv \
  --output-tsv demo/ranked_promoter_variants.tsv \
  --summary-html demo/summary.html \
  --strength-distribution-png demo/strength_distribution.png \
  --conservation-vs-strength-png demo/conservation_vs_strength.png \
  --log demo/ranker.log
```

## Galaxy Installation

Install the shared Python package into the Python environment that will run
Galaxy jobs:

```bash
python -m pip install -e /data/tools/galaxy_bio/tools/promoter-design-tools
```

The repository-level Galaxy config can load these tools directly from
`tools/promoter-design-tools`. If installing into a different Galaxy instance,
add this section to a Galaxy tool config:

```xml
<section id="promoter_design_tools" name="Promoter Design Tools">
  <tool file="promoter-design-tools/tools/promoter_window_normalizer/promoter_window_normalizer.xml" />
  <tool file="promoter-design-tools/tools/promoter_motif_scanner/promoter_motif_scanner.xml" />
  <tool file="promoter-design-tools/tools/promoter_motif_pair_score/promoter_motif_pair_score.xml" />
  <tool file="promoter-design-tools/tools/promoter_random_mutagenesis/promoter_random_mutagenesis.xml" />
  <tool file="promoter-design-tools/tools/promoter_strength_predictor/promoter_strength_predictor.xml" />
  <tool file="promoter-design-tools/tools/promoter_variant_ranker/promoter_variant_ranker.xml" />
</section>
```

Restart Galaxy after editing the tool config. Import
`workflows/promoter_conservation_mutagenesis_strength.ga` from the Galaxy
workflow menu.

## Planemo

```bash
planemo lint tools/promoter_window_normalizer/promoter_window_normalizer.xml
planemo test tools/promoter_window_normalizer/promoter_window_normalizer.xml

for xml in tools/promoter_*/*.xml; do
  planemo lint "$xml"
done
```

## MEME Suite and FIMO

`promoter_motif_scanner` works without external tools in `simple` mode. FIMO
mode requires MEME Suite and a MEME motif file, for example
`data/sigma70_default_motifs/sigma70_combined.meme`.

## Known Limitations

- -35/-10 conservation is not equivalent to real promoter strength.
- The heuristic backend is a relative screening score, not an absolute
  transcription initiation rate.
- Promoter Calculator support is an adapter; users must install and configure
  the external calculator themselves.
- Non-sigma70 promoters should use appropriate consensus or PWM files.
