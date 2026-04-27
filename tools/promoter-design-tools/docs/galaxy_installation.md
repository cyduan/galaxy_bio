# Galaxy Installation

## Local Install

1. Create the dependency environment:

   ```bash
   conda env create -f environment.yml
   conda activate promoter-design-tools
   python -m pip install -e /data/tools/galaxy_bio/tools/promoter-design-tools
   ```

2. If Galaxy jobs run in an existing Galaxy virtualenv instead of this conda
   environment, install the package there as well:

   ```bash
   /data/tools/galaxy_bio/.venv/bin/python -m pip install -e /data/tools/galaxy_bio/tools/promoter-design-tools
   ```

3. Add a section to a Galaxy tool config. In this repository,
   `config/tool_conf.promoter_design.xml` already contains:

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

4. Restart Galaxy.

## Planemo Testing

Run lint and tests from this repository root:

```bash
for xml in tools/promoter_*/*.xml; do
  planemo lint "$xml"
  planemo test "$xml"
done
```

## Optional FIMO

The simple scanner requires no external motif software. To use `mode=fimo`,
install MEME Suite and ensure `fimo` is available in the Galaxy job
environment. Use `data/sigma70_default_motifs/sigma70_combined.meme` as a
starting motif file.
