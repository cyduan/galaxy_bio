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

Recommended server environment:

```bash
conda create -p /data/conda_envs/promoter_tools -y \
  -c conda-forge -c bioconda \
  python=3.11 biopython pandas numpy matplotlib scipy meme promotercalculator
```

`config/job_conf.yml` maps all promoter tools to `promoter_design_local`, which
uses `/data/conda_envs/promoter_tools/bin` in `PATH`.

Check FIMO:

```bash
/data/conda_envs/promoter_tools/bin/fimo --version

mkdir -p /data/test/promoter_tools/fimo
cat > /data/test/promoter_tools/promoters.fasta <<'FASTA'
>promoter_A
ACGTACGTACGTACGTACGTTGACAGCGTATGCGTATGCGTATATAATCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGT
FASTA

/data/conda_envs/promoter_tools/bin/fimo \
  --oc /data/test/promoter_tools/fimo/out \
  /data/tools/galaxy_bio/tools/promoter-design-tools/data/sigma70_default_motifs/sigma70_combined.meme \
  /data/test/promoter_tools/promoters.fasta

head /data/test/promoter_tools/fimo/out/fimo.tsv
```

Check the Galaxy scanner wrapper in FIMO mode:

```bash
cd /data/tools/galaxy_bio
PATH=/data/conda_envs/promoter_tools/bin:$PATH \
/data/conda_envs/promoter_tools/bin/python \
  tools/promoter-design-tools/tools/promoter_motif_scanner/promoter_motif_scanner.py \
  --input-fasta /data/test/promoter_tools/promoters.fasta \
  --mode fimo \
  --meme-file tools/promoter-design-tools/data/sigma70_default_motifs/sigma70_combined.meme \
  --fimo-binary /data/conda_envs/promoter_tools/bin/fimo \
  --output-tsv /data/test/promoter_tools/fimo/motif_hits.tsv \
  --log /data/test/promoter_tools/fimo/scanner.log
```

## Optional Promoter Calculator Backend

The heuristic backend does not need Promoter Calculator. To enable
`backend=promoter_calculator`, install `promotercalculator` in the
`promoter_tools` environment and check that the CLI is available:

```bash
/data/conda_envs/promoter_tools/bin/promoter-calculator --help
```

The Galaxy wrapper calls the command template from `PROMOTER_CALCULATOR_COMMAND`
in `config/job_conf.yml`:

```text
/data/conda_envs/promoter_tools/bin/promoter-calculator -i {input_fasta} -o {output_tsv} --threads {threads}
```

The command must write a TSV/CSV containing `sequence_id` and
`predicted_strength`, or aliases such as `id`/`strength`,
`name`/`expression`, or `promoter`/`tx_rate`.

Wrapper smoke test:

```bash
mkdir -p /data/test/promoter_tools/promoter_calculator

PROMOTER_CALCULATOR_COMMAND="/data/conda_envs/promoter_tools/bin/promoter-calculator -i {input_fasta} -o {output_tsv} --threads {threads}" \
PATH=/data/conda_envs/promoter_tools/bin:$PATH \
/data/conda_envs/promoter_tools/bin/python \
  /data/tools/galaxy_bio/tools/promoter-design-tools/tools/promoter_strength_predictor/promoter_strength_predictor.py \
  --input-fasta /data/test/promoter_tools/promoters.fasta \
  --backend promoter_calculator \
  --output-tsv /data/test/promoter_tools/promoter_calculator/promoter_strength.tsv \
  --log /data/test/promoter_tools/promoter_calculator/strength.log

cat /data/test/promoter_tools/promoter_calculator/promoter_strength.tsv
```
