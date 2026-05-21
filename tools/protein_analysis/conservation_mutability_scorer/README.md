# Conservation / Mutability Scorer

Tool 4 for the HotSpot Wizard-like protein engineering workflow.

It reads an MSA FASTA and writes the workflow-facing table:

```text
residue_conservation.tsv
```

The table contains per-residue conservation, mutability, entropy, consensus, and
accepted amino-acid information. The wrapper can call Rate4Site and also
preserves the raw Rate4Site output.

## Install Rate4Site

First install the shared HotSpot Wizard command-line environment if it does not
exist. Then install Rate4Site. Depending on the server image, either a system
package or a conda package may be available.

System package option:

```bash
sudo apt-get update
sudo apt-get install -y rate4site
which rate4site
rate4site -h
```

Conda option, if available on your configured channels:

```bash
conda install -p /data/conda_envs/hotspot_wizard -y -c bioconda -c conda-forge rate4site
/data/conda_envs/hotspot_wizard/bin/rate4site -h
```

If Rate4Site is installed system-wide, set this in `config/job_conf.yml`:

```yaml
RATE4SITE_BINARY: /usr/bin/rate4site
```

If it is installed in the conda environment:

```yaml
RATE4SITE_BINARY: /data/conda_envs/hotspot_wizard/bin/rate4site
```

## Smoke Test

Entropy-only smoke test:

```bash
/data/conda_envs/hotspot_wizard/bin/python \
  /data/tools/galaxy_bio/tools/protein_analysis/conservation_mutability_scorer/run_conservation_mutability_scorer.py \
  --msa-fasta /data/tools/galaxy_bio/tools/protein_analysis/conservation_mutability_scorer/test-data/mini_msa.fasta \
  --output-tsv /tmp/residue_conservation.tsv \
  --rate4site-output /tmp/rate4site_raw_output.txt \
  --run-log /tmp/conservation_run_log.txt \
  --backend entropy_only \
  --reference-mode first_record \
  --chain A
```

Rate4Site smoke test:

```bash
/data/conda_envs/hotspot_wizard/bin/python \
  /data/tools/galaxy_bio/tools/protein_analysis/conservation_mutability_scorer/run_conservation_mutability_scorer.py \
  --msa-fasta /data/tools/galaxy_bio/tools/protein_analysis/conservation_mutability_scorer/test-data/mini_msa.fasta \
  --output-tsv /tmp/residue_conservation.tsv \
  --rate4site-output /tmp/rate4site_raw_output.txt \
  --run-log /tmp/conservation_run_log.txt \
  --backend rate4site \
  --reference-mode first_record \
  --chain A \
  --rate4site-binary /usr/bin/rate4site
```

## Galaxy

The tool is listed under `Protein Analysis` as `Conservation / Mutability
Scorer` and is mapped to `hotspot_wizard_local` in `config/job_conf.yml`.

Restart Galaxy after editing config:

```bash
cd /data/tools/galaxy_bio
./run.sh --stop-daemon
./run.sh --daemon
```
