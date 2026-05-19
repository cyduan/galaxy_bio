# Homolog Search and MSA

Tool 3 for the HotSpot Wizard-like protein engineering workflow.

It performs:

1. Homolog search with BLAST+ `blastp` or MMseqs2 `easy-search`.
2. Search-hit filtering by E-value, identity, coverage, and alignment length.
3. Redundancy filtering with CD-HIT, MMseqs2 `easy-cluster`, or exact internal deduplication.
4. MSA construction with MAFFT or MUSCLE.

## Install

Install the required command-line tools into the shared HotSpot Wizard environment:

```bash
conda install -p /data/conda_envs/hotspot_wizard -y \
  -c conda-forge -c bioconda \
  blast cd-hit mafft muscle mmseqs2
```

The Galaxy job destination is configured in `config/job_conf.yml` with binary
paths such as:

```yaml
BLASTP_BINARY: /data/conda_envs/hotspot_wizard/bin/blastp
MAKEBLASTDB_BINARY: /data/conda_envs/hotspot_wizard/bin/makeblastdb
CDHIT_BINARY: /data/conda_envs/hotspot_wizard/bin/cd-hit
MAFFT_BINARY: /data/conda_envs/hotspot_wizard/bin/mafft
MUSCLE_BINARY: /data/conda_envs/hotspot_wizard/bin/muscle
MMSEQS_BINARY: /data/conda_envs/hotspot_wizard/bin/mmseqs
```

## Environment Checks

```bash
/data/conda_envs/hotspot_wizard/bin/blastp -version
/data/conda_envs/hotspot_wizard/bin/makeblastdb -version
/data/conda_envs/hotspot_wizard/bin/cd-hit -h
/data/conda_envs/hotspot_wizard/bin/mafft --version
/data/conda_envs/hotspot_wizard/bin/muscle -version
/data/conda_envs/hotspot_wizard/bin/mmseqs version
```

## Smoke Test Without External Binaries

This uses the built-in curated-input path and exact deduplication:

```bash
python /data/tools/galaxy_bio/tools/protein_analysis/homolog_search_msa/run_homolog_search_msa.py \
  --query-fasta /data/tools/galaxy_bio/tools/protein_analysis/homolog_search_msa/test-data/homolog_query.fasta \
  --subject-fasta /data/tools/galaxy_bio/tools/protein_analysis/homolog_search_msa/test-data/homolog_subjects.fasta \
  --homologs-fasta /tmp/homologs.fasta \
  --filtered-homologs-fasta /tmp/filtered_homologs.fasta \
  --msa-fasta /tmp/msa.fasta \
  --summary-tsv /tmp/homolog_summary.tsv \
  --raw-search-output /tmp/raw_search.tsv \
  --cluster-report /tmp/cluster_report.txt \
  --msa-log /tmp/msa_log.txt \
  --run-log /tmp/run_log.txt \
  --search-backend all_subjects \
  --filter-backend internal \
  --msa-backend copy \
  --include-query-in-msa
```

## Production Example

```bash
/data/conda_envs/hotspot_wizard/bin/python \
  /data/tools/galaxy_bio/tools/protein_analysis/homolog_search_msa/run_homolog_search_msa.py \
  --query-fasta target.fasta \
  --subject-fasta local_family_or_database.fasta \
  --homologs-fasta homologs.fasta \
  --filtered-homologs-fasta filtered_homologs.fasta \
  --msa-fasta msa.fasta \
  --summary-tsv homolog_summary.tsv \
  --raw-search-output raw_search.tsv \
  --cluster-report cluster_report.txt \
  --msa-log msa_log.txt \
  --run-log run_log.txt \
  --search-backend blastp \
  --filter-backend cd-hit \
  --msa-backend mafft \
  --evalue 1e-5 \
  --min-identity 20 \
  --min-query-coverage 40 \
  --identity-threshold 0.9 \
  --threads 4 \
  --include-query-in-msa
```

After editing Galaxy config, restart Galaxy:

```bash
cd /data/tools/galaxy_bio
./run.sh --stop-daemon
./run.sh --daemon
```
