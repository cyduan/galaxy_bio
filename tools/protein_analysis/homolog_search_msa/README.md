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
HOTSPOT_SWISSPROT_FASTA: /data/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta
HOTSPOT_SWISSPROT_BLASTDB: /data/databases/uniprot/blastdb/uniprot_sprot
HOTSPOT_UNIREF90_FASTA: /data/databases/uniprot/current_release/uniref/uniref90/uniref90.fasta
HOTSPOT_UNIREF90_BLASTDB: /data/databases/uniprot/blastdb/uniref90
HOTSPOT_UNIREF50_FASTA: /data/databases/uniprot/current_release/uniref/uniref50/uniref50.fasta
HOTSPOT_UNIREF50_BLASTDB: /data/databases/uniprot/blastdb/uniref50
```

## Download UniProt Databases

Recommended server locations:

```text
/data/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta
/data/databases/uniprot/current_release/uniref/uniref90/uniref90.fasta
/data/databases/uniprot/current_release/uniref/uniref50/uniref50.fasta
/data/databases/uniprot/blastdb/
```

Download FASTA files:

```bash
mkdir -p /data/databases/uniprot/current_release/knowledgebase/complete
mkdir -p /data/databases/uniprot/current_release/uniref/uniref90
mkdir -p /data/databases/uniprot/current_release/uniref/uniref50
mkdir -p /data/databases/uniprot/blastdb

cd /data/databases/uniprot/current_release/knowledgebase/complete
wget -c https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
wget -c https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz.md5
md5sum -c uniprot_sprot.fasta.gz.md5
gunzip -kf uniprot_sprot.fasta.gz

cd /data/databases/uniprot/current_release/uniref/uniref90
wget -c https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/uniref90.fasta.gz
wget -c https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/uniref90.fasta.gz.md5
md5sum -c uniref90.fasta.gz.md5
gunzip -kf uniref90.fasta.gz

cd /data/databases/uniprot/current_release/uniref/uniref50
wget -c https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref50/uniref50.fasta.gz
wget -c https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref50/uniref50.fasta.gz.md5
md5sum -c uniref50.fasta.gz.md5
gunzip -kf uniref50.fasta.gz
```

Build BLAST databases once. This is important for UniRef90/UniRef50 because
temporary per-job `makeblastdb` would be very slow and disk-heavy:

```bash
/data/conda_envs/hotspot_wizard/bin/makeblastdb \
  -in /data/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta \
  -dbtype prot \
  -parse_seqids \
  -out /data/databases/uniprot/blastdb/uniprot_sprot

/data/conda_envs/hotspot_wizard/bin/makeblastdb \
  -in /data/databases/uniprot/current_release/uniref/uniref90/uniref90.fasta \
  -dbtype prot \
  -parse_seqids \
  -out /data/databases/uniprot/blastdb/uniref90

/data/conda_envs/hotspot_wizard/bin/makeblastdb \
  -in /data/databases/uniprot/current_release/uniref/uniref50/uniref50.fasta \
  -dbtype prot \
  -parse_seqids \
  -out /data/databases/uniprot/blastdb/uniref50
```

Check the files:

```bash
ls -lh /data/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta
ls -lh /data/databases/uniprot/current_release/uniref/uniref90/uniref90.fasta
ls -lh /data/databases/uniprot/current_release/uniref/uniref50/uniref50.fasta
ls -lh /data/databases/uniprot/blastdb/uniprot_sprot.*
ls -lh /data/databases/uniprot/blastdb/uniref90.*
ls -lh /data/databases/uniprot/blastdb/uniref50.*
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
