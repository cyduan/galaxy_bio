# BLAST Search Galaxy Tool

This directory contains the first BLAST+ integration for this repository.

## Included programs

- `blastn`
- `blastp`
- `makeblastdb` (used internally to build a temporary subject database)

## Runtime expectation

The wrapper expects NCBI BLAST+ command line tools to be available in a
dedicated environment. In this repository, the recommended setup is:

- `/home/ubuntu/miniconda3/envs/blastenv/bin/makeblastdb`
- `/home/ubuntu/miniconda3/envs/blastenv/bin/blastn`
- `/home/ubuntu/miniconda3/envs/blastenv/bin/blastp`

These paths are wired into `config/job_conf.yml`.

## Recommended installation

For this server layout, a dedicated Conda environment keeps BLAST isolated from
Galaxy and ESMFold:

```bash
conda create -n blastenv -c conda-forge -c bioconda blast=2.16.0 -y
conda activate blastenv
makeblastdb -version
blastn -version
blastp -version
```

## Quick functional test

```bash
mkdir -p /tmp/blast_selftest
cat > /tmp/blast_selftest/query_nt.fa <<'EOF'
>q1
ACGTACGT
EOF
cat > /tmp/blast_selftest/subject_nt.fa <<'EOF'
>s1
ACGTACGT
EOF
makeblastdb -in /tmp/blast_selftest/subject_nt.fa -dbtype nucl -out /tmp/blast_selftest/subject_nt_db
blastn -query /tmp/blast_selftest/query_nt.fa -db /tmp/blast_selftest/subject_nt_db -outfmt 6
```

The wrapper itself can be smoke-tested without a real BLAST installation by
using the bundled mock commands through Galaxy tool tests.
