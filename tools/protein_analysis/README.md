# Protein Analysis Tools

This directory contains Galaxy tools for sequence-level and structure-level
protein analysis:

- `ProtParam`: computes physicochemical properties from FASTA records.
- `Protein-Sol`: reports solubility-related sequence features and can delegate
  to a configured local Protein-Sol command.
- `Structure Quality Annotator`: calls DSSP/mkdssp and optionally FreeSASA to
  report residue-level secondary structure, solvent exposure, and quality flags.
- `Homolog Search and MSA`: finds homologous proteins, removes redundancy, and
  builds an MSA for conservation/back-to-consensus analysis.
- `Conservation / Mutability Scorer`: calculates per-residue conservation,
  entropy, consensus amino acids, and mutability scores from an MSA.
- `Pocket and Tunnel Finder`: detects fpocket pockets/cavities and annotates
  pocket residues for hotspot ranking.
- `Hotspot Residue Ranker`: integrates structure, conservation, pocket/tunnel,
  and optional active-site evidence into ranked hotspot residues.

## ProtParam install and checks

The ProtParam wrapper uses Biopython `ProteinAnalysis`.

```bash
/data/tools/galaxy_bio/.venv/bin/python -m pip install biopython
/data/tools/galaxy_bio/.venv/bin/python - <<'PY'
from Bio.SeqUtils.ProtParam import ProteinAnalysis
print(ProteinAnalysis("MKTAYIAK").molecular_weight())
PY
```

## Protein-Sol install and checks

The University of Manchester Protein-Sol service provides a downloadable
software package from:

```text
https://protein-sol.manchester.ac.uk/software
```

Because the local package format can vary, this Galaxy wrapper supports a command
template through `PROTEINSOL_COMMAND`. Create a small launcher that accepts an
input FASTA and output directory, then set the template in `config/job_conf.yml`.

Example launcher path:

```text
/data/tools/protein-sol/run_protein_sol.sh
```

Expected template:

```bash
/data/tools/protein-sol/run_protein_sol.sh {input_fasta} {output_dir}
```

Check the launcher:

```bash
mkdir -p /data/test/protein_sol_selftest
cp /data/tools/galaxy_bio/tools/protein_analysis/test-data/protein_analysis_input.fasta /data/test/protein_sol_selftest/input.fasta
/data/tools/protein-sol/run_protein_sol.sh /data/test/protein_sol_selftest/input.fasta /data/test/protein_sol_selftest/out
ls -lh /data/test/protein_sol_selftest/out
```

If no local official Protein-Sol package is configured, the Galaxy tool still
returns an offline feature table. That table is useful for triage, but it is not
the official QuerySol score.

## Structure Quality Annotator install and checks

Install DSSP and FreeSASA into the shared HotSpot Wizard environment:

```bash
conda install -p /data/conda_envs/hotspot_wizard -y -c conda-forge dssp freesasa-c biopython
```

Check DSSP:

```bash
/data/conda_envs/hotspot_wizard/bin/mkdssp --help
```

Check FreeSASA. Some conda builds do not support `--output-depth=residue`, so
the Galaxy wrapper uses RSA output as a compatible fallback:

```bash
/data/conda_envs/hotspot_wizard/bin/freesasa \
  --format=rsa \
  /data/test/dssp_selftest/input.pdb \
  > /data/test/dssp_selftest/freesasa.rsa
```

Run the full wrapper:

```bash
/data/conda_envs/hotspot_wizard/bin/python \
  /data/tools/galaxy_bio/tools/protein_analysis/structure_quality_annotator/run_structure_quality_annotator.py \
  --input-structure /data/test/dssp_selftest/input.pdb \
  --residue-features /data/test/dssp_selftest/residue_structure_features.tsv \
  --quality-report /data/test/dssp_selftest/quality_report.html \
  --dssp-output /data/test/dssp_selftest/output.dssp \
  --dssp-binary /data/conda_envs/hotspot_wizard/bin/mkdssp \
  --run-freesasa \
  --freesasa-output /data/test/dssp_selftest/freesasa.rsa \
  --freesasa-binary /data/conda_envs/hotspot_wizard/bin/freesasa
```

## Homolog Search and MSA install and checks

Install the homology-search and MSA tools into the shared HotSpot Wizard
environment:

```bash
conda install -p /data/conda_envs/hotspot_wizard -y \
  -c conda-forge -c bioconda \
  blast cd-hit mafft muscle mmseqs2
```

Check the binaries:

```bash
/data/conda_envs/hotspot_wizard/bin/blastp -version
/data/conda_envs/hotspot_wizard/bin/makeblastdb -version
/data/conda_envs/hotspot_wizard/bin/cd-hit -h
/data/conda_envs/hotspot_wizard/bin/mafft --version
/data/conda_envs/hotspot_wizard/bin/muscle -version
/data/conda_envs/hotspot_wizard/bin/mmseqs version
```

The Homolog Search and MSA wrapper can use uploaded FASTA datasets or server
databases configured in `config/job_conf.yml`:

```text
/data/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta
/data/databases/uniprot/current_release/uniref/uniref90/uniref90.fasta
/data/databases/uniprot/current_release/uniref/uniref50/uniref50.fasta
/data/databases/uniprot/blastdb/
```

See `tools/protein_analysis/homolog_search_msa/README.md` for download and
`makeblastdb` commands.

## Conservation / Mutability Scorer install and checks

Install Rate4Site either system-wide or in the shared HotSpot Wizard
environment:

```bash
sudo apt-get update
sudo apt-get install -y rate4site
rate4site -h
```

If the binary is installed elsewhere, update `RATE4SITE_BINARY` in
`config/job_conf.yml`.

## Pocket and Tunnel Finder install and checks

Install fpocket, Java, and configure P2Rank/CAVER:

```bash
conda install -p /data/conda_envs/hotspot_wizard -y -c conda-forge fpocket
/data/conda_envs/hotspot_wizard/bin/fpocket -h

conda install -p /data/conda_envs/hotspot_wizard -y -c conda-forge openjdk
/data/conda_envs/hotspot_wizard/bin/java -version
```

P2Rank is expected at `/data/tools/p2rank/prank`, and CAVER at
`/data/tools/caver_3.0/caver/caver.jar`. If installed elsewhere, update
`FPOCKET_BINARY`, `P2RANK_BINARY`, `JAVA_BINARY`, `CAVER_HOME`, and `CAVER_JAR`
in `config/job_conf.yml`.
