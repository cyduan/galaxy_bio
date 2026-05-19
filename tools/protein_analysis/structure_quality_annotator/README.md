# Structure Quality Annotator

Galaxy protein-analysis tool for residue-level structural annotation.

The tool calls `mkdssp` and, by default, `freesasa` to annotate each residue with:

- DSSP secondary structure
- coarse secondary-structure class
- DSSP solvent-accessible surface area
- FreeSASA total/relative SASA
- FreeSASA main-chain/side-chain and polar/apolar SASA
- relative ASA and exposure class
- phi/psi angles
- average PDB B-factor when available
- simple quality flags for downstream hotspot ranking

## Server Environment

Install DSSP and FreeSASA into the shared HotSpot Wizard environment:

```bash
conda create -p /data/conda_envs/hotspot_wizard -y \
  -c conda-forge \
  python=3.11 dssp freesasa-c biopython pandas numpy matplotlib
```

Or add them to an existing environment:

```bash
conda install -p /data/conda_envs/hotspot_wizard -y -c conda-forge dssp freesasa-c biopython
```

The Galaxy job destination is configured in `config/job_conf.yml` with:

```yaml
DSSP_BINARY: /data/conda_envs/hotspot_wizard/bin/mkdssp
FREESASA_BINARY: /data/conda_envs/hotspot_wizard/bin/freesasa
```

## Command-Line Smoke Test

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

To test FreeSASA alone on older CLI builds, use RSA output:

```bash
/data/conda_envs/hotspot_wizard/bin/freesasa \
  --format=rsa \
  /data/test/dssp_selftest/input.pdb \
  > /data/test/dssp_selftest/freesasa.rsa
```

## Galaxy

After installing DSSP/FreeSASA, restart Galaxy:

```bash
cd /data/tools/galaxy_bio
./run.sh --stop-daemon
./run.sh --daemon
```

The tool appears under `Protein Analysis` as `Structure Quality Annotator`.
FreeSASA is enabled by the `Run FreeSASA residue SASA calculation` checkbox.
