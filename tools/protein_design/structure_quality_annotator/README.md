# Structure Quality Annotator

Galaxy Tool 2 for a HotSpot Wizard-like semi-rational protein engineering workflow.

The tool calls `mkdssp` to annotate each residue with:

- DSSP secondary structure
- coarse secondary-structure class
- solvent-accessible surface area
- relative ASA and exposure class
- phi/psi angles
- average PDB B-factor when available
- simple quality flags for downstream hotspot ranking

## Server Environment

Install DSSP into the shared HotSpot Wizard environment:

```bash
conda create -p /data/conda_envs/hotspot_wizard -y \
  -c conda-forge \
  python=3.11 dssp biopython pandas numpy matplotlib
```

Or add DSSP to an existing environment:

```bash
conda install -p /data/conda_envs/hotspot_wizard -y -c conda-forge dssp
```

The Galaxy job destination is configured in `config/job_conf.yml` with:

```yaml
DSSP_BINARY: /data/conda_envs/hotspot_wizard/bin/mkdssp
```

## Command-Line Smoke Test

```bash
/data/conda_envs/hotspot_wizard/bin/python \
  /data/tools/galaxy_bio/tools/protein_design/structure_quality_annotator/run_structure_quality_annotator.py \
  --input-structure /data/test/dssp_selftest/input.pdb \
  --residue-features /data/test/dssp_selftest/residue_structure_features.tsv \
  --quality-report /data/test/dssp_selftest/quality_report.html \
  --dssp-output /data/test/dssp_selftest/output.dssp \
  --dssp-binary /data/conda_envs/hotspot_wizard/bin/mkdssp
```

## Galaxy

After installing DSSP, restart Galaxy:

```bash
cd /data/tools/galaxy_bio
./run.sh --stop-daemon
./run.sh --daemon
```

The tool appears under `Protein Design` as `Structure Quality Annotator`.
