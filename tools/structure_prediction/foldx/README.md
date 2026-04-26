# FoldX Galaxy Tools

This directory contains Galaxy wrappers for FoldX:

- `FoldX Stability`: optional `RepairPDB`, then `Stability`, with the primary
  output converted to CSV using descriptive FoldX energy-term headers
- `FoldX BuildModel`: optional `RepairPDB`, then mutation `BuildModel`, with the
  primary energy table converted to CSV

FoldX is licensed software and must be downloaded from the FoldX Suite site by
an authorized user. Do not commit the FoldX executable or license materials into
this repository.

Recommended server layout:

```bash
mkdir -p /data/tools/foldx
# Place the authorized FoldX executable in /data/tools/foldx/foldx
chmod +x /data/tools/foldx/foldx
```

Verify the install:

```bash
/data/tools/foldx/foldx --help
/data/tools/foldx/foldx --version
```

Minimal CLI checks:

```bash
mkdir -p /data/test/foldx_selftest
cp /data/tools/galaxy_bio/tools/structure_prediction/foldx/test-data/foldx_input.pdb /data/test/foldx_selftest/input.pdb
cd /data/test/foldx_selftest
/data/tools/foldx/foldx --command=RepairPDB --pdb=input.pdb
/data/tools/foldx/foldx --command=Stability --pdb=input_Repair.pdb --output-file=stability
cat > individual_list.txt <<'EOF'
GA1V;
EOF
/data/tools/foldx/foldx --command=BuildModel --pdb=input_Repair.pdb --mutant-file=individual_list.txt --numberOfRuns=5 --output-file=buildmodel
ls -lh
```

Then set `FOLDX_BINARY` in `config/job_conf.yml`:

```yaml
    foldx_local:
      runner: local
      env:
        - name: FOLDX_BINARY
          value: /data/tools/foldx/foldx
```
