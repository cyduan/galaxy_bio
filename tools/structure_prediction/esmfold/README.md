# ESMFold Galaxy Tool

This directory contains the first Galaxy wrapper for ESMFold in this repository.

## What it does

- accepts one FASTA record per dataset
- runs the `esm-fold` CLI
- writes a `pdb` structure model
- writes a small `json` summary for downstream workflows

## Runtime expectation

The wrapper currently assumes the job environment already provides the `esm-fold` executable described in the official Meta ESM repository:

- https://github.com/facebookresearch/esm

That keeps the wrapper neutral with respect to how you want to deploy it:

- directly in the Galaxy runner environment
- inside a custom Docker image
- through a GPU-specific Slurm or Kubernetes destination later

## Local verification

You can verify the wrapper logic without a real ESMFold installation by using the built-in mock mode:

```bash
python tools/structure_prediction/esmfold/run_esmfold.py \
  --input-fasta tools/structure_prediction/esmfold/test-data/esmfold_single.fasta \
  --output-pdb result.pdb \
  --summary-json summary.json \
  --esmfold-command mock-esmfold
```

For a real deployment, keep the tool default at `esm-fold` and make that executable available in the Galaxy job environment, or export `ESMFOLD_BINARY` to point at a custom wrapper binary.
