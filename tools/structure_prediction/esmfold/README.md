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

When ESMFold is installed through an OpenFold-based environment, make sure the
runtime dependencies pulled in by OpenFold are complete as well. In particular,
if Galaxy reports ``No module named 'modelcif'`` while loading the model, install
that package into the same ESMFold environment:

```bash
/home/ubuntu/miniconda3/envs/esmfold39/bin/python -m pip install modelcif
```

If Galaxy instead reports ``No module named 'torch._six'`` while importing
DeepSpeed, the ESMFold environment has drifted away from Meta's archived ESMFold
stack. The official repository instead recommends either building the provided
``environment.yml`` or installing the pinned OpenFold commit from the README:

```bash
/home/ubuntu/miniconda3/envs/esmfold39/bin/python -m pip install \
  'openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307'
```

If Galaxy reports missing keys such as
``trunk.structure_module.ipa.linear_kv_points.linear.*`` or
``trunk.structure_module.ipa.linear_q_points.linear.*``, the installed OpenFold
implementation does not match the ESMFold weights. This is commonly caused by a
newer editable OpenFold checkout shadowing the pinned ESMFold dependency. In
that case, remove the conflicting OpenFold install and reinstall the pinned
commit above.

For the most reproducible setup, rebuilding a clean environment from the
archived Meta ``environment.yml`` is safer than mixing newer PyTorch/OpenFold
packages into an existing environment.

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
