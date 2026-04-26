# ESMFold Galaxy Tool

This directory contains the first Galaxy wrapper for ESMFold in this repository.

## What it does

- accepts one FASTA record per dataset with `esmfold.xml`
- accepts multi-FASTA batch jobs with `esmfold_batch.xml`
- runs the `esm-fold` CLI
- writes `pdb` structure models
- writes small `json` summaries for downstream workflows
- optionally writes a complete ZIP archive of the raw `esm-fold` output directory

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
/data/conda_envs/esmfold_official/bin/python -m pip install modelcif
```

If Galaxy instead reports ``No module named 'torch._six'`` while importing
DeepSpeed, the ESMFold environment has drifted away from Meta's archived ESMFold
stack. The official repository instead recommends either building the provided
``environment.yml`` or installing the pinned OpenFold commit from the README:

```bash
/data/conda_envs/esmfold_official/bin/python -m pip install \
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

If installing the pinned OpenFold commit fails because the server cannot reach
GitHub, use a local source checkout or a transferred source archive instead of
letting `pip` clone the repository directly. For example, if an existing
OpenFold checkout is available on the server and already contains the pinned
commit:

```bash
git clone --no-hardlinks /data/tools/Repo/openfold_1 /data/tools/Repo/openfold_esmfold_pinned
cd /data/tools/Repo/openfold_esmfold_pinned
git checkout 4b41059694619831a7db195b7e0988fc4ff3a307
/data/conda_envs/esmfold_official/bin/python -m pip install -e .
```

If the local checkout does not contain that commit, download the source archive
on a machine that can access GitHub, transfer it to the server, and install the
archive with:

```bash
/data/conda_envs/esmfold_official/bin/python -m pip install /data/tools/sources/openfold-4b41059694619831a7db195b7e0988fc4ff3a307.tar.gz
```

If that local archive install fails while importing PyTorch with
``undefined symbol: iJIT_NotifyEvent``, repair the environment's MKL packages
first and then retry the OpenFold install:

```bash
conda install -p /data/conda_envs/esmfold_official --dry-run \
  --override-channels --no-channel-priority \
  -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge \
  -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/pytorch \
  "mkl=2024.0.0"

conda install -p /data/conda_envs/esmfold_official -y \
  --override-channels --no-channel-priority \
  -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge \
  -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/pytorch \
  "mkl=2024.0.0"
/data/conda_envs/esmfold_official/bin/python -c "import torch; print(torch.__version__)"
```

If Conda still cannot resolve the environment cleanly, keep the existing CUDA
libraries in place and replace only the MKL runtime with the exact conda-forge
package:

```bash
conda install -p /data/conda_envs/esmfold_official -y --no-deps \
  https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/linux-64/mkl-2024.0.0-ha957f24_49657.conda
/data/conda_envs/esmfold_official/bin/python -c "import torch; print(torch.__version__)"
```

After installing OpenFold, verify that the ESM package includes the official
`esm-fold` CLI:

```bash
/data/conda_envs/esmfold_official/bin/python -c "import importlib.util; print(importlib.util.find_spec('esm.scripts.fold'))"
/data/conda_envs/esmfold_official/bin/esm-fold -h
```

If `esm.scripts.fold` is missing, reinstall `fair-esm` from the official source
archive or repository instead of relying on a partial wheel:

```bash
/data/conda_envs/esmfold_official/bin/python -m pip uninstall -y fair-esm esm
/data/conda_envs/esmfold_official/bin/python -m pip install /data/tools/sources/esm-main.zip
```

If `esm-fold` fails while importing `einops` with a `SyntaxError` at a function
argument containing `/`, the environment has installed a newer `einops` release
that no longer supports Python 3.7. Pin `einops` to the last compatible line:

```bash
/data/conda_envs/esmfold_official/bin/python -m pip install "einops==0.6.1"
```

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

## Single versus batch prediction

The single-sequence tool intentionally keeps the previous behavior: one FASTA
record in, one PDB model and one JSON summary out. ESMFold generally produces
one model per input sequence, not a set of ranked alternative models like
AlphaFold. If `Save complete esm-fold output archive` is enabled, Galaxy also
returns a ZIP containing the raw `esm-fold` output folder.

Use `ESMFold Batch` for multi-FASTA inputs. It runs `esm-fold` once, then
organizes outputs into:

- a collection of PDB structures
- a collection of per-record JSON summaries
- one combined JSON summary
- one plain-text run log
- an optional ZIP archive containing raw and organized outputs

The batch wrapper runs records one at a time. A failed record is written as a
failed per-record summary and logged, but later records continue. Within one
active Galaxy job, if a structure and summary already exist for a record, that
record is skipped. If the whole Galaxy job is killed and resubmitted, Galaxy
normally starts with a fresh job directory, so cross-job resume requires a
persistent external working directory rather than the default job sandbox.
