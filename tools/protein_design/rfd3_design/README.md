# RFDiffusion3 Galaxy Tool

This wrapper runs the Foundry `rfd3 design` command from Galaxy.

Recommended server layout:

```bash
conda create -p /data/conda_envs/foundry python=3.12 -y
conda activate /data/conda_envs/foundry
python -m pip install --upgrade pip
python -m pip install "rc-foundry[all]"
foundry install base-models --checkpoint-dir /data/models/foundry/checkpoints
```

Clone the production branch as well if you want the official tutorial input
files locally:

```bash
git clone --branch production https://github.com/RosettaCommons/foundry.git /data/tools/Repo/foundry
```

Verify the environment before starting Galaxy:

```bash
/data/conda_envs/foundry/bin/python - <<'PY'
import importlib.util, sys
print("python:", sys.executable)
print("rfd3:", importlib.util.find_spec("rfd3"))
print("mpnn:", importlib.util.find_spec("mpnn"))
PY
/data/conda_envs/foundry/bin/rfd3 --help
```

For a minimal real RFD3 smoke test, use a tiny design and one batch after the
checkpoint install has completed:

```bash
mkdir -p /data/test/rfd3_selftest
cat > /data/test/rfd3_selftest/input.json <<'EOF'
{
  "galaxy_smoke_test": {
    "length": "30"
  }
}
EOF
TORCH_HOME=/data/cache/torch \
FOUNDRY_CHECKPOINT_DIRS=/data/models/foundry/checkpoints \
/data/conda_envs/foundry/bin/rfd3 design \
  out_dir=/data/test/rfd3_selftest/out \
  inputs=/data/test/rfd3_selftest/input.json \
  n_batches=1 diffusion_batch_size=1 num_timesteps=1
```

The Galaxy wrapper collects generated `.cif.gz`, `.cif`, and `.pdb` files into a
structure collection and preserves the full Foundry output directory as a ZIP
archive.
