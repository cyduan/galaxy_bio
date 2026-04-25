# Foundry ProteinMPNN Galaxy Tool

This wrapper delegates to the Foundry `mpnn` command and collects FASTA files
plus the full output directory.

Recommended server layout is shared with RFDiffusion3:

```bash
conda create -p /data/conda_envs/foundry python=3.12 -y
conda activate /data/conda_envs/foundry
python -m pip install --upgrade pip
python -m pip install "rc-foundry[all]"
foundry install base-models --checkpoint-dir /data/models/foundry/checkpoints
```

Clone the production branch if you want the official examples and tutorial
inputs locally:

```bash
git clone --branch production https://github.com/RosettaCommons/foundry.git /data/tools/Repo/foundry
```

Verify the environment:

```bash
/data/conda_envs/foundry/bin/python - <<'PY'
import importlib.util, sys
print("python:", sys.executable)
print("mpnn:", importlib.util.find_spec("mpnn"))
PY
/data/conda_envs/foundry/bin/mpnn --help
```

The Foundry ProteinMPNN README currently marks detailed command-line and JSON
inference documentation as forthcoming, so the Galaxy tool exposes a CLI
argument template instead of hard-coding unstable option names. Use
`__STRUCTURE__`, `__OUTPUT_DIR__`, `__CONFIG__`, and `__WEIGHTS__` placeholders
inside the template after checking `mpnn --help` in the installed environment.
