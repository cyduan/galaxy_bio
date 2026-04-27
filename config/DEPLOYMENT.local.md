# Galaxy Local Deployment Notes

This repository now includes a deployable local configuration:

- `config/galaxy.yml`: Galaxy web UI listening on `0.0.0.0:2718`
- `config/job_conf.yml`: local runner plus Docker-backed Interactive Tools
- `config/tool_conf.interactivetools.xml`: enables a focused set of built-in Interactive Tools

## Start Galaxy

```bash
./run.sh
```

Then open:

```text
http://<current-host-ip>:2718
```

## Interactive Tools

The current default listens on `0.0.0.0:2718` and sets `galaxy_infrastructure_url` to `http://${HOST_IP}:2718`, so the repository no longer hard-codes a specific IP address.

For stable access from other computers, replace `${HOST_IP}` and `galaxy.example.com` in `config/galaxy.yml` with a DNS name or DDNS name that keeps following your server even when its IP changes. If Docker-backed Interactive Tools need a host alias inside containers, mirror that same hostname in `config/job_conf.yml` by uncommenting the matching `docker_run_extra_arguments` example.

## ESMFold

The repository now includes a first structure-prediction toolbox entry:

- `config/tool_conf.structure_prediction.xml`
- `tools/structure_prediction/esmfold/esmfold.xml`

By default, `config/job_conf.yml` routes the `esmfold` tool to `esmfold_local`, which expects the real `esm-fold` executable to be available in the Galaxy job environment.

To enable real inference:

1. Install the ESMFold runtime where Galaxy jobs can see it, or set `ESMFOLD_BINARY` in `config/job_conf.yml` to your own launcher path.
2. Restart Galaxy.
3. Upload a single-record FASTA and run `Structure Prediction > ESMFold`.

If the job reaches `Loading model` and then fails with `No module named 'modelcif'`,
the ESMFold/OpenFold runtime is still missing a Python dependency. Install it into
the same environment that provides `esm-fold`, for example:

```bash
/data/conda_envs/esmfold_official/bin/python -m pip install modelcif
```

If the traceback instead ends with `No module named 'torch._six'`, the current
DeepSpeed package in the ESMFold environment is mismatched with the rest of the
stack. For the official ESMFold CLI, prefer Meta's archived environment pins or
the pinned OpenFold commit from the README instead of mixing in newer OpenFold
and PyTorch packages:

```bash
/data/conda_envs/esmfold_official/bin/python -m pip install \
  'openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307'
```

If the traceback instead reports missing keys like
`trunk.structure_module.ipa.linear_kv_points.linear.*` and
`trunk.structure_module.ipa.linear_q_points.linear.*`, the ESMFold weights do
not match the installed OpenFold code. This usually means a newer editable
OpenFold checkout is shadowing the pinned ESMFold dependency. Remove the
conflicting OpenFold install and rebuild the ESMFold environment against the
official pinned commit.

If the pinned OpenFold install fails while cloning from GitHub, the server is
usually hitting a network timeout rather than an ESMFold package error. Reuse a
local checkout if it already contains the pinned commit:

```bash
git clone --no-hardlinks /data/tools/Repo/openfold_1 /data/tools/Repo/openfold_esmfold_pinned
cd /data/tools/Repo/openfold_esmfold_pinned
git checkout 4b41059694619831a7db195b7e0988fc4ff3a307
/data/conda_envs/esmfold_official/bin/python -m pip install -e .
```

If that commit is not available locally, download the OpenFold source archive on
a machine with GitHub access, transfer it to the server, and install the archive
with:

```bash
/data/conda_envs/esmfold_official/bin/python -m pip install /data/tools/sources/openfold-4b41059694619831a7db195b7e0988fc4ff3a307.tar.gz
```

If the archive install fails during `import torch` with
`undefined symbol: iJIT_NotifyEvent`, downgrade the MKL runtime in the ESMFold
environment and retry:

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

If the solver still cannot resolve the environment cleanly, keep the existing
CUDA libraries and only replace the MKL runtime with the exact conda-forge
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
argument containing `/`, pin `einops` to a Python-3.7 compatible release:

```bash
/data/conda_envs/esmfold_official/bin/python -m pip install "einops==0.6.1"
```

If you later move ESMFold into a GPU container, uncomment the `esmfold_gpu` example in `config/job_conf.yml`, set a real image name, and remap the `esmfold` tool to that environment.

The ESMFold single-sequence tool keeps the original output behavior: one input
FASTA record produces one PDB structure and one JSON summary. Enable `Save
complete esm-fold output archive` only when you also want the raw `esm-fold`
output directory as a ZIP. For multi-FASTA jobs, use `ESMFold Batch`; it produces
a PDB collection, a per-record summary JSON collection, one combined summary, a
plain-text run log, and an optional complete output archive. Batch records are
processed one at a time, so a failed sequence is reported in the summaries/log
and later records continue. Full resume after a killed Galaxy job requires a
persistent external working directory; the default Galaxy job directory is fresh
on resubmission.

## Foundry RFDiffusion3 and ProteinMPNN

The repository includes Galaxy wrappers for Foundry RFDiffusion3 and
ProteinMPNN:

- `tools/protein_design/rfd3_design/rfd3_design.xml`
- `tools/protein_design/foundry_mpnn/foundry_mpnn.xml`

Both tools are routed to `foundry_local` in `config/job_conf.yml`, which expects:

```text
/data/conda_envs/foundry/bin/rfd3
/data/conda_envs/foundry/bin/mpnn
/data/models/foundry/checkpoints
```

Recommended install:

```bash
conda create -p /data/conda_envs/foundry python=3.12 -y
conda activate /data/conda_envs/foundry
python -m pip install --upgrade pip
python -m pip install "rc-foundry[all]"
foundry install base-models --checkpoint-dir /data/models/foundry/checkpoints
```

Clone the production branch separately if you want the official tutorial inputs
and example files:

```bash
git clone --branch production https://github.com/RosettaCommons/foundry.git /data/tools/Repo/foundry
```

Environment checks:

```bash
/data/conda_envs/foundry/bin/python - <<'PY'
import importlib.util, sys
print("python:", sys.executable)
print("rfd3:", importlib.util.find_spec("rfd3"))
print("mpnn:", importlib.util.find_spec("mpnn"))
PY
/data/conda_envs/foundry/bin/rfd3 --help
/data/conda_envs/foundry/bin/mpnn --help
```

Minimal RFD3 smoke test after checkpoints are installed:

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

ProteinMPNN should be checked with `mpnn --help` first. Its wrapper exposes a
CLI argument template because the upstream Foundry README still marks detailed
command-line and JSON inference documentation as forthcoming.

## FoldX

The repository includes FoldX wrappers for stability and mutation-energy
calculations:

- `tools/protein_stability/foldx/foldx_stability.xml`
- `tools/protein_stability/foldx/foldx_buildmodel.xml`

FoldX is licensed software. Download the Linux executable from the official
FoldX Suite site with an authorized account and place it outside the repository,
for example:

```bash
mkdir -p /data/tools/foldx
# Put the authorized FoldX executable at /data/tools/foldx/foldx
chmod +x /data/tools/foldx/foldx
```

The Galaxy job destination `foldx_local` in `config/job_conf.yml` expects:

```text
/data/tools/foldx/foldx
```

Environment checks:

```bash
/data/tools/foldx/foldx --help
/data/tools/foldx/foldx --version
```

Minimal FoldX checks:

```bash
mkdir -p /data/test/foldx_selftest
cp /data/tools/galaxy_bio/tools/protein_stability/foldx/test-data/foldx_input.pdb /data/test/foldx_selftest/input.pdb
cd /data/test/foldx_selftest
/data/tools/foldx/foldx --command=RepairPDB --pdb=input.pdb
/data/tools/foldx/foldx --command=Stability --pdb=input_Repair.pdb --output-file=stability
cat > individual_list.txt <<'EOF'
GA1V;
EOF
/data/tools/foldx/foldx --command=BuildModel --pdb=input_Repair.pdb --mutant-file=individual_list.txt --numberOfRuns=5 --output-file=buildmodel
ls -lh
```

In Galaxy:

- `FoldX Stability` optionally runs `RepairPDB` and then `Stability`; the main
  output is converted to CSV with descriptive FoldX energy-term headers.
- `FoldX BuildModel` optionally runs `RepairPDB`, accepts typed mutations or an
  `individual_list.txt`, then returns mutation-energy CSV, mutant
  structures, logs, summary JSON, and a full ZIP archive.

## ProtParam and Protein-Sol

The repository includes sequence-level protein analysis tools:

- `tools/protein_analysis/protparam.xml`
- `tools/protein_analysis/protein_sol.xml`

### ProtParam

ProtParam uses Biopython's `ProteinAnalysis` implementation. Check that
Biopython is available in the Galaxy Python environment:

```bash
/data/tools/galaxy_bio/.venv/bin/python - <<'PY'
from Bio.SeqUtils.ProtParam import ProteinAnalysis
print(ProteinAnalysis("MKTAYIAK").molecular_weight())
PY
```

If that fails:

```bash
/data/tools/galaxy_bio/.venv/bin/python -m pip install biopython
```

### Protein-Sol

The University of Manchester Protein-Sol service provides downloadable software
from:

```text
https://protein-sol.manchester.ac.uk/software
```

Because the local package invocation can vary, this Galaxy wrapper supports a
command template through `PROTEINSOL_COMMAND` in `config/job_conf.yml`.

Expected launcher shape:

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

The Galaxy `Protein-Sol` tool also has a default offline sequence-feature mode.
That mode is useful for screening and workflow wiring, but it is not the
official QuerySol score. Use external mode with the configured local package for
official Protein-Sol predictions.

## BLAST

The repository now includes a first local BLAST+ integration:

- `tools/similarity_search/blast_search/blast_search.xml`
- `tools/similarity_search/blast_search/blast_search.py`

It currently exposes:

- `blastn`
- `blastp`

The wrapper accepts query and subject FASTA datasets from Galaxy history,
builds a temporary local database with `makeblastdb`, and then runs BLAST.

### Recommended environment install

NCBI's official BLAST+ manual documents the supported installer formats
(installer, RPM, tarball, or source build). In this repository's Conda-based
server layout, a dedicated Conda environment is the simplest operational
choice:

```bash
conda create -n blastenv -c conda-forge -c bioconda blast=2.16.0 -y
conda activate blastenv
makeblastdb -version
blastn -version
blastp -version
```

Then keep the default paths in `config/job_conf.yml`:

```yaml
    blast_local:
      runner: local
      env:
        - name: MAKEBLASTDB_BINARY
          value: /home/ubuntu/miniconda3/envs/blastenv/bin/makeblastdb
        - name: BLASTN_BINARY
          value: /home/ubuntu/miniconda3/envs/blastenv/bin/blastn
        - name: BLASTP_BINARY
          value: /home/ubuntu/miniconda3/envs/blastenv/bin/blastp
```

### How to test the environment

1. Confirm the binaries are present:

```bash
/home/ubuntu/miniconda3/envs/blastenv/bin/makeblastdb -version
/home/ubuntu/miniconda3/envs/blastenv/bin/blastn -version
/home/ubuntu/miniconda3/envs/blastenv/bin/blastp -version
```

2. Run a tiny nucleotide self-test:

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
/home/ubuntu/miniconda3/envs/blastenv/bin/makeblastdb -in /tmp/blast_selftest/subject_nt.fa -dbtype nucl -out /tmp/blast_selftest/subject_nt_db
/home/ubuntu/miniconda3/envs/blastenv/bin/blastn -query /tmp/blast_selftest/query_nt.fa -db /tmp/blast_selftest/subject_nt_db -outfmt 6
```

3. Run a tiny protein self-test:

```bash
cat > /tmp/blast_selftest/query_aa.fa <<'EOF'
>q1
MKTAYIAK
EOF
cat > /tmp/blast_selftest/subject_aa.fa <<'EOF'
>s1
MKTAYIAK
EOF
/home/ubuntu/miniconda3/envs/blastenv/bin/makeblastdb -in /tmp/blast_selftest/subject_aa.fa -dbtype prot -out /tmp/blast_selftest/subject_aa_db
/home/ubuntu/miniconda3/envs/blastenv/bin/blastp -query /tmp/blast_selftest/query_aa.fa -db /tmp/blast_selftest/subject_aa_db -outfmt 6
```

If these commands print version strings and tabular alignments, the BLAST
environment is ready for Galaxy.

## Reverse Proxy

This configuration uses direct Interactive Tools proxy mode:

- Galaxy UI: `http://${HOST_IP}:2718`
- gx-it-proxy host placeholder: `galaxy.example.com:4002`

If you later place Galaxy behind nginx or another reverse proxy, switch `interactivetools_upstream_proxy` in `config/galaxy.yml` to `true` and route wildcard/subpath Interactive Tools traffic through the upstream proxy instead of exposing `4002` directly.

## Notes

- `${HOST_IP}` is resolved by Galaxy at startup from the current host network configuration. It is convenient for lab or single-host setups, but a DDNS name or FQDN is still the better long-term public address.
- `interactivetools_proxy_host` does not auto-expand `${HOST_IP}` in the same way, so set it explicitly to your stable host or DDNS name before relying on Interactive Tools from other machines.
- For a longer-lived deployment, prefer a domain name plus nginx/Caddy and HTTPS over exposing a raw IP and port directly.
