# Pocket and Tunnel Finder

Tool 5 for the HotSpot Wizard-like protein engineering workflow.

Current implementation: fpocket only.

It writes the HotSpot-facing outputs:

- `pockets.tsv`
- `tunnels.tsv` placeholder
- `residue_pocket_tunnel_annotation.tsv`
- `pocket_structure.pdb`

It can also preserve raw fpocket files as a ZIP archive.

## Install fpocket

Preferred conda installation:

```bash
conda install -p /data/conda_envs/hotspot_wizard -y -c conda-forge fpocket
/data/conda_envs/hotspot_wizard/bin/fpocket -h
```

Alternative source installation:

```bash
cd /data/tools
git clone https://github.com/Discngine/fpocket.git
cd fpocket
make
./bin/fpocket -h
```

If using the source build, update `FPOCKET_BINARY` in `config/job_conf.yml`.

## Smoke Test

```bash
mkdir -p /data/test/fpocket_selftest
cp /data/tools/galaxy_bio/tools/protein_analysis/pocket_tunnel_finder/test-data/mini_pocket_test.pdb \
  /data/test/fpocket_selftest/input.pdb

cd /data/test/fpocket_selftest
/data/conda_envs/hotspot_wizard/bin/fpocket -f input.pdb
find input_out -maxdepth 2 -type f | head
```

Wrapper test:

```bash
/data/conda_envs/hotspot_wizard/bin/python \
  /data/tools/galaxy_bio/tools/protein_analysis/pocket_tunnel_finder/run_pocket_tunnel_finder.py \
  --input-pdb /data/test/fpocket_selftest/input.pdb \
  --pockets-tsv /data/test/fpocket_selftest/pockets.tsv \
  --tunnels-tsv /data/test/fpocket_selftest/tunnels.tsv \
  --residue-annotation-tsv /data/test/fpocket_selftest/residue_pocket_tunnel_annotation.tsv \
  --pocket-structure-pdb /data/test/fpocket_selftest/pocket_structure.pdb \
  --fpocket-info /data/test/fpocket_selftest/fpocket_info.txt \
  --raw-archive /data/test/fpocket_selftest/fpocket_raw_outputs.zip \
  --run-log /data/test/fpocket_selftest/run_log.txt \
  --fpocket-binary /data/conda_envs/hotspot_wizard/bin/fpocket
```

## Galaxy

The tool appears under `Protein Analysis` as `Pocket and Tunnel Finder`.

Restart Galaxy after installing fpocket or editing config:

```bash
cd /data/tools/galaxy_bio
./run.sh --stop-daemon
./run.sh --daemon
```
