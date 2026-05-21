# Pocket and Tunnel Finder

Tool 5 for the HotSpot Wizard-like protein engineering workflow.

Current backend integrations:

- `fpocket` for geometric pockets/cavities.
- `P2Rank` for ligand-binding site prediction.
- `CAVER` for substrate/product tunnel analysis.

Workflow-facing outputs:

- `pockets.tsv`
- `tunnels.tsv`
- `residue_pocket_tunnel_annotation.tsv`
- `pocket_structure.pdb`

Optional raw outputs can be retained through Galaxy checkboxes.

## Install fpocket

```bash
conda install -p /data/conda_envs/hotspot_wizard -y -c conda-forge fpocket
/data/conda_envs/hotspot_wizard/bin/fpocket -h
```

## Install Java

P2Rank and CAVER are Java-based tools. Install OpenJDK into the shared
environment:

```bash
conda install -p /data/conda_envs/hotspot_wizard -y -c conda-forge openjdk
/data/conda_envs/hotspot_wizard/bin/java -version
```

## Install P2Rank

Download a P2Rank release under `/data/tools/p2rank`.

```bash
mkdir -p /data/tools/p2rank
cd /data/tools/p2rank

# Example. Replace version if a newer release is preferred.
wget -c https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz
tar -xzf p2rank_2.5.1.tar.gz --strip-components=1
chmod +x prank

/data/tools/p2rank/prank --help
/data/tools/p2rank/prank predict --help
```

`config/job_conf.yml` should contain:

```yaml
P2RANK_BINARY: /data/tools/p2rank/prank
```

## Test P2Rank

```bash
mkdir -p /data/test/p2rank_selftest
cp /data/tools/galaxy_bio/tools/protein_analysis/pocket_tunnel_finder/test-data/mini_pocket_test.pdb \
  /data/test/p2rank_selftest/input.pdb

/data/tools/p2rank/prank predict \
  -f /data/test/p2rank_selftest/input.pdb \
  -o /data/test/p2rank_selftest/p2rank_out \
  -threads 1

find /data/test/p2rank_selftest/p2rank_out -maxdepth 2 -type f | head
```

Expected useful files usually include `*_predictions.csv` and
`*_residues.csv`.

## Install CAVER

CAVER distributions vary. Put the extracted CAVER directory at:

```text
/data/tools/caver_3.0/caver
```

Expected important files:

```text
/data/tools/caver_3.0/caver/caver.jar
/data/tools/caver_3.0/caver/lib/
```

Example layout check:

```bash
ls -lh /data/tools/caver_3.0/caver/caver.jar
ls -lh /data/tools/caver_3.0/caver/lib | head
```

`config/job_conf.yml` should contain:

```yaml
JAVA_BINARY: /data/conda_envs/hotspot_wizard/bin/java
CAVER_HOME: /data/tools/caver_3.0/caver
CAVER_JAR: /data/tools/caver_3.0/caver/caver.jar
CAVER_JAVA_MEM: 4000m
```

If your CAVER build uses a different command line, set `CAVER_COMMAND` with
placeholders:

```yaml
CAVER_COMMAND: "{java} -Xmx{java_mem} -jar {caver_jar} -home {caver_home} -pdb {pdb_dir} -conf {config} -out {out_dir}"
```

Available placeholders:

- `{java}`
- `{java_mem}`
- `{caver_home}`
- `{caver_jar}`
- `{pdb_dir}`
- `{config}`
- `{out_dir}`

## Test CAVER Through The Wrapper

CAVER needs a starting point. The wrapper can use the geometric centroid for a
smoke test, but real enzyme engineering should use active-site or ligand-near
coordinates.

```bash
mkdir -p /data/test/caver_selftest
cp /data/tools/galaxy_bio/tools/protein_analysis/pocket_tunnel_finder/test-data/mini_pocket_test.pdb \
  /data/test/caver_selftest/input.pdb
```

Run the wrapper with P2Rank/CAVER enabled and fpocket disabled:

```bash
/data/conda_envs/hotspot_wizard/bin/python \
  /data/tools/galaxy_bio/tools/protein_analysis/pocket_tunnel_finder/run_pocket_tunnel_finder.py \
  --input-pdb /data/test/caver_selftest/input.pdb \
  --pockets-tsv /data/test/caver_selftest/pockets.tsv \
  --tunnels-tsv /data/test/caver_selftest/tunnels.tsv \
  --residue-annotation-tsv /data/test/caver_selftest/residue_pocket_tunnel_annotation.tsv \
  --pocket-structure-pdb /data/test/caver_selftest/pocket_structure.pdb \
  --fpocket-info /data/test/caver_selftest/fpocket_info.txt \
  --p2rank-predictions /data/test/caver_selftest/p2rank_predictions.csv \
  --p2rank-residues /data/test/caver_selftest/p2rank_residues.csv \
  --caver-config /data/test/caver_selftest/caver_config.txt \
  --raw-archive /data/test/caver_selftest/pocket_tunnel_raw_outputs.zip \
  --run-log /data/test/caver_selftest/run_log.txt \
  --run-p2rank \
  --run-caver \
  --p2rank-binary /data/tools/p2rank/prank \
  --java-binary /data/conda_envs/hotspot_wizard/bin/java \
  --caver-home /data/tools/caver_3.0/caver \
  --caver-jar /data/tools/caver_3.0/caver/caver.jar \
  --caver-start-mode auto_centroid
```

For meaningful CAVER runs, use manual coordinates:

```bash
--caver-start-mode manual \
--caver-start-x 10.5 \
--caver-start-y 22.1 \
--caver-start-z 8.3
```

## Galaxy

The tool appears under `Protein Analysis` as `Pocket and Tunnel Finder`.

Restart Galaxy after installing external tools or editing config:

```bash
cd /data/tools/galaxy_bio
./run.sh --stop-daemon
./run.sh --daemon
```
