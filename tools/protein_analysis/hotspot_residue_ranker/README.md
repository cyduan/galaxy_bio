# Hotspot Residue Ranker

Tool 6 for the HotSpot Wizard-like protein engineering workflow.

It integrates:

- `residue_structure_features.tsv`
- `residue_conservation.tsv`
- `pockets.tsv`
- `tunnels.tsv`
- optional `optional_active_site.tsv`

Outputs:

- `hotspot_residue_ranked.tsv`
- `hotspot_report.html`
- `hotspot_viewer.html`

## Active-site TSV format

Recommended columns:

```text
site_id	chain	residue_number	residue_name	role
cat1	A	2	VAL	catalytic
```

Optional coordinate columns are accepted for forward compatibility:

```text
x	y	z
```

The current ranker applies exact residue penalties when active-site residue IDs
match. Future versions can use residue centroids for geometric distance scoring.

## Smoke Test

```bash
python /data/tools/galaxy_bio/tools/protein_analysis/hotspot_residue_ranker/run_hotspot_residue_ranker.py \
  --structure-features /data/tools/galaxy_bio/tools/protein_analysis/hotspot_residue_ranker/test-data/residue_structure_features.tsv \
  --conservation-tsv /data/tools/galaxy_bio/tools/protein_analysis/hotspot_residue_ranker/test-data/residue_conservation.tsv \
  --pockets-tsv /data/tools/galaxy_bio/tools/protein_analysis/hotspot_residue_ranker/test-data/pockets.tsv \
  --tunnels-tsv /data/tools/galaxy_bio/tools/protein_analysis/hotspot_residue_ranker/test-data/tunnels.tsv \
  --active-site-tsv /data/tools/galaxy_bio/tools/protein_analysis/hotspot_residue_ranker/test-data/optional_active_site.tsv \
  --output-tsv /tmp/hotspot_residue_ranked.tsv \
  --report-html /tmp/hotspot_report.html \
  --viewer-html /tmp/hotspot_viewer.html \
  --run-log /tmp/hotspot_ranker_run_log.txt
```

## Galaxy

The tool appears under `Protein Analysis` as `Hotspot Residue Ranker`.

Restart Galaxy after editing config:

```bash
cd /data/tools/galaxy_bio
./run.sh --stop-daemon
./run.sh --daemon
```
