# Protein Analysis Tools

This directory contains two Galaxy tools:

- `ProtParam`: computes physicochemical properties from FASTA records.
- `Protein-Sol`: reports solubility-related sequence features and can delegate
  to a configured local Protein-Sol command.

## ProtParam install and checks

The ProtParam wrapper uses Biopython `ProteinAnalysis`.

```bash
/data/tools/galaxy_bio/.venv/bin/python -m pip install biopython
/data/tools/galaxy_bio/.venv/bin/python - <<'PY'
from Bio.SeqUtils.ProtParam import ProteinAnalysis
print(ProteinAnalysis("MKTAYIAK").molecular_weight())
PY
```

## Protein-Sol install and checks

The University of Manchester Protein-Sol service provides a downloadable
software package from:

```text
https://protein-sol.manchester.ac.uk/software
```

Because the local package format can vary, this Galaxy wrapper supports a command
template through `PROTEINSOL_COMMAND`. Create a small launcher that accepts an
input FASTA and output directory, then set the template in `config/job_conf.yml`.

Example launcher path:

```text
/data/tools/protein-sol/run_protein_sol.sh
```

Expected template:

```bash
/data/tools/protein-sol/run_protein_sol.sh {input_fasta} {output_dir}
```

Check the launcher:

```bash
mkdir -p /data/test/protein_sol_selftest
cp /data/tools/galaxy_bio/tools/structure_prediction/protein_analysis/test-data/protein_analysis_input.fasta /data/test/protein_sol_selftest/input.fasta
/data/tools/protein-sol/run_protein_sol.sh /data/test/protein_sol_selftest/input.fasta /data/test/protein_sol_selftest/out
ls -lh /data/test/protein_sol_selftest/out
```

If no local official Protein-Sol package is configured, the Galaxy tool still
returns an offline feature table. That table is useful for triage, but it is not
the official QuerySol score.
