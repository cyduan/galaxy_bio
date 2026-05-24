# Protein Analysis Tools

This directory contains Galaxy tools for sequence-level and structure-level
protein analysis:

- `ProtParam`: computes physicochemical properties from FASTA records.
- `Protein-Sol`: reports solubility-related sequence features and can delegate
  to a configured local Protein-Sol command.
- `Structure Quality Annotator`: calls DSSP/mkdssp and optionally FreeSASA to
  report residue-level secondary structure, solvent exposure, and quality flags.
- `Functional Residue Annotator`: marks catalytic, ligand-contact,
  metal-binding, essential, and user/database-provided functional residues.
- `Residue Numbering Validator`: checks that PDB numbering, target FASTA
  positions, MSA columns, structure features, and conservation rows are aligned.
- `HotSpot Workflow Preset Advisor`: chooses objective-specific parameters for
  the structure-to-smart-library workflow without hiding intermediate tools.
- `Homolog Search and MSA`: finds homologous proteins, removes redundancy, and
  builds an MSA for conservation/back-to-consensus analysis.
- `Conservation / Mutability Scorer`: calculates per-residue conservation,
  entropy, consensus amino acids, and mutability scores from an MSA.
- `Pocket and Tunnel Finder`: detects fpocket pockets/cavities and annotates
  pocket residues for hotspot ranking.
- `Hotspot Residue Ranker`: integrates structure, conservation, pocket/tunnel,
  and optional active-site evidence into ranked hotspot residues.
- `Mutation Designer`: recommends hotspot substitutions from MSA variation,
  consensus residues, conservative replacements, alanine scanning, and reduced
  amino-acid alphabets; also writes smart-library and degenerate-codon tables.
- `Stability Predictor`: converts candidate mutations into FoldX
  `individual_list.txt`, runs RepairPDB/BuildModel, and formats ΔΔG results for
  HotSpot filtering while optionally preserving raw FoldX outputs.
- `Smart Library Report Generator`: merges ranked hotspots, recommended
  mutations, FoldX ΔΔG results, smart-library design, and optional structure
  context into final TSV/HTML reports for experimental planning.

## HotSpot Structure-to-Smart-Library Workflow

The recommended user-facing workflow is:

```text
tools/protein_analysis/workflows/hotspot_structure_to_smart_library.ga
```

It connects:

```text
Protein structure PDB + protein FASTA
  -> Structure Quality Annotator
  -> Functional Residue Annotator
  -> Residue Numbering Validator
  -> Homolog Search and MSA
  -> Conservation / Mutability Scorer
  -> Pocket and Tunnel Finder
  -> Hotspot Residue Ranker
  -> Mutation Designer
  -> Stability Predictor
  -> Smart Library Report Generator
```

Use `HotSpot Workflow Preset Advisor` before running the workflow. It produces:

```text
recommended_parameters.tsv
workflow_preset.json
workflow_guide.html
```

The advisor does not replace the workflow. It gives objective-specific defaults
for activity/substrate specificity, stability improvement, conservative small
libraries, exploratory large libraries, or balanced design.

Import the workflow in Galaxy:

1. Go to `Workflow` -> `Import`.
2. Upload `tools/protein_analysis/workflows/hotspot_structure_to_smart_library.ga`.
3. Open the imported workflow editor.
4. Apply values from `recommended_parameters.tsv` to the corresponding steps.
5. Run with one PDB structure dataset and one protein FASTA dataset.

To make it a public workflow inside a running Galaxy instance, import it as an
admin user and publish/share it from the Galaxy workflow menu. Repository files
alone do not automatically create database-backed public workflows in Galaxy.

### HotSpot Tool Input/Output Quick Reference

| Step | Galaxy tool | Main input type | Main input file | Main output files | Downstream use |
| --- | --- | --- | --- | --- | --- |
| Advisor | HotSpot Workflow Preset Advisor | form parameters only | no dataset required | `recommended_parameters.tsv`, `workflow_preset.json`, `workflow_guide.html` | choose objective-specific workflow parameters |
| Tool 2 | Structure Quality Annotator | `pdb`, `cif`, `mmcif` | protein structure with atom coordinates | `residue_structure_features.tsv`, `quality_report.html` | structure features for hotspot ranking |
| Tool 2b | Functional Residue Annotator | `pdb`, optional `tabular` | protein PDB plus optional curated functional residue TSV | `residue_functional_annotation.tsv`, `active_site_candidates.tsv`, `ligand_contacts.tsv`, `functional_annotation_report.html` | protects catalytic/essential residues and feeds active-site penalties to Tool 6 |
| QC | Residue Numbering Validator | `pdb`, `fasta`, optional aligned `fasta`, optional `tabular` | PDB, target FASTA, MSA, Tool 2/Tool 4 TSVs | `residue_numbering_validation.tsv`, `numbering_summary.tsv`, `numbering_report.html` | detects chain/numbering/mapping problems before trusting hotspot ranking |
| Tool 3 | Homolog Search and MSA | `fasta` | target protein FASTA plus server/uploaded homolog database | `homologs.fasta`, `filtered_homologs.fasta`, `msa.fasta`, `homolog_summary.tsv` | MSA for conservation and mutation design |
| Tool 4 | Conservation / Mutability Scorer | aligned `fasta` | `msa.fasta` from Tool 3 | `residue_conservation.tsv`, raw Rate4Site output, log | conservation, mutability, consensus, accepted AAs |
| Tool 5 | Pocket and Tunnel Finder | `pdb` | protein structure PDB | `pockets.tsv`, `tunnels.tsv`, `residue_pocket_tunnel_annotation.tsv`, `pocket_structure.pdb` | functional-region evidence for hotspot ranking |
| Tool 6 | Hotspot Residue Ranker | `tabular` plus optional `pdb` | Tool 2, Tool 4, Tool 5 TSVs | `hotspot_residue_ranked.tsv`, `hotspot_report.html`, `hotspot_viewer.html` | ranked hotspot residues for mutation design |
| Tool 7 | Mutation Designer | `tabular`, aligned `fasta`, sequence `fasta` | Tool 6 hotspots, Tool 3 MSA, target protein FASTA | `candidate_mutations.tsv`, `smart_library.tsv`, `degenerate_codons.tsv` | candidate substitutions and codon/library design |
| Tool 8 | Stability Predictor | `tabular`, `pdb` | Tool 7 candidates plus FoldX-compatible PDB | `mutation_ddg.tsv`, `stability_filtered_mutations.tsv`, optional FoldX raw outputs/models | FoldX stability filtering |
| Tool 9 | Smart Library Report Generator | `tabular` plus optional `pdb` | Tool 6/7/8 TSVs plus optional structure | `ranked_hotspots.tsv`, `ranked_mutations.tsv`, `smart_library_design.tsv`, `summary.html`, `structure_viewer.html` | final experimental design report |

For routine users, the most important final files are:

```text
summary.html
ranked_mutations.tsv
smart_library_design.tsv
structure_viewer.html
```

Intermediate files are intentionally preserved for auditability and expert
debugging, but they do not all need to be shown to experimental collaborators.

## History Output Naming Convention

HotSpot workflow tools use a consistent Galaxy History label pattern:

```text
Tool Name on input-dataset | 01 output_file.ext
```

For example:

```text
Hotspot Residue Ranker on dataset 20 | 01 hotspot_residue_ranked.tsv
Hotspot Residue Ranker on dataset 20 | 02 hotspot_report.html
Hotspot Residue Ranker on dataset 20 | 03 hotspot_viewer.html
```

This does not modify Galaxy's internal History ID numbering. Galaxy History IDs
are global database IDs and should not be repurposed as per-tool batch numbers.
The numbered suffix keeps outputs from the same tool run easy to search, sort,
and review without touching Galaxy core code.

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
cp /data/tools/galaxy_bio/tools/protein_analysis/test-data/protein_analysis_input.fasta /data/test/protein_sol_selftest/input.fasta
/data/tools/protein-sol/run_protein_sol.sh /data/test/protein_sol_selftest/input.fasta /data/test/protein_sol_selftest/out
ls -lh /data/test/protein_sol_selftest/out
```

If no local official Protein-Sol package is configured, the Galaxy tool still
returns an offline feature table. That table is useful for triage, but it is not
the official QuerySol score.

## Structure Quality Annotator install and checks

Install DSSP and FreeSASA into the shared HotSpot Wizard environment:

```bash
conda install -p /data/conda_envs/hotspot_wizard -y -c conda-forge dssp freesasa-c biopython
```

Check DSSP:

```bash
/data/conda_envs/hotspot_wizard/bin/mkdssp --help
```

Check FreeSASA. Some conda builds do not support `--output-depth=residue`, so
the Galaxy wrapper uses RSA output as a compatible fallback:

```bash
/data/conda_envs/hotspot_wizard/bin/freesasa \
  --format=rsa \
  /data/test/dssp_selftest/input.pdb \
  > /data/test/dssp_selftest/freesasa.rsa
```

Run the full wrapper:

```bash
/data/conda_envs/hotspot_wizard/bin/python \
  /data/tools/galaxy_bio/tools/protein_analysis/structure_quality_annotator/run_structure_quality_annotator.py \
  --input-structure /data/test/dssp_selftest/input.pdb \
  --residue-features /data/test/dssp_selftest/residue_structure_features.tsv \
  --quality-report /data/test/dssp_selftest/quality_report.html \
  --dssp-output /data/test/dssp_selftest/output.dssp \
  --dssp-binary /data/conda_envs/hotspot_wizard/bin/mkdssp \
  --run-freesasa \
  --freesasa-output /data/test/dssp_selftest/freesasa.rsa \
  --freesasa-binary /data/conda_envs/hotspot_wizard/bin/freesasa
```

## Functional Residue Annotator checks

Functional Residue Annotator is a pure Python Galaxy tool. It does not require
an external binary. It can run in two modes:

```text
PDB-only mode:
  infer ligand-contact and metal-binding residues from HETATM distance.

PDB + annotation TSV mode:
  additionally load curated catalytic/essential/active-site residues from a TSV.
```

Recommended optional annotation TSV columns:

```text
chain
residue_number
residue_name
role
source
evidence
protect_from_mutation
```

Example sources for a curated TSV are UniProt feature annotations, M-CSA/CSA,
BioLiP, literature curation, or active-site residues from a project notebook.

Run a local smoke test:

```bash
mkdir -p /data/test/functional_residue_selftest
/data/conda_envs/hotspot_wizard/bin/python \
  /data/tools/galaxy_bio/tools/protein_analysis/functional_residue_annotator/run_functional_residue_annotator.py \
  --input-pdb /data/tools/galaxy_bio/tools/protein_analysis/functional_residue_annotator/test-data/mini_functional.pdb \
  --functional-annotation-tsv /data/tools/galaxy_bio/tools/protein_analysis/functional_residue_annotator/test-data/functional_annotations.tsv \
  --residue-annotation-tsv /data/test/functional_residue_selftest/residue_functional_annotation.tsv \
  --active-site-tsv /data/test/functional_residue_selftest/active_site_candidates.tsv \
  --ligand-contacts-tsv /data/test/functional_residue_selftest/ligand_contacts.tsv \
  --report-html /data/test/functional_residue_selftest/functional_annotation_report.html \
  --run-log /data/test/functional_residue_selftest/run_log.txt
```

Connect `active_site_candidates.tsv` to Hotspot Residue Ranker's optional
active-site input to protect catalytic cores, metal-binding residues, and other
essential sites from being over-prioritized as mutation targets.

## Residue Numbering Validator checks

Residue Numbering Validator is also a pure Python tool. It is recommended before
reviewing hotspot results, especially for PDB files whose residue numbers do
not start at 1 or contain insertion codes.

Run a local smoke test:

```bash
mkdir -p /data/test/residue_numbering_selftest
/data/conda_envs/hotspot_wizard/bin/python \
  /data/tools/galaxy_bio/tools/protein_analysis/residue_numbering_validator/run_residue_numbering_validator.py \
  --input-pdb /data/tools/galaxy_bio/tools/protein_analysis/residue_numbering_validator/test-data/mini_numbering.pdb \
  --protein-fasta /data/tools/galaxy_bio/tools/protein_analysis/residue_numbering_validator/test-data/protein_sequence.fasta \
  --msa-fasta /data/tools/galaxy_bio/tools/protein_analysis/residue_numbering_validator/test-data/msa.fasta \
  --structure-features /data/tools/galaxy_bio/tools/protein_analysis/residue_numbering_validator/test-data/residue_structure_features.tsv \
  --conservation-tsv /data/tools/galaxy_bio/tools/protein_analysis/residue_numbering_validator/test-data/residue_conservation.tsv \
  --chain-id A \
  --validation-tsv /data/test/residue_numbering_selftest/residue_numbering_validation.tsv \
  --summary-tsv /data/test/residue_numbering_selftest/numbering_summary.tsv \
  --report-html /data/test/residue_numbering_selftest/numbering_report.html \
  --run-log /data/test/residue_numbering_selftest/run_log.txt
```

The key summary field is `hotspot_merge_risk`:

```text
low       PDB numbering and FASTA positions are directly mergeable.
moderate  Residue identities match, but PDB numbering differs from sequence positions.
high      Length or amino-acid mismatches should be fixed before hotspot ranking.
```

## Homolog Search and MSA install and checks

Install the homology-search and MSA tools into the shared HotSpot Wizard
environment:

```bash
conda install -p /data/conda_envs/hotspot_wizard -y \
  -c conda-forge -c bioconda \
  blast cd-hit mafft muscle mmseqs2
```

Check the binaries:

```bash
/data/conda_envs/hotspot_wizard/bin/blastp -version
/data/conda_envs/hotspot_wizard/bin/makeblastdb -version
/data/conda_envs/hotspot_wizard/bin/cd-hit -h
/data/conda_envs/hotspot_wizard/bin/mafft --version
/data/conda_envs/hotspot_wizard/bin/muscle -version
/data/conda_envs/hotspot_wizard/bin/mmseqs version
```

The Homolog Search and MSA wrapper can use uploaded FASTA datasets or server
databases configured in `config/job_conf.yml`:

```text
/data/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta
/data/databases/uniprot/current_release/uniref/uniref90/uniref90.fasta
/data/databases/uniprot/current_release/uniref/uniref50/uniref50.fasta
/data/databases/uniprot/blastdb/
```

See `tools/protein_analysis/homolog_search_msa/README.md` for download and
`makeblastdb` commands.

## Conservation / Mutability Scorer install and checks

Install Rate4Site either system-wide or in the shared HotSpot Wizard
environment:

```bash
sudo apt-get update
sudo apt-get install -y rate4site
rate4site -h
```

If the binary is installed elsewhere, update `RATE4SITE_BINARY` in
`config/job_conf.yml`.

## Pocket and Tunnel Finder install and checks

Install fpocket, Java, and configure P2Rank/CAVER:

```bash
conda install -p /data/conda_envs/hotspot_wizard -y -c conda-forge fpocket
/data/conda_envs/hotspot_wizard/bin/fpocket -h

conda install -p /data/conda_envs/hotspot_wizard -y -c conda-forge openjdk
/data/conda_envs/hotspot_wizard/bin/java -version
```

P2Rank is expected at `/data/tools/p2rank/prank`, and CAVER at
`/data/tools/caver_3.0/caver/caver.jar`. If installed elsewhere, update
`FPOCKET_BINARY`, `P2RANK_BINARY`, `JAVA_BINARY`, `CAVER_HOME`, and `CAVER_JAR`
in `config/job_conf.yml`.

## Mutation Designer checks

Mutation Designer is a pure Python wrapper and does not require additional
external binaries. It uses the ranked hotspot table from Tool 6, the MSA from
Tool 3, and the target protein FASTA to propose candidate substitutions.

Run a local smoke test:

```bash
/data/conda_envs/hotspot_wizard/bin/python \
  /data/tools/galaxy_bio/tools/protein_analysis/mutation_designer/run_mutation_designer.py \
  --hotspots-tsv /data/tools/galaxy_bio/tools/protein_analysis/mutation_designer/test-data/hotspot_residue_ranked.tsv \
  --msa-fasta /data/tools/galaxy_bio/tools/protein_analysis/mutation_designer/test-data/msa.fasta \
  --protein-sequence /data/tools/galaxy_bio/tools/protein_analysis/mutation_designer/test-data/protein_sequence.fasta \
  --candidate-mutations /data/test/mutation_designer_selftest/candidate_mutations.tsv \
  --smart-library /data/test/mutation_designer_selftest/smart_library.tsv \
  --degenerate-codons /data/test/mutation_designer_selftest/degenerate_codons.tsv \
  --run-log /data/test/mutation_designer_selftest/run_log.txt
```

Expected main outputs:

```text
candidate_mutations.tsv
smart_library.tsv
degenerate_codons.tsv
```

## Stability Predictor / FoldX checks

FoldX is licensed software and must be installed separately by an authorized
user. This repository does not include FoldX binaries or license files. The
configured server path is:

```text
/data/tools/foldx/foldx
```

Check the authorized executable:

```bash
/data/tools/foldx/foldx --help
/data/tools/foldx/foldx --version
```

Run a local Tool 8 smoke test:

```bash
mkdir -p /data/test/stability_predictor_selftest
cd /data/tools/galaxy_bio
python tools/protein_analysis/stability_predictor/run_stability_predictor.py \
  --candidate-mutations tools/protein_analysis/stability_predictor/test-data/candidate_mutations.tsv \
  --pdb tools/protein_analysis/stability_predictor/test-data/foldx_input.pdb \
  --mutation-ddg /data/test/stability_predictor_selftest/mutation_ddg.tsv \
  --filtered-mutations /data/test/stability_predictor_selftest/stability_filtered_mutations.tsv \
  --raw-ddg-output /data/test/stability_predictor_selftest/foldx_raw_ddg.fxout \
  --repaired-pdb /data/test/stability_predictor_selftest/repaired.pdb \
  --mutant-models-dir /data/test/stability_predictor_selftest/mutant_models \
  --archive /data/test/stability_predictor_selftest/foldx_outputs.zip \
  --individual-list /data/test/stability_predictor_selftest/individual_list.txt \
  --summary-json /data/test/stability_predictor_selftest/summary.json \
  --run-log /data/test/stability_predictor_selftest/run_log.txt \
  --foldx-command /data/tools/foldx/foldx
```

For wrapper-only testing without a licensed FoldX binary, use the built-in mock:

```bash
python tools/protein_analysis/stability_predictor/run_stability_predictor.py \
  --candidate-mutations tools/protein_analysis/stability_predictor/test-data/candidate_mutations.tsv \
  --pdb tools/protein_analysis/stability_predictor/test-data/foldx_input.pdb \
  --mutation-ddg /tmp/mutation_ddg.tsv \
  --filtered-mutations /tmp/stability_filtered_mutations.tsv \
  --raw-ddg-output /tmp/foldx_raw_ddg.fxout \
  --repaired-pdb /tmp/repaired.pdb \
  --mutant-models-dir /tmp/mutant_models \
  --archive /tmp/foldx_outputs.zip \
  --individual-list /tmp/individual_list.txt \
  --summary-json /tmp/summary.json \
  --run-log /tmp/run_log.txt \
  --foldx-command mock-foldx
```

The HotSpot filtering categories are:

```text
ΔΔG <= 0 kcal/mol      priority
0 < ΔΔG <= 1 kcal/mol acceptable
1 < ΔΔG <= 2 kcal/mol caution
ΔΔG > 2 kcal/mol      high_risk
```

## Smart Library Report Generator checks

The final report generator is a pure Python summarization wrapper. It consumes
Tool 6, Tool 7, and Tool 8 outputs and writes three final TSV tables plus two
HTML reports.

Run a local smoke test:

```bash
mkdir -p /data/test/smart_library_report_selftest
cd /data/tools/galaxy_bio
/data/conda_envs/hotspot_wizard/bin/python \
  tools/protein_analysis/smart_library_report_generator/run_smart_library_report_generator.py \
  --hotspots-tsv tools/protein_analysis/smart_library_report_generator/test-data/hotspot_residue_ranked.tsv \
  --candidate-mutations tools/protein_analysis/smart_library_report_generator/test-data/candidate_mutations.tsv \
  --smart-library tools/protein_analysis/smart_library_report_generator/test-data/smart_library.tsv \
  --degenerate-codons tools/protein_analysis/smart_library_report_generator/test-data/degenerate_codons.tsv \
  --mutation-ddg tools/protein_analysis/smart_library_report_generator/test-data/mutation_ddg.tsv \
  --structure-pdb tools/protein_analysis/smart_library_report_generator/test-data/foldx_input.pdb \
  --ranked-hotspots /data/test/smart_library_report_selftest/ranked_hotspots.tsv \
  --ranked-mutations /data/test/smart_library_report_selftest/ranked_mutations.tsv \
  --smart-library-design /data/test/smart_library_report_selftest/smart_library_design.tsv \
  --summary-html /data/test/smart_library_report_selftest/summary.html \
  --structure-viewer-html /data/test/smart_library_report_selftest/structure_viewer.html \
  --run-log /data/test/smart_library_report_selftest/run_log.txt
```

Expected final outputs:

```text
ranked_hotspots.tsv
ranked_mutations.tsv
smart_library_design.tsv
summary.html
structure_viewer.html
```
