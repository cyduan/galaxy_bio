# Workflow Design

The workflow in `workflows/promoter_conservation_mutagenesis_strength.ga`
contains nine logical steps:

1. Normalize WT promoter windows.
2. Scan WT promoters for -35 and -10 candidates.
3. Score WT motif pairs.
4. Generate a reproducible random mutant library.
5. Scan mutant promoters.
6. Score mutant motif pairs.
7. Predict WT strength.
8. Predict mutant strength.
9. Merge, rank, and visualize variants.

The workflow keeps intermediate TSV and FASTA datasets. This is deliberate:
users can inspect motif calls, reject poor windows, or swap the mutagenesis
region before ranking variants.

Recommended first run:

- `mutants_per_sequence`: 10-100
- `mutation_mode`: `snv`
- `region_mode`: `non_core` for conservative tuning, `full` for broader search
- `seed`: fixed integer for reproducibility

