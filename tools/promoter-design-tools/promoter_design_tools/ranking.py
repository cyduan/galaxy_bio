from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd


def aggregate_mutations(mutations: pd.DataFrame) -> pd.DataFrame:
    if mutations.empty:
        return pd.DataFrame(columns=["mutant_id", "parent_id", "mutation_summary", "mutated_regions"])
    rows = []
    for mutant_id, group in mutations.groupby("mutant_id", sort=False):
        parent_id = str(group["parent_id"].iloc[0])
        summaries = []
        for _, row in group.iterrows():
            if row["mutation_type"] == "snv":
                summaries.append(f"{row['ref']}{row['position_1based']}{row['alt']}")
            else:
                summaries.append(f"{row['mutation_type']}@{row['position_1based']}:{row['ref']}->{row['alt']}")
        rows.append(
            {
                "mutant_id": mutant_id,
                "parent_id": parent_id,
                "mutation_summary": ";".join(summaries),
                "mutated_regions": ";".join(sorted(set(str(value) for value in group["region"]))),
            }
        )
    return pd.DataFrame(rows)


def rank_variants(
    wt_conservation: Path,
    mutant_conservation: Path,
    wt_strength: Path,
    mutant_strength: Path,
    mutations_tsv: Path,
    rank_by: str = "predicted_strength",
    min_conservation_score: float | None = None,
    max_conservation_score: float | None = None,
    min_fold_change: float | None = None,
    max_fold_change: float | None = None,
    keep_core_motif: bool = False,
) -> pd.DataFrame:
    wt_cons = pd.read_csv(wt_conservation, sep="\t")
    mut_cons = pd.read_csv(mutant_conservation, sep="\t")
    wt_str = pd.read_csv(wt_strength, sep="\t")
    mut_str = pd.read_csv(mutant_strength, sep="\t")
    mutations = pd.read_csv(mutations_tsv, sep="\t")
    mut_agg = aggregate_mutations(mutations)

    wt_cons = wt_cons.rename(columns={"sequence_id": "parent_id", "conservation_score": "parent_conservation_score"})
    mut_cons = mut_cons.rename(columns={"sequence_id": "mutant_id", "conservation_score": "mutant_conservation_score"})
    wt_str = wt_str.rename(columns={"sequence_id": "parent_id", "predicted_strength": "parent_strength"})
    mut_str = mut_str.rename(columns={"sequence_id": "mutant_id", "predicted_strength": "mutant_strength"})

    df = mut_agg.merge(mut_cons[["mutant_id", "mutant_conservation_score"]], on="mutant_id", how="left")
    df = df.merge(mut_str[["mutant_id", "mutant_strength"]], on="mutant_id", how="left")
    df = df.merge(wt_cons[["parent_id", "parent_conservation_score"]], on="parent_id", how="left")
    df = df.merge(wt_str[["parent_id", "parent_strength"]], on="parent_id", how="left")

    for column in ["parent_conservation_score", "mutant_conservation_score", "parent_strength", "mutant_strength"]:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    df["fold_change_vs_parent"] = df["mutant_strength"] / df["parent_strength"].replace(0, pd.NA)
    df["delta_conservation_score"] = df["mutant_conservation_score"] - df["parent_conservation_score"]

    if min_conservation_score is not None:
        df = df[df["mutant_conservation_score"] >= min_conservation_score]
    if max_conservation_score is not None:
        df = df[df["mutant_conservation_score"] <= max_conservation_score]
    if min_fold_change is not None:
        df = df[df["fold_change_vs_parent"] >= min_fold_change]
    if max_fold_change is not None:
        df = df[df["fold_change_vs_parent"] <= max_fold_change]
    if keep_core_motif:
        df = df[~df["mutated_regions"].str.contains("minus35|minus10", na=False)]

    if rank_by == "predicted_strength":
        sort_column = "mutant_strength"
    elif rank_by == "fold_change":
        sort_column = "fold_change_vs_parent"
    elif rank_by == "conservation_score":
        sort_column = "mutant_conservation_score"
    elif rank_by == "custom":
        df["custom_score"] = (
            percentile(df["mutant_strength"]) * 0.5
            + percentile(df["fold_change_vs_parent"]) * 0.3
            + percentile(df["mutant_conservation_score"]) * 0.2
        )
        sort_column = "custom_score"
    else:
        raise ValueError(f"Unsupported rank_by value: {rank_by}")

    df = df.sort_values(sort_column, ascending=False).reset_index(drop=True)
    df["rank"] = range(1, len(df) + 1)
    df["recommendation"] = df.apply(recommendation, axis=1)
    columns = [
        "mutant_id",
        "parent_id",
        "mutation_summary",
        "mutated_regions",
        "parent_conservation_score",
        "mutant_conservation_score",
        "delta_conservation_score",
        "parent_strength",
        "mutant_strength",
        "fold_change_vs_parent",
        "rank",
        "recommendation",
    ]
    return df[columns]


def percentile(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.nunique(dropna=True) <= 1:
        return pd.Series([1.0] * len(series), index=series.index)
    return numeric.rank(pct=True).fillna(0.0)


def recommendation(row: pd.Series) -> str:
    if row["fold_change_vs_parent"] >= 1.5 and row["mutant_conservation_score"] >= 0.65:
        return "high_priority"
    if row["fold_change_vs_parent"] >= 1.0 and row["mutant_conservation_score"] >= 0.45:
        return "review"
    return "low_priority"


def write_summary_html(df: pd.DataFrame, output_html: Path) -> None:
    top = df.head(20).copy()
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Promoter Variant Ranking Summary</title>
  <style>
    body {{ font-family: sans-serif; margin: 2rem; color: #1f2937; }}
    table {{ border-collapse: collapse; width: 100%; font-size: 0.9rem; }}
    th, td {{ border: 1px solid #d1d5db; padding: 0.4rem; }}
    th {{ background: #e5e7eb; }}
    .metric {{ display: inline-block; margin-right: 1rem; padding: 0.7rem 1rem; background: #f3f4f6; border-radius: 0.5rem; }}
  </style>
</head>
<body>
  <h1>Promoter Variant Ranking Summary</h1>
  <p class="metric"><strong>Variants ranked:</strong> {len(df)}</p>
  <p class="metric"><strong>Median fold-change:</strong> {df['fold_change_vs_parent'].median():.3f}</p>
  <p class="metric"><strong>Best predicted strength:</strong> {df['mutant_strength'].max():.3f}</p>
  <h2>Top candidates</h2>
  {top.to_html(index=False, escape=True)}
  <p><em>Heuristic strength is a relative score, not an absolute transcription rate.</em></p>
</body>
</html>
"""
    output_html.write_text(html, encoding="utf-8")


def write_plots(df: pd.DataFrame, strength_png: Path, scatter_png: Path) -> None:
    plt.figure(figsize=(7, 4.5))
    plt.hist(df["mutant_strength"].dropna(), bins=20, color="#2563eb", alpha=0.85)
    plt.xlabel("Predicted relative strength")
    plt.ylabel("Variant count")
    plt.title("Promoter library strength distribution")
    plt.tight_layout()
    plt.savefig(strength_png, dpi=160)
    plt.close()

    plt.figure(figsize=(6, 5))
    plt.scatter(df["mutant_conservation_score"], df["mutant_strength"], c=df["fold_change_vs_parent"], cmap="viridis", edgecolor="black", linewidth=0.2)
    plt.xlabel("Mutant conservation score")
    plt.ylabel("Predicted relative strength")
    plt.title("Conservation vs strength")
    cbar = plt.colorbar()
    cbar.set_label("Fold-change vs parent")
    plt.tight_layout()
    plt.savefig(scatter_png, dpi=160)
    plt.close()

