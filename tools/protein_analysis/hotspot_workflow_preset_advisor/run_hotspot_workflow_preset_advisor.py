#!/usr/bin/env python3
"""Create objective-specific parameter presets for the HotSpot workflow."""

from __future__ import annotations

import argparse
import csv
import html
import json
from pathlib import Path


PRESETS = {
    "balanced": {
        "label": "Balanced semi-rational design",
        "summary": "General-purpose default for activity and stability exploration with moderate library size.",
        "hotspot_ranker": {
            "pocket_or_tunnel_weight": 0.25,
            "mutability_weight": 0.25,
            "surface_accessibility_weight": 0.15,
            "flexibility_weight": 0.10,
            "active_site_distance_weight": 0.10,
            "catalytic_core_penalty": 0.25,
            "highly_conserved_penalty": 0.20,
            "buried_core_penalty": 0.15,
        },
        "mutation_designer": {
            "top_n_hotspots": 20,
            "min_hotspot_score": 0.0,
            "min_msa_frequency": 0.05,
            "max_suggestions_per_site": 8,
            "include_natural": True,
            "include_consensus": True,
            "include_conservative": True,
            "include_alanine": True,
            "include_reduced_alphabet": True,
            "reduced_alphabet": "ADKSTVFY",
        },
        "stability_predictor": {"max_ddg_for_filter": 1.0, "include_caution": False},
    },
    "activity_substrate_specificity": {
        "label": "Activity / substrate specificity",
        "summary": "Prioritizes pocket and tunnel residues that are not absolutely conserved.",
        "hotspot_ranker": {
            "pocket_or_tunnel_weight": 0.38,
            "mutability_weight": 0.24,
            "surface_accessibility_weight": 0.12,
            "flexibility_weight": 0.08,
            "active_site_distance_weight": 0.12,
            "catalytic_core_penalty": 0.30,
            "highly_conserved_penalty": 0.18,
            "buried_core_penalty": 0.10,
        },
        "mutation_designer": {
            "top_n_hotspots": 30,
            "min_hotspot_score": 0.0,
            "min_msa_frequency": 0.03,
            "max_suggestions_per_site": 10,
            "include_natural": True,
            "include_consensus": True,
            "include_conservative": True,
            "include_alanine": True,
            "include_reduced_alphabet": True,
            "reduced_alphabet": "ADKSTVFY",
        },
        "stability_predictor": {"max_ddg_for_filter": 1.0, "include_caution": True},
    },
    "stability": {
        "label": "Stability improvement",
        "summary": "Prioritizes consensus/back-to-consensus mutations and FoldX-friendly substitutions.",
        "hotspot_ranker": {
            "pocket_or_tunnel_weight": 0.12,
            "mutability_weight": 0.22,
            "surface_accessibility_weight": 0.20,
            "flexibility_weight": 0.18,
            "active_site_distance_weight": 0.04,
            "catalytic_core_penalty": 0.30,
            "highly_conserved_penalty": 0.15,
            "buried_core_penalty": 0.18,
        },
        "mutation_designer": {
            "top_n_hotspots": 25,
            "min_hotspot_score": 0.0,
            "min_msa_frequency": 0.05,
            "max_suggestions_per_site": 6,
            "include_natural": True,
            "include_consensus": True,
            "include_conservative": True,
            "include_alanine": False,
            "include_reduced_alphabet": False,
            "reduced_alphabet": "AVILMFYSTNQ",
        },
        "stability_predictor": {"max_ddg_for_filter": 0.0, "include_caution": False},
    },
    "conservative_library": {
        "label": "Conservative small library",
        "summary": "Keeps library size small by using natural/consensus/conservative substitutions only.",
        "hotspot_ranker": {
            "pocket_or_tunnel_weight": 0.25,
            "mutability_weight": 0.30,
            "surface_accessibility_weight": 0.15,
            "flexibility_weight": 0.08,
            "active_site_distance_weight": 0.08,
            "catalytic_core_penalty": 0.30,
            "highly_conserved_penalty": 0.25,
            "buried_core_penalty": 0.15,
        },
        "mutation_designer": {
            "top_n_hotspots": 12,
            "min_hotspot_score": 0.05,
            "min_msa_frequency": 0.08,
            "max_suggestions_per_site": 4,
            "include_natural": True,
            "include_consensus": True,
            "include_conservative": True,
            "include_alanine": False,
            "include_reduced_alphabet": False,
            "reduced_alphabet": "AVILMFYSTNQ",
        },
        "stability_predictor": {"max_ddg_for_filter": 1.0, "include_caution": False},
    },
    "exploratory_library": {
        "label": "Exploratory large library",
        "summary": "Allows broader reduced-alphabet and alanine-scanning suggestions for high-throughput screening.",
        "hotspot_ranker": {
            "pocket_or_tunnel_weight": 0.32,
            "mutability_weight": 0.20,
            "surface_accessibility_weight": 0.14,
            "flexibility_weight": 0.12,
            "active_site_distance_weight": 0.12,
            "catalytic_core_penalty": 0.22,
            "highly_conserved_penalty": 0.15,
            "buried_core_penalty": 0.10,
        },
        "mutation_designer": {
            "top_n_hotspots": 40,
            "min_hotspot_score": 0.0,
            "min_msa_frequency": 0.02,
            "max_suggestions_per_site": 12,
            "include_natural": True,
            "include_consensus": True,
            "include_conservative": True,
            "include_alanine": True,
            "include_reduced_alphabet": True,
            "reduced_alphabet": "ADKSTVFY",
        },
        "stability_predictor": {"max_ddg_for_filter": 2.0, "include_caution": True},
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate HotSpot workflow parameter presets.")
    parser.add_argument("--objective", choices=sorted(PRESETS), default="balanced")
    parser.add_argument("--screening-capacity", choices=["small", "medium", "large"], default="medium")
    parser.add_argument("--has-active-site", choices=["yes", "no"], default="no")
    parser.add_argument("--run-caver", choices=["auto", "yes", "no"], default="auto")
    parser.add_argument("--parameters-tsv", required=True)
    parser.add_argument("--preset-json", required=True)
    parser.add_argument("--workflow-guide-html", required=True)
    return parser.parse_args()


def adjust_for_capacity(preset: dict, capacity: str) -> dict:
    adjusted = json.loads(json.dumps(preset))
    designer = adjusted["mutation_designer"]
    if capacity == "small":
        designer["top_n_hotspots"] = min(designer["top_n_hotspots"], 10)
        designer["max_suggestions_per_site"] = min(designer["max_suggestions_per_site"], 4)
    elif capacity == "large":
        designer["top_n_hotspots"] = max(designer["top_n_hotspots"], 30)
        designer["max_suggestions_per_site"] = max(designer["max_suggestions_per_site"], 10)
    return adjusted


def flatten_params(preset: dict, args: argparse.Namespace) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for tool, params in preset.items():
        if not isinstance(params, dict):
            continue
        for name, value in params.items():
            rows.append({"tool": tool, "parameter": name, "recommended_value": str(value)})
    caver = args.run_caver
    if caver == "auto":
        caver = "yes" if args.has_active_site == "yes" else "no"
    rows.append({"tool": "pocket_tunnel_finder", "parameter": "run_caver", "recommended_value": caver})
    rows.append({"tool": "pocket_tunnel_finder", "parameter": "run_fpocket", "recommended_value": "true"})
    rows.append({"tool": "pocket_tunnel_finder", "parameter": "run_p2rank", "recommended_value": "true"})
    rows.append({"tool": "homolog_search_msa", "parameter": "database", "recommended_value": "swissprot first; uniref90 for broader coverage"})
    rows.append({"tool": "smart_library_report_generator", "parameter": "final_outputs", "recommended_value": "ranked_hotspots.tsv, ranked_mutations.tsv, smart_library_design.tsv, summary.html, structure_viewer.html"})
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["tool", "parameter", "recommended_value"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_html(path: Path, preset: dict, rows: list[dict[str, str]], args: argparse.Namespace) -> None:
    table = "\n".join(
        f"<tr><td>{esc(row['tool'])}</td><td>{esc(row['parameter'])}</td><td>{esc(row['recommended_value'])}</td></tr>"
        for row in rows
    )
    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>HotSpot Workflow Preset Guide</title>
  <style>
    body {{ font-family: ui-sans-serif, system-ui, -apple-system, Segoe UI, sans-serif; background:#f3f6fb; color:#172033; margin:0; }}
    main {{ max-width:1100px; margin:0 auto; padding:32px; }}
    .hero {{ background:white; border:1px solid #dbe3ef; border-radius:18px; padding:24px; box-shadow:0 12px 30px rgba(20,32,52,.07); }}
    h1 {{ margin:0 0 8px; }}
    table {{ width:100%; border-collapse:collapse; background:white; margin-top:18px; }}
    th, td {{ border-bottom:1px solid #dbe3ef; padding:8px 10px; text-align:left; }}
    th {{ background:#172033; color:white; }}
    code {{ background:#eef2ff; padding:2px 5px; border-radius:5px; }}
  </style>
</head>
<body><main>
  <section class="hero">
    <h1>{esc(preset['label'])}</h1>
    <p>{esc(preset['summary'])}</p>
    <p><strong>Screening capacity:</strong> {esc(args.screening_capacity)} |
       <strong>Active-site coordinates available:</strong> {esc(args.has_active_site)}</p>
  </section>
  <h2>Recommended Workflow Order</h2>
  <ol>
    <li>Protein structure PDB + protein FASTA input</li>
    <li>Structure Quality Annotator</li>
    <li>Homolog Search and MSA</li>
    <li>Conservation / Mutability Scorer</li>
    <li>Pocket and Tunnel Finder</li>
    <li>Hotspot Residue Ranker</li>
    <li>Mutation Designer</li>
    <li>Stability Predictor</li>
    <li>Smart Library Report Generator</li>
  </ol>
  <h2>Preset Parameters</h2>
  <table><thead><tr><th>Tool</th><th>Parameter</th><th>Recommended value</th></tr></thead><tbody>{table}</tbody></table>
  <h2>How To Use</h2>
  <p>Import <code>tools/protein_analysis/workflows/hotspot_structure_to_smart_library.ga</code> into Galaxy.
  Then apply the parameter values in this table to the corresponding workflow steps before running.</p>
</main></body></html>
"""
    path.write_text(html_text, encoding="utf-8")


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def main() -> int:
    args = parse_args()
    preset = adjust_for_capacity(PRESETS[args.objective], args.screening_capacity)
    payload = {
        "objective": args.objective,
        "screening_capacity": args.screening_capacity,
        "has_active_site": args.has_active_site,
        "preset": preset,
    }
    rows = flatten_params(preset, args)
    write_tsv(Path(args.parameters_tsv), rows)
    Path(args.preset_json).write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_html(Path(args.workflow_guide_html), preset, rows, args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
