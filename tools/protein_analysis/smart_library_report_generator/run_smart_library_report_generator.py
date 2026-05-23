#!/usr/bin/env python3
"""Generate final HotSpot smart-library design reports."""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path


RANKED_HOTSPOT_COLUMNS = [
    "rank",
    "residue_id",
    "chain",
    "residue_number",
    "wt_aa",
    "hotspot_score",
    "hotspot_class",
    "strategy",
    "recommendation",
    "in_pocket",
    "pocket_ids",
    "in_tunnel",
    "tunnel_ids",
    "exposure_class",
    "structural_context",
    "secondary_structure_class",
    "relative_asa",
    "conservation_score",
    "mutability_score",
    "consensus_aa",
    "accepted_aas",
    "best_mutation",
    "best_ddg_kcal_per_mol",
    "best_stability_class",
    "library_aas",
    "degenerate_codon",
    "site_variant_count",
    "notes",
]

RANKED_MUTATION_COLUMNS = [
    "rank",
    "mutation_id",
    "position",
    "residue_id",
    "chain",
    "wt_aa",
    "suggested_aa",
    "foldx_mutation",
    "strategy",
    "hotspot_class",
    "structural_context",
    "conservation_score",
    "mutability_score",
    "hotspot_score",
    "ddg_kcal_per_mol",
    "stability_class",
    "keep_for_hotspot",
    "mutation_type",
    "candidate_reason",
    "msa_frequency",
    "source",
    "degenerate_codon",
    "library_id",
    "recommendation",
    "notes",
]

SMART_LIBRARY_COLUMNS = [
    "library_id",
    "design_type",
    "positions",
    "residue_ids",
    "wt_aas",
    "designed_aas",
    "recommended_degenerate_codons",
    "theoretical_variant_count",
    "coverage_fraction",
    "hotspot_classes",
    "best_ddg_kcal_per_mol",
    "mean_ddg_kcal_per_mol",
    "max_ddg_kcal_per_mol",
    "kept_mutations",
    "caution_mutations",
    "high_risk_mutations",
    "rationale",
    "notes",
]


class ToolError(RuntimeError):
    """User-facing CLI error."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge HotSpot ranking, mutation design, FoldX stability, and library design into final reports."
    )
    parser.add_argument("--hotspots-tsv", required=True, help="Tool 6 hotspot_residue_ranked.tsv.")
    parser.add_argument("--candidate-mutations", required=True, help="Tool 7 candidate_mutations.tsv.")
    parser.add_argument("--smart-library", required=True, help="Tool 7 smart_library.tsv.")
    parser.add_argument("--degenerate-codons", required=True, help="Tool 7 degenerate_codons.tsv.")
    parser.add_argument("--mutation-ddg", required=True, help="Tool 8 mutation_ddg.tsv.")
    parser.add_argument("--ranked-hotspots", required=True, help="Output ranked_hotspots.tsv.")
    parser.add_argument("--ranked-mutations", required=True, help="Output ranked_mutations.tsv.")
    parser.add_argument("--smart-library-design", required=True, help="Output smart_library_design.tsv.")
    parser.add_argument("--summary-html", required=True, help="Output summary.html.")
    parser.add_argument("--structure-viewer-html", required=True, help="Output structure_viewer.html.")
    parser.add_argument("--run-log", required=True, help="Run log.")
    parser.add_argument("--structure-pdb", default="", help="Optional PDB structure to embed in the viewer HTML.")
    parser.add_argument("--project-title", default="HotSpot Smart Library Design Report")
    parser.add_argument("--top-n-hotspots", type=int, default=50)
    parser.add_argument("--top-n-mutations", type=int, default=100)
    parser.add_argument("--max-viewer-residues", type=int, default=80)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ToolError(f"{path} is empty or does not contain a TSV header.")
        return [{key: (value or "") for key, value in row.items()} for row in reader]


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def parse_float(value: str) -> float | None:
    try:
        if value in {"", "NA", "None", "nan"}:
            return None
        result = float(value)
        if math.isnan(result):
            return None
        return result
    except (TypeError, ValueError):
        return None


def parse_int(value: str) -> int | None:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def format_float(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.6g}"


def get_position(row: dict[str, str]) -> str:
    return row.get("residue_number") or row.get("position") or parse_position_from_residue_id(row.get("residue_id", ""))


def parse_position_from_residue_id(value: str) -> str:
    match = re.search(r"(-?\d+)", value or "")
    return match.group(1) if match else ""


def row_key(row: dict[str, str]) -> str:
    chain = row.get("chain") or parse_chain_from_residue_id(row.get("residue_id", "")) or "A"
    position = get_position(row)
    return f"{chain}:{position}" if position else row.get("residue_id", "")


def parse_chain_from_residue_id(value: str) -> str:
    match = re.match(r"\s*([A-Za-z0-9])\s*:", value or "")
    return match.group(1) if match else ""


def normalize_bool(value: str) -> bool:
    return str(value).strip().lower() in {"yes", "true", "1", "y"}


def structural_context(row: dict[str, str]) -> str:
    contexts: list[str] = []
    if normalize_bool(row.get("in_pocket", "")):
        contexts.append("pocket")
    if normalize_bool(row.get("in_tunnel", "")):
        contexts.append("tunnel")
    exposure = row.get("exposure_class") or row.get("freesasa_exposure_class") or ""
    if exposure:
        contexts.append(exposure)
    if not contexts and parse_float(row.get("relative_asa", "")) is not None:
        asa = parse_float(row.get("relative_asa", "")) or 0.0
        contexts.append("core" if asa < 0.09 else "surface")
    return ",".join(dict.fromkeys(contexts)) or "unknown"


def infer_strategy(hotspot: dict[str, str], mutation: dict[str, str] | None = None) -> str:
    text = " ".join(
        [
            hotspot.get("hotspot_class", ""),
            hotspot.get("recommendation", ""),
            mutation.get("candidate_reason", "") if mutation else "",
            mutation.get("reason", "") if mutation else "",
            mutation.get("source", "") if mutation else "",
        ]
    ).lower()
    if "coevolution" in text or "co-evolution" in text:
        return "coevolution"
    if "consensus" in text:
        return "consensus"
    if "stability" in text:
        return "stability_hotspot"
    if "functional" in text or "function" in text or "pocket" in text or "tunnel" in text:
        return "functional_hotspot"
    if "risky" in text:
        return "risky_hotspot"
    return "exploratory"


def ddg_sort_key(row: dict[str, str]) -> tuple[int, float, int]:
    keep_rank = 0 if normalize_bool(row.get("keep_for_hotspot", "")) else 1
    ddg = parse_float(row.get("ddg_kcal_per_mol", ""))
    ddg_value = ddg if ddg is not None else 9999.0
    risk_rank = {"priority": 0, "acceptable": 1, "caution": 2, "high_risk": 3, "unknown": 4}.get(
        row.get("stability_class", ""), 4
    )
    return (keep_rank, ddg_value, risk_rank)


def build_mutation_maps(
    candidates: list[dict[str, str]], ddg_rows: list[dict[str, str]]
) -> tuple[dict[str, dict[str, str]], dict[str, list[dict[str, str]]]]:
    candidate_by_mutation: dict[str, dict[str, str]] = {}
    candidate_by_site_mut: dict[str, dict[str, str]] = {}
    for row in candidates:
        position = row.get("position", "")
        chain = row.get("chain") or parse_chain_from_residue_id(row.get("residue_id", "")) or "A"
        wt = row.get("wt_aa", "")
        alt = row.get("suggested_aa", "")
        mutation = row.get("foldx_mutation") or f"{wt}{chain}{position}{alt}"
        candidate_by_mutation[mutation] = row
        candidate_by_site_mut[f"{chain}:{position}:{alt}"] = row

    merged_by_position: dict[str, list[dict[str, str]]] = defaultdict(list)
    mutation_by_foldx: dict[str, dict[str, str]] = {}
    for ddg in ddg_rows:
        mutation = ddg.get("foldx_mutation", "")
        candidate = candidate_by_mutation.get(mutation)
        if candidate is None:
            key = f"{ddg.get('chain', 'A')}:{ddg.get('position', '')}:{ddg.get('mutant_aa', '')}"
            candidate = candidate_by_site_mut.get(key, {})
        merged = {**candidate, **ddg}
        if "suggested_aa" not in merged or not merged.get("suggested_aa"):
            merged["suggested_aa"] = ddg.get("mutant_aa", "")
        if "reason" in candidate and "candidate_reason" not in merged:
            merged["candidate_reason"] = candidate["reason"]
        position_key = f"{merged.get('chain', 'A')}:{merged.get('position', '')}"
        merged_by_position[position_key].append(merged)
        mutation_by_foldx[mutation] = merged

    for rows in merged_by_position.values():
        rows.sort(key=ddg_sort_key)
    return mutation_by_foldx, merged_by_position


def build_degenerate_map(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        key = row_key(row)
        result[key] = row
    return result


def build_library_position_map(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    result: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        for position in split_positions(row.get("positions", "")):
            chains = split_positions(row.get("residue_ids", ""))
            chain = "A"
            for residue_id in chains:
                parsed_pos = parse_position_from_residue_id(residue_id)
                if parsed_pos == position:
                    chain = parse_chain_from_residue_id(residue_id) or chain
                    break
            result[f"{chain}:{position}"].append(row)
    return result


def split_positions(value: str) -> list[str]:
    return [part.strip() for part in re.split(r"[;,]", value or "") if part.strip()]


def build_ranked_hotspots(
    hotspots: list[dict[str, str]],
    mutations_by_position: dict[str, list[dict[str, str]]],
    degenerate_by_position: dict[str, dict[str, str]],
    library_by_position: dict[str, list[dict[str, str]]],
    top_n: int,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    sorted_hotspots = sorted(
        hotspots,
        key=lambda row: (
            parse_int(row.get("rank", "")) or 999999,
            -(parse_float(row.get("hotspot_score", "")) or 0.0),
        ),
    )
    for hotspot in sorted_hotspots[:top_n]:
        key = row_key(hotspot)
        best_mutation = mutations_by_position.get(key, [{}])[0]
        degenerate = degenerate_by_position.get(key, {})
        libraries = library_by_position.get(key, [])
        site_library = first_site_library(libraries)
        best_alt = best_mutation.get("suggested_aa") or best_mutation.get("mutant_aa", "")
        wt = hotspot.get("wt_aa", "")
        best_label = f"{wt}{hotspot.get('residue_number', get_position(hotspot))}{best_alt}" if best_alt else ""
        row = {
            "rank": hotspot.get("rank", ""),
            "residue_id": hotspot.get("residue_id", key),
            "chain": hotspot.get("chain", key.split(":")[0] if ":" in key else ""),
            "residue_number": get_position(hotspot),
            "wt_aa": wt,
            "hotspot_score": hotspot.get("hotspot_score", ""),
            "hotspot_class": hotspot.get("hotspot_class", ""),
            "strategy": infer_strategy(hotspot),
            "recommendation": hotspot.get("recommendation", ""),
            "in_pocket": hotspot.get("in_pocket", ""),
            "pocket_ids": hotspot.get("pocket_ids", ""),
            "in_tunnel": hotspot.get("in_tunnel", ""),
            "tunnel_ids": hotspot.get("tunnel_ids", ""),
            "exposure_class": hotspot.get("exposure_class", ""),
            "structural_context": structural_context(hotspot),
            "secondary_structure_class": hotspot.get("secondary_structure_class", ""),
            "relative_asa": hotspot.get("relative_asa", "") or hotspot.get("freesasa_relative_asa", ""),
            "conservation_score": hotspot.get("conservation_score", ""),
            "mutability_score": hotspot.get("mutability_score", ""),
            "consensus_aa": hotspot.get("consensus_aa", ""),
            "accepted_aas": hotspot.get("accepted_aas", ""),
            "best_mutation": best_label,
            "best_ddg_kcal_per_mol": best_mutation.get("ddg_kcal_per_mol", ""),
            "best_stability_class": best_mutation.get("stability_class", ""),
            "library_aas": degenerate.get("designed_aas", ""),
            "degenerate_codon": degenerate.get("degenerate_codon", ""),
            "site_variant_count": site_library.get("theoretical_variant_count", ""),
            "notes": hotspot.get("notes", ""),
        }
        rows.append(row)
    return rows


def first_site_library(libraries: list[dict[str, str]]) -> dict[str, str]:
    for library in libraries:
        if library.get("design_type") == "single_site":
            return library
    return libraries[0] if libraries else {}


def build_ranked_mutations(
    hotspots_by_position: dict[str, dict[str, str]],
    mutations_by_position: dict[str, list[dict[str, str]]],
    degenerate_by_position: dict[str, dict[str, str]],
    library_by_position: dict[str, list[dict[str, str]]],
    top_n: int,
) -> list[dict[str, str]]:
    all_mutations: list[dict[str, str]] = []
    for rows in mutations_by_position.values():
        all_mutations.extend(rows)
    all_mutations.sort(
        key=lambda row: (
            0 if normalize_bool(row.get("keep_for_hotspot", "")) else 1,
            parse_float(row.get("ddg_kcal_per_mol", "")) if parse_float(row.get("ddg_kcal_per_mol", "")) is not None else 9999,
            parse_int(row.get("hotspot_rank", "")) or 999999,
            row.get("foldx_mutation", ""),
        )
    )

    output: list[dict[str, str]] = []
    for rank, mutation in enumerate(all_mutations[:top_n], start=1):
        position_key = f"{mutation.get('chain', 'A')}:{mutation.get('position', '')}"
        hotspot = hotspots_by_position.get(position_key, {})
        degenerate = degenerate_by_position.get(position_key, {})
        libraries = library_by_position.get(position_key, [])
        row = {
            "rank": str(rank),
            "mutation_id": mutation.get("mutation_id", ""),
            "position": mutation.get("position", ""),
            "residue_id": hotspot.get("residue_id", position_key),
            "chain": mutation.get("chain", hotspot.get("chain", "")),
            "wt_aa": mutation.get("wt_aa", hotspot.get("wt_aa", "")),
            "suggested_aa": mutation.get("suggested_aa", "") or mutation.get("mutant_aa", ""),
            "foldx_mutation": mutation.get("foldx_mutation", ""),
            "strategy": infer_strategy(hotspot, mutation),
            "hotspot_class": mutation.get("hotspot_class", hotspot.get("hotspot_class", "")),
            "structural_context": structural_context(hotspot),
            "conservation_score": hotspot.get("conservation_score", ""),
            "mutability_score": hotspot.get("mutability_score", ""),
            "hotspot_score": mutation.get("hotspot_score", hotspot.get("hotspot_score", "")),
            "ddg_kcal_per_mol": mutation.get("ddg_kcal_per_mol", ""),
            "stability_class": mutation.get("stability_class", ""),
            "keep_for_hotspot": mutation.get("keep_for_hotspot", ""),
            "mutation_type": mutation.get("mutation_type", ""),
            "candidate_reason": mutation.get("candidate_reason", "") or mutation.get("reason", ""),
            "msa_frequency": mutation.get("msa_frequency", ""),
            "source": mutation.get("source", ""),
            "degenerate_codon": degenerate.get("degenerate_codon", ""),
            "library_id": ",".join(lib.get("library_id", "") for lib in libraries if lib.get("design_type") == "single_site"),
            "recommendation": mutation_recommendation(mutation),
            "notes": mutation.get("warning", "") or mutation.get("notes", ""),
        }
        output.append(row)
    return output


def mutation_recommendation(row: dict[str, str]) -> str:
    stability = row.get("stability_class", "")
    keep = normalize_bool(row.get("keep_for_hotspot", ""))
    if stability == "priority":
        return "prioritize"
    if stability == "acceptable" and keep:
        return "include"
    if stability == "caution":
        return "consider_if_functionally_interesting"
    if stability == "high_risk":
        return "deprioritize_unless_functional_gain"
    return "review"


def build_smart_library_design(
    smart_library: list[dict[str, str]],
    mutations_by_position: dict[str, list[dict[str, str]]],
) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    for row in smart_library:
        positions = split_positions(row.get("positions", ""))
        ddgs: list[float] = []
        kept = caution = high_risk = 0
        for position in positions:
            for key, mutations in mutations_by_position.items():
                if key.endswith(f":{position}"):
                    for mutation in mutations:
                        ddg = parse_float(mutation.get("ddg_kcal_per_mol", ""))
                        if ddg is not None:
                            ddgs.append(ddg)
                        if normalize_bool(mutation.get("keep_for_hotspot", "")):
                            kept += 1
                        if mutation.get("stability_class") == "caution":
                            caution += 1
                        if mutation.get("stability_class") == "high_risk":
                            high_risk += 1
        output.append(
            {
                "library_id": row.get("library_id", ""),
                "design_type": row.get("design_type", ""),
                "positions": row.get("positions", ""),
                "residue_ids": row.get("residue_ids", ""),
                "wt_aas": row.get("wt_aas", ""),
                "designed_aas": row.get("designed_aas", ""),
                "recommended_degenerate_codons": row.get("recommended_degenerate_codons", ""),
                "theoretical_variant_count": row.get("theoretical_variant_count", ""),
                "coverage_fraction": row.get("coverage_fraction", ""),
                "hotspot_classes": row.get("hotspot_classes", ""),
                "best_ddg_kcal_per_mol": format_float(min(ddgs) if ddgs else None),
                "mean_ddg_kcal_per_mol": format_float(sum(ddgs) / len(ddgs) if ddgs else None),
                "max_ddg_kcal_per_mol": format_float(max(ddgs) if ddgs else None),
                "kept_mutations": str(kept),
                "caution_mutations": str(caution),
                "high_risk_mutations": str(high_risk),
                "rationale": row.get("rationale", ""),
                "notes": row.get("notes", ""),
            }
        )
    return output


def html_table(rows: list[dict[str, str]], columns: list[str], limit: int | None = None) -> str:
    limited = rows[:limit] if limit else rows
    header = "".join(f"<th>{escape(column)}</th>" for column in columns)
    body_parts = []
    for row in limited:
        cells = "".join(f"<td>{format_cell(column, row.get(column, ''))}</td>" for column in columns)
        body_parts.append(f"<tr>{cells}</tr>")
    if not body_parts:
        body_parts.append(f"<tr><td colspan='{len(columns)}'>No rows available.</td></tr>")
    return f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(body_parts)}</tbody></table>"


def format_cell(column: str, value: str) -> str:
    if column in {"stability_class", "hotspot_class", "strategy", "keep_for_hotspot", "recommendation"}:
        css = re.sub(r"[^a-zA-Z0-9_-]+", "_", value.lower() or "empty")
        return f"<span class='badge {css}'>{escape(value)}</span>"
    return escape(value)


def escape(value: object) -> str:
    return html.escape(str(value), quote=True)


def render_summary_html(
    path: Path,
    title: str,
    hotspots: list[dict[str, str]],
    mutations: list[dict[str, str]],
    libraries: list[dict[str, str]],
    viewer_name: str,
    log_summary: dict[str, object],
) -> None:
    priority = sum(1 for row in mutations if row.get("stability_class") == "priority")
    acceptable = sum(1 for row in mutations if row.get("stability_class") == "acceptable")
    caution = sum(1 for row in mutations if row.get("stability_class") == "caution")
    high_risk = sum(1 for row in mutations if row.get("stability_class") == "high_risk")
    kept = sum(1 for row in mutations if normalize_bool(row.get("keep_for_hotspot", "")))
    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{escape(title)}</title>
  <style>{base_css()}</style>
</head>
<body>
  <main>
    <h1>{escape(title)}</h1>
    <p class="lede">Final HotSpot smart-library report for experimental design. Tables are ranked for practical screening: stable or minimally destabilizing mutations first, then caution and high-risk variants.</p>
    <section class="cards">
      {metric_card("Hotspot sites", len(hotspots))}
      {metric_card("Candidate mutations", len(mutations))}
      {metric_card("Kept after stability filter", kept)}
      {metric_card("Library designs", len(libraries))}
    </section>
    <section class="cards">
      {metric_card("Priority ddG <= 0", priority)}
      {metric_card("Acceptable 0-1", acceptable)}
      {metric_card("Caution 1-2", caution)}
      {metric_card("High risk > 2", high_risk)}
    </section>

    <section>
      <h2>1. Recommended Hotspot Sites</h2>
      <p>Includes structural context, pocket/tunnel/surface/core labels, conservation, mutability, best FoldX-supported mutation, and codon/library information.</p>
      {html_table(hotspots, RANKED_HOTSPOT_COLUMNS[:25], limit=30)}
    </section>

    <section>
      <h2>2. Ranked Mutations</h2>
      <p>Mutations are sorted by keep flag, FoldX ΔΔG, hotspot rank, and mutation identifier.</p>
      {html_table(mutations, RANKED_MUTATION_COLUMNS[:24], limit=50)}
    </section>

    <section>
      <h2>3. Smart Library Design</h2>
      <p>Variant counts are theoretical amino-acid/codon coverage estimates. Split pooled libraries if the total size is too large for the screening capacity.</p>
      {html_table(libraries, SMART_LIBRARY_COLUMNS, limit=30)}
    </section>

    <section>
      <h2>4. Interpretation Rules</h2>
      <ul>
        <li><strong>priority</strong>: ΔΔG &le; 0 kcal/mol, likely stabilizing or non-destabilizing.</li>
        <li><strong>acceptable</strong>: 0 &lt; ΔΔG &le; 1 kcal/mol, generally acceptable.</li>
        <li><strong>caution</strong>: 1 &lt; ΔΔG &le; 2 kcal/mol, keep only with strong functional rationale.</li>
        <li><strong>high_risk</strong>: ΔΔG &gt; 2 kcal/mol, usually deprioritize unless functional gain is expected.</li>
      </ul>
      <p><a href="{escape(viewer_name)}">Open structure viewer / residue selection helper</a></p>
    </section>

    <section>
      <h2>Run Summary</h2>
      <pre>{escape(json.dumps(log_summary, indent=2, sort_keys=True))}</pre>
    </section>
  </main>
</body>
</html>
"""
    path.write_text(html_text, encoding="utf-8")


def metric_card(label: str, value: object) -> str:
    return f"<div class='card'><div class='metric'>{escape(value)}</div><div class='label'>{escape(label)}</div></div>"


def base_css() -> str:
    return """
:root { color-scheme: light; --ink:#18212f; --muted:#5c677a; --line:#d7dee8; --panel:#f7f9fc; --accent:#1f6feb; }
body { margin:0; font-family: ui-sans-serif, system-ui, -apple-system, Segoe UI, sans-serif; color:var(--ink); background:#eef3f8; }
main { max-width: 1320px; margin: 0 auto; padding: 32px; }
h1 { margin: 0 0 8px; font-size: 34px; letter-spacing:-0.03em; }
h2 { margin-top: 34px; padding-top: 10px; border-top: 2px solid var(--line); }
.lede { color:var(--muted); max-width: 960px; font-size: 16px; }
.cards { display:grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr)); gap: 12px; margin: 18px 0; }
.card { background:white; border:1px solid var(--line); border-radius: 14px; padding: 16px; box-shadow:0 8px 24px rgba(18,31,53,0.06); }
.metric { font-size: 28px; font-weight: 750; color:var(--accent); }
.label { color:var(--muted); margin-top: 4px; }
table { width:100%; border-collapse: collapse; background:white; border:1px solid var(--line); font-size: 13px; }
th, td { border-bottom:1px solid var(--line); padding: 7px 8px; text-align:left; vertical-align:top; }
th { background:#172033; color:white; position:sticky; top:0; }
tr:nth-child(even) td { background:#f8fafc; }
pre { background:#111827; color:#e5e7eb; padding:16px; border-radius:12px; overflow:auto; }
.badge { display:inline-block; padding:2px 7px; border-radius:999px; background:#e9eef7; color:#1f2937; font-weight:650; }
.priority, .acceptable, .yes, .include, .prioritize { background:#dcfce7; color:#14532d; }
.caution, .consider_if_functionally_interesting { background:#fef3c7; color:#7c2d12; }
.high_risk, .no, .deprioritize_unless_functional_gain { background:#fee2e2; color:#7f1d1d; }
.functional_hotspot, .functional_hotspot { background:#dbeafe; color:#1e3a8a; }
.stability_hotspot { background:#ede9fe; color:#4c1d95; }
.consensus { background:#cffafe; color:#164e63; }
"""


def render_viewer_html(
    path: Path,
    title: str,
    hotspots: list[dict[str, str]],
    mutations: list[dict[str, str]],
    structure_pdb: str,
    max_residues: int,
) -> None:
    selected_hotspots = hotspots[:max_residues]
    selections = make_structure_selections(selected_hotspots)
    pdb_text = ""
    if structure_pdb:
        try:
            pdb_text = Path(structure_pdb).read_text(encoding="utf-8", errors="replace")
        except OSError:
            pdb_text = ""
    data_script = f"<script type='text/plain' id='pdb-data'>{escape(pdb_text)}</script>" if pdb_text else ""
    mol_script = viewer_script() if pdb_text else ""
    viewer_panel = "<div id='viewer'></div>" if pdb_text else "<div class='notice'>No PDB was provided. Use the selection commands below with your preferred structure viewer.</div>"
    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{escape(title)} - structure viewer</title>
  <style>{base_css()}{viewer_css()}</style>
  <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
</head>
<body>
  <main>
    <h1>{escape(title)} - Structure Viewer</h1>
    <p class="lede">This page embeds hotspot residue selections and optional PDB content. If the browser cannot load 3Dmol.js, the tables and PyMOL/ChimeraX commands remain usable.</p>
    {viewer_panel}
    {data_script}
    <section>
      <h2>Hotspot Selection Commands</h2>
      <div class="grid2">
        <div>
          <h3>PyMOL</h3>
          <pre>{escape(selections["pymol"])}</pre>
        </div>
        <div>
          <h3>ChimeraX</h3>
          <pre>{escape(selections["chimerax"])}</pre>
        </div>
      </div>
    </section>
    <section>
      <h2>Hotspot Residues</h2>
      {html_table(selected_hotspots, ["rank", "residue_id", "chain", "residue_number", "wt_aa", "hotspot_class", "structural_context", "best_mutation", "best_ddg_kcal_per_mol"], limit=max_residues)}
    </section>
    <section>
      <h2>Top Mutations</h2>
      {html_table(mutations, ["rank", "foldx_mutation", "strategy", "stability_class", "ddg_kcal_per_mol", "recommendation"], limit=50)}
    </section>
  </main>
  {mol_script}
</body>
</html>
"""
    path.write_text(html_text, encoding="utf-8")


def viewer_css() -> str:
    return """
#viewer { width:100%; height:560px; border-radius:16px; border:1px solid var(--line); background:#0f172a; margin:18px 0; overflow:hidden; }
.notice { padding:16px; background:#fff7ed; border:1px solid #fed7aa; border-radius:14px; color:#7c2d12; }
.grid2 { display:grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap:16px; }
"""


def make_structure_selections(hotspots: list[dict[str, str]]) -> dict[str, str]:
    by_chain: dict[str, list[str]] = defaultdict(list)
    for row in hotspots:
        chain = row.get("chain", "A") or "A"
        residue = row.get("residue_number") or parse_position_from_residue_id(row.get("residue_id", ""))
        if residue:
            by_chain[chain].append(residue)
    pymol_lines = ["hide everything", "show cartoon", "color gray80"]
    chimerax_lines = ["hide atoms", "show cartoons", "color gray"]
    for chain, residues in by_chain.items():
        unique_residues = ",".join(dict.fromkeys(residues))
        pymol_lines.append(f"select hotspots_{chain}, chain {chain} and resi {unique_residues}")
        pymol_lines.append(f"show sticks, hotspots_{chain}")
        pymol_lines.append(f"color orange, hotspots_{chain}")
        chimerax_lines.append(f"select /{chain}:{unique_residues}")
        chimerax_lines.append("style sel stick")
        chimerax_lines.append("color sel orange")
    return {"pymol": "\n".join(pymol_lines), "chimerax": "\n".join(chimerax_lines)}


def viewer_script() -> str:
    return """
<script>
(function () {
  const element = document.getElementById("viewer");
  const pdbNode = document.getElementById("pdb-data");
  if (!element || !pdbNode || typeof $3Dmol === "undefined") {
    if (element) {
      element.innerHTML = "<div style='color:white;padding:24px'>3Dmol.js could not be loaded. Use the selection commands below.</div>";
    }
    return;
  }
  const viewer = $3Dmol.createViewer(element, { backgroundColor: "#0f172a" });
  viewer.addModel(pdbNode.textContent, "pdb");
  viewer.setStyle({}, { cartoon: { color: "lightgray" } });
  const rows = Array.from(document.querySelectorAll("table tbody tr"));
  rows.slice(0, 80).forEach(function (row) {
    const cells = row.querySelectorAll("td");
    if (cells.length < 4) return;
    const chain = cells[2].textContent.trim();
    const resi = parseInt(cells[3].textContent.trim(), 10);
    if (!chain || !resi) return;
    viewer.addStyle({ chain: chain, resi: resi }, { stick: { colorscheme: "orangeCarbon", radius: 0.22 } });
  });
  viewer.zoomTo();
  viewer.render();
})();
</script>
"""


def main() -> int:
    args = parse_args()
    log_lines: list[str] = []
    try:
        hotspots = read_tsv(Path(args.hotspots_tsv))
        candidates = read_tsv(Path(args.candidate_mutations))
        smart_library = read_tsv(Path(args.smart_library))
        degenerate_codons = read_tsv(Path(args.degenerate_codons))
        mutation_ddg = read_tsv(Path(args.mutation_ddg))

        _, mutations_by_position = build_mutation_maps(candidates, mutation_ddg)
        hotspots_by_position = {row_key(row): row for row in hotspots}
        degenerate_by_position = build_degenerate_map(degenerate_codons)
        library_by_position = build_library_position_map(smart_library)

        ranked_hotspots = build_ranked_hotspots(
            hotspots,
            mutations_by_position,
            degenerate_by_position,
            library_by_position,
            args.top_n_hotspots,
        )
        ranked_mutations = build_ranked_mutations(
            hotspots_by_position,
            mutations_by_position,
            degenerate_by_position,
            library_by_position,
            args.top_n_mutations,
        )
        smart_library_design = build_smart_library_design(smart_library, mutations_by_position)

        write_tsv(Path(args.ranked_hotspots), RANKED_HOTSPOT_COLUMNS, ranked_hotspots)
        write_tsv(Path(args.ranked_mutations), RANKED_MUTATION_COLUMNS, ranked_mutations)
        write_tsv(Path(args.smart_library_design), SMART_LIBRARY_COLUMNS, smart_library_design)

        summary = {
            "hotspot_rows": len(hotspots),
            "ranked_hotspots": len(ranked_hotspots),
            "candidate_mutation_rows": len(candidates),
            "mutation_ddg_rows": len(mutation_ddg),
            "ranked_mutations": len(ranked_mutations),
            "smart_library_rows": len(smart_library_design),
            "structure_pdb_embedded": bool(args.structure_pdb),
        }
        render_summary_html(
            Path(args.summary_html),
            args.project_title,
            ranked_hotspots,
            ranked_mutations,
            smart_library_design,
            Path(args.structure_viewer_html).name,
            summary,
        )
        render_viewer_html(
            Path(args.structure_viewer_html),
            args.project_title,
            ranked_hotspots,
            ranked_mutations,
            args.structure_pdb,
            args.max_viewer_residues,
        )

        log_lines.extend(
            [
                "Smart Library Report Generator completed successfully.",
                json.dumps(summary, sort_keys=True),
            ]
        )
        Path(args.run_log).write_text("\n".join(log_lines) + "\n", encoding="utf-8")
        return 0
    except Exception as exc:
        log_lines.append(f"ERROR: {exc}")
        Path(args.run_log).write_text("\n".join(log_lines) + "\n", encoding="utf-8")
        print(f"smart_library_report_generator: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
