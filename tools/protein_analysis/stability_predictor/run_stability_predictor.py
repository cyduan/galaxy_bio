#!/usr/bin/env python3
"""Run FoldX BuildModel for HotSpot candidate mutations and classify ddG."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path


AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

MUTATION_DDG_COLUMNS = [
    "mutation_id",
    "position",
    "chain",
    "wt_aa",
    "mutant_aa",
    "foldx_mutation",
    "ddg_kcal_per_mol",
    "stability_class",
    "keep_for_hotspot",
    "filter_reason",
    "hotspot_rank",
    "hotspot_score",
    "hotspot_class",
    "mutation_type",
    "candidate_reason",
    "msa_frequency",
    "source",
    "foldx_source_row",
    "warning",
]


@dataclass
class CandidateMutation:
    mutation_id: str
    position: int
    chain: str
    wt_aa: str
    mutant_aa: str
    foldx_mutation: str
    hotspot_rank: str = ""
    hotspot_score: str = ""
    hotspot_class: str = ""
    mutation_type: str = ""
    candidate_reason: str = ""
    msa_frequency: str = ""
    source: str = ""


@dataclass
class DdgResult:
    foldx_mutation: str
    ddg: float | None
    source_row: str
    warning: str = ""


class ToolError(RuntimeError):
    """User-facing CLI error."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Predict mutation stability impact with FoldX BuildModel and format outputs for HotSpot workflows."
    )
    parser.add_argument("--candidate-mutations", required=True, help="Tool 7 candidate_mutations.tsv.")
    parser.add_argument("--pdb", required=True, help="Input PDB structure.")
    parser.add_argument("--mutation-ddg", required=True, help="HotSpot mutation_ddg.tsv output.")
    parser.add_argument("--filtered-mutations", required=True, help="stability_filtered_mutations.tsv output.")
    parser.add_argument("--raw-ddg-output", required=True, help="Preserved raw FoldX ddG output.")
    parser.add_argument("--repaired-pdb", required=True, help="Optional repaired PDB output path.")
    parser.add_argument("--mutant-models-dir", required=True, help="Directory for mutant PDB models.")
    parser.add_argument("--archive", required=True, help="Complete FoldX outputs ZIP.")
    parser.add_argument("--individual-list", required=True, help="FoldX individual_list.txt output.")
    parser.add_argument("--summary-json", required=True, help="Run summary JSON.")
    parser.add_argument("--run-log", required=True, help="Run log.")
    parser.add_argument("--foldx-command", default=os.environ.get("FOLDX_BINARY", "foldx"))
    parser.add_argument("--repair", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--number-of-runs", type=int, default=5)
    parser.add_argument("--pH", default="7.0")
    parser.add_argument("--temperature", default="298")
    parser.add_argument("--ionic-strength", default="0.05")
    parser.add_argument("--max-ddg-for-filter", type=float, default=1.0)
    parser.add_argument("--include-caution", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--max-candidates", type=int, default=0, help="0 means use all candidate rows.")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ToolError(f"{path} is empty or does not contain a TSV header.")
        return [{key: (value or "") for key, value in row.items()} for row in reader]


def parse_candidate_mutations(path: Path, max_candidates: int) -> list[CandidateMutation]:
    rows = read_tsv(path)
    candidates: list[CandidateMutation] = []
    seen: set[str] = set()
    for row in rows:
        position = parse_int(row.get("position", ""))
        chain = (row.get("chain") or parse_chain_from_residue_id(row.get("residue_id", "")) or "A").strip()
        wt_aa = normalize_aa(row.get("wt_aa", ""))
        mutant_aa = normalize_aa(row.get("suggested_aa", "") or row.get("mutant_aa", ""))
        if position is None or not wt_aa or not mutant_aa or wt_aa == mutant_aa:
            continue
        foldx_mutation = f"{wt_aa}{chain}{position}{mutant_aa}"
        if foldx_mutation in seen:
            continue
        seen.add(foldx_mutation)
        mutation_id = row.get("mutation_id", "") or f"{chain}{position}{wt_aa}_to_{mutant_aa}"
        candidates.append(
            CandidateMutation(
                mutation_id=mutation_id,
                position=position,
                chain=chain,
                wt_aa=wt_aa,
                mutant_aa=mutant_aa,
                foldx_mutation=foldx_mutation,
                hotspot_rank=row.get("hotspot_rank", ""),
                hotspot_score=row.get("hotspot_score", ""),
                hotspot_class=row.get("hotspot_class", ""),
                mutation_type=row.get("mutation_type", ""),
                candidate_reason=row.get("reason", ""),
                msa_frequency=row.get("msa_frequency", ""),
                source=row.get("source", ""),
            )
        )
        if max_candidates > 0 and len(candidates) >= max_candidates:
            break
    if not candidates:
        raise ToolError("No usable candidate mutation rows were found.")
    return candidates


def parse_int(value: str) -> int | None:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def normalize_aa(value: str) -> str:
    value = (value or "").strip().upper()
    if len(value) == 1 and value in AMINO_ACIDS:
        return value
    return ""


def parse_chain_from_residue_id(value: str) -> str:
    match = re.match(r"\s*([A-Za-z0-9])\s*:", value or "")
    return match.group(1) if match else ""


def resolve_command(command: str) -> list[str]:
    if command == "mock-foldx":
        return ["mock-foldx"]
    parts = shlex.split(command, posix=(os.name != "nt"))
    if not parts:
        raise ToolError("FoldX command is empty.")
    return parts


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("._") or "foldx_input"


def copy_input_pdb(source: Path, work_dir: Path) -> Path:
    target = work_dir / f"{safe_name(source.stem)}.pdb"
    shutil.copyfile(source, target)
    return target


def build_base_options(args: argparse.Namespace) -> list[str]:
    return [
        f"--pH={args.pH}",
        f"--temperature={args.temperature}",
        f"--ionStrength={args.ionic_strength}",
    ]


def run_command(command: list[str], cwd: Path, log_lines: list[str]) -> subprocess.CompletedProcess[str]:
    log_lines.append("$ " + " ".join(command))
    if command[0] == "mock-foldx":
        completed = run_mock(command, cwd)
    else:
        try:
            completed = subprocess.run(command, cwd=cwd, capture_output=True, text=True)
        except FileNotFoundError as exc:
            raise ToolError(
                f"FoldX executable not found: {command[0]}. Set FOLDX_BINARY or --foldx-command."
            ) from exc
    log_lines.append(f"[returncode] {completed.returncode}")
    if completed.stdout:
        log_lines.append("[stdout]\n" + completed.stdout.rstrip())
    if completed.stderr:
        log_lines.append("[stderr]\n" + completed.stderr.rstrip())
    return completed


def run_mock(command: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    command_text = " ".join(command)
    pdb_name = next((item.split("=", 1)[1] for item in command if item.startswith("--pdb=")), "input.pdb")
    pdb_path = cwd / pdb_name
    pdb_stem = Path(pdb_name).stem
    if "--command=RepairPDB" in command_text:
        repaired = cwd / f"{pdb_stem}_Repair.pdb"
        shutil.copyfile(pdb_path, repaired)
        return subprocess.CompletedProcess(command, 0, f"RepairPDB wrote {repaired.name}\n", "")
    if "--command=BuildModel" in command_text:
        mutation_file = next(
            (cwd / item.split("=", 1)[1] for item in command if item.startswith("--mutant-file=")),
            cwd / "individual_list.txt",
        )
        mutation_lines = [line.strip().rstrip(";") for line in mutation_file.read_text(encoding="utf-8").splitlines()]
        rows = ["Pdb\tMutations\tddG"]
        for index, mutation in enumerate(mutation_lines, start=1):
            if not mutation:
                continue
            ddg = mock_ddg_for_mutation(mutation)
            rows.append(f"{pdb_name}\t{mutation}\t{ddg:.3f}")
            (cwd / f"{pdb_stem}_{index}.pdb").write_text(
                "\n".join(
                    [
                        "HEADER    MOCK FOLDX MUTANT",
                        f"REMARK    MUTATION {mutation}",
                        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N",
                        "TER",
                        "END",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
        (cwd / "Dif_buildmodel.fxout").write_text("\n".join(rows) + "\n", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "BuildModel completed\n", "")
    return subprocess.CompletedProcess(command, 1, "", "Unknown mock FoldX command\n")


def mock_ddg_for_mutation(mutation: str) -> float:
    mutant = mutation[-1]
    if mutant in {"A", "V", "I", "L"}:
        return 0.5
    if mutant in {"D", "E", "K", "R"}:
        return 1.5
    if mutant == "P":
        return 2.5
    return -0.2


def run_repair(args: argparse.Namespace, foldx: list[str], input_pdb: Path, work_dir: Path, log_lines: list[str]) -> Path:
    command = [
        *foldx,
        "--command=RepairPDB",
        f"--pdb={input_pdb.name}",
        *build_base_options(args),
    ]
    completed = run_command(command, work_dir, log_lines)
    if completed.returncode != 0:
        raise ToolError(foldx_failure_message("RepairPDB", completed))
    repaired = find_repaired_pdb(work_dir, input_pdb)
    if repaired is None:
        raise ToolError("FoldX RepairPDB completed but no *_Repair.pdb file was found.")
    return repaired


def run_buildmodel(
    args: argparse.Namespace,
    foldx: list[str],
    active_pdb: Path,
    mutation_file: Path,
    work_dir: Path,
    log_lines: list[str],
) -> Path:
    command = [
        *foldx,
        "--command=BuildModel",
        f"--pdb={active_pdb.name}",
        f"--mutant-file={mutation_file.name}",
        f"--numberOfRuns={args.number_of_runs}",
        "--output-file=buildmodel",
        *build_base_options(args),
    ]
    completed = run_command(command, work_dir, log_lines)
    if completed.returncode != 0:
        raise ToolError(foldx_failure_message("BuildModel", completed))
    result = first_matching_file(work_dir, ["Dif_*.fxout", "*Dif*.fxout", "Average_*.fxout", "*.fxout"])
    if result is None:
        raise ToolError("FoldX BuildModel completed but no .fxout file was found.")
    return result


def foldx_failure_message(step: str, completed: subprocess.CompletedProcess[str]) -> str:
    parts = [f"FoldX {step} failed with status {completed.returncode}."]
    if completed.stdout:
        parts.append("FoldX stdout:\n" + completed.stdout.strip())
    if completed.stderr:
        parts.append("FoldX stderr:\n" + completed.stderr.strip())
    return "\n\n".join(parts)


def find_repaired_pdb(work_dir: Path, input_pdb: Path) -> Path | None:
    expected = work_dir / f"{input_pdb.stem}_Repair.pdb"
    if expected.exists():
        return expected
    matches = sorted(work_dir.glob("*_Repair.pdb"))
    return matches[0] if matches else None


def first_matching_file(work_dir: Path, patterns: list[str]) -> Path | None:
    matches: list[Path] = []
    for pattern in patterns:
        matches.extend(sorted(work_dir.glob(pattern)))
    return matches[0] if matches else None


def write_individual_list(candidates: list[CandidateMutation], path: Path) -> None:
    path.write_text(";\n".join(candidate.foldx_mutation for candidate in candidates) + ";\n", encoding="utf-8")


def parse_foldx_ddg(path: Path, candidates: list[CandidateMutation]) -> dict[str, DdgResult]:
    rows = read_foldx_rows(path)
    if not rows:
        raise ToolError(f"FoldX ddG output is empty: {path}")

    header = rows[0]
    data_rows = rows[1:] if looks_like_header(header) else rows
    colmap = {normalize_header(name): index for index, name in enumerate(header)} if rows[1:] else {}
    mutation_col = first_existing_col(colmap, ["mutations", "mutation", "mutant", "mutationlist"])
    ddg_col = first_existing_col(colmap, ["ddg", "delta_delta_g", "deltaenergy", "difference", "totalenergy"])

    results: dict[str, DdgResult] = {}
    candidate_by_index = {index: candidate for index, candidate in enumerate(candidates)}
    candidate_by_mutation = {candidate.foldx_mutation: candidate for candidate in candidates}

    for index, row in enumerate(data_rows):
        mutation = ""
        if mutation_col is not None and mutation_col < len(row):
            mutation = normalize_foldx_mutation(row[mutation_col])
        if not mutation and index in candidate_by_index:
            mutation = candidate_by_index[index].foldx_mutation
        if mutation not in candidate_by_mutation:
            extracted = extract_mutation_from_row(row, candidate_by_mutation)
            if extracted:
                mutation = extracted
        if not mutation:
            continue
        ddg = extract_ddg(row, ddg_col)
        results[mutation] = DdgResult(
            foldx_mutation=mutation,
            ddg=ddg,
            source_row="\t".join(row),
            warning="" if ddg is not None else "Could not identify numeric ddG in FoldX row.",
        )
    return results


def read_foldx_rows(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            rows.append(re.split(r"\t+|\s{2,}|\s+", line))
    return rows


def looks_like_header(row: list[str]) -> bool:
    return any(not looks_numeric(value) for value in row[1:])


def normalize_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def first_existing_col(colmap: dict[str, int], names: list[str]) -> int | None:
    for name in names:
        normalized = normalize_header(name)
        if normalized in colmap:
            return colmap[normalized]
    return None


def normalize_foldx_mutation(value: str) -> str:
    return value.strip().rstrip(";").replace(",", "")


def extract_mutation_from_row(row: list[str], candidate_by_mutation: dict[str, CandidateMutation]) -> str:
    row_text = " ".join(row)
    for mutation in candidate_by_mutation:
        if mutation in row_text:
            return mutation
    return ""


def extract_ddg(row: list[str], ddg_col: int | None) -> float | None:
    if ddg_col is not None and ddg_col < len(row) and looks_numeric(row[ddg_col]):
        return float(row[ddg_col])
    numeric_values = [float(value) for value in row if looks_numeric(value)]
    if not numeric_values:
        return None
    return numeric_values[-1]


def looks_numeric(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


def classify_ddg(ddg: float | None) -> tuple[str, str]:
    if ddg is None:
        return "unknown", "no_ddg"
    if ddg <= 0:
        return "priority", "ddg_le_0"
    if ddg <= 1:
        return "acceptable", "ddg_0_to_1"
    if ddg <= 2:
        return "caution", "ddg_1_to_2"
    return "high_risk", "ddg_gt_2"


def write_outputs(
    candidates: list[CandidateMutation],
    results: dict[str, DdgResult],
    args: argparse.Namespace,
) -> tuple[int, int]:
    mutation_rows: list[dict[str, str]] = []
    filtered_rows: list[dict[str, str]] = []
    for candidate in candidates:
        result = results.get(candidate.foldx_mutation)
        ddg = result.ddg if result else None
        stability_class, reason = classify_ddg(ddg)
        keep = should_keep(ddg, stability_class, args)
        warning = result.warning if result else "FoldX result row was not found for this mutation."
        row = {
            "mutation_id": candidate.mutation_id,
            "position": str(candidate.position),
            "chain": candidate.chain,
            "wt_aa": candidate.wt_aa,
            "mutant_aa": candidate.mutant_aa,
            "foldx_mutation": candidate.foldx_mutation,
            "ddg_kcal_per_mol": format_float(ddg),
            "stability_class": stability_class,
            "keep_for_hotspot": "yes" if keep else "no",
            "filter_reason": reason,
            "hotspot_rank": candidate.hotspot_rank,
            "hotspot_score": candidate.hotspot_score,
            "hotspot_class": candidate.hotspot_class,
            "mutation_type": candidate.mutation_type,
            "candidate_reason": candidate.candidate_reason,
            "msa_frequency": candidate.msa_frequency,
            "source": candidate.source,
            "foldx_source_row": result.source_row if result else "",
            "warning": warning,
        }
        mutation_rows.append(row)
        if keep:
            filtered_rows.append(row)
    write_tsv(Path(args.mutation_ddg), MUTATION_DDG_COLUMNS, mutation_rows)
    write_tsv(Path(args.filtered_mutations), MUTATION_DDG_COLUMNS, filtered_rows)
    return len(mutation_rows), len(filtered_rows)


def should_keep(ddg: float | None, stability_class: str, args: argparse.Namespace) -> bool:
    if ddg is None:
        return False
    if ddg <= args.max_ddg_for_filter:
        return True
    return args.include_caution and stability_class == "caution"


def format_float(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.6g}"


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def collect_mutant_models(work_dir: Path, destination: Path, excluded_stems: set[str]) -> list[str]:
    destination.mkdir(parents=True, exist_ok=True)
    copied: list[str] = []
    for path in sorted(work_dir.glob("*.pdb")):
        if path.name.endswith("_Repair.pdb"):
            continue
        if path.stem in excluded_stems:
            continue
        target = destination / path.name
        shutil.copyfile(path, target)
        copied.append(str(target))
    return copied


def make_archive(work_dir: Path, archive: Path) -> None:
    archive.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as zip_handle:
        for path in sorted(work_dir.rglob("*")):
            if path.is_file():
                zip_handle.write(path, path.relative_to(work_dir))


def main() -> int:
    args = parse_args()
    work_dir = Path.cwd() / "stability_predictor_foldx_work"
    log_lines = [
        "HotSpot Stability Predictor run log",
        f"candidate_mutations={args.candidate_mutations}",
        f"input_pdb={args.pdb}",
        "",
    ]
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True)

    try:
        candidates = parse_candidate_mutations(Path(args.candidate_mutations), args.max_candidates)
        foldx = resolve_command(args.foldx_command)
        copied_pdb = copy_input_pdb(Path(args.pdb), work_dir)
        active_pdb = copied_pdb
        repaired_pdb = ""
        mutation_file = work_dir / "individual_list.txt"
        write_individual_list(candidates, mutation_file)
        shutil.copyfile(mutation_file, args.individual_list)

        if args.repair:
            active_pdb = run_repair(args, foldx, copied_pdb, work_dir, log_lines)
            repaired_pdb = str(active_pdb)
            shutil.copyfile(active_pdb, args.repaired_pdb)
        else:
            Path(args.repaired_pdb).write_text("", encoding="utf-8")

        raw_ddg = run_buildmodel(args, foldx, active_pdb, mutation_file, work_dir, log_lines)
        shutil.copyfile(raw_ddg, args.raw_ddg_output)
        results = parse_foldx_ddg(raw_ddg, candidates)
        total_rows, filtered_rows = write_outputs(candidates, results, args)
        mutant_models = collect_mutant_models(
            work_dir,
            Path(args.mutant_models_dir),
            excluded_stems={copied_pdb.stem, active_pdb.stem},
        )
        make_archive(work_dir, Path(args.archive))

        summary = {
            "tool": "HotSpot Stability Predictor",
            "foldx_command": args.foldx_command,
            "status": "success",
            "repair": args.repair,
            "input_pdb": args.pdb,
            "active_pdb": str(active_pdb),
            "repaired_pdb": repaired_pdb,
            "candidate_count": len(candidates),
            "mutation_ddg_rows": total_rows,
            "filtered_mutation_rows": filtered_rows,
            "mutant_models": mutant_models,
            "ddg_rules": {
                "priority": "ddG <= 0 kcal/mol",
                "acceptable": "0 < ddG <= 1 kcal/mol",
                "caution": "1 < ddG <= 2 kcal/mol",
                "high_risk": "ddG > 2 kcal/mol",
            },
            "max_ddg_for_filter": args.max_ddg_for_filter,
            "include_caution": args.include_caution,
            "number_of_runs": args.number_of_runs,
            "pH": args.pH,
            "temperature": args.temperature,
            "ionic_strength": args.ionic_strength,
        }
        Path(args.summary_json).write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
        log_lines.extend(
            [
                "",
                "Status: success",
                f"candidate_count={len(candidates)}",
                f"mutation_ddg_rows={total_rows}",
                f"filtered_mutation_rows={filtered_rows}",
                f"mutant_models={len(mutant_models)}",
            ]
        )
        Path(args.run_log).write_text("\n".join(log_lines) + "\n", encoding="utf-8")
        return 0
    except Exception as exc:
        log_lines.extend(["", "Status: failed", str(exc)])
        Path(args.run_log).write_text("\n".join(log_lines) + "\n", encoding="utf-8")
        print(f"stability_predictor: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
