#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import html
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import uuid
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


MAX_ASA = {
    "A": 129.0,
    "R": 274.0,
    "N": 195.0,
    "D": 193.0,
    "C": 167.0,
    "Q": 225.0,
    "E": 223.0,
    "G": 104.0,
    "H": 224.0,
    "I": 197.0,
    "L": 201.0,
    "K": 236.0,
    "M": 224.0,
    "F": 240.0,
    "P": 159.0,
    "S": 155.0,
    "T": 172.0,
    "W": 285.0,
    "Y": 263.0,
    "V": 174.0,
}


SS_CLASS = {
    "H": "helix",
    "G": "helix",
    "I": "helix",
    "E": "strand",
    "B": "strand",
    "T": "turn",
    "S": "bend",
    " ": "coil",
    "": "coil",
}


@dataclass
class ResidueFeature:
    structure_id: str
    chain_id: str
    residue_number: str
    insertion_code: str
    residue_id: str
    amino_acid: str
    dssp_secondary_structure: str
    secondary_structure_class: str
    asa: str
    relative_asa: str
    exposure_class: str
    freesasa_total_asa: str
    freesasa_relative_asa: str
    freesasa_exposure_class: str
    freesasa_main_chain_asa: str
    freesasa_side_chain_asa: str
    freesasa_polar_asa: str
    freesasa_apolar_asa: str
    phi: str
    psi: str
    avg_b_factor: str
    quality_flags: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate structure quality, secondary structure, and residue exposure with DSSP/mkdssp."
    )
    parser.add_argument("--input-structure", required=True, help="Input PDB or mmCIF structure.")
    parser.add_argument("--residue-features", required=True, help="Output residue-level TSV.")
    parser.add_argument("--quality-report", required=True, help="Output HTML quality report.")
    parser.add_argument("--dssp-output", default="structure.dssp", help="Intermediate raw DSSP output.")
    parser.add_argument("--dssp-binary", default="", help="Path or command for mkdssp.")
    parser.add_argument("--run-freesasa", action="store_true", help="Also calculate residue SASA with FreeSASA.")
    parser.add_argument("--freesasa-output", default="freesasa.rsa", help="Intermediate raw FreeSASA output.")
    parser.add_argument("--freesasa-binary", default="", help="Path or command for freesasa.")
    parser.add_argument("--buried-threshold", type=float, default=0.09)
    parser.add_argument("--exposed-threshold", type=float, default=0.36)
    parser.add_argument("--b-factor-warning", type=float, default=70.0)
    return parser.parse_args()


def command_parts(command: str) -> list[str]:
    return shlex.split(command, posix=(os.name != "nt"))


def find_dssp_binary(explicit: str = "") -> str:
    for candidate in [
        explicit,
        os.environ.get("DSSP_BINARY", ""),
        os.environ.get("MKDSSP_BINARY", ""),
        shutil.which("mkdssp") or "",
        shutil.which("dssp") or "",
    ]:
        if candidate and str(candidate).strip():
            return str(candidate).strip()
    raise RuntimeError("Could not find mkdssp. Set DSSP_BINARY or install dssp in the Galaxy job environment.")


def find_freesasa_binary(explicit: str = "") -> str:
    for candidate in [
        explicit,
        os.environ.get("FREESASA_BINARY", ""),
        os.environ.get("FREE_SASA_BINARY", ""),
        shutil.which("freesasa") or "",
    ]:
        if candidate and str(candidate).strip():
            return str(candidate).strip()
    raise RuntimeError("Could not find freesasa. Set FREESASA_BINARY or install freesasa in the Galaxy job environment.")


def detect_structure_format(path: Path) -> str:
    """Return pdb, mmcif, or unknown using file content rather than Galaxy's .dat suffix."""
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines()[:200]:
        clean_line = raw_line.lstrip("\ufeff").lstrip()
        line = clean_line.strip()
        if not line:
            continue
        if line.startswith("data_") or line.startswith("loop_") or line.startswith("_"):
            return "mmcif"
        if clean_line.startswith(("HEADER", "TITLE ", "CRYST1", "MODEL ", "ATOM  ", "HETATM", "TER", "END")):
            return "pdb"
    return "unknown"


def sanitize_pdb_text(text: str) -> str:
    """Normalize common upload artifacts while preserving PDB fixed-column records."""
    known_records = (
        "HEADER",
        "TITLE ",
        "COMPND",
        "SOURCE",
        "KEYWDS",
        "EXPDTA",
        "AUTHOR",
        "REMARK",
        "DBREF ",
        "SEQRES",
        "MODRES",
        "HET   ",
        "HETNAM",
        "FORMUL",
        "HELIX ",
        "SHEET ",
        "SSBOND",
        "LINK  ",
        "SITE  ",
        "CRYST1",
        "ORIGX1",
        "ORIGX2",
        "ORIGX3",
        "SCALE1",
        "SCALE2",
        "SCALE3",
        "MODEL ",
        "ATOM  ",
        "HETATM",
        "ANISOU",
        "TER   ",
        "ENDMDL",
        "END   ",
    )
    sanitized: list[str] = []
    for raw_line in text.splitlines():
        line = raw_line.lstrip("\ufeff").rstrip("\r\n")
        shifted = line.lstrip()
        if any(shifted.startswith(record) for record in known_records):
            line = shifted
        if line.strip():
            sanitized.append(line)
    return "\n".join(sanitized) + "\n"


def pdb_text_with_default_cryst1(path: Path) -> str:
    text = sanitize_pdb_text(path.read_text(encoding="utf-8", errors="replace"))
    lines = text.splitlines()
    if any(line.startswith("CRYST1") for line in lines[:200]):
        return text
    # mkdssp is stricter than many visualization tools. A placeholder CRYST1
    # keeps simple/model PDB files readable without changing atom coordinates.
    cryst1 = "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n"
    insert_after_records = ("HEADER", "TITLE ", "COMPND", "SOURCE", "KEYWDS", "EXPDTA", "AUTHOR", "REMARK")
    insert_at = 0
    while insert_at < len(lines) and lines[insert_at].startswith(insert_after_records):
        insert_at += 1
    return "\n".join(lines[:insert_at] + [cryst1.rstrip("\n")] + lines[insert_at:]) + "\n"


def prepare_structure_for_dssp(input_structure: Path, work_dir: Path) -> tuple[Path, str]:
    structure_format = detect_structure_format(input_structure)
    if structure_format == "mmcif":
        prepared = work_dir / "input_structure.cif"
        shutil.copyfile(input_structure, prepared)
        return prepared, structure_format
    if structure_format == "pdb":
        prepared = work_dir / "input_structure.pdb"
        prepared.write_text(pdb_text_with_default_cryst1(input_structure), encoding="utf-8")
        return prepared, structure_format
    # If the format is ambiguous, use a .pdb suffix first. Most Galaxy uploads
    # for this tool are PDB files stored internally as extensionless .dat files.
    prepared = work_dir / "input_structure.pdb"
    prepared.write_text(pdb_text_with_default_cryst1(input_structure), encoding="utf-8")
    return prepared, structure_format


def convert_mmcif_to_pdb(mmcif_path: Path, pdb_path: Path) -> Path:
    try:
        from Bio.PDB import MMCIFParser, PDBIO
    except ImportError as exc:
        raise RuntimeError(
            "DSSP could not read this mmCIF directly, and Biopython is not available for mmCIF-to-PDB fallback."
        ) from exc

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(mmcif_path.stem, str(mmcif_path))
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(pdb_path))
    pdb_path.write_text(pdb_text_with_default_cryst1(pdb_path), encoding="utf-8")
    return pdb_path


def convert_pdb_to_clean_pdb(input_pdb: Path, output_pdb: Path) -> Path:
    """Rewrite a PDB through Biopython, falling back to coordinate-line extraction.

    Some RCSB PDB downloads contain large or unusual headers that are valid
    enough for viewers but rejected by mkdssp. DSSP only needs coordinates, so
    this fallback intentionally drops most metadata.
    """
    try:
        from Bio.PDB import PDBIO, PDBParser

        parser = PDBParser(QUIET=True, PERMISSIVE=True)
        structure = parser.get_structure(input_pdb.stem, str(input_pdb))
        writer = PDBIO()
        writer.set_structure(structure)
        writer.save(str(output_pdb))
    except Exception:
        coordinate_lines = []
        for raw_line in input_pdb.read_text(encoding="utf-8", errors="replace").splitlines():
            line = raw_line.lstrip("\ufeff").lstrip()
            if line.startswith(("ATOM  ", "HETATM", "TER", "END")):
                coordinate_lines.append(line)
        if not any(line.startswith(("ATOM  ", "HETATM")) for line in coordinate_lines):
            raise RuntimeError("Could not extract any ATOM/HETATM coordinate records from the PDB file.")
        output_pdb.write_text("\n".join(coordinate_lines) + "\n", encoding="utf-8")
    output_pdb.write_text(pdb_text_with_default_cryst1(output_pdb), encoding="utf-8")
    return output_pdb


def run_dssp(dssp_binary: str, input_structure: Path, dssp_output: Path) -> tuple[list[str], str]:
    base_command = command_parts(dssp_binary)
    attempts = [
        base_command + ["--output-format", "dssp", str(input_structure), str(dssp_output)],
        base_command + [str(input_structure), str(dssp_output)],
    ]
    errors: list[str] = []
    for command in attempts:
        completed = subprocess.run(command, capture_output=True, text=True)
        if completed.returncode == 0 and dssp_output.exists() and dssp_output.stat().st_size > 0:
            return command, completed.stderr.strip()
        errors.append(
            f"Command: {' '.join(command)}\n"
            f"Exit code: {completed.returncode}\n"
            f"STDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )
    raise RuntimeError("DSSP/mkdssp failed.\n" + "\n\n".join(errors))


def run_freesasa(freesasa_binary: str, input_structure: Path, freesasa_output: Path) -> tuple[list[str], str]:
    base_command = command_parts(freesasa_binary)
    attempts = [
        # Newer FreeSASA builds can emit residue-level JSON directly.
        base_command + ["--format=json", "--output-depth=residue", str(input_structure)],
        # Older conda builds commonly support NACCESS/RSA-style residue output.
        base_command + ["--format=rsa", str(input_structure)],
        base_command + ["-f", "rsa", str(input_structure)],
    ]
    errors: list[str] = []
    for command in attempts:
        completed = subprocess.run(command, capture_output=True, text=True)
        if completed.returncode == 0 and completed.stdout.strip():
            freesasa_output.write_text(completed.stdout, encoding="utf-8")
            return command, completed.stderr.strip()
        errors.append(
            f"Command: {' '.join(command)}\n"
            f"Exit code: {completed.returncode}\n"
            f"STDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )
    raise RuntimeError("FreeSASA failed.\n" + "\n\n".join(errors))


def parse_float(text: str) -> str:
    text = text.strip()
    if not text:
        return ""
    try:
        return str(round(float(text), 3))
    except ValueError:
        return ""


def parse_int(text: str) -> str:
    text = text.strip()
    if not text:
        return ""
    try:
        return str(int(text))
    except ValueError:
        return ""


def normalize_amino_acid(aa: str) -> str:
    aa = (aa or "").strip()
    if not aa:
        return "X"
    if aa == "!":
        return "!"
    return aa.upper()


def exposure_class(relative_asa: float | None, buried_threshold: float, exposed_threshold: float) -> str:
    if relative_asa is None or math.isnan(relative_asa):
        return "unknown"
    if relative_asa < buried_threshold:
        return "buried"
    if relative_asa >= exposed_threshold:
        return "exposed"
    return "intermediate"


def split_residue_number(value: object) -> tuple[str, str]:
    text = str(value).strip()
    match = re.match(r"^(-?\d+)([A-Za-z]?)$", text)
    if match:
        return match.group(1), match.group(2)
    return text, ""


def as_float(value: object) -> float | None:
    try:
        if value is None:
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def get_nested_area(mapping: dict, *keys: str) -> float | None:
    current: object = mapping
    for key in keys:
        if not isinstance(current, dict) or key not in current:
            return None
        current = current[key]
    return as_float(current)


def parse_freesasa_json(path: Path) -> dict[tuple[str, str, str], dict[str, str]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    residue_features: dict[tuple[str, str, str], dict[str, str]] = {}
    result_blocks = data.get("results", [])
    for result in result_blocks:
        for structure in result.get("structures", []):
            for chain in structure.get("chains", []):
                chain_id = str(chain.get("label", "") or chain.get("name", "") or ".")
                for residue in chain.get("residues", []):
                    residue_number, insertion_code = split_residue_number(residue.get("number", ""))
                    area = residue.get("area", {}) if isinstance(residue.get("area", {}), dict) else {}
                    relative_area = residue.get("relative-area", {}) if isinstance(residue.get("relative-area", {}), dict) else {}
                    total = get_nested_area(area, "total")
                    relative_total_percent = get_nested_area(relative_area, "total")
                    relative_total = relative_total_percent / 100.0 if relative_total_percent is not None else None
                    main_chain = get_nested_area(area, "main-chain")
                    side_chain = get_nested_area(area, "side-chain")
                    polar = get_nested_area(area, "polar")
                    apolar = get_nested_area(area, "apolar")
                    residue_features[(chain_id, residue_number, insertion_code)] = {
                        "freesasa_total_asa": str(round(total, 3)) if total is not None else "",
                        "freesasa_relative_asa": str(round(relative_total, 4)) if relative_total is not None else "",
                        "freesasa_main_chain_asa": str(round(main_chain, 3)) if main_chain is not None else "",
                        "freesasa_side_chain_asa": str(round(side_chain, 3)) if side_chain is not None else "",
                        "freesasa_polar_asa": str(round(polar, 3)) if polar is not None else "",
                        "freesasa_apolar_asa": str(round(apolar, 3)) if apolar is not None else "",
                    }
    return residue_features


def parse_freesasa_rsa(path: Path) -> dict[tuple[str, str, str], dict[str, str]]:
    """Parse NACCESS/RSA-style FreeSASA residue output.

    Expected RES columns are:
    RES, resname, chain, resnum, total_abs, total_rel,
    side_abs, side_rel, main_abs, main_rel, apolar_abs, apolar_rel,
    polar_abs, polar_rel.
    Relative values are percentages in RSA output and are converted to 0-1.
    """
    residue_features: dict[tuple[str, str, str], dict[str, str]] = {}
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if not line.startswith("RES"):
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        chain_id = parts[2] if len(parts) > 2 else "."
        residue_number, insertion_code = split_residue_number(parts[3] if len(parts) > 3 else "")
        total = as_float(parts[4] if len(parts) > 4 else None)
        relative_total_percent = as_float(parts[5] if len(parts) > 5 else None)
        side_chain = as_float(parts[6] if len(parts) > 6 else None)
        main_chain = as_float(parts[8] if len(parts) > 8 else None)
        apolar = as_float(parts[10] if len(parts) > 10 else None)
        polar = as_float(parts[12] if len(parts) > 12 else None)
        relative_total = relative_total_percent / 100.0 if relative_total_percent is not None else None
        residue_features[(chain_id, residue_number, insertion_code)] = {
            "freesasa_total_asa": str(round(total, 3)) if total is not None else "",
            "freesasa_relative_asa": str(round(relative_total, 4)) if relative_total is not None else "",
            "freesasa_main_chain_asa": str(round(main_chain, 3)) if main_chain is not None else "",
            "freesasa_side_chain_asa": str(round(side_chain, 3)) if side_chain is not None else "",
            "freesasa_polar_asa": str(round(polar, 3)) if polar is not None else "",
            "freesasa_apolar_asa": str(round(apolar, 3)) if apolar is not None else "",
        }
    return residue_features


def parse_freesasa_output(path: Path) -> dict[tuple[str, str, str], dict[str, str]]:
    text = path.read_text(encoding="utf-8", errors="replace").lstrip()
    if text.startswith("{"):
        return parse_freesasa_json(path)
    return parse_freesasa_rsa(path)


def apply_freesasa_features(
    features: list[ResidueFeature],
    freesasa_features: dict[tuple[str, str, str], dict[str, str]],
    buried_threshold: float,
    exposed_threshold: float,
) -> int:
    matched = 0
    for feature in features:
        values = freesasa_features.get((feature.chain_id, feature.residue_number, feature.insertion_code))
        if values is None:
            continue
        matched += 1
        feature.freesasa_total_asa = values.get("freesasa_total_asa", "")
        feature.freesasa_relative_asa = values.get("freesasa_relative_asa", "")
        feature.freesasa_main_chain_asa = values.get("freesasa_main_chain_asa", "")
        feature.freesasa_side_chain_asa = values.get("freesasa_side_chain_asa", "")
        feature.freesasa_polar_asa = values.get("freesasa_polar_asa", "")
        feature.freesasa_apolar_asa = values.get("freesasa_apolar_asa", "")
        relative = as_float(feature.freesasa_relative_asa)
        feature.freesasa_exposure_class = exposure_class(relative, buried_threshold, exposed_threshold)
    return matched


def parse_dssp(
    dssp_path: Path,
    structure_id: str,
    b_factors: dict[tuple[str, str, str], float],
    buried_threshold: float,
    exposed_threshold: float,
    b_factor_warning: float,
) -> list[ResidueFeature]:
    features: list[ResidueFeature] = []
    in_table = False
    for raw_line in dssp_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if raw_line.lstrip().startswith("#") and "RESIDUE" in raw_line and "AA" in raw_line:
            in_table = True
            continue
        if not in_table or len(raw_line) < 17:
            continue
        if raw_line[13:14] == "!":
            continue

        residue_number = raw_line[5:10].strip()
        insertion_code = raw_line[10:11].strip()
        chain_id = raw_line[11:12].strip() or "."
        amino_acid = normalize_amino_acid(raw_line[13:14])
        if not residue_number or amino_acid == "!":
            continue

        dssp_ss = raw_line[16:17].strip() or "C"
        ss_class = SS_CLASS.get(dssp_ss, "coil")
        asa_text = parse_int(raw_line[34:38] if len(raw_line) >= 38 else "")
        asa_value: float | None = float(asa_text) if asa_text else None
        max_asa = MAX_ASA.get(amino_acid)
        if asa_value is not None and max_asa:
            rel_value = min(1.5, max(0.0, asa_value / max_asa))
            rel_text = str(round(rel_value, 4))
        else:
            rel_value = None
            rel_text = ""
        phi = parse_float(raw_line[103:109] if len(raw_line) >= 109 else "")
        psi = parse_float(raw_line[109:115] if len(raw_line) >= 115 else "")

        b_key = (chain_id, residue_number, insertion_code)
        avg_b = b_factors.get(b_key)
        flags: list[str] = []
        if not asa_text:
            flags.append("missing_asa")
        if amino_acid not in MAX_ASA:
            flags.append("noncanonical_or_unknown_residue")
        if avg_b is not None and avg_b >= b_factor_warning:
            flags.append("high_b_factor")

        residue_id = f"{chain_id}:{residue_number}{insertion_code}"
        features.append(
            ResidueFeature(
                structure_id=structure_id,
                chain_id=chain_id,
                residue_number=residue_number,
                insertion_code=insertion_code,
                residue_id=residue_id,
                amino_acid=amino_acid,
                dssp_secondary_structure=dssp_ss,
                secondary_structure_class=ss_class,
                asa=asa_text,
                relative_asa=rel_text,
                exposure_class=exposure_class(rel_value, buried_threshold, exposed_threshold),
                freesasa_total_asa="",
                freesasa_relative_asa="",
                freesasa_exposure_class="not_run",
                freesasa_main_chain_asa="",
                freesasa_side_chain_asa="",
                freesasa_polar_asa="",
                freesasa_apolar_asa="",
                phi=phi,
                psi=psi,
                avg_b_factor=str(round(avg_b, 3)) if avg_b is not None else "",
                quality_flags=";".join(flags) if flags else "ok",
            )
        )
    if not features:
        raise RuntimeError("DSSP output was produced but no residue rows could be parsed.")
    return features


def parse_pdb_b_factors(path: Path) -> dict[tuple[str, str, str], float]:
    values: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except UnicodeDecodeError:
        return {}
    for line in lines:
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 66:
            continue
        chain_id = line[21:22].strip() or "."
        residue_number = line[22:26].strip()
        insertion_code = line[26:27].strip()
        try:
            b_factor = float(line[60:66])
        except ValueError:
            continue
        values[(chain_id, residue_number, insertion_code)].append(b_factor)
    return {key: sum(v) / len(v) for key, v in values.items() if v}


def write_features(features: list[ResidueFeature], path: Path) -> None:
    fieldnames = list(ResidueFeature.__dataclass_fields__.keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for feature in features:
            writer.writerow(feature.__dict__)


def pct(value: int, total: int) -> str:
    return f"{(100.0 * value / total):.1f}%" if total else "0.0%"


def write_report(
    features: list[ResidueFeature],
    path: Path,
    input_structure: Path,
    dssp_command: list[str],
    dssp_stderr: str,
    freesasa_command: list[str] | None = None,
    freesasa_stderr: str = "",
    freesasa_matched: int = 0,
) -> None:
    total = len(features)
    ss_counts = Counter(feature.secondary_structure_class for feature in features)
    exposure_counts = Counter(feature.exposure_class for feature in features)
    freesasa_exposure_counts = Counter(feature.freesasa_exposure_class for feature in features)
    flag_counts = Counter()
    for feature in features:
        for flag in feature.quality_flags.split(";"):
            flag_counts[flag] += 1

    rows = "\n".join(
        f"<tr><td>{html.escape(feature.residue_id)}</td>"
        f"<td>{html.escape(feature.amino_acid)}</td>"
        f"<td>{html.escape(feature.secondary_structure_class)}</td>"
        f"<td>{html.escape(feature.relative_asa)}</td>"
        f"<td>{html.escape(feature.exposure_class)}</td>"
        f"<td>{html.escape(feature.freesasa_relative_asa)}</td>"
        f"<td>{html.escape(feature.freesasa_exposure_class)}</td>"
        f"<td>{html.escape(feature.avg_b_factor)}</td>"
        f"<td>{html.escape(feature.quality_flags)}</td></tr>"
        for feature in features[:50]
    )
    ss_items = "".join(
        f"<li>{html.escape(name)}: {count} ({pct(count, total)})</li>" for name, count in sorted(ss_counts.items())
    )
    exposure_items = "".join(
        f"<li>{html.escape(name)}: {count} ({pct(count, total)})</li>"
        for name, count in sorted(exposure_counts.items())
    )
    freesasa_exposure_items = "".join(
        f"<li>{html.escape(name)}: {count} ({pct(count, total)})</li>"
        for name, count in sorted(freesasa_exposure_counts.items())
    )
    flag_items = "".join(
        f"<li>{html.escape(name)}: {count} ({pct(count, total)})</li>" for name, count in sorted(flag_counts.items())
    )
    dssp_stderr_html = f"<pre>{html.escape(dssp_stderr)}</pre>" if dssp_stderr else "<p>No DSSP warnings reported.</p>"
    if freesasa_command:
        freesasa_html = (
            f"<p><strong>FreeSASA command:</strong> <code>{html.escape(' '.join(freesasa_command))}</code></p>"
            f"<p><strong>Matched residues:</strong> {freesasa_matched} / {total}</p>"
            + (f"<pre>{html.escape(freesasa_stderr)}</pre>" if freesasa_stderr else "<p>No FreeSASA warnings reported.</p>")
        )
    else:
        freesasa_html = "<p>FreeSASA was not run for this job.</p>"

    path.write_text(
        f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>Structure Quality Annotator Report</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2rem; color: #1f2937; }}
    .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 1rem; }}
    .card {{ border: 1px solid #d1d5db; border-radius: 10px; padding: 1rem; background: #f9fafb; }}
    table {{ border-collapse: collapse; width: 100%; margin-top: 1rem; }}
    th, td {{ border: 1px solid #d1d5db; padding: 0.35rem 0.5rem; text-align: left; }}
    th {{ background: #e5e7eb; }}
    code, pre {{ background: #f3f4f6; padding: 0.2rem 0.35rem; border-radius: 4px; }}
  </style>
</head>
<body>
  <h1>Structure Quality Annotator Report</h1>
  <p><strong>Input:</strong> {html.escape(str(input_structure))}</p>
  <p><strong>DSSP command:</strong> <code>{html.escape(" ".join(dssp_command))}</code></p>
  <div class="grid">
    <div class="card"><h2>Residues</h2><p>{total}</p></div>
    <div class="card"><h2>Secondary Structure</h2><ul>{ss_items}</ul></div>
    <div class="card"><h2>Exposure</h2><ul>{exposure_items}</ul></div>
    <div class="card"><h2>FreeSASA Exposure</h2><ul>{freesasa_exposure_items}</ul></div>
    <div class="card"><h2>Quality Flags</h2><ul>{flag_items}</ul></div>
  </div>
  <h2>DSSP Messages</h2>
  {dssp_stderr_html}
  <h2>FreeSASA</h2>
  {freesasa_html}
  <h2>First 50 Residues</h2>
  <table>
    <thead><tr><th>Residue</th><th>AA</th><th>SS class</th><th>DSSP relative ASA</th><th>DSSP exposure</th><th>FreeSASA relative ASA</th><th>FreeSASA exposure</th><th>Avg B-factor</th><th>Flags</th></tr></thead>
    <tbody>{rows}</tbody>
  </table>
  <p><em>Use the TSV output for downstream hotspot ranking. DSSP-derived exposure is useful for avoiding buried destabilizing mutations unless stability design is intended.</em></p>
</body>
</html>
""",
        encoding="utf-8",
    )


def main() -> int:
    args = parse_args()
    work_dir = Path.cwd() / f"structure_quality_tmp_{uuid.uuid4().hex}"
    try:
        input_structure = Path(args.input_structure)
        work_dir.mkdir(parents=True, exist_ok=False)
        dssp_input, detected_format = prepare_structure_for_dssp(input_structure, work_dir)
        dssp_output = Path(args.dssp_output)
        dssp_binary = find_dssp_binary(args.dssp_binary)
        try:
            dssp_command, dssp_stderr = run_dssp(dssp_binary, dssp_input, dssp_output)
        except RuntimeError as first_error:
            fallback_pdb = work_dir / "input_structure_coordinates_only.pdb"
            try:
                if detected_format == "mmcif":
                    convert_mmcif_to_pdb(dssp_input, fallback_pdb)
                    fallback_note = "DSSP failed on the original mmCIF, so the wrapper converted coordinates to PDB and retried."
                else:
                    convert_pdb_to_clean_pdb(dssp_input, fallback_pdb)
                    fallback_note = "DSSP failed on the original PDB, so the wrapper rewrote coordinate records to a clean PDB and retried."
                dssp_command, dssp_stderr = run_dssp(dssp_binary, fallback_pdb, dssp_output)
                dssp_stderr = (
                    f"{fallback_note}\n"
                    f"Original DSSP error:\n{first_error}\n\n{dssp_stderr}"
                ).strip()
                dssp_input = fallback_pdb
            except Exception as fallback_error:
                raise RuntimeError(
                    f"{first_error}\n\nCoordinate-only PDB fallback also failed:\n{fallback_error}"
                ) from fallback_error
        b_factors = parse_pdb_b_factors(dssp_input)
        features = parse_dssp(
            dssp_output,
            structure_id=input_structure.stem,
            b_factors=b_factors,
            buried_threshold=args.buried_threshold,
            exposed_threshold=args.exposed_threshold,
            b_factor_warning=args.b_factor_warning,
        )
        freesasa_command = None
        freesasa_stderr = ""
        freesasa_matched = 0
        if args.run_freesasa:
            freesasa_output = Path(args.freesasa_output)
            freesasa_binary = find_freesasa_binary(args.freesasa_binary)
            freesasa_command, freesasa_stderr = run_freesasa(freesasa_binary, dssp_input, freesasa_output)
            freesasa_features = parse_freesasa_output(freesasa_output)
            freesasa_matched = apply_freesasa_features(
                features,
                freesasa_features,
                buried_threshold=args.buried_threshold,
                exposed_threshold=args.exposed_threshold,
            )
            if freesasa_matched == 0:
                raise RuntimeError("FreeSASA completed, but no residue-level SASA rows matched DSSP residues.")
        write_features(features, Path(args.residue_features))
        write_report(
            features,
            Path(args.quality_report),
            input_structure,
            dssp_command,
            dssp_stderr,
            freesasa_command=freesasa_command,
            freesasa_stderr=freesasa_stderr,
            freesasa_matched=freesasa_matched,
        )
        return 0
    except Exception as exc:
        print(f"structure_quality_annotator: {exc}", file=sys.stderr)
        return 1
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
