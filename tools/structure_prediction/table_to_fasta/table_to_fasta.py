#!/usr/bin/env python

from __future__ import annotations

import argparse
import csv
import re
import sys
import zipfile
from pathlib import Path
from typing import Iterable, NoReturn
from xml.etree import ElementTree


XML_NS = {
    "main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "rel": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "pkgrel": "http://schemas.openxmlformats.org/package/2006/relationships",
}


def stop_err(message: str) -> NoReturn:
    raise SystemExit(f"ERROR: {message}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert CSV/XLSX sequence tables into FASTA.")
    parser.add_argument("--input", required=True, dest="input_path")
    parser.add_argument("--input-format", required=True)
    parser.add_argument("--name-column", default="name")
    parser.add_argument("--sequence-column", default="sequence")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def iter_non_empty_rows(rows: Iterable[list[str]]) -> Iterable[list[str]]:
    for row in rows:
        if any(cell.strip() for cell in row):
            yield row


def load_delimited_rows(path: Path) -> list[list[str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        except csv.Error:
            dialect = csv.excel
        reader = csv.reader(handle, dialect)
        return [[cell.strip() for cell in row] for row in reader]


def load_shared_strings(archive: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in archive.namelist():
        return []
    root = ElementTree.fromstring(archive.read("xl/sharedStrings.xml"))
    shared_strings: list[str] = []
    for item in root.findall("main:si", XML_NS):
        text_parts = [node.text or "" for node in item.findall(".//main:t", XML_NS)]
        shared_strings.append("".join(text_parts))
    return shared_strings


def column_letters_to_index(letters: str) -> int:
    value = 0
    for char in letters.upper():
        if not ("A" <= char <= "Z"):
            stop_err(f"Invalid Excel column reference '{letters}'.")
        value = value * 26 + (ord(char) - ord("A") + 1)
    return value - 1


def cell_reference_to_index(reference: str) -> int:
    match = re.match(r"([A-Z]+)", reference.upper())
    if not match:
        stop_err(f"Could not parse XLSX cell reference '{reference}'.")
    return column_letters_to_index(match.group(1))


def get_first_sheet_path(archive: zipfile.ZipFile) -> str:
    workbook = ElementTree.fromstring(archive.read("xl/workbook.xml"))
    first_sheet = workbook.find("main:sheets/main:sheet", XML_NS)
    if first_sheet is None:
        stop_err("The XLSX workbook does not contain any worksheets.")
    relation_id = first_sheet.attrib.get(f"{{{XML_NS['rel']}}}id")
    if not relation_id:
        stop_err("Could not locate the first worksheet relationship in the XLSX workbook.")
    relationships = ElementTree.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    for relationship in relationships.findall("pkgrel:Relationship", XML_NS):
        if relationship.attrib.get("Id") == relation_id:
            target = relationship.attrib["Target"]
            if target.startswith("/"):
                return target.lstrip("/")
            return f"xl/{target.removeprefix('xl/')}"
    stop_err("Could not resolve the first worksheet inside the XLSX workbook.")


def extract_xlsx_cell_value(cell: ElementTree.Element, shared_strings: list[str]) -> str:
    cell_type = cell.attrib.get("t")
    if cell_type == "inlineStr":
        return "".join(node.text or "" for node in cell.findall(".//main:t", XML_NS)).strip()
    value_node = cell.find("main:v", XML_NS)
    if value_node is None or value_node.text is None:
        return ""
    raw_value = value_node.text
    if cell_type == "s":
        shared_index = int(raw_value)
        return shared_strings[shared_index].strip()
    if cell_type == "b":
        return "TRUE" if raw_value == "1" else "FALSE"
    return raw_value.strip()


def load_xlsx_rows(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    with zipfile.ZipFile(path) as archive:
        shared_strings = load_shared_strings(archive)
        sheet_path = get_first_sheet_path(archive)
        root = ElementTree.fromstring(archive.read(sheet_path))
        for row in root.findall(".//main:sheetData/main:row", XML_NS):
            values: list[str] = []
            current_index = 0
            for cell in row.findall("main:c", XML_NS):
                reference = cell.attrib.get("r")
                target_index = current_index if reference is None else cell_reference_to_index(reference)
                while len(values) < target_index:
                    values.append("")
                values.append(extract_xlsx_cell_value(cell, shared_strings))
                current_index = target_index + 1
            rows.append(values)
    return rows


def normalize_input_format(input_format: str) -> str:
    return input_format.lower().strip()


def load_rows(path: Path, input_format: str) -> list[list[str]]:
    if normalize_input_format(input_format) == "xlsx":
        return load_xlsx_rows(path)
    return load_delimited_rows(path)


def normalize_headers(header_row: list[str]) -> list[str]:
    return [header.strip() for header in header_row]


def resolve_column(column_spec: str, headers: list[str], column_role: str) -> int | None:
    spec = (column_spec or "").strip()
    if not spec:
        return None
    if spec.isdigit():
        index = int(spec) - 1
        if index < 0 or index >= len(headers):
            stop_err(f"{column_role} column '{column_spec}' is out of range for {len(headers)} columns.")
        return index
    if re.fullmatch(r"[A-Za-z]+", spec):
        index = column_letters_to_index(spec)
        if index < len(headers):
            return index
    lower_headers = {header.lower(): idx for idx, header in enumerate(headers)}
    if spec.lower() in lower_headers:
        return lower_headers[spec.lower()]
    stop_err(
        f"Could not find the {column_role} column '{column_spec}'. "
        f"Available headers: {', '.join(headers) if headers else '(none)'}"
    )


def wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index : index + width] for index in range(0, len(sequence), width))


def extract_records(rows: list[list[str]], name_column: str, sequence_column: str) -> list[tuple[str, str]]:
    non_empty_rows = list(iter_non_empty_rows(rows))
    if not non_empty_rows:
        stop_err("The input table is empty.")

    headers = normalize_headers(non_empty_rows[0])
    if not headers:
        stop_err("The input table does not contain a readable header row.")

    sequence_index = resolve_column(sequence_column, headers, "sequence")
    assert sequence_index is not None

    try:
        name_index = resolve_column(name_column, headers, "name")
    except SystemExit:
        if len(headers) == 1:
            name_index = None
        else:
            raise

    records: list[tuple[str, str]] = []
    auto_name_counter = 1
    for offset, row in enumerate(non_empty_rows[1:], start=2):
        padded_row = row + [""] * max(0, len(headers) - len(row))
        sequence_value = padded_row[sequence_index].strip() if sequence_index < len(padded_row) else ""
        sequence = re.sub(r"\s+", "", sequence_value)
        if not sequence:
            if any(cell.strip() for cell in padded_row):
                stop_err(f"Row {offset} does not contain a value in the selected sequence column.")
            continue
        if name_index is None or name_index >= len(padded_row):
            name = str(auto_name_counter)
        else:
            name = padded_row[name_index].strip() or str(auto_name_counter)
        records.append((name, sequence))
        auto_name_counter += 1

    if not records:
        stop_err("No FASTA records were produced from the input table.")
    return records


def write_fasta(records: list[tuple[str, str]], output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n{wrap_sequence(sequence)}\n")


def main() -> int:
    args = parse_args()
    rows = load_rows(Path(args.input_path), args.input_format)
    records = extract_records(rows, args.name_column, args.sequence_column)
    write_fasta(records, Path(args.output))
    sys.stdout.write(f"Wrote {len(records)} FASTA records.\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
