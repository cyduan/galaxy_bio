#!/usr/bin/env python

from __future__ import annotations

import argparse
import csv
import sys
import zipfile
from pathlib import Path
from typing import NoReturn
from xml.sax.saxutils import escape


FIXED_ZIP_TIME = (2026, 1, 1, 0, 0, 0)


def stop_err(message: str) -> NoReturn:
    raise SystemExit(f"ERROR: {message}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert FASTA files into CSV or XLSX tables.")
    parser.add_argument("--input", required=True, dest="input_path")
    parser.add_argument("--output-format", required=True, choices=["csv", "xlsx"])
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    identifier: str | None = None
    sequence_lines: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    sequence = "".join(sequence_lines).strip()
                    if not sequence:
                        stop_err(f"FASTA record '{identifier}' is empty.")
                    records.append((identifier, sequence))
                identifier = line[1:].strip() or str(len(records) + 1)
                sequence_lines = []
            else:
                if identifier is None:
                    stop_err(f"Line {line_number} is not part of a valid FASTA record.")
                sequence_lines.append(line)
    if identifier is not None:
        sequence = "".join(sequence_lines).strip()
        if not sequence:
            stop_err(f"FASTA record '{identifier}' is empty.")
        records.append((identifier, sequence))
    if not records:
        stop_err("No FASTA records were found in the input dataset.")
    return records


def write_csv(records: list[tuple[str, str]], output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["name", "sequence"])
        writer.writerows(records)


def excel_column_name(index: int) -> str:
    result = []
    current = index + 1
    while current > 0:
        current, remainder = divmod(current - 1, 26)
        result.append(chr(ord("A") + remainder))
    return "".join(reversed(result))


def cell_xml(column_index: int, row_index: int, value: str) -> str:
    cell_ref = f"{excel_column_name(column_index)}{row_index}"
    safe_value = escape(value)
    return f'<c r="{cell_ref}" t="inlineStr"><is><t>{safe_value}</t></is></c>'


def build_sheet_xml(records: list[tuple[str, str]]) -> str:
    rows = []
    all_rows = [("name", "sequence"), *records]
    for row_index, row_values in enumerate(all_rows, start=1):
        cells = "".join(cell_xml(column_index, row_index, value) for column_index, value in enumerate(row_values))
        rows.append(f'<row r="{row_index}">{cells}</row>')
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f"<sheetData>{''.join(rows)}</sheetData>"
        "</worksheet>"
    )


def write_zip_text(archive: zipfile.ZipFile, member_name: str, content: str) -> None:
    info = zipfile.ZipInfo(member_name)
    info.date_time = FIXED_ZIP_TIME
    info.compress_type = zipfile.ZIP_DEFLATED
    archive.writestr(info, content)


def write_xlsx(records: list[tuple[str, str]], output_path: Path) -> None:
    content_types = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
  <Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
  <Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>
  <Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
  <Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
</Types>
"""
    rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
  <Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>
"""
    workbook = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <sheets>
    <sheet name="Sequences" sheetId="1" r:id="rId1"/>
  </sheets>
</workbook>
"""
    workbook_rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>
</Relationships>
"""
    styles = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <fonts count="1"><font><sz val="11"/><name val="Calibri"/></font></fonts>
  <fills count="1"><fill><patternFill patternType="none"/></fill></fills>
  <borders count="1"><border/></borders>
  <cellStyleXfs count="1"><xf/></cellStyleXfs>
  <cellXfs count="1"><xf xfId="0"/></cellXfs>
  <cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>
</styleSheet>
"""
    core = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:dcmitype="http://purl.org/dc/dcmitype/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <dc:creator>Codex</dc:creator>
  <cp:lastModifiedBy>Codex</cp:lastModifiedBy>
  <dcterms:created xsi:type="dcterms:W3CDTF">2026-01-01T00:00:00Z</dcterms:created>
  <dcterms:modified xsi:type="dcterms:W3CDTF">2026-01-01T00:00:00Z</dcterms:modified>
</cp:coreProperties>
"""
    app = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">
  <Application>Codex</Application>
</Properties>
"""
    with zipfile.ZipFile(output_path, "w") as archive:
        write_zip_text(archive, "[Content_Types].xml", content_types)
        write_zip_text(archive, "_rels/.rels", rels)
        write_zip_text(archive, "xl/workbook.xml", workbook)
        write_zip_text(archive, "xl/_rels/workbook.xml.rels", workbook_rels)
        write_zip_text(archive, "xl/worksheets/sheet1.xml", build_sheet_xml(records))
        write_zip_text(archive, "xl/styles.xml", styles)
        write_zip_text(archive, "docProps/core.xml", core)
        write_zip_text(archive, "docProps/app.xml", app)


def main() -> int:
    args = parse_args()
    records = parse_fasta(Path(args.input_path))
    output_path = Path(args.output)
    if args.output_format == "csv":
        write_csv(records, output_path)
    else:
        write_xlsx(records, output_path)
    sys.stdout.write(f"Wrote {len(records)} sequence records to {args.output_format.upper()}.\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
