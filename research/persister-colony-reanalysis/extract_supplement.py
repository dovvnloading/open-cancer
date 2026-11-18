#!/usr/bin/env python3
"""Re-extract the selected published cells from their pinned source workbook.

Nonclinical provenance audit, not a model fit. No treatment advice. Standard
library only. The original CC BY 4.0 workbook is not redistributed here.
Figure annotations were visually transcribed; this script cannot verify them
from pixels or substitute for the inaccessible Mendeley raw-data archive.
"""

import argparse
import csv
import hashlib
import math
from pathlib import Path, PurePosixPath
import re
import xml.etree.ElementTree as ET
import zipfile


HERE = Path(__file__).resolve().parent
WORKBOOK_SHA256 = "12a201b7ef42678daca80c57e972558a383dba9fedf48366a1aa7177cc775b4f"
MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
NS = {"m": MAIN_NS}


def number_from_text(value):
    """Interpret only explicit decimal/scientific notation, never guess exponents."""
    compact = value.replace(" ", "").replace("\u2212", "-")
    match = re.fullmatch(r"10\^(-?\d+)", compact)
    if match:
        result = 10.0 ** int(match.group(1))
    else:
        match = re.fullmatch(r"([0-9.]+)[×∗]10\^(-?\d+)", compact)
        if match:
            result = float(match.group(1)) * 10.0 ** int(match.group(2))
        else:
            result = float(compact)
    if not math.isfinite(result):
        raise ValueError("Source value must be finite")
    return result


def string_value(node):
    """Keep superscripts/subscripts explicit, so 10^5 never becomes 105."""
    plain, rich = [], []
    for item in node:
        if item.tag == f"{{{MAIN_NS}}}t":
            plain.append(item.text or "")
            rich.append(item.text or "")
        elif item.tag == f"{{{MAIN_NS}}}r":
            text_node = item.find("m:t", NS)
            text = "" if text_node is None else text_node.text or ""
            plain.append(text)
            vertical = item.find("m:rPr/m:vertAlign", NS)
            if vertical is not None:
                if vertical.get("val") == "superscript":
                    text = "^" + text
                elif vertical.get("val") == "subscript":
                    text = "_" + text
            rich.append(text)
    return "".join(plain), "".join(rich)


def workbook_cells(path):
    """Read the exact pinned XLSX bytes, without evaluating workbook formulas."""
    content = Path(path).read_bytes()
    if hashlib.sha256(content).hexdigest() != WORKBOOK_SHA256:
        raise ValueError("Source workbook SHA-256 does not match the inspected version")
    with zipfile.ZipFile(path) as archive:
        strings = [string_value(node) for node in ET.fromstring(
            archive.read("xl/sharedStrings.xml"))]
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {node.get("Id"): node.get("Target") for node in relationships
                   if node.get("TargetMode") != "External"}
        result = {}
        for sheet in workbook.find("m:sheets", NS):
            target = targets[sheet.get(f"{{{REL_NS}}}id")]
            path_parts = PurePosixPath(target)
            if path_parts.is_absolute() or ".." in path_parts.parts:
                raise ValueError("Unexpected workbook relationship path")
            cells = {}
            for cell in ET.fromstring(archive.read("xl/" + target)).findall(".//m:c", NS):
                value = cell.find("m:v", NS)
                if value is None:
                    continue
                text = value.text or ""
                cells[cell.get("r")] = strings[int(text)] if cell.get("t") == "s" else (text, text)
            result[sheet.get("name")] = cells
        return result


def verify_selected(workbook, aggregates=HERE / "data" / "published_parameters.csv"):
    """Every committed selected number must match its exact sheet/cell locator."""
    cells = workbook_cells(workbook)
    with Path(aggregates).open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    seen = set()
    for row in rows:
        key = (row["sheet"], row["cell"])
        if key in seen:
            raise ValueError(f"Repeated source cell: {key}")
        seen.add(key)
        plain, rich = cells[key[0]][key[1]]
        if (plain, rich) != (row["raw_cell_value"], row["rich_text_value"]):
            raise ValueError(f"Source transcription mismatch at {key}")
        if number_from_text(rich) != float(row["value"]):
            raise ValueError(f"Source numeric conversion mismatch at {key}")
        if row["locator"] != "!".join(key):
            raise ValueError(f"Inconsistent source locator: {key}")
    return len(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", type=Path, required=True)
    args = parser.parse_args()
    count = verify_selected(args.workbook)
    print(f"Re-extracted and matched {count} published numeric cells; raw archive remains unaudited.")


if __name__ == "__main__":
    main()
