#!/usr/bin/env python3
"""Select frozen NCI aggregate cells from already acquired official CSVs.

Nonclinical, unreviewed computational research. No medical advice.
Standard library only; this script never accesses the network. Source CSVs
and pdftotext -layout outputs remain outside the repository. For example:

    python extract_sources.py --source-dir /path/to/sources \
        --output-dir data --pdf-text-dir /path/to/sources --check

The optional PDF audit compares every selected printed string with a second
representation of the same NCI source. Agreement is not independent biological
validation, source authentication, or confirmation of the original statistics.
"""

import argparse
import csv
import io
import json
from pathlib import Path
import re


BASE_URL = "https://seer.cancer.gov/archive/csr/1975_2017/"
SECTION = "Age-Specific SEER Incidence Rates, 2013-2017"
PDF_SECTION = "Age-Specific Rates, 2013-2017"
ACCESS_DATE = "2026-11-18"
RACES = ("All Races", "Whites", "Blacks")
SEXES = ("Both Sexes", "Males", "Females")
AGES = tuple(f"{start}-{start + 4}" for start in range(30, 85, 5))
SPECS = (
    {
        "site": "female_breast", "table": "4.11", "section": "4",
        "stem": "sect_04_table.11", "title": "Cancer of the Breast (Invasive)",
        "sexes": ("Females",),
    },
    {
        "site": "colon_and_rectum", "table": "6.10", "section": "6",
        "stem": "sect_06_table.10", "title": "Cancer of the Colon and Rectum (Invasive)",
        "sexes": SEXES,
    },
    {
        "site": "lung_and_bronchus", "table": "15.10", "section": "15",
        "stem": "sect_15_table.10", "title": "Cancer of the Lung and Bronchus (Invasive)",
        "sexes": SEXES,
    },
    {
        "site": "pancreas", "table": "22.7", "section": "22",
        "stem": "sect_22_table.07", "title": "Cancer of the Pancreas (Invasive)",
        "sexes": SEXES,
    },
)
FIELDS = (
    "group_id", "cancer_site", "behavior", "race", "sex", "registry_set",
    "period_start", "period_end", "age_group", "age_start", "age_end_exclusive",
    "rate_per_100000", "rate_units", "reported_decimal_places", "case_count",
    "person_years", "source_name", "source_version", "source_url",
    "source_pdf_url", "source_html_url", "source_table", "source_section",
    "source_column", "source_page", "access_date", "reuse_basis", "reuse_terms_url",
)


def normalized_text(path):
    """Read UTF-8 with universal newline normalization, preserving other text."""
    return Path(path).read_text(encoding="utf-8")


def _column_header(value):
    return " / ".join(part.strip() for part in value.split(","))


def _slug(value):
    return value.lower().replace(" ", "_")


def extract_table(path, spec):
    """Return the prespecified long-format rows, checking the source schema."""
    raw_rows = list(csv.reader(io.StringIO(normalized_text(path))))
    rows = [[cell.strip() for cell in row] for row in raw_rows]
    if rows[:2] != [[f"Table {spec['table']}"], [spec["title"]]]:
        raise ValueError(f"Wrong source table/title: {Path(path).as_posix()}")
    locations = [i for i, row in enumerate(rows) if row == [SECTION]]
    if len(locations) != 1:
        raise ValueError("Expected exactly one age-specific incidence subsection")
    start = locations[0]
    header = rows[start + 1]
    expected = ["Age at Diagnosis"] + [
        f"{race} / {sex}" for race in RACES for sex in SEXES
    ]
    if [_column_header(cell) for cell in header] != expected:
        raise ValueError("Unexpected incidence column order or labels")

    by_age = {}
    for row in rows[start + 2:]:
        if not row:
            break
        if len(row) != len(header):
            raise ValueError("Malformed age-specific incidence row")
        if row[0] in AGES:
            if row[0] in by_age:
                raise ValueError(f"Duplicate source age: {row[0]}")
            by_age[row[0]] = row
    if set(by_age) != set(AGES):
        raise ValueError("Source does not contain every frozen age band")

    selected = []
    for race in RACES:
        for sex in spec["sexes"]:
            column = 1 + RACES.index(race) * len(SEXES) + SEXES.index(sex)
            for age in AGES:
                rate = by_age[age][column]
                if rate != "-" and re.fullmatch(r"\d+\.\d", rate) is None:
                    raise ValueError(f"Unexpected printed rate: {rate!r}")
                age_start = int(age.split("-")[0])
                selected.append({
                    "group_id": "__".join((spec["site"], _slug(race), _slug(sex))),
                    "cancer_site": spec["site"],
                    "behavior": "invasive",
                    "race": race,
                    "sex": sex,
                    "registry_set": "SEER 21 areas",
                    "period_start": "2013",
                    "period_end": "2017",
                    "age_group": age,
                    "age_start": str(age_start),
                    "age_end_exclusive": str(age_start + 5),
                    "rate_per_100000": rate,
                    "rate_units": "cases per 100000 population per year",
                    "reported_decimal_places": "1",
                    "case_count": "",
                    "person_years": "",
                    "source_name": "National Cancer Institute, SEER Cancer Statistics Review 1975-2017",
                    "source_version": "November 2019 SEER submission; released 2020-04-15",
                    "source_url": f"{BASE_URL}download_csr_datafile.php/{spec['stem']}.csv",
                    "source_pdf_url": f"{BASE_URL}results_single/{spec['stem']}.pdf",
                    "source_html_url": (
                        f"{BASE_URL}browse_csr.php?pageSEL={spec['stem']}"
                        f"&sectionSEL={spec['section']}"
                    ),
                    "source_table": spec["table"],
                    "source_section": SECTION,
                    "source_column": f"{race} / {sex}",
                    "source_page": "1",
                    "access_date": ACCESS_DATE,
                    "reuse_basis": (
                        "NCI report is public domain; reproduction or copying permitted "
                        "without permission; source citation requested"
                    ),
                    "reuse_terms_url": BASE_URL,
                })
    return selected


def audit_pdf_text(rows, text, spec):
    """Compare every selected cell with independently extracted PDF text.

    Input is the corresponding one-page PDF's `pdftotext -layout` output.
    Use the incidence section before U.S. Mortality and verify the displayed
    race/sex order. Only whitespace is normalized; printed rate strings must
    agree exactly. Return the number of matched cells or raise ValueError.
    """
    lines = [" ".join(line.split()) for line in text.splitlines()]
    if f"Table {spec['table']}" not in lines or spec["title"] not in lines:
        raise ValueError("PDF table/title mismatch")
    if "All Races Whites Blacks" not in lines:
        raise ValueError("PDF race headings mismatch")
    if "Total Males Females Total Males Females Total Males Females" not in lines:
        raise ValueError("PDF sex headings mismatch")
    if "SEER Incidence" not in lines or "Age at Diagnosis" not in lines:
        raise ValueError("PDF incidence heading missing")
    if not any("SEER 21 areas." in line and "per 100,000" in line for line in lines):
        raise ValueError("PDF registry/units footnote missing")
    incidence_start = lines.index("SEER Incidence")
    mortality_start = lines.index("U.S. Mortality", incidence_start + 1)
    incidence = lines[incidence_start:mortality_start]
    if incidence.count(PDF_SECTION) != 1:
        raise ValueError("PDF age-specific incidence subsection mismatch")
    body = incidence[incidence.index(PDF_SECTION) + 1:]
    by_age = {}
    for line in body:
        cells = line.split()
        if cells and cells[0] in AGES:
            if len(cells) != 10 or cells[0] in by_age:
                raise ValueError("Malformed or duplicate PDF age row")
            by_age[cells[0]] = cells
    if set(by_age) != set(AGES):
        raise ValueError("Missing selected PDF age")

    checked = 0
    for row in rows:
        column = 1 + RACES.index(row["race"]) * len(SEXES) + SEXES.index(row["sex"])
        expected = by_age[row["age_group"]][column]
        if expected != row["rate_per_100000"]:
            raise ValueError(
                f"PDF/CSV mismatch: {row['group_id']} {row['age_group']} "
                f"PDF={expected!r} CSV={row['rate_per_100000']!r}"
            )
        checked += 1
    expected_count = len(RACES) * len(spec["sexes"]) * len(AGES)
    if checked != expected_count:
        raise ValueError("Selected row count does not match the protocol")
    return checked


def csv_text(rows):
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--pdf-text-dir", type=Path)
    parser.add_argument("--check", action="store_true", help="Check stored CSVs without writing")
    args = parser.parse_args()
    generated = []
    counts = {}
    for spec in SPECS:
        selected = extract_table(args.source_dir / f"{spec['stem']}.csv", spec)
        audited = 0
        if args.pdf_text_dir is not None:
            audited = audit_pdf_text(
                selected, normalized_text(args.pdf_text_dir / f"{spec['stem']}.txt"), spec
            )
        output = args.output_dir / f"{spec['site']}.csv"
        generated.append((output, csv_text(selected)))
        counts[output.as_posix()] = {"selected_cells": len(selected), "pdf_matches": audited}
    # Validate every source before writing any output.
    for output, text in generated:
        if args.check:
            if normalized_text(output) != text:
                raise ValueError(f"Stored selection differs: {output.as_posix()}")
        else:
            output.parent.mkdir(parents=True, exist_ok=True)
            with output.open("w", encoding="utf-8", newline="") as stream:
                stream.write(text)
    print(json.dumps(counts, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
