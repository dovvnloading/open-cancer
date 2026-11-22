"""Offline, nonclinical GALAXY denominator audit. No occupancy model is run."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REQUIRED = frozenset({
    "surveillance_patients", "surveillance_samples", "min_tests_per_person",
    "max_tests_per_person", "ever_positive_patients", "serial_negative_patients",
    "ever_positive_recurrences", "serial_negative_recurrences",
    "reported_patient_specificity_numerator", "reported_patient_specificity_denominator",
    "reported_sample_specificity_numerator", "reported_sample_specificity_denominator",
})
GATE_NAMES = frozenset({"G1", "G2", "G3", "G4"})


def normalized_sha256(data: bytes) -> str:
    """The same text checksum on LF, CRLF and bare-CR checkouts."""
    return hashlib.sha256(data.replace(b"\r\n", b"\n").replace(b"\r", b"\n")).hexdigest()


def hashes(root: Path = HERE) -> dict[str, str]:
    paths = (root / "PROTOCOL.md", root / "data" / "printed_counts.csv",
             root / "data" / "denominator_gate.json")
    return {path.relative_to(root).as_posix(): normalized_sha256(path.read_bytes())
            for path in paths}


def counts(path: Path) -> dict[str, int]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ["key", "value", "unit", "locator", "source_url",
                                 "doi", "accessed", "licence"]:
            raise ValueError("unexpected printed-counts schema")
        items: dict[str, int] = {}
        for row in reader:
            key = row["key"]
            if key not in REQUIRED or key in items:
                raise ValueError("missing/duplicate/unexpected count key")
            if (not row["value"].isascii() or not row["value"].isdecimal()
                    or not row["unit"] or not row["locator"]
                    or row["source_url"] != "https://www.nature.com/articles/s41698-026-01277-5"
                    or row["doi"] != "10.1038/s41698-026-01277-5"
                    or row["accessed"] != "2026-11-18"
                    or row["licence"] != "CC BY 4.0"):
                raise ValueError("unattributed or malformed printed count")
            items[key] = int(row["value"])
    if items.keys() != REQUIRED:
        raise ValueError("incomplete printed counts")
    return items


def audit(root: Path = HERE) -> dict[str, object]:
    """Check only printed aggregate consistency and the blocking source gate."""
    c = counts(root / "data" / "printed_counts.csv")
    d = json.loads((root / "data" / "denominator_gate.json").read_text(encoding="utf-8"))
    if (d.get("article_html") != "https://www.nature.com/articles/s41698-026-01277-5"
            or d.get("doi") != "10.1038/s41698-026-01277-5"
            or d.get("licence") != "CC BY 4.0 (article HTML rights statement)"
            or set(d.get("gates", {})) != GATE_NAMES):
        raise ValueError("missing denominator evidence or provenance")
    for gate in d["gates"].values():
        if gate.get("status") not in {"supported", "not_documented", "contradicted"}:
            raise ValueError("unknown gate status")
        if not gate.get("locator") or not gate.get("reason"):
            raise ValueError("undocumented source-gate rationale")
    if not (c["ever_positive_patients"] + c["serial_negative_patients"]
            == c["surveillance_patients"] == 190
            and c["ever_positive_recurrences"] == 54
            and c["serial_negative_recurrences"] == 10
            and c["min_tests_per_person"] == 1
            and c["max_tests_per_person"] == 8
            and c["surveillance_samples"] == 816):
        raise ValueError("reported surveillance aggregates no longer reconcile")
    positive_no_recurrence = (c["ever_positive_patients"]
                              - c["ever_positive_recurrences"])
    negative_no_recurrence = (c["serial_negative_patients"]
                              - c["serial_negative_recurrences"])
    total_no_recurrence = positive_no_recurrence + negative_no_recurrence
    if not (c["reported_patient_specificity_numerator"] == negative_no_recurrence
            and c["reported_patient_specificity_denominator"] == total_no_recurrence
            and 0 <= c["reported_sample_specificity_numerator"]
                <= c["reported_sample_specificity_denominator"]
                <= c["surveillance_samples"]):
        raise ValueError("printed specificity fractions conflict with count identities")
    blocked = sorted(g for g, item in d["gates"].items()
                     if item["status"] != "supported")
    # This audit deliberately cannot be promoted to an occupancy model by
    # changing flags. A new protocol, source linkage and analysis are required.
    if not blocked:
        raise ValueError("all gates marked supported, but occupancy analysis is not implemented")
    return {
        "verdict": "UNTESTABLE", "reason": "same-person sample denominator not documented",
        "blocked_gates": blocked, "occupancy_executed": False,
        "occupancy_p_value": None, "within_person_concentration": None,
        "patient_table": {
            "ever_positive_recurred": c["ever_positive_recurrences"],
            "ever_positive_nonrecurred": positive_no_recurrence,
            "serial_negative_recurred": c["serial_negative_recurrences"],
            "serial_negative_nonrecurred": negative_no_recurrence,
            "nonrecurrent_total": total_no_recurrence,
        },
        "patient_specificity_exact": str(Fraction(
            c["reported_patient_specificity_numerator"],
            c["reported_patient_specificity_denominator"])),
        "sample_specificity_exact": str(Fraction(
            c["reported_sample_specificity_numerator"],
            c["reported_sample_specificity_denominator"])),
        "sample_remainder_unassigned": (c["reported_sample_specificity_denominator"]
                                        - c["reported_sample_specificity_numerator"]),
        "source_hashes_sha256_lf": hashes(root),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true", help="run the stop-gate audit")
    group.add_argument("--json", action="store_true", help="print stable JSON audit")
    args = parser.parse_args()
    result = audit()
    if args.json:
        print(json.dumps(result, sort_keys=True, indent=2))
    else:
        print("UNTESTABLE: G2-G4 have no published exact sample-to-person bridge; "
              "no occupancy test was executed.")


if __name__ == "__main__":
    main()

