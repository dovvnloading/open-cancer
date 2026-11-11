#!/usr/bin/env python3
"""Compare frozen EGFR/LUAD evidence IDs with new PUBLIC metadata only.

This is a nonclinical data-drift audit queue, not a scientific review. No
mutation calls, case identifiers, sample IDs, or patient data are downloaded.
"""
from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import tempfile

try:
    from .egfr_luad_pilot import fetch, SOURCE_URLS, REQUIRED_CIVIC, GENE_PATTERN
    from .validate import ROOT, validate_repository
except ImportError:
    from egfr_luad_pilot import fetch, SOURCE_URLS, REQUIRED_CIVIC, GENE_PATTERN
    from validate import ROOT, validate_repository

RESULT = ROOT / "results/egfr-luad-evidence-coverage/pilot-20261110"
ORIGINAL_MANIFEST_SHA = "9fce9beaaf873b5f0b0a14e338c3b372d8aa5d19c03be11b911ea1eb2b92931e"
SELECTED_FIELDS = (
    "record_id", "original_molecular_profile", "current_molecular_profile",
    "original_disease", "current_disease", "original_doid", "current_doid",
    "original_evidence_type", "current_evidence_type", "current_status",
    "current_strict_luad", "changed_metadata_fields", "record_url",
    "source_record_independently_assessed",
)
BROADER_FIELDS = ("record_id", "molecular_profile", "disease", "doid",
                  "evidence_type", "record_url", "reason",
                  "source_record_independently_assessed")


def _sha(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _csv(rows: list[dict], fields: tuple[str, ...]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for record in rows:
        safe = {}
        for field in fields:
            value = str(record.get(field, ""))
            safe[field] = "'" + value if value.startswith(("=", "+", "-", "@")) else value
        writer.writerow(safe)
    return stream.getvalue().encode("utf-8")


def canonical_doid(raw: object) -> str:
    """Normalize optional DOID namespace prefixes without erasing identifiers."""
    value = str(raw or "").strip()
    if value.upper().startswith("DOID:"):
        value = value[5:].strip()
    return value


def is_strict_luad(row: dict) -> bool:
    doid = canonical_doid(row.get("doid"))
    disease = str(row.get("disease") or "").strip().casefold()
    return doid == "3910" or (not doid and disease in {
        "lung adenocarcinoma", "adenocarcinoma of the lung", "adenocarcinoma of lung"
    })


def parse_civic_snapshot(raw: bytes) -> dict:
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), delimiter="\t")
    if not REQUIRED_CIVIC.issubset(set(reader.fieldnames or ())):
        raise ValueError("CIViC nightly export columns changed or are unavailable")
    all_by_id = {}
    strict = {}
    broader = {}
    total = 0
    for row in reader:
        total += 1
        evidence_id = str(row.get("evidence_id") or "").strip()
        if not evidence_id.isdigit() or evidence_id in all_by_id:
            raise ValueError("CIViC duplicate or missing evidence ID; refusing incomplete audit")
        all_by_id[evidence_id] = row
        if not GENE_PATTERN.search(row.get("molecular_profile") or ""):
            continue
        if str(row.get("evidence_status") or "").strip().casefold() != "accepted":
            continue
        if is_strict_luad(row):
            strict[evidence_id] = row
        elif "lung" in str(row.get("disease") or "").casefold():
            broader[evidence_id] = row
    if not total or not strict:
        raise ValueError("No usable CIViC source or no exact-LUAD evidence; never report false zero")
    return {"all_by_id": all_by_id, "strict": strict, "broader": broader,
            "total": total}


def compare_civic(original_rows: list[dict], current: dict):
    seen = set()
    selected = []
    for original in sorted(original_rows, key=lambda row: int(row["record_id"])):
        evidence_id = original["record_id"]
        if not evidence_id.isdigit() or evidence_id in seen:
            raise ValueError("Original CIViC ID artifact has invalid/duplicate IDs")
        seen.add(evidence_id)
        row = current["all_by_id"].get(evidence_id)
        fields = ("molecular_profile", "disease", "doid", "evidence_type")
        changed = [
            field for field in fields
            if row is None or (
                canonical_doid(row.get(field)) != canonical_doid(original.get(field))
                if field == "doid"
                else str(row.get(field) or "").strip() != str(original[field]).strip()
            )
        ]
        status = str(row.get("evidence_status") or "").strip() if row else "not-present"
        if not row or status.casefold() != "accepted" or not is_strict_luad(row):
            changed.append("current_inclusion")
        record = {
            "record_id": evidence_id,
            "current_status": status,
            "current_strict_luad": "yes" if evidence_id in current["strict"] else "no",
            "changed_metadata_fields": ";".join(changed) if changed else "none",
            "record_url": "https://civicdb.org/links/evidence_items/" + evidence_id,
            "source_record_independently_assessed": "no",
        }
        for field in fields:
            record["original_" + field] = original[field]
            record["current_" + field] = row.get(field, "") if row else ""
        selected.append(record)
    if not selected:
        raise ValueError("No original selected CIViC records")
    broader = []
    for evidence_id, row in sorted(current["broader"].items(), key=lambda item: int(item[0])):
        broader.append({
            "record_id": evidence_id,
            "molecular_profile": row.get("molecular_profile", ""),
            "disease": row.get("disease", ""),
            "doid": row.get("doid", ""),
            "evidence_type": row.get("evidence_type", ""),
            "record_url": "https://civicdb.org/links/evidence_items/" + evidence_id,
            "reason": "live-accepted-EGFR-broader-lung-label-excluded-by-strict-LUAD-filter",
            "source_record_independently_assessed": "no",
        })
    return selected, broader


def make_live_report(*, civic_raw: bytes, cbio_study_raw: bytes, gdc_project_raw: bytes,
                     original_csv: bytes, manifest_bytes: bytes, coverage: dict,
                     checked_at: str) -> dict[str, bytes]:
    if _sha(manifest_bytes) != ORIGINAL_MANIFEST_SHA:
        raise ValueError("Original version-bound scientific manifest has changed")
    snapshot = parse_civic_snapshot(civic_raw)
    original = list(csv.DictReader(io.StringIO(original_csv.decode("utf-8-sig"))))
    selected, broader = compare_civic(original, snapshot)
    if len(original) != coverage["civic"]["egfr_accepted_evidence_exact_luad"]:
        raise ValueError("Original selected CIViC IDs no longer match frozen count")
    cbio = json.loads(cbio_study_raw)
    if not isinstance(cbio, dict) or cbio.get("studyId") != "luad_tcga_pan_can_atlas_2018":
        raise ValueError("Unexpected cBioPortal study metadata")
    gdc_payload = json.loads(gdc_project_raw)
    gdc = gdc_payload.get("data") if isinstance(gdc_payload, dict) else None
    if not isinstance(gdc, dict) or gdc.get("project_id") != "TCGA-LUAD":
        raise ValueError("Unexpected GDC project metadata")
    checked = datetime.fromisoformat(checked_at.replace("Z", "+00:00"))
    if checked.tzinfo is None:
        raise ValueError("Source audit requires timezone-aware UTC timestamp")
    summary = {
        "schema_version": 1,
        "status": "public-metadata-drift-not-independent-scientific-review",
        "checked_at": checked_at,
        "original_manifest_sha256": ORIGINAL_MANIFEST_SHA,
        "original_civic_response_sha256": coverage["source_content_sha256"]["civic"],
        "live_source_sha256": {
            "civic": _sha(civic_raw),
            "cbio_study_metadata": _sha(cbio_study_raw),
            "gdc_project_metadata": _sha(gdc_project_raw),
        },
        "source_urls": {
            "civic": SOURCE_URLS["civic"],
            "cbio_study_metadata": SOURCE_URLS["cbio_study"],
            "gdc_project_metadata": SOURCE_URLS["gdc"],
        },
        "civic": {
            "original_selected": len(original),
            "live_strict_selected": len(snapshot["strict"]),
            "original_broader_lung_excluded_count_only": coverage["civic"]["egfr_accepted_other_lung_labels_excluded"],
            "live_broader_lung_candidates": len(broader),
            "selected_original_IDs_no_longer_strict": sum(row["current_strict_luad"] == "no" for row in selected),
            "selected_original_IDs_with_metadata_or_inclusion_changes": sum(
                row["changed_metadata_fields"] != "none" for row in selected),
            "civic_snapshot_bytes_changed": _sha(civic_raw) != coverage["source_content_sha256"]["civic"],
            "source_record_review": "not-performed",
        },
        "cbioportal": {
            "study_id": cbio["studyId"],
            "study_name": cbio.get("name"),
            "study_level_sequenced_sample_count_if_provided": cbio.get("sequencedSampleCount")
                if type(cbio.get("sequencedSampleCount")) is int else None,
            "egfr_mutation_profile_denominator": None,
            "source_record_review": "not-performed",
        },
        "gdc": {
            "project_id": gdc["project_id"],
            "project_level_case_count_if_provided": (
                (gdc.get("summary") or {}).get("case_count")
                if type((gdc.get("summary") or {}).get("case_count")) is int else None),
            "egfr_profiled_cohort_denominator": None,
            "source_record_review": "not-performed",
        },
        "cross_source_prevalence_calculable": False,
        "scientific_interpretation_assessed": False,
        "clinical_use_assessed": False,
        "limitations": [
            "This compares public metadata only; no individual evidence statement, therapy outcome or primary study has been evaluated.",
            "The original excluded broader-lung IDs were not retained; live broader-lung candidates cannot be equated to original excluded records.",
            "CIViC sources are mutable and can change labels, statuses, and evidence details without preserving old versions.",
            "cBioPortal study-level sequenced samples are not necessarily an EGFR-mutation-profile assay denominator.",
            "GDC project cases are not a denominator for filtered gene-specific SSM occurrences.",
            "No patient/sample/case-level records or mutation API rows were fetched for this report.",
            "No biomedical efficacy, prevalence or clinical conclusion can be calculated or endorsed.",
        ],
    }
    outputs = {
        "civic-selected-record-drift.csv": _csv(selected, SELECTED_FIELDS),
        "civic-broader-lung-candidates.csv": _csv(broader, BROADER_FIELDS),
        "public-source-metadata-summary.json":
            (json.dumps(summary, sort_keys=True, indent=2, ensure_ascii=False) + "\n").encode("utf-8"),
    }
    outputs["SHA256SUMS"] = "".join(
        f"{_sha(payload)}  {name}\n" for name, payload in sorted(outputs.items())
    ).encode("ascii")
    return outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if args.output_dir.exists():
        parser.error("Choose a new output directory; old snapshots must not be overwritten")
    errors = validate_repository(ROOT)
    if errors:
        parser.exit(1, "Repository integrity validation failed:\n" + "\n".join(errors))
    manifest_bytes = (RESULT / "manifest.json").read_bytes()
    original_csv = (RESULT / "civic_evidence_ids.csv").read_bytes()
    coverage = json.loads((RESULT / "coverage.json").read_bytes())
    # Only accepted public CC0 CIViC export and project/study-level metadata.
    # Never fetch mutation-call rows, GDC case data, restricted files or tokens.
    civic_raw = fetch(SOURCE_URLS["civic"])
    cbio_raw = fetch(SOURCE_URLS["cbio_study"])
    gdc_raw = fetch(SOURCE_URLS["gdc"])
    now = datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")
    results = make_live_report(
        civic_raw=civic_raw, cbio_study_raw=cbio_raw, gdc_project_raw=gdc_raw,
        original_csv=original_csv, manifest_bytes=manifest_bytes,
        coverage=coverage, checked_at=now,
    )
    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=args.output_dir.parent) as tmp:
        staging = Path(tmp) / "report"
        staging.mkdir()
        for filename, payload in results.items():
            (staging / filename).write_bytes(payload)
        staging.rename(args.output_dir)
    summary = json.loads(results["public-source-metadata-summary.json"])
    print(json.dumps({
        "checked_at": summary["checked_at"],
        "original_selected": summary["civic"]["original_selected"],
        "live_strict_selected": summary["civic"]["live_strict_selected"],
        "live_broader_lung_candidates": summary["civic"]["live_broader_lung_candidates"],
        "selected_original_IDs_with_metadata_or_inclusion_changes":
            summary["civic"]["selected_original_IDs_with_metadata_or_inclusion_changes"],
        "public_metadata_only": True,
        "human_scientific_review": False,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
