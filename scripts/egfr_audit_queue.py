#!/usr/bin/env python3
"""Offline, nonclinical triage aids for EGFR/LUAD source-level scientific QA.

Only the reviewed pilot's existing CC0 CIViC evidence IDs, profile labels,
aggregate counts and source fingerprints are used. This NEVER recontacts
upstream sources, adjudicates scientific correctness, or processes patient IDs.
"""
from __future__ import annotations

from collections import Counter
import csv
import hashlib
import io
import json
import re

REVIEW_FIELDS = (
    "record_id", "molecular_profile", "disease", "doid", "evidence_type",
    "record_url", "text_only_flags", "upstream_record_checked",
    "scientific_interpretation_checked",
)
PROFILE_GENE = re.compile(r"\bEGFR\b", re.IGNORECASE)
COMPOUND = re.compile(r"\b(?:AND|OR)\b", re.IGNORECASE)
GENERIC = re.compile(r"^EGFR\s+(?:Rare\s+)?Mutation$", re.IGNORECASE)


def text_only_flags(profile: str) -> list[str]:
    """Return *review priorities*, not biological variant classifications."""
    if not isinstance(profile, str) or not PROFILE_GENE.search(profile):
        raise ValueError("Expected an EGFR-containing public CIViC profile label")
    flags = []
    if COMPOUND.search(profile):
        flags.append("multiple-terms-AND-OR")
    if "::" in profile or re.search(r"\bfusion\b", profile, re.IGNORECASE):
        flags.append("fusion-or-partner-profile")
    if GENERIC.fullmatch(profile.strip()):
        flags.append("categorical-unspecified-mutation")
    return flags


def _safe_csv(value) -> str:
    s = str(value)
    # Safe for manual spreadsheet inspection; never evaluate source text.
    return "'" + s if s.startswith(("=", "+", "-", "@")) else s


def make_civic_triage(rows: list[dict]) -> tuple[bytes, dict]:
    """Inventory all selected public CIViC IDs without silently adjudicating them."""
    if not rows:
        raise ValueError("No selected evidence: cannot create an empty review queue")
    ids = set()
    counts = Counter()
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(REVIEW_FIELDS),
                            lineterminator="\n")
    writer.writeheader()
    for row in sorted(rows, key=lambda row: int(row["record_id"])):
        id_ = str(row["record_id"])
        if not id_.isdigit() or id_ in ids:
            raise ValueError("Missing, duplicated or invalid CIViC evidence ID")
        ids.add(id_)
        url = row.get("record_url")
        if url != "https://civicdb.org/links/evidence_items/" + id_:
            raise ValueError("CIViC record link does not bind to the evidence ID")
        if row.get("doid") != "DOID:3910":
            raise ValueError("Original strict LUAD filter was changed")
        flags = text_only_flags(row["molecular_profile"])
        counts[row["evidence_type"]] += 1
        record = {
            "record_id": id_,
            "molecular_profile": row["molecular_profile"],
            "disease": row["disease"],
            "doid": row["doid"],
            "evidence_type": row["evidence_type"],
            "record_url": url,
            "text_only_flags": ";".join(flags) if flags else "none-detected",
            # An offline label inventory is NEVER a live source-record audit.
            "upstream_record_checked": "no",
            "scientific_interpretation_checked": "no",
        }
        writer.writerow({key: _safe_csv(record[key]) for key in REVIEW_FIELDS})
    return stream.getvalue().encode("utf-8"), dict(sorted(counts.items()))


def make_cross_source_units(coverage: dict, manifest_bytes: bytes, civic_count: int,
                            cbio_aggregate_rows: int) -> bytes:
    """Compare published aggregates with their own units; no cohort pooling."""
    civic = coverage["civic"]
    cbio = coverage["cbioportal"]
    gdc = coverage["gdc"]["egfr_occurrences"]
    if civic_count != civic["egfr_accepted_evidence_exact_luad"]:
        raise ValueError("Published CIViC selected count differs from saved IDs")
    if cbio_aggregate_rows != cbio["egfr_mutation_call_rows"]:
        raise ValueError("cBioPortal aggregate counts do not sum to saved mutation-call rows")
    if gdc["status"] == "retrieved":
        occurrences = gdc["somatic_mutation_occurrence_rows"]
        cases = gdc["distinct_case_ids_with_recorded_occurrence"]
        if (type(occurrences) is not int or type(cases) is not int or
                not 0 < cases <= occurrences <= gdc["limit"]):
            raise ValueError("Inconsistent GDC occurrence and case-ID summary")
    else:
        occurrences = gdc["somatic_mutation_occurrence_rows"]
        cases = gdc["distinct_case_ids_with_recorded_occurrence"]
        if occurrences is not None or cases is not None:
            raise ValueError("Unavailable GDC occurrences must remain null, never zero")
    summary = {
        "schema_version": 1,
        "status": "offline-consistency-only-not-independent-scientific-review",
        "original_manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
        "original_snapshot_created_at": coverage["created_at"],
        "source_response_sha256": coverage["source_content_sha256"],
        "sources": {
            "civic": {
                "source_url": coverage["source_urls"]["civic"],
                "count": civic_count,
                "unit": "editor-accepted CIViC evidence-item row selected by exact LUAD/EGFR text filter",
                "excluded_other_lung_labels_reported": civic["egfr_accepted_other_lung_labels_excluded"],
                "excluded_record_ids_available_in_original_artifact": False,
                "upstream_item_verification": "not-performed",
            },
            "cbioportal": {
                "source_url": coverage["source_urls"]["cbio_mutations"],
                "count": cbio_aggregate_rows,
                "unit": "EGFR mutation-call API row from study/profile/sample list",
                "molecular_profiled_denominator": None,
                "unique_patients_reported": None,
                "sample_cohort_harmonization": "not-performed",
            },
            "gdc": {
                "source_url": coverage["source_urls"]["gdc_occurrences"],
                "occurrence_rows": occurrences,
                "distinct_case_ids_aggregated_in_memory": cases,
                "unit": "simple somatic mutation occurrence versus distinct matched case identifier",
                "project_profiled_denominator": None,
                "mutation_consequence_and_driver_filter": "not-applied",
                "sample_cohort_harmonization": "not-performed",
            },
        },
        "cross_source_prevalence_calculable": False,
        "clinical_relevance_assessed": False,
        "independent_domain_review_completed": False,
        "limitations": [
            "Original source response hashes are fingerprints, not an archived immutable upstream release.",
            "Different observation units, assays and selection criteria prohibit pooled frequencies.",
            "No original raw CIViC rows for excluded broader lung labels were retained; false-negative assessment requires a separately documented live source retrieval.",
            "No individual CIViC source records were independently reopened by this offline check.",
            "No source-specific profiled denominator, treatment response or clinical efficacy was established.",
            "No patient or case identifiers are included in this output.",
        ],
    }
    return (json.dumps(summary, indent=2, sort_keys=True, ensure_ascii=False) + "\n").encode("utf-8")
