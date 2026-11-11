#!/usr/bin/env python3
"""An exploratory, NON-CLINICAL EGFR/LUAD evidence coverage snapshot.

No patient-level records are written. CIViC evidence is CC0; cBioPortal
mutation records are aggregated entirely in memory. See README in run dir.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import re
import sys
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

CIVIC = "https://civicdb.org/downloads/nightly/nightly-ClinicalEvidenceSummaries.tsv"
CBIO_STUDY_ID = "luad_tcga_pan_can_atlas_2018"
CBIO_STUDY = f"https://www.cbioportal.org/api/studies/{CBIO_STUDY_ID}"
CBIO_MUTATIONS = (
    f"https://www.cbioportal.org/api/molecular-profiles/"
    f"{CBIO_STUDY_ID}_mutations/mutations"
    f"?sampleListId={CBIO_STUDY_ID}_all&entrezGeneId=1956"
)
GDC = "https://api.gdc.cancer.gov/projects/TCGA-LUAD?format=json"
# Do not use /analysis/top_cases_counts_by_genes: its aggregate semantics
# differed from independently requested occurrence-level data on 2026-11-10.
GDC_OCCURRENCE_FILTER = {
    "op": "and", "content": [
        {"op": "in", "content": {
            "field": "ssm.consequence.transcript.gene.gene_id",
            "value": ["ENSG00000146648"]}},
        {"op": "in", "content": {
            "field": "case.project.project_id", "value": ["TCGA-LUAD"]}},
    ]
}
MAX_GDC_OCCURRENCES = 1000
GDC_OCCURRENCES = "https://api.gdc.cancer.gov/ssm_occurrences?" + urlencode({
    "filters": json.dumps(GDC_OCCURRENCE_FILTER, separators=(",", ":"), sort_keys=True),
    "fields": "case.case_id,ssm.ssm_id",
    "size": MAX_GDC_OCCURRENCES, "from": 0, "format": "json",
})
SOURCE_URLS = {"civic": CIVIC, "cbio_study": CBIO_STUDY,
               "cbio_mutations": CBIO_MUTATIONS, "gdc": GDC,
               "gdc_occurrences": GDC_OCCURRENCES}
MANDATORY_SOURCES = frozenset({"civic", "cbio_study", "cbio_mutations", "gdc"})
REQUIRED_CIVIC = {"molecular_profile", "evidence_id", "evidence_status",
                  "disease", "doid", "evidence_type", "evidence_civic_url"}
GENE_PATTERN = re.compile(r"\bEGFR\b", flags=re.IGNORECASE)
MAX_SOURCE_BYTES = 8_000_000


def fetch(url: str, retries: int = 2) -> bytes:
    """Bound downloads and fail closed if a source is unavailable."""
    for attempt in range(retries + 1):
        try:
            req = Request(url, headers={"User-Agent": "open-cancer-pilot/0.2",
                                        "Accept": "application/json"})
            with urlopen(req, timeout=40) as response:
                data = response.read(MAX_SOURCE_BYTES + 1)
                if len(data) > MAX_SOURCE_BYTES:
                    raise ValueError("Source exceeds 8 MB safety limit")
                return data
        except (HTTPError, URLError, TimeoutError) as exc:
            if attempt >= retries:
                raise RuntimeError(f"Upstream request failed for {url}: {exc}") from exc
            time.sleep(2 ** attempt)
    raise AssertionError("unreachable")


def sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def as_csv(rows: list[dict], fields: list[str]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        # Reduce the chance of spreadsheet formula injection by source data.
        clean = {}
        for key in fields:
            val = str(row.get(key, ""))
            clean[key] = "'" + val if val.startswith(("=", "+", "-", "@")) else val
        writer.writerow(clean)
    return buffer.getvalue()


def civic_analysis(raw: bytes):
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), delimiter="\t")
    if not REQUIRED_CIVIC.issubset(set(reader.fieldnames or [])):
        raise ValueError(f"CIViC columns changed: {reader.fieldnames}")
    total_rows = 0
    egfr_accepted = 0
    broader_lung_not_luad = 0
    included = []
    status_counts = Counter()
    for row in reader:
        total_rows += 1
        if not GENE_PATTERN.search(row["molecular_profile"] or ""):
            continue
        status = (row["evidence_status"] or "").strip().casefold()
        status_counts[status] += 1
        if status != "accepted":
            continue
        egfr_accepted += 1
        doid = (row.get("doid") or "").upper().replace("DOID:", "").strip()
        disease = (row.get("disease") or "").strip().casefold()
        is_luad = doid == "3910" or (not doid and disease in {
            "lung adenocarcinoma", "adenocarcinoma of the lung",
            "adenocarcinoma of lung"})
        if not is_luad:
            if "lung" in disease:
                broader_lung_not_luad += 1
            continue
        evidence_id = str(row["evidence_id"] or "").strip()
        if not evidence_id.isdigit():
            raise ValueError(f"Unexpected CIViC evidence ID: {evidence_id!r}")
        evidence_url = (row.get("evidence_civic_url") or "").strip()
        if not evidence_url.startswith("https://civicdb.org/"):
            evidence_url = f"https://civicdb.org/evidence/{evidence_id}/summary"
        included.append({
            "source": "civic",
            "record_id": evidence_id,
            "molecular_profile": row["molecular_profile"].strip(),
            "disease": row["disease"].strip(),
            "doid": "DOID:3910",
            "evidence_type": (row.get("evidence_type") or "").strip(),
            "record_url": evidence_url,
        })
    if not total_rows or not egfr_accepted:
        raise ValueError("No usable CIViC EGFR evidence; refusing an apparent zero result")
    ids = [item["record_id"] for item in included]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate CIViC EIDs in selected records")
    included.sort(key=lambda item: int(item["record_id"]))
    counts = {
        "all_export_rows": total_rows,
        "egfr_accepted_rows_all_diseases": egfr_accepted,
        "egfr_accepted_evidence_exact_luad": len(included),
        "egfr_accepted_other_lung_labels_excluded": broader_lung_not_luad,
        "egfr_evidence_status_counts_all_diseases": dict(sorted(status_counts.items())),
        "evidence_type_counts_exact_luad": dict(sorted(Counter(
            r["evidence_type"] for r in included).items())),
        "filter": "EGFR token in molecular_profile; accepted status; DOID:3910 or exact LUAD label when DOID missing",
        "unit": "CIViC evidence item; not a patient, variant prevalence estimate, or independent publication",
    }
    return counts, included


def cbio_analysis(study_bytes: bytes, mutations_bytes: bytes):
    study = json.loads(study_bytes)
    if not isinstance(study, dict) or study.get("studyId") != CBIO_STUDY_ID:
        raise ValueError("Unexpected cBioPortal study response")
    records = json.loads(mutations_bytes)
    if not isinstance(records, list):
        raise ValueError("cBioPortal mutations endpoint did not return a list")
    changes = Counter()
    types = Counter()
    for record in records:
        if not isinstance(record, dict):
            raise ValueError("Non-object cBioPortal mutation record")
        if str(record.get("entrezGeneId")) != "1956":
            raise ValueError("Non-EGFR cBioPortal mutation record")
        if record.get("studyId") != CBIO_STUDY_ID:
            raise ValueError("Cross-study mutation record")
        change = (record.get("proteinChange") or "").strip() or "<missing>"
        changes[change] += 1
        types[(record.get("mutationType") or "<missing>").strip()] += 1
    if not records:
        raise ValueError("Zero mutation rows; cannot silently infer absence")
    # Retain only aggregate mutation-call counts, never sample/patient IDs.
    top = sorted(changes.items(), key=lambda item: (-item[1], item[0]))[:12]
    top_total = sum(n for _, n in top)
    rows = [{"protein_change_as_reported": name, "mutation_call_rows": n}
            for name, n in top]
    if top_total < len(records):
        rows.append({"protein_change_as_reported": "<other>",
                     "mutation_call_rows": len(records) - top_total})
    summary = {
        "study_id": CBIO_STUDY_ID,
        "study_name": study.get("name"),
        "egfr_mutation_call_rows": len(records),
        "protein_changes_reported": len(changes),
        "reported_mutation_type_counts": dict(sorted(types.items())),
        "top_protein_changes_limit": 12,
        "unit": "mutation-call API row; not a distinct patient or a clinical efficacy measurement",
        "selection": "entrezGeneId=1956, mutation profile and all-samples list from one LUAD study",
    }
    return summary, rows


def gdc_analysis(raw: bytes):
    obj = json.loads(raw)
    data = obj.get("data") if isinstance(obj, dict) else None
    if not isinstance(data, dict) or data.get("project_id") != "TCGA-LUAD":
        raise ValueError("Unexpected GDC project response; requested JSON")
    summary = data.get("summary") or {}
    return {
        "project_id": data["project_id"],
        "name": data.get("name"),
        "case_count_reported": summary.get("case_count"),
        "note": "Project metadata only. EGFR counts come from the separately filtered GDC /ssm_occurrences endpoint; case IDs are deduplicated only in memory, without a profiled denominator.",
        "unit": "GDC project, not a gene-level observation",
    }


def gdc_occurrences_analysis(raw: bytes) -> dict:
    """Count distinct public GDC case IDs for EGFR SSM occurrences in memory.

    No case-level identifiers, observed mutations, or raw response are written.
    This counts distinct cases IN THE MUTATION-OBSERVED COHORT only, not the
    full project population or clinical driver mutation prevalence.
    """
    obj = json.loads(raw)
    if not isinstance(obj, dict) or obj.get("warnings"):
        raise ValueError("GDC occurrence endpoint returned warnings or a non-object")
    data = obj.get("data")
    if not isinstance(data, dict):
        raise ValueError("GDC occurrence endpoint missing data")
    pagination = data.get("pagination")
    hits = data.get("hits")
    if not isinstance(pagination, dict) or not isinstance(hits, list):
        raise ValueError("GDC occurrence response lacks pagination or hits")
    total = pagination.get("total")
    if type(total) is not int or total <= 0 or total > MAX_GDC_OCCURRENCES:
        raise ValueError("GDC occurrence count is zero, invalid, or exceeds the safe single-page limit")
    if pagination.get("from") != 0 or len(hits) != total:
        raise ValueError("GDC occurrence response is incomplete or truncated")
    case_ids = set()
    occurrence_ids = set()
    for hit in hits:
        if not isinstance(hit, dict) or not isinstance(hit.get("case"), dict):
            raise ValueError("GDC occurrence record missing case")
        case_id = hit["case"].get("case_id")
        occurrence_id = hit.get("id")
        if not isinstance(case_id, str) or not case_id:
            raise ValueError("GDC occurrence missing case_id")
        if not isinstance(occurrence_id, str) or not occurrence_id:
            raise ValueError("GDC occurrence missing id")
        if occurrence_id in occurrence_ids:
            raise ValueError("Duplicate GDC occurrence ID; response cannot be trusted")
        case_ids.add(case_id)
        occurrence_ids.add(occurrence_id)
    return {
        "status": "retrieved",
        "project_id": "TCGA-LUAD",
        "gene_id": "ENSG00000146648",
        "somatic_mutation_occurrence_rows": len(occurrence_ids),
        "distinct_case_ids_with_recorded_occurrence": len(case_ids),
        "unit": "Distinct GDC public case identifiers associated with an EGFR simple somatic mutation occurrence",
        "interpretation": "A descriptive, project-filtered count among indexed mutation occurrences; neither driver status nor full-cohort prevalence.",
        "privacy": "Case identifiers are deduplicated only in memory; none are written to the repository or logs.",
        "limit": MAX_GDC_OCCURRENCES,
    }


def run(downloads: dict[str, bytes], outdir: Path, run_id: str, accessed_on: str,
        captured_at: str, emit_review_payload: bool = False,
        optional_gdc_error: str | None = None) -> dict:
    received = set(downloads)
    if not MANDATORY_SOURCES.issubset(received) or received - set(SOURCE_URLS):
        raise ValueError("Missing mandatory upstream source or unexpected source key")
    civic, rows = civic_analysis(downloads["civic"])
    cbio, changes = cbio_analysis(downloads["cbio_study"], downloads["cbio_mutations"])
    gdc = gdc_analysis(downloads["gdc"])
    if "gdc_occurrences" in downloads:
        gdc["egfr_occurrences"] = gdc_occurrences_analysis(downloads["gdc_occurrences"])
    else:
        gdc["egfr_occurrences"] = {
            "status": "unavailable",
            "somatic_mutation_occurrence_rows": None,
            "distinct_case_ids_with_recorded_occurrence": None,
            "unit": "GDC public simple somatic mutation occurrence records",
            "error": optional_gdc_error or "No GDC mutation occurrences retrieved",
            "interpretation": "Null means an unavailable upstream response, NOT zero mutation occurrences or patients.",
        }
    source_hashes = {key: sha256(value) for key, value in sorted(downloads.items())}
    outdir.mkdir(parents=True, exist_ok=True)
    coverage = {
        "schema_version": 1,
        "status": "unreviewed",
        "created_at": captured_at,
        "scope": "Descriptive EGFR/LUAD data coverage; not biomedical or therapeutic inference",
        "gene": {"symbol": "EGFR", "entrez_id": "1956", "ensembl_id": "ENSG00000146648"},
        "disease": {"name": "lung adenocarcinoma", "doid": "DOID:3910"},
        "source_content_sha256": source_hashes,
        "source_urls": SOURCE_URLS,
        "source_availability": {key: ("retrieved" if key in downloads else "unavailable")
                                for key in SOURCE_URLS},
        "civic": civic, "cbioportal": cbio, "gdc": gdc,
        "limitations": [
            "CIViC evidence items, cBioPortal mutation-call rows and GDC project case identifiers are different units and study selections; their counts are not directly comparable.",
            "CIViC DOID:3910 uses an exact ontology filter; broader non-small cell lung cancer evidence is excluded even when it might overlap LUAD.",
            "CIViC molecular profiles can be compound; matching the EGFR token does not imply an isolated EGFR alteration.",
            "cBioPortal counts reflect this specific study and assay/profiling coverage, not all lung cancers or prevalence.",
            "GDC case IDs are deduplicated from filtered SSM occurrences in memory. No variant-consequence, driver, or sequenced-denominator filter is applied; this is not a prevalence estimate. A missing response is null, never zero.",
            "Nightly CIViC files and cBioPortal API results are mutable: match source hashes and dates before comparing reruns.",
            "No patient-level records are saved; result has not undergone scientific or clinical review.",
        ],
    }
    artifact_text = {
        "coverage.json": json.dumps(coverage, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        "civic_evidence_ids.csv": as_csv(rows, [
            "source", "record_id", "molecular_profile", "disease", "doid",
            "evidence_type", "record_url"]),
        "cbio_protein_changes.csv": as_csv(changes, [
            "protein_change_as_reported", "mutation_call_rows"]),
    }
    inputs = [
        {"source_id": "civic", "accessed_on": accessed_on,
         "query_or_filter": "Nightly full clinical evidence TSV; accepted; EGFR token; DOID:3910",
         "url": CIVIC, "source_release": "sha256:" + source_hashes["civic"]},
        {"source_id": "cbioportal", "accessed_on": accessed_on,
         "query_or_filter": f"study={CBIO_STUDY_ID}, mutation profile, sampleListId={CBIO_STUDY_ID}_all, entrezGeneId=1956",
         "url": CBIO_MUTATIONS, "source_release": "sha256:" + source_hashes["cbio_mutations"]},
        {"source_id": "gdc", "accessed_on": accessed_on,
         "query_or_filter": "Project TCGA-LUAD metadata; format=json",
         "url": GDC, "source_release": "sha256:" + source_hashes["gdc"]},
    ]
    if "gdc_occurrences" in downloads:
        inputs.append({"source_id": "gdc", "accessed_on": accessed_on,
            "query_or_filter": "SSM occurrences: case.project.project_id=TCGA-LUAD and ssm.consequence.transcript.gene.gene_id=ENSG00000146648; single-page bound 1000; count distinct case IDs in memory",
            "url": GDC_OCCURRENCES,
            "source_release": "sha256:" + source_hashes["gdc_occurrences"]})
    manifest = {
        "schema_version": 1,
        "task_id": "egfr-luad-evidence-coverage",
        "run_id": run_id,
        "submitted_by": {"type": "agent", "name": "open-cancer-research-pilot", "model_or_tool": "Python standard library"},
        "created_at": captured_at,
        "inputs": inputs,
        "execution": {"code_path": "analyze.py",
                      "command": f"python {outdir.as_posix()}/analyze.py --out-dir {outdir.as_posix()} --run-id {run_id}",
                      "environment": "Python 3.12+ standard library with outbound HTTPS; repository root as current directory"},
        "artifacts": [{"path": name, "sha256": sha256(content.encode("utf-8")),
                       "description": {
                           "coverage.json": "Source hashes, aggregate counts and interpretation limitations",
                           "civic_evidence_ids.csv": "CIViC CC0 evidence identifiers and public record links",
                           "cbio_protein_changes.csv": "Aggregated EGFR mutation-call rows; no sample identifiers",
                       }[name]} for name, content in sorted(artifact_text.items())],
        "findings": [],
        "limitations": "Unreviewed source-coverage snapshot; GDC occurrence-level case counts are not full project prevalence, clinical relevance, or comparable to cBioPortal mutation-call counts.",
    }
    artifact_text["manifest.json"] = json.dumps(manifest, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    for name, content in artifact_text.items():
        (outdir / name).write_text(content, encoding="utf-8", newline="\n")
    if emit_review_payload:
        for name, content in sorted(artifact_text.items()):
            # Logs contain only CC0 evidence IDs, source links, and aggregate summaries.
            print("OPEN_CANCER_ARTIFACT::" + json.dumps({
                "path": (outdir / name).as_posix(), "content": content}, ensure_ascii=True))
    print("RESULT " + json.dumps({
        "civic_exact_luad_evidence_items": civic["egfr_accepted_evidence_exact_luad"],
        "cbio_egfr_mutation_calls": cbio["egfr_mutation_call_rows"],
        "gdc_project": gdc["project_id"],
        "gdc_egfr_occurrences": gdc["egfr_occurrences"]["somatic_mutation_occurrence_rows"],
        "gdc_egfr_distinct_cases": gdc["egfr_occurrences"]["distinct_case_ids_with_recorded_occurrence"],
        "source_hashes": source_hashes}))
    return coverage


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--emit-review-payload", action="store_true")
    args = parser.parse_args()
    today = datetime.now(timezone.utc)
    downloads = {key: fetch(SOURCE_URLS[key]) for key in sorted(MANDATORY_SOURCES)}
    optional_gdc_error = None
    try:
        downloads["gdc_occurrences"] = fetch(GDC_OCCURRENCES, retries=3)
        # Do not silently turn GDC warnings, schema changes, or truncation into results.
        gdc_occurrences_analysis(downloads["gdc_occurrences"])
    except (RuntimeError, ValueError, json.JSONDecodeError) as exc:
        downloads.pop("gdc_occurrences", None)
        optional_gdc_error = str(exc)
        print("WARNING: Optional GDC SSM occurrences unavailable; recording null, not zero.", flush=True)
    run(downloads, args.out_dir, args.run_id,
        today.date().isoformat(), today.isoformat().replace("+00:00", "Z"),
        args.emit_review_payload, optional_gdc_error)
    return 0


if __name__ == "__main__":
    sys.exit(main())
