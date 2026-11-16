"""Bounded EGFR/LUAD job packets from the validated Agent 1 CIViC snapshot."""
from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
from hashlib import sha256
import json

if __package__:
    from . import civic_evidence
    from .compute_fanout import ROOT, SHA, digest, require
    from .egfr_live_source_audit import canonical_doid, parse_civic_snapshot
    from .egfr_luad_pilot import CIVIC, MAX_SOURCE_BYTES, fetch
else:
    import civic_evidence
    from compute_fanout import ROOT, SHA, digest, require
    from egfr_live_source_audit import canonical_doid, parse_civic_snapshot
    from egfr_luad_pilot import CIVIC, MAX_SOURCE_BYTES, fetch

FIXTURE_PATH = "tests/fixtures/civic/egfr_luad_page.json"
TASK = civic_evidence.TASK
ALTERATIONS = ("EGFR L858R", "EGFR exon 19 deletion", "EGFR exon 20 insertion")
LIMIT = 8
FIELDS = ("molecular_profile", "disease", "doid", "evidence_status", "evidence_type",
          "evidence_direction", "evidence_level", "significance", "therapies",
          "therapy_interaction_type", "source_type", "citation_id", "is_flagged")
SAMPLING = (" Up to 8 items per alteration, taking one item in turn from each reported "
            "direction/significance group in lexical order, lowest numeric EID first "
            "within each group. Compound profiles, aliases and broader lung terms are "
            "excluded. This selection is not a representative sample.")
GRAPHQL_SELECTION = ("Agent 1 CIViC GraphQL pages. Accepted, unflagged evidence with a "
                     "nondeprecated exact profile and nondeprecated Lung Adenocarcinoma "
                     "with DOID 3910. Counts refer only to retrieved pages." + SAMPLING)
NIGHTLY_SELECTION = ("CIViC nightly TSV following unavailable GraphQL transport. Accepted, "
                     "unflagged evidence, exact Lung Adenocarcinoma and DOID 3910, exact "
                     "alteration label ignoring case and repeated whitespace. The TSV has no "
                     "profile/disease deprecation flags; that exclusion is unassessed." + SAMPLING)
GAP_REASONS = {
    "gdc": "This pilot has no GDC occurrence connector. Assay denominator unknown.",
    "cbioportal": "This pilot has no cBioPortal mutation-call connector. Profiled denominator unknown.",
}
FALLBACK_ERRORS = {"network-error", "http-error"}


class EvidenceUnavailable(ValueError):
    """A safe retrieval failure with inspectable Agent 1 request lineage."""

    def __init__(self, bundle, fallback_attempted=False):
        self.sources = {"civic": {"status": "unavailable",
                                 "agent1": {k: v for k, v in bundle.items() if k != "records"},
                                 "nightly_fallback_attempted": fallback_attempted}, **_source_gaps()}
        message = "Agent 1 snapshot unavailable: " + bundle["failure"]
        if fallback_attempted:
            message += ". Nightly fallback retrieval or projection failed."
        super().__init__(message)


def _label(value: str) -> str:
    return " ".join(value.casefold().split())


def _source_catalog() -> dict:
    return json.loads((ROOT / "catalog/sources/civic.json").read_text(encoding="utf-8"))


def _source_gaps() -> dict:
    return {name: {"status": "not-retrieved", "reason": reason} for name, reason in GAP_REASONS.items()}


def _source(bundle, mode: str) -> dict:
    catalog = _source_catalog()
    return {"status": "fixture" if mode == "fixture" else "retrieved",
            "unit": "CIViC evidence item", "license": "MIT" if mode == "fixture" else bundle["license_id"],
            "terms_url": "LICENSE" if mode == "fixture" else bundle["terms_url"],
            "attribution": "Open Cancer synthetic fixture" if mode == "fixture" else catalog["publisher"],
            "source_catalog_sha256": digest(catalog),
            "agent1": {k: v for k, v in bundle.items() if k != "records"}}


def _snapshot(mode, bundle, source, rows, selection):
    counts = {}
    selected = []
    for alt in ALTERATIONS:
        matches = [r for r in rows if _label(r["molecular_profile"]) == _label(alt)]
        counts[alt] = len(matches)
        groups = defaultdict(list)
        for row in sorted(matches, key=lambda r: int(r["record_id"].removeprefix("EID"))):
            groups[(row["evidence_direction"], row["significance"])].append(row)
        indexed = [(i, group, r) for group, values in groups.items() for i, r in enumerate(values)]
        selected.extend(r for _, _, r in sorted(indexed, key=lambda t: (t[0], t[1]))[:LIMIT])
    source["matching_by_alteration"] = counts
    body = {"schema_version": 2, "mode": mode, "synthetic": mode == "fixture", "selection": selection,
            "sources": {"civic": source, **_source_gaps()}, "records": selected, "agent1_bundle": bundle}
    return {**body, "snapshot_sha256": digest(body)}


def from_agent1_bundle(bundle: dict, mode: str) -> dict:
    """Import native evidence without network access. Preserve pages, IDs, hashes and truncation."""
    require(mode in {"fixture", "real"}, "Unknown evidence mode")
    civic_evidence.validate_bundle(bundle, _source_catalog())
    require(bundle["synthetic"] is (mode == "fixture"), "Agent 1 bundle mode mismatch")
    if bundle["status"] == "unavailable":
        raise EvidenceUnavailable(bundle)
    metadata = _source(bundle, mode)
    metadata.update(connector="agent1-graphql", retrieval_status=bundle["status"],
                    sha256=bundle["sha256"], url=bundle["endpoint"], retrieved_at=bundle["retrieved_at"],
                    release="Agent 1 query, response pages and normalized records pinned in native bundle")
    rows = []
    for alt in ALTERATIONS:
        # Agent 1 accepts single-source alteration jobs. Scope it to CIViC here, then
        # restore the complete workload source list and its missing-source reasons.
        source_job = {"task_id": TASK, "source_ids": ["civic"], "dimensions": {"alteration": alt}}
        for item in civic_evidence.records_for_job(bundle, source_job):
            rid = item["record_id"]
            row = {k: item[k] for k in ("source_id", "record_id", "url", "molecular_profile", "disease",
                                         "doid", *civic_evidence.FIELD_VALUES)}
            row.update(evidence_status=item["status"], is_flagged="false",
                       record_sha256=item["sha256"], page=item["page"])
            if mode == "fixture":
                row["url"] = "https://fixture.invalid/civic/" + rid
            rows.append(row)
    return _snapshot(mode, bundle, metadata, rows, GRAPHQL_SELECTION)


def _nightly_policy() -> dict:
    dataset = json.loads((ROOT / "catalog/datasets/civic-nightly-clinical-evidence.json").read_text())
    require(dataset["source_id"] == "civic" and dataset["retrieval_url"] == CIVIC and
            dataset["access"]["level"] == "open" and dataset["reuse"]["status"] == "verified-permissive" and
            dataset["reuse"]["license_id"] == "CC0-1.0" and
            dataset["reuse"]["terms_url"] == "https://civicdb.org/releases/licensing",
            "CIViC nightly access, reuse or retrieval allowlist changed")
    return {"connector": "nightly-tsv-fallback", "retrieval_status": "complete",
            "url": CIVIC, "release": dataset["release"]["label"], "terms_url": dataset["reuse"]["terms_url"],
            "dataset_catalog_sha256": digest(dataset)}


def _nightly(bundle):
    metadata = {**_source(bundle, "real"), **_nightly_policy()}
    raw = fetch(CIVIC, retries=0)
    parsed = parse_civic_snapshot(raw)
    if not set(FIELDS).issubset(next(iter(parsed["strict"].values()))):
        raise ValueError("CIViC projected columns changed")
    rows = []
    profiles = {_label(a) for a in ALTERATIONS}
    for r in parsed["strict"].values():
        if (_label(r["molecular_profile"]) not in profiles or _label(r["disease"]) != "lung adenocarcinoma" or
                canonical_doid(r["doid"]) != "3910" or _label(r["is_flagged"] or "") != "false"):
            continue
        eid = r["evidence_id"].strip()
        url = civic_evidence.LINK + eid
        require((r["evidence_civic_url"] or "").strip() in {url, f"https://civicdb.org/evidence/{eid}/summary"},
                "CIViC source record URL differs from evidence ID")
        row = {"source_id": "civic", "record_id": "EID" + eid, "url": url}
        row.update({k: (r[k] or "").strip() for k in FIELDS})
        rows.append(row)
    metadata.update(sha256=sha256(raw).hexdigest(), bytes=len(raw), retrieved_at=datetime.now(timezone.utc).isoformat(),
                    export_rows=parsed["total"], strict_rows=len(parsed["strict"]))
    snapshot = _snapshot("real", bundle, metadata, rows, NIGHTLY_SELECTION)
    validate_snapshot(snapshot, "real")
    return snapshot


def prepare_snapshot(mode: str) -> dict:
    """Prefer Agent 1. Retry eligible transport unavailability once through the nightly export."""
    require(mode in {"real", "fixture"}, "Unknown evidence mode")
    catalog = _source_catalog()
    if mode == "fixture":
        return from_agent1_bundle(civic_evidence.retrieve_fixture(catalog), mode)
    bundle = civic_evidence.retrieve(catalog)
    civic_evidence.validate_bundle(bundle, catalog)
    require(bundle["synthetic"] is False, "Agent 1 bundle mode mismatch")
    if bundle["status"] == "unavailable" and bundle["failure"] in FALLBACK_ERRORS:
        try:
            return _nightly(bundle)
        except (ValueError, OSError, RuntimeError, KeyError, TypeError) as exc:
            raise EvidenceUnavailable(bundle, fallback_attempted=True) from exc
    return from_agent1_bundle(bundle, mode)


def validate_snapshot(snapshot: dict, mode: str) -> None:
    require(mode in {"fixture", "real"} and isinstance(snapshot, dict) and
            set(snapshot) == {"schema_version", "mode", "synthetic", "selection", "sources", "records",
                              "snapshot_sha256", "agent1_bundle"} and snapshot["schema_version"] == 2 and
            snapshot["mode"] == mode and snapshot["synthetic"] is (mode == "fixture"), "Evidence snapshot mode differs")
    require(snapshot["snapshot_sha256"] == digest({k: v for k, v in snapshot.items() if k != "snapshot_sha256"}),
            "Evidence snapshot fingerprint mismatch")
    bundle = snapshot["agent1_bundle"]
    sources = snapshot["sources"]
    civic = sources["civic"]
    if civic["connector"] == "agent1-graphql":
        require(snapshot == from_agent1_bundle(bundle, mode), "Snapshot differs from Agent 1 projection")
        return
    civic_evidence.validate_bundle(bundle, _source_catalog())
    require(mode == "real" and not bundle["synthetic"] and bundle["status"] == "unavailable" and
            bundle["failure"] in FALLBACK_ERRORS and snapshot["selection"] == NIGHTLY_SELECTION and
            sources == {"civic": civic, **_source_gaps()}, "Invalid fallback source or coverage")
    policy = {**_source(bundle, mode), **_nightly_policy()}
    require(set(civic) == {*policy, "sha256", "bytes", "retrieved_at", "export_rows", "strict_rows",
                           "matching_by_alteration"} and all(civic[k] == v for k, v in policy.items()) and
            isinstance(civic["sha256"], str) and SHA.fullmatch(civic["sha256"]) and
            type(civic["bytes"]) is int and 0 < civic["bytes"] <= MAX_SOURCE_BYTES and
            isinstance(civic["retrieved_at"], str) and datetime.fromisoformat(civic["retrieved_at"]).tzinfo,
            "CIViC nightly provenance mismatch")
    counts = civic["matching_by_alteration"]
    require(isinstance(counts, dict) and set(counts) == set(ALTERATIONS) and
            all(type(n) is int and n >= 0 for n in counts.values()) and
            all(type(civic[k]) is int for k in ("export_rows", "strict_rows")) and
            0 < civic["strict_rows"] <= civic["export_rows"] and sum(counts.values()) <= civic["strict_rows"],
            "Invalid nightly selection counts")
    records = snapshot["records"]
    require(isinstance(records, list) and len(records) == sum(min(LIMIT, n) for n in counts.values()),
            "Invalid number of projected records")
    ids = set()
    included = defaultdict(int)
    for row in records:
        require(isinstance(row, dict) and set(row) == {"source_id", "record_id", "url", *FIELDS} and
                all(isinstance(v, str) and len(v) <= 512 for v in row.values()), "Invalid projected record")
        rid = row["record_id"]
        eid = rid.removeprefix("EID")
        require(rid == "EID" + eid and eid.isascii() and eid.isdigit() and int(eid) > 0 and rid not in ids and
                row["source_id"] == "civic" and row["url"] == civic_evidence.LINK + eid and
                _label(row["disease"]) == "lung adenocarcinoma" and canonical_doid(row["doid"]) == "3910" and
                _label(row["is_flagged"]) == "false" and _label(row["evidence_status"]) == "accepted",
                "Invalid evidence identifier or scope")
        ids.add(rid)
        included[_label(row["molecular_profile"])] += 1
    require(all(included[_label(k)] == min(LIMIT, n) for k, n in counts.items()), "Alteration scope mismatch")


def evidence_for_job(snapshot: dict, job: dict) -> dict:
    """Supply only this alteration's projected records and source-level provenance."""
    alteration = job["dimensions"]["alteration"]
    require(job["task_id"] == TASK and set(job["source_ids"]) == set(snapshot["sources"]) and
            alteration in ALTERATIONS, "Job outside EGFR/LUAD pilot scope")
    records = [row for row in snapshot["records"] if _label(row["molecular_profile"]) == _label(alteration)]
    civic = dict(snapshot["sources"]["civic"])
    matched = civic.pop("matching_by_alteration")[alteration]
    reason = f"{matched} exact-label matches in retrieved input. {len(records)} projected evidence items supplied."
    if civic["retrieval_status"] == "truncated":
        reason += " Pagination was truncated. Additional matches and their annotations are unknown."
    civic.update(matching_records=matched, included_records=len(records), reason=reason)
    return {"snapshot_sha256": snapshot["snapshot_sha256"], "synthetic": snapshot["synthetic"],
            "alteration": alteration, "selection": snapshot["selection"], "records": records,
            "sources": {**snapshot["sources"], "civic": civic}}
