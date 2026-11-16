"""Anonymous, bounded retrieval of accepted CIViC annotations for the EGFR/LUAD task."""
from __future__ import annotations

import base64
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
import hashlib
from http.client import HTTPException
import json
from pathlib import Path
import re
import time
from urllib.error import HTTPError, URLError
from urllib.request import HTTPRedirectHandler, Request, build_opener

ROOT = Path(__file__).resolve().parents[1]
CATALOG = ROOT / "catalog/sources/civic.json"
FIXTURE = ROOT / "tests/fixtures/civic/egfr_luad_page.json"
ENDPOINT = "https://civicdb.org/api/graphql"
TERMS = "https://docs.civicdb.org/en/latest/about.html"
LINK = "https://civicdb.org/links/evidence_items/"
TASK = "egfr-luad-evidence-coverage"
PROFILES = {name.casefold(): name for name in (
    "EGFR L858R", "EGFR Exon 19 Deletion", "EGFR Exon 20 Insertion")}
DISEASE = "Lung Adenocarcinoma"
QUERY = """query OpenCancerEgfrLuad($first: Int!, $after: String) {
  evidenceItems(first: $first, after: $after, status: ACCEPTED,
                molecularProfileName: "EGFR", diseaseName: "Lung Adenocarcinoma",
                sortBy: {column: ID, direction: ASC}) {
    totalCount
    nodes {
      id status flagged evidenceType evidenceLevel evidenceDirection significance
      molecularProfile { name deprecated }
      disease { name doid deprecated }
    }
    pageInfo { hasNextPage endCursor }
  }
}"""

FIELD_VALUES = {
    "evidence_level": frozenset("ABCDE"),
    "evidence_type": frozenset("PREDICTIVE DIAGNOSTIC PROGNOSTIC PREDISPOSING FUNCTIONAL ONCOGENIC".split()),
    "evidence_direction": frozenset(("SUPPORTS", "DOES_NOT_SUPPORT", "NA")),
    "significance": frozenset("""
        SENSITIVITYRESPONSE RESISTANCE BETTER_OUTCOME POOR_OUTCOME POSITIVE NEGATIVE NA
        ADVERSE_RESPONSE PATHOGENIC LIKELY_PATHOGENIC BENIGN LIKELY_BENIGN UNCERTAIN_SIGNIFICANCE
        REDUCED_SENSITIVITY GAIN_OF_FUNCTION LOSS_OF_FUNCTION UNALTERED_FUNCTION NEOMORPHIC
        UNKNOWN DOMINANT_NEGATIVE PREDISPOSITION PROTECTIVENESS ONCOGENICITY LIKELY_ONCOGENIC
    """.split()),
}
ANNOTATIONS = {"evidence_level": "evidenceLevel", "evidence_type": "evidenceType",
               "evidence_direction": "evidenceDirection", "significance": "significance"}
MAX_PAGES, MAX_PAGE_SIZE = 8, 25
MAX_PAGE_BYTES = 128 * 1024
MAX_ATTEMPTS, TIMEOUT, MAX_BACKOFF = 3, 8, 10
FAILURES = {"rate-limit", "network-error", "http-error", "malformed-response", "source-changed"}
HASH = re.compile(r"[0-9a-f]{64}")
TIMESTAMP = re.compile(r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z")


def _digest(value):
    data = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)
    return hashlib.sha256(data.encode("utf-8")).hexdigest()


def _check(condition, message="Invalid CIViC snapshot"):
    if not condition:
        raise ValueError(message)


def _integer(value, low, high):
    return type(value) is int and low <= value <= high


def _keys(value, names):
    return isinstance(value, dict) and value.keys() == set(names)


def _permission(catalog):
    if catalog is None:
        catalog = json.loads(CATALOG.read_text(encoding="utf-8"))
    _check(isinstance(catalog, dict) and catalog.get("id") == "civic" and
           catalog.get("access") == "open" and catalog.get("api_url") == ENDPOINT and
           isinstance(catalog.get("reuse"), dict) and
           catalog["reuse"].get("status") == "verified-permissive" and
           catalog["reuse"].get("license_id") == "CC0-1.0" and
           catalog["reuse"].get("terms_url") == TERMS, "CIViC catalog permission or endpoint changed")


def supports_job(job):
    dims, ids = job.get("dimensions"), job.get("source_ids")
    return (job.get("task_id") == TASK and ids == ["civic"] and
            isinstance(dims, dict) and dims.keys() == {"alteration"} and
            isinstance(dims.get("alteration"), str) and
            dims["alteration"].casefold() in PROFILES)


def records_for_job(bundle, job):
    """Select only exact single-alteration records from an already validated snapshot."""
    if not supports_job(job):
        return []
    profile = PROFILES[job["dimensions"]["alteration"].casefold()]
    return [r for r in bundle["records"] if r["molecular_profile"] == profile]


class _Unavailable(Exception):
    """Closed reason code only. Upstream messages, headers and bodies are never logged."""


class _NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def _retry_delay(header, attempt):
    if header is None:
        return 2 ** attempt
    try:
        if re.fullmatch(r"[0-9]{1,9}", header.strip()):
            return int(header)
        dt = parsedate_to_datetime(header)
        return max(0, (dt - datetime.now(timezone.utc)).total_seconds())
    except (ValueError, TypeError, OverflowError):
        # A malformed hint is not permission to retry early.
        return MAX_BACKOFF + 1


def _fetch_page(variables):
    request = Request(ENDPOINT, data=json.dumps({"query": QUERY, "variables": variables}).encode(),
                      headers={"Content-Type": "application/json", "Accept": "application/json",
                               "User-Agent": "open-cancer-evidence/1.0"})
    opener = build_opener(_NoRedirect())
    for attempt in range(MAX_ATTEMPTS):
        delay = 2 ** attempt
        try:
            deadline = time.monotonic() + TIMEOUT
            with opener.open(request, timeout=TIMEOUT) as response:
                if response.status != 200 or response.geturl() != ENDPOINT:
                    raise _Unavailable("http-error")
                if response.headers.get_content_type() not in {
                        "application/json", "application/graphql-response+json"}:
                    raise _Unavailable("malformed-response")
                body = bytearray()
                while len(body) <= MAX_PAGE_BYTES:
                    if time.monotonic() >= deadline:
                        raise TimeoutError
                    chunk = response.read1(min(16_384, MAX_PAGE_BYTES + 1 - len(body)))
                    if not chunk:
                        return bytes(body)
                    body.extend(chunk)
                raise _Unavailable("malformed-response")
        except HTTPError as exc:
            with exc:
                reason = "rate-limit" if exc.code == 429 else "http-error"
                if exc.code not in {429, 502, 503, 504}:
                    raise _Unavailable(reason) from None
                delay = max(delay, _retry_delay(exc.headers.get("Retry-After"), attempt))
        except (URLError, OSError, HTTPException):
            reason = "network-error"
        if attempt + 1 == MAX_ATTEMPTS or delay > MAX_BACKOFF:
            raise _Unavailable(reason)
        time.sleep(delay)


def _cursor(count):
    # CIViC uses graphql-ruby offset cursors (base64url of an ASCII integer).
    # Persisting opaque arbitrary upstream strings could carry secrets or free text.
    return base64.urlsafe_b64encode(str(count).encode("ascii")).decode("ascii").rstrip("=")


def _parse_page(raw, first, offset):
    if len(raw) > MAX_PAGE_BYTES:
        raise _Unavailable("malformed-response")

    def unique_object(pairs):
        obj = dict(pairs)
        _check(len(obj) == len(pairs))
        return obj

    page = json.loads(raw.decode("utf-8"), object_pairs_hook=unique_object)
    _check(isinstance(page, dict) and "errors" not in page)
    data = page["data"]["evidenceItems"]
    nodes, info, total = data["nodes"], data["pageInfo"], data["totalCount"]
    _check(_integer(total, 0, 2**31 - 1))
    _check(isinstance(nodes, list) and len(nodes) <= first)
    _check(type(info["hasNextPage"]) is bool)
    if not nodes:
        _check(not info["hasNextPage"] and info["endCursor"] is None)
    else:
        _check(info["endCursor"] == _cursor(offset + len(nodes)))
    # The collection can change in between Relay offset queries.
    if ((not nodes and offset) or offset + len(nodes) > total or
            info["hasNextPage"] != (offset + len(nodes) < total)):
        raise _Unavailable("source-changed")
    return nodes, info["hasNextPage"], total


def _project(node, index):
    """Only controlled annotation values and constructed public identifiers can leave here."""
    eid = node["id"]
    _check(_integer(eid, 1, 2**31 - 1))
    _check(node["status"] in ("ACCEPTED", "SUBMITTED", "REJECTED"))
    _check(type(node["flagged"]) is bool)
    mp, disease = node["molecularProfile"], node["disease"]
    _check(isinstance(mp["name"], str) and type(mp["deprecated"]) is bool)
    if disease is None:
        return None
    _check(isinstance(disease["name"], str) and type(disease["deprecated"]) is bool and
           (disease["doid"] is None or isinstance(disease["doid"], str)))
    profile = PROFILES.get(mp["name"].casefold())
    if (node["status"] != "ACCEPTED" or node["flagged"] or mp["deprecated"] or
            disease["deprecated"] or profile is None or disease["doid"] not in ("3910", "DOID:3910")
            or disease["name"].casefold() != DISEASE.casefold()):
        return None
    record = {"source_id": "civic", "record_id": f"EID{eid}", "url": LINK + str(eid),
              "status": "ACCEPTED", "molecular_profile": profile, "disease": DISEASE, "doid": "DOID:3910"}
    for key, upstream in ANNOTATIONS.items():
        value = node[upstream]
        _check(isinstance(value, str) and value in FIELD_VALUES[key])
        record[key] = value
    return {**record, "sha256": _digest(record), "page": index}


def _snapshot(catalog, max_pages, page_size, fetch, retrieved_at, synthetic):
    _permission(catalog)
    _check(_integer(max_pages, 1, MAX_PAGES) and _integer(page_size, 1, MAX_PAGE_SIZE),
           "Invalid CIViC retrieval limits")
    bundle = {"schema_version": 1, "source_id": "civic", "license_id": "CC0-1.0", "terms_url": TERMS,
              "endpoint": ENDPOINT, "query": QUERY, "page_size": page_size, "max_pages": max_pages,
              "retrieved_at": retrieved_at, "synthetic": synthetic, "status": "complete",
              "failure": None, "failed_variables": None, "pages": [], "records": []}
    offset = previous_id = 0
    try:
        for index in range(max_pages):
            if index and not synthetic:
                time.sleep(0.35)  # Anonymous limit is 3 requests/s, measured over five minutes.
            variables = {"first": page_size, "after": _cursor(offset) if index else None}
            raw = fetch(variables)
            nodes, more, total = _parse_page(raw, page_size, offset)
            if index and total != bundle["pages"][0]["total_count"]:
                raise _Unavailable("source-changed")
            records = []
            for node in nodes:
                record = _project(node, index)
                if node["id"] <= previous_id:
                    raise _Unavailable("source-changed")
                previous_id = node["id"]
                if record is not None:
                    records.append(record)
            bundle["pages"].append({"variables": variables, "node_count": len(nodes),
                                    "has_next": more, "total_count": total,
                                    "response_sha256": hashlib.sha256(raw).hexdigest()})
            bundle["records"].extend(records)
            offset += len(nodes)
            if not more:
                break
        else:
            bundle["status"] = "truncated"
    except (_Unavailable, ValueError, TypeError, KeyError, RecursionError) as exc:
        bundle["status"] = "unavailable"
        bundle["failed_variables"] = variables
        bundle["failure"] = str(exc) if isinstance(exc, _Unavailable) else "malformed-response"
        bundle["records"].clear()  # A late failure cannot masquerade as a usable subset.
    bundle["sha256"] = _digest(bundle)
    return validate_bundle(bundle, catalog)


def retrieve(catalog, *, max_pages=4, page_size=20):
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    return _snapshot(catalog, max_pages, page_size, _fetch_page, stamp, False)


def retrieve_fixture(catalog):
    """Exercise the same decoder and filters with an explicitly synthetic offline response."""
    return _snapshot(catalog, 1, 20, lambda _: FIXTURE.read_bytes(), "2026-01-01T00:00:00Z", True)


def validate_bundle(bundle, catalog=None):
    """Recheck local rights, safe schema, query lineage and checksums on every replay."""
    _permission(catalog)
    fields = ("schema_version source_id license_id terms_url endpoint query page_size max_pages "
              "retrieved_at synthetic status failure failed_variables pages records sha256").split()
    _check(_keys(bundle, fields))
    _check(_integer(bundle["schema_version"], 1, 1) and bundle["source_id"] == "civic" and
           bundle["license_id"] == "CC0-1.0" and bundle["terms_url"] == TERMS and
           bundle["endpoint"] == ENDPOINT and bundle["query"] == QUERY)
    _check(isinstance(bundle["retrieved_at"], str) and TIMESTAMP.fullmatch(bundle["retrieved_at"]))
    datetime.strptime(bundle["retrieved_at"], "%Y-%m-%dT%H:%M:%SZ")
    _check(type(bundle["synthetic"]) is bool)
    size, limit, pages, records = (bundle[k] for k in ("page_size", "max_pages", "pages", "records"))
    _check(_integer(size, 1, MAX_PAGE_SIZE) and _integer(limit, 1, MAX_PAGES) and
           isinstance(pages, list) and len(pages) <= limit and isinstance(records, list) and
           len(records) <= len(pages) * size)
    offset = 0
    for i, page in enumerate(pages):
        _check(_keys(page, ("variables", "node_count", "has_next", "total_count", "response_sha256")))
        variables = page["variables"]
        _check(_keys(variables, ("first", "after")) and _integer(variables["first"], size, size) and
               variables["after"] == (_cursor(offset) if i else None))
        _check(_integer(page["node_count"], 0, size) and type(page["has_next"]) is bool and
               isinstance(page["response_sha256"], str) and HASH.fullmatch(page["response_sha256"]))
        _check(not page["has_next"] or page["node_count"] > 0)
        _check(i == len(pages) - 1 or page["has_next"])
        offset += page["node_count"]
        total = page["total_count"]
        _check(_integer(total, offset, 2**31 - 1) and total == pages[0]["total_count"] and
               page["has_next"] == (offset < total) and (not i or page["node_count"] > 0))
    status = bundle["status"]
    if status == "unavailable":
        failed = bundle["failed_variables"]
        _check(_keys(failed, ("first", "after")) and _integer(failed["first"], size, size) and
               failed["after"] == (_cursor(offset) if pages else None))
        _check(not records and isinstance(bundle["failure"], str) and bundle["failure"] in FAILURES and
               len(pages) < limit and (not pages or pages[-1]["has_next"]))
    else:
        _check(bundle["failure"] is None and bundle["failed_variables"] is None and pages and
               ((status == "complete" and not pages[-1]["has_next"]) or
                (status == "truncated" and len(pages) == limit and pages[-1]["has_next"])))
    last_id, last_page = 0, 0
    counts = [0] * len(pages)
    for record in records:
        _check(_keys(record, {"source_id", "record_id", "url", "sha256", "page", "molecular_profile",
                              "disease", "doid", "status", *FIELD_VALUES}))
        _check(isinstance(record["record_id"], str) and re.fullmatch(r"EID[1-9][0-9]{0,9}", record["record_id"]))
        eid = int(record["record_id"][3:])
        _check(last_id < eid < 2**31 and record["url"] == LINK + str(eid) and
               record["source_id"] == "civic" and record["status"] == "ACCEPTED" and
               record["molecular_profile"] in PROFILES.values() and record["disease"] == DISEASE and
               record["doid"] == "DOID:3910" and _integer(record["page"], last_page, len(pages) - 1))
        for key, values in FIELD_VALUES.items():
            _check(isinstance(record[key], str) and record[key] in values)
        _check(record["sha256"] == _digest({k: v for k, v in record.items() if k not in ("page", "sha256")}))
        last_page, last_id = record["page"], eid
        counts[last_page] += 1
        _check(counts[last_page] <= pages[last_page]["node_count"])
    _check(bundle["sha256"] == _digest({k: v for k, v in bundle.items() if k != "sha256"}))
    return bundle
