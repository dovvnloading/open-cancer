#!/usr/bin/env python3
"""Monitor upstream public metadata without editing curated catalog records."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from http.client import HTTPException
import json
from pathlib import Path
import socket
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

try:
    from . import refresh_catalog
except ImportError:
    import refresh_catalog

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ("gdc", "cbioportal", "civic")
# Only compare connector-observed metadata with the same connector's *previous*
# observations. Curated descriptions, ontology mappings and rights statements
# are editorial and must not be treated as source-native upstream changes.
OBSERVED_FIELDS = (
    "accession", "title", "disease", "modalities", "release",
    "citation_url", "retrieval_url", "landing_url",
)


def observation(record):
    return {field: record[field] for field in OBSERVED_FIELDS if field in record}


def managed(source, key):
    if source == "gdc":
        return key.startswith("gdc-tcga-") and key.endswith("-project")
    if source == "cbioportal":
        return key.startswith("cbioportal-") and key.endswith("-tcga-pancan-atlas-2018")
    return source == "civic" and key.startswith("civic-nightly-")


def baseline(directory):
    result = {source: [] for source in SOURCES}
    for path in sorted(directory.glob("*.json")):
        row = json.loads(path.read_text(encoding="utf-8"))
        source = row.get("source_id")
        if source in result and managed(source, row.get("id", "")):
            result[source].append(row)
    return result


def compare(source, old_rows, new_rows, *, previous_observations=None):
    if not new_rows:
        raise ValueError(f"{source}: zero upstream records")
    ids = [row.get("id") for row in new_rows]
    if any(not isinstance(key, str) for key in ids) or len(ids) != len(set(ids)):
        raise ValueError(f"{source}: duplicate or invalid upstream identifiers")
    if any(row.get("source_id") != source or not managed(source, row["id"]) for row in new_rows):
        raise ValueError(f"{source}: unexpected upstream source or identifier")
    committed_ids = {row["id"] for row in old_rows}
    observed = {row["id"]: observation(row) for row in new_rows}
    if previous_observations is not None and not isinstance(previous_observations, dict):
        raise ValueError(f"{source}: corrupt previous observation ledger")
    changes = {}
    for key in sorted(observed.keys() & (previous_observations or {}).keys()):
        old = previous_observations[key]
        if not isinstance(old, dict):
            raise ValueError(f"{source}: invalid previous observation for {key}")
        fields = [name for name in OBSERVED_FIELDS if old.get(name) != observed[key].get(name)]
        if fields:
            changes[key] = fields
    return {
        "committed_count": len(committed_ids),
        "upstream_count": len(observed),
        "unexpected_shrinkage": len(observed) < len(committed_ids),
        "added": sorted(observed.keys() - committed_ids),
        "removed": sorted(committed_ids - observed.keys()),
        "changed_fields": changes,
        "baseline_established": previous_observations is not None,
        "observations": observed,
    }


def retry(operation, attempts=3, sleep=time.sleep):
    for attempt in range(attempts):
        try:
            return operation()
        except (HTTPError, URLError, TimeoutError, socket.timeout, HTTPException) as exc:
            if isinstance(exc, HTTPError) and exc.code not in (408, 429, 500, 502, 503, 504):
                raise
            if attempt == attempts - 1:
                raise
            sleep(min(2 ** attempt, 4))


def probe_terms(url):
    if not url.startswith("https://"):
        raise ValueError("source terms must use HTTPS")
    headers = {"User-Agent": "open-cancer-catalog-health/1"}
    try:
        response = urlopen(Request(url, headers=headers, method="HEAD"), timeout=12)
    except HTTPError as exc:
        if exc.code != 405:
            raise
        response = urlopen(Request(url, headers=headers | {"Range": "bytes=0-0"}), timeout=12)
    with response:
        return {"url": url, "http_status": response.status,
                "etag": response.headers.get("ETag"),
                "last_modified": response.headers.get("Last-Modified")}


def scan(*, current, discover, terms_urls=None, probe=probe_terms, checked_at=None,
         previous_report=None, sleep=time.sleep):
    now = checked_at or datetime.now(timezone.utc).isoformat(timespec="seconds")
    report = {"schema_version": 1, "checked_at": now, "sources": {}}
    for source in SOURCES:
        row = {"last_checked": now, "ok": False}
        previous = ((previous_report or {}).get("sources") or {}).get(source, {})
        old_observations = previous.get("observations") if isinstance(previous, dict) else None
        try:
            items = retry(lambda: discover(source), sleep=sleep)
            if isinstance(items, dict):
                items = list(items.values())
            diff = compare(source, current[source], items,
                           previous_observations=old_observations)
            row["observations"] = diff.pop("observations")
            row["drift"] = diff
            row["ok"] = True
        except (HTTPError, URLError, TimeoutError, socket.timeout, HTTPException,
                ValueError, KeyError, OSError, TypeError) as exc:
            row["error"] = f"{type(exc).__name__}: {exc}"[:400]
        if terms_urls and source in terms_urls:
            try:
                row["terms_link"] = retry(lambda: probe(terms_urls[source]), sleep=sleep)
            except (HTTPError, URLError, TimeoutError, socket.timeout, HTTPException,
                    ValueError, KeyError, OSError, TypeError) as exc:
                row["terms_link_error"] = f"{type(exc).__name__}: {exc}"[:400]
        report["sources"][source] = row
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--previous-report", type=Path,
                        help="Last durable health ledger report; never an upstream data source")
    args = parser.parse_args()
    previous = None
    if args.previous_report and args.previous_report.is_file():
        previous = json.loads(args.previous_report.read_text(encoding="utf-8"))
    terms = {
        name: json.loads((ROOT / "catalog" / "sources" / f"{name}.json").read_text())["reuse"]["terms_url"]
        for name in SOURCES
    }
    report = scan(current=baseline(ROOT / "catalog" / "datasets"),
                  discover=refresh_catalog.discover, terms_urls=terms,
                  previous_report=previous)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    for name, row in report["sources"].items():
        d = row.get("drift", {})
        print(f"{name}: {'available' if row['ok'] else 'FAILED'}; "
              f"{d.get('committed_count', '?')} committed / {d.get('upstream_count', '?')} upstream; "
              f"{len(d.get('added', []))} additions, {len(d.get('removed', []))} removals, "
              f"{len(d.get('changed_fields', {}))} changed records")
    print("Results describe metadata health only, not scientific or licensing verification.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
