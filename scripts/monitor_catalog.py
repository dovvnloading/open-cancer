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
FIELDS = (
    "accession", "source_id", "title", "description", "disease", "modalities",
    "cohort_definition", "unit_of_observation", "release", "reuse", "access",
    "retrieval_url", "landing_url", "documentation_url",
)


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


def compare(source, old_rows, new_rows):
    if not new_rows:
        raise ValueError(f"{source}: zero upstream records")
    ids = [row.get("id") for row in new_rows]
    if any(not isinstance(key, str) for key in ids) or len(ids) != len(set(ids)):
        raise ValueError(f"{source}: duplicate or invalid upstream identifiers")
    if any(row.get("source_id") != source or not managed(source, row["id"]) for row in new_rows):
        raise ValueError(f"{source}: unexpected upstream source or identifier")
    old = {row["id"]: row for row in old_rows}
    new = {row["id"]: row for row in new_rows}
    changed = {key: [field for field in FIELDS if old[key].get(field) != new[key].get(field)]
               for key in sorted(old.keys() & new.keys())}
    return {
        "committed_count": len(old),
        "upstream_count": len(new),
        "unexpected_shrinkage": len(new) < len(old),
        "added": sorted(new.keys() - old.keys()),
        "removed": sorted(old.keys() - new.keys()),
        "changed_fields": {key: fields for key, fields in changed.items() if fields},
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


def scan(*, current, discover, terms_urls=None, probe=probe_terms, checked_at=None, sleep=time.sleep):
    now = checked_at or datetime.now(timezone.utc).isoformat(timespec="seconds")
    report = {"schema_version": 1, "checked_at": now, "sources": {}}
    for source in SOURCES:
        row = {"last_checked": now, "ok": False}
        try:
            items = retry(lambda: discover(source), sleep=sleep)
            if isinstance(items, dict):
                items = list(items.values())
            row["drift"] = compare(source, current[source], items)
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
    args = parser.parse_args()
    terms = {
        name: json.loads((ROOT / "catalog" / "sources" / f"{name}.json").read_text())["reuse"]["terms_url"]
        for name in SOURCES
    }
    report = scan(current=baseline(ROOT / "catalog" / "datasets"),
                  discover=refresh_catalog.discover, terms_urls=terms)
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
