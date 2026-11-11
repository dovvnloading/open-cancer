#!/usr/bin/env python3
"""Smoke-check a built or deployed Explorer snapshot, without scientific inference."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
from urllib.parse import urljoin, urlsplit
from urllib.request import Request, urlopen

MAX_RESPONSE = 4 * 1024 * 1024
EXPECTED_DATASET_IDS = {
    "gdc-tcga-luad-project",
    "gdc-tcga-luad-egfr-ssm-occurrences",
    "cbioportal-luad-tcga-pancan-atlas-2018",
    "civic-nightly-clinical-evidence",
}
REPO = "dovvnloading/open-cancer"


def read_json(fetch, name):
    data = fetch("api/v1/" + name + ".json")
    if len(data) > MAX_RESPONSE:
        raise ValueError(name + " response exceeds size limit")
    return json.loads(data)


def check_snapshot(fetch, revision, *, minimum_datasets=69, require_build_info=False):
    """Inspect full API roundtrips and selected UI links, not biology or clinical claims."""
    errors = []
    if not re.fullmatch(r"[0-9a-f]{40}", revision):
        return ["Invalid expected Git revision"]
    try:
        index = read_json(fetch, "index")
        if not isinstance(index, dict):
            return ["Explorer index is not a JSON object"]
        if index.get("version") != 1 or index.get("repository") != REPO:
            errors.append("Unexpected index version or repository identity")
        if index.get("revision") != revision:
            errors.append("Deployed Explorer revision differs from requested commit")
        collections = ("datasets", "sources", "questions", "reviews", "pilots",
                       "experimental_results")
        for name in collections:
            payload = read_json(fetch, name)
            if not isinstance(payload, list) or payload != index.get(name):
                errors.append(name + " endpoint differs from index or is not an array")
        datasets = index.get("datasets", [])
        if not isinstance(datasets, list) or len(datasets) < minimum_datasets:
            errors.append("Dataset collection smaller than expected")
            datasets = []
        records = {row.get("id"): row for row in datasets if isinstance(row, dict)}
        if not EXPECTED_DATASET_IDS.issubset(records):
            errors.append("Missing linked EGFR/LUAD GDC, cBioPortal or CIViC datasets")
        for id_ in EXPECTED_DATASET_IDS & records.keys():
            row = records[id_]
            links = row.get("research_context", [])
            if not any(item.get("task_id") == "egfr-luad-evidence-coverage"
                       and item.get("gene_symbol") == "EGFR"
                       for item in links if isinstance(item, dict)):
                errors.append("Missing EGFR research crosswalk for " + id_)
            source_url = row.get("record_url", "")
            if f"/blob/{revision}/catalog/datasets/" not in source_url:
                errors.append("Non-version-bound dataset link for " + id_)
        pilots = index.get("pilots", [])
        pilot = next((row for row in pilots if isinstance(row, dict)
                      and row.get("task_id") == "egfr-luad-evidence-coverage"
                      and row.get("run_id") == "pilot-20261110"), None)
        if pilot is None:
            errors.append("Featured EGFR/LUAD pilot missing")
        else:
            if pilot.get("coverage", {}).get("status") != "unreviewed":
                errors.append("Original EGFR/LUAD snapshot not marked unreviewed")
            if pilot.get("clinical_validation_status") != "not-assessed":
                errors.append("Clinical nonvalidation disclaimer missing")
            if pilot.get("independent_reviewer_credentials_verified") is not False:
                errors.append("Independent reviewer credential caveat missing")
            if pilot.get("publication_status") != "legacy-featured":
                errors.append("Historical pilot publication status changed")
            if f"/blob/{revision}/results/" not in pilot.get("manifest_url", ""):
                errors.append("Pilot manifest link is not version bound")
            reviews = index.get("reviews", [])
            for review in reviews:
                if (review.get("task_id"), review.get("run_id")) == (
                        pilot["task_id"], pilot["run_id"]):
                    if review.get("reviewed_manifest_sha256") != pilot.get("manifest_sha256"):
                        errors.append("Review points at the wrong pilot manifest")
                    if f"/blob/{revision}/reviews/" not in review.get("record_url", ""):
                        errors.append("Review record link is not version bound")
        html = fetch("index.html").decode("utf-8")
        js = fetch("app.js").decode("utf-8")
        if "Open Cancer Explorer" not in html or "api/v1/index.json" not in js:
            errors.append("Missing Explorer page or API loading code")
        if require_build_info:
            build = read_json(fetch, "release")
            if build.get("revision") != revision or build.get("repository") != REPO:
                errors.append("Build metadata does not match deployed revision")
            value = build.get("built_at")
            try:
                timestamp = datetime.fromisoformat(value.replace("Z", "+00:00"))
                if timestamp.tzinfo is None or timestamp > datetime.now(timezone.utc):
                    raise ValueError("timestamp missing timezone or in future")
            except (AttributeError, TypeError, ValueError):
                errors.append("Missing/invalid UTC build timestamp")
    except (OSError, ValueError, UnicodeError, KeyError, TypeError) as exc:
        errors.append(f"Missing or malformed published asset: {exc}")
    return errors


def local_fetch(site):
    root = site.resolve()
    def fetch(path):
        candidate = (root / path).resolve()
        if not candidate.is_relative_to(root) or not candidate.is_file():
            raise ValueError("Unsafe or missing local release asset: " + path)
        if candidate.stat().st_size > MAX_RESPONSE:
            raise ValueError("Local release asset exceeds limit: " + path)
        return candidate.read_bytes()
    return fetch


def remote_fetch(site_url, opener=urlopen):
    url = urlsplit(site_url)
    if (url.scheme != "https" or not url.netloc or url.username or url.password
            or url.query or url.fragment):
        raise ValueError("Expected a public HTTPS site root without query or credentials")
    base = site_url.rstrip("/") + "/"
    def fetch(path):
        request = Request(urljoin(base, path), headers={
            "Accept": "application/json, text/html, text/javascript",
            "User-Agent": "open-cancer-explorer-smoke/1"})
        with opener(request, timeout=20) as response:
            return response.read(MAX_RESPONSE + 1)
    return fetch


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    target = parser.add_mutually_exclusive_group(required=True)
    target.add_argument("--site", type=Path, help="Local built site directory")
    target.add_argument("--url", help="Deployed HTTPS site root")
    parser.add_argument("--expected-revision", required=True)
    parser.add_argument("--minimum-datasets", type=int, default=69)
    parser.add_argument("--require-build-info", action="store_true")
    args = parser.parse_args()
    if args.minimum_datasets < 1:
        parser.error("Expected a positive dataset count")
    try:
        fetch = local_fetch(args.site) if args.site else remote_fetch(args.url)
        errors = check_snapshot(fetch, args.expected_revision,
                                minimum_datasets=args.minimum_datasets,
                                require_build_info=args.require_build_info)
    except ValueError as exc:
        errors = [str(exc)]
    if errors:
        for error in errors:
            print("FAIL: " + error)
        return 1
    print("OK: Explorer version, API collections, research crosswalk, pilot and review links, "
          "and optional build metadata are coherent. No scientific validation implied.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
