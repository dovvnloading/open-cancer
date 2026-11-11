#!/usr/bin/env python3
"""Fail closed when a public Explorer API is stale, malformed or unreachable.

This is a release-health check, NOT a biomedical or clinical validation.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

MAX_INDEX_BYTES = 4 * 1024 * 1024


def validate_index(index: object, *, revision: str, repository: str, minimum_datasets: int = 1) -> list[str]:
    errors = []
    if not isinstance(index, dict):
        return ["Expected an Explorer JSON index object"]
    if index.get("version") != 1:
        errors.append("Wrong Explorer JSON API version")
    if index.get("repository") != repository:
        errors.append("Unexpected Explorer repository identity")
    observed = index.get("revision")
    if observed != revision:
        errors.append(f"STALE Explorer: deployed revision {observed!r}, expected {revision!r}")
    for key in ("datasets", "sources", "questions", "pilots", "reviews"):
        if not isinstance(index.get(key), list):
            errors.append(f"Missing or malformed collection: {key}")
    if isinstance(index.get("datasets"), list) and len(index["datasets"]) < minimum_datasets:
        errors.append("Explorer catalog has fewer records than expected")
    return errors


def check_remote(url: str, *, revision: str, repository: str, minimum_datasets: int = 1,
                 opener=urlopen) -> list[str]:
    parsed = urlsplit(url)
    if parsed.scheme != "https" or not parsed.netloc or parsed.username or parsed.password or parsed.fragment:
        return ["Expected a public HTTPS JSON endpoint with no embedded credentials or fragment"]
    try:
        with opener(Request(url, headers={"Accept": "application/json",
                                          "User-Agent": "open-cancer-release-qa/1"}), timeout=20) as resp:
            payload = resp.read(MAX_INDEX_BYTES + 1)
        if len(payload) > MAX_INDEX_BYTES:
            return ["Explorer index exceeds bounded response size"]
        data = json.loads(payload)
    except (OSError, ValueError, UnicodeError) as exc:
        return [f"Could not fetch and parse public Explorer index: {exc}"]
    return validate_index(data, revision=revision, repository=repository, minimum_datasets=minimum_datasets)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", required=True, help="Live public index.json HTTPS URL")
    parser.add_argument("--expected-revision", help="Full Git commit SHA; defaults to current checkout")
    parser.add_argument("--repository", default="dovvnloading/open-cancer")
    parser.add_argument("--minimum-datasets", type=int, default=69)
    args = parser.parse_args()
    try:
        revision = args.expected_revision or subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f"Could not determine source revision: {exc}\n")
    if not re.fullmatch(r"[0-9a-f]{40}", revision) or args.minimum_datasets < 1:
        parser.exit(2, "Expected full SHA-1 Git revision and positive dataset minimum\n")
    errors = check_remote(args.url, revision=revision, repository=args.repository,
                          minimum_datasets=args.minimum_datasets)
    for error in errors:
        print(f"FAIL: {error}")
    if errors:
        return 1
    print(f"OK: public Explorer index corresponds to {args.repository}@{revision} "
          f"and contains at least {args.minimum_datasets} dataset records.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
