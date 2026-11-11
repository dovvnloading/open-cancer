#!/usr/bin/env python3
"""Build the public Explorer and its versioned JSON API from the local checkout."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from urllib.parse import quote

if __package__:
    from .validate import ROOT, validate_repository
else:
    from validate import ROOT, validate_repository


SITE_ASSETS = ("index.html", "style.css", "app.js", "catalog.mjs")
RECORDS = {
    "datasets": "catalog/datasets/*.json",
    "sources": "catalog/sources/*.json",
    "questions": "tasks/*.json",
    "reviews": "reviews/*/*/*.json",
}


def make_index(root: Path, repository: str, revision: str) -> dict:
    if not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", repository):
        raise ValueError("Expected a GitHub repository in owner/name form.")
    if not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise ValueError("Expected a full, lowercase Git commit SHA.")
    errors = validate_repository(root)
    if errors:
        raise ValueError("\n".join(errors))

    def source_url(path: Path) -> str:
        relative = quote(path.relative_to(root).as_posix())
        return f"https://github.com/{repository}/blob/{revision}/{relative}"

    index = {"version": 1, "repository": repository, "revision": revision}
    for collection, pattern in RECORDS.items():
        records = []
        for path in sorted(root.glob(pattern)):
            record = json.loads(path.read_text(encoding="utf-8"))
            record["record_url"] = source_url(path)
            records.append(record)
        index[collection] = records

    pilots = []
    for path in sorted(root.glob("results/*/*/manifest.json")):
        data = path.read_bytes()
        pilot = json.loads(data)
        pilot.update(manifest_sha256=hashlib.sha256(data).hexdigest(),
                     manifest_url=source_url(path),
                     code_url=source_url(path.parent / pilot["execution"]["code_path"]),
                     coverage=None)
        for artifact in pilot["artifacts"]:
            artifact_path = path.parent / artifact["path"]
            artifact["url"] = source_url(artifact_path)
            if (pilot["task_id"] == "egfr-luad-evidence-coverage"
                    and artifact["path"] == "coverage.json"):
                coverage = json.loads(artifact_path.read_text(encoding="utf-8"))
                if isinstance(coverage, dict) and coverage.get("schema_version") == 1:
                    pilot["coverage"] = coverage
        pilots.append(pilot)
    index["pilots"] = pilots
    return index


def build(root: Path, output: Path, repository: str, revision: str) -> None:
    output = output.resolve()
    assets = root / "explorer"
    if output.is_relative_to(assets.resolve()):
        raise ValueError("Choose an output directory outside explorer/.")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError(f"Output must be a new or empty directory: {output}")
    if assets.is_symlink():
        raise ValueError("Explorer asset directory may not be a symlink.")
    for name in SITE_ASSETS:
        path = assets / name
        if path.is_symlink() or not path.is_file():
            raise ValueError(f"Expected regular public asset: explorer/{name}")
    index = make_index(root, repository, revision)

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as staging:
        site = Path(staging) / "site"
        site.mkdir()
        for name in SITE_ASSETS:
            shutil.copyfile(assets / name, site / name)
        api = site / "api" / "v1"
        api.mkdir(parents=True)
        exports = {"index": index}
        exports.update({name: index[name] for name in (*RECORDS, "pilots")})
        for name, payload in exports.items():
            content = json.dumps(payload, sort_keys=True, ensure_ascii=False,
                                 indent=2, allow_nan=False)
            (api / f"{name}.json").write_text(content + "\n", encoding="utf-8")
        site.replace(output)


def source_revision(root: Path, requested: str | None) -> str:
    if not (root / ".git").exists():
        if requested:
            return requested
        raise ValueError("A source archive needs --revision with its source commit SHA.")

    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    if requested and requested != head:
        raise ValueError("Requested revision differs from the checked-out commit.")
    paths = ["catalog", "tasks", "results", "reviews", "schemas", "explorer",
             "scripts/validate.py", "scripts/build_explorer.py"]
    status = subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=all",
                                      "--", *paths], cwd=root, text=True)
    if status:
        raise ValueError("Explorer inputs differ from HEAD. Commit the changes or build "
                         "from a clean checkout so provenance links match the published files.")
    return head


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "_site",
                        help="New or empty output directory (default: _site)")
    parser.add_argument("--repository", default=os.environ.get(
        "GITHUB_REPOSITORY", "dovvnloading/open-cancer"))
    parser.add_argument("--revision", help="Source commit SHA (default: checked-out HEAD)")
    args = parser.parse_args()
    try:
        revision = source_revision(ROOT, args.revision)
        build(ROOT, args.output, args.repository, revision)
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f"Explorer build failed: {exc}\n")
    print(f"Built Explorer at {args.output} from {args.repository}@{revision}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
