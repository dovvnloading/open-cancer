"""Deterministic/no-network release smoke tests, including failure modes."""
from datetime import datetime, timezone
import json
from pathlib import Path
import tempfile
import unittest

from scripts.smoke_explorer_release import (
    EXPECTED_DATASET_IDS, check_snapshot, local_fetch, remote_fetch,
)


REV = "a" * 40
MANIFEST = "b" * 64


def assets():
    rows = []
    for key in sorted(EXPECTED_DATASET_IDS):
        rows.append({
            "id": key,
            "research_context": [{"task_id": "egfr-luad-evidence-coverage",
                                  "gene_symbol": "EGFR"}],
            "record_url": "https://github.com/dovvnloading/open-cancer/blob/" +
                          REV + "/catalog/datasets/" + key + ".json",
        })
    pilot = {
        "task_id": "egfr-luad-evidence-coverage", "run_id": "pilot-20261110",
        "manifest_sha256": MANIFEST, "coverage": {"status": "unreviewed"},
        "publication_status": "legacy-featured",
        "clinical_validation_status": "not-assessed",
        "independent_reviewer_credentials_verified": False,
        "manifest_url": "https://github.com/dovvnloading/open-cancer/blob/" +
                        REV + "/results/egfr-luad-evidence-coverage/pilot-20261110/manifest.json",
    }
    review = {
        "task_id": pilot["task_id"], "run_id": pilot["run_id"],
        "reviewed_manifest_sha256": MANIFEST,
        "record_url": "https://github.com/dovvnloading/open-cancer/blob/" +
                      REV + "/reviews/egfr-luad-evidence-coverage/pilot-20261110/review.json",
    }
    index = {
        "version": 1, "repository": "dovvnloading/open-cancer", "revision": REV,
        "datasets": rows, "sources": [], "questions": [], "pilots": [pilot],
        "reviews": [review], "experimental_results": [],
    }
    all_files = {"api/v1/index.json": json.dumps(index).encode(),
                 "index.html": b"<!doctype html><title>Open Cancer Explorer</title>",
                 "app.js": b'const response = await fetch("./api/v1/index.json");'}
    for name in ("datasets", "sources", "questions", "pilots", "reviews", "experimental_results"):
        all_files["api/v1/" + name + ".json"] = json.dumps(index[name]).encode()
    all_files["api/v1/release.json"] = json.dumps({
        "version": 1, "repository": index["repository"],
        "revision": REV, "built_at": "2026-11-10T10:00:00Z",
        "workflow_run_id": "1234",
    }).encode()
    return all_files


def reader(files):
    def read(name):
        return files[name]
    return read


class ExplorerSmokeTests(unittest.TestCase):
    def test_good_local_snapshot_and_build_info(self):
        self.assertEqual([], check_snapshot(reader(assets()), REV,
                                            minimum_datasets=4, require_build_info=True))

    def test_stale_revision_and_partial_collection_fail(self):
        files = assets()
        errors = check_snapshot(reader(files), "b" * 40, minimum_datasets=4)
        self.assertTrue(any("revision" in err for err in errors))
        files["api/v1/reviews.json"] = b"[]"
        errors = check_snapshot(reader(files), REV, minimum_datasets=4)
        self.assertTrue(any("reviews endpoint" in err for err in errors))

    def test_cancer_crosswalk_and_version_bound_links_required(self):
        files = assets()
        index = json.loads(files["api/v1/index.json"])
        index["datasets"][0]["research_context"] = []
        index["pilots"][0]["manifest_url"] = "https://example.invalid/unpinned"
        files["api/v1/index.json"] = json.dumps(index).encode()
        files["api/v1/datasets.json"] = json.dumps(index["datasets"]).encode()
        files["api/v1/pilots.json"] = json.dumps(index["pilots"]).encode()
        errors = check_snapshot(reader(files), REV, minimum_datasets=4)
        self.assertTrue(any("crosswalk" in err for err in errors))
        self.assertTrue(any("not version bound" in err for err in errors))

    def test_bad_reviewer_link_and_clinical_claim_fail(self):
        files = assets()
        index = json.loads(files["api/v1/index.json"])
        index["reviews"][0]["reviewed_manifest_sha256"] = "c" * 64
        index["pilots"][0]["clinical_validation_status"] = "approved"
        files["api/v1/index.json"] = json.dumps(index).encode()
        files["api/v1/reviews.json"] = json.dumps(index["reviews"]).encode()
        files["api/v1/pilots.json"] = json.dumps(index["pilots"]).encode()
        errors = check_snapshot(reader(files), REV, minimum_datasets=4)
        self.assertTrue(any("Clinical" in err for err in errors))
        self.assertTrue(any("wrong pilot" in err for err in errors))

    def test_build_metadata_missing_stale_or_invalid_fails_closed(self):
        files = assets()
        del files["api/v1/release.json"]
        errors = check_snapshot(reader(files), REV, minimum_datasets=4, require_build_info=True)
        self.assertTrue(any("Missing or malformed" in err for err in errors))
        files = assets()
        meta = json.loads(files["api/v1/release.json"])
        meta["revision"] = "c" * 40
        meta["built_at"] = "3000-01-01T00:00:00Z"
        files["api/v1/release.json"] = json.dumps(meta).encode()
        errors = check_snapshot(reader(files), REV, minimum_datasets=4, require_build_info=True)
        self.assertTrue(any("Build metadata" in err for err in errors))
        self.assertTrue(any("timestamp" in err for err in errors))

    def test_fetcher_rejects_unsafe_urls_and_missing_files(self):
        for url in ("http://example.invalid/", "https://user:pw@example.invalid/",
                    "https://example.invalid/?q=secret"):
            with self.assertRaisesRegex(ValueError, "HTTPS site root"):
                remote_fetch(url)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "index.html").write_text("test")
            self.assertEqual(local_fetch(root)("index.html"), b"test")
            with self.assertRaisesRegex(ValueError, "missing local"):
                local_fetch(root)("../outside.txt")


if __name__ == "__main__":
    unittest.main()
