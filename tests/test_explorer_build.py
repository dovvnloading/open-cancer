"""Offline checks of Explorer JSON contract, version links, and publication boundary."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from scripts.build_explorer import ROOT, SITE_ASSETS, build, make_index, source_revision


REPOSITORY = "example/open-cancer"
REVISION = "a" * 40


class ExplorerBuildTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.base = Path(directory.name)
        self.root = self.base / "repo"
        for name in ("schemas", "catalog", "tasks", "results", "examples"):
            shutil.copytree(ROOT / name, self.root / name)
        (self.root / "explorer").mkdir()
        for name in SITE_ASSETS:
            (self.root / "explorer" / name).write_text("Fixture: " + name)

    def test_release_links_require_committed_inputs(self):
        subprocess.run(["git", "init", "-q"], cwd=self.root, check=True)
        subprocess.run(["git", "add", "."], cwd=self.root, check=True)
        subprocess.run(["git", "-c", "user.name=Test", "-c", "user.email=test@example.invalid",
                        "commit", "-qm", "Fixture"], cwd=self.root, check=True)
        revision = source_revision(self.root, None)
        self.assertEqual(revision, source_revision(self.root, revision))
        with self.assertRaisesRegex(ValueError, "differs from the checked-out commit"):
            source_revision(self.root, REVISION)
        (self.root / "catalog/datasets/uncommitted.json").write_text("{}")
        with self.assertRaisesRegex(ValueError, "inputs differ from HEAD"):
            source_revision(self.root, None)

    def test_only_regular_public_assets_are_copied(self):
        assets = self.root / "explorer"
        (assets / ".gitignore").write_text("local-only.txt\n")
        (assets / "local-only.txt").write_text("development fixture")
        (assets / "stray-link").symlink_to(assets / "local-only.txt")
        output = self.base / "site"
        build(self.root, output, REPOSITORY, REVISION)
        self.assertEqual({"index.html", "app.js", "style.css", "catalog.mjs", "api"},
                         {path.name for path in output.iterdir()})
        (assets / "app.js").unlink()
        (assets / "app.js").symlink_to(assets / "local-only.txt")
        with self.assertRaisesRegex(ValueError, "Expected regular public asset: explorer/app.js"):
            build(self.root, self.base / "site2", REPOSITORY, REVISION)
        self.assertFalse((self.base / "site2").exists())

    def test_site_api_roundtrip_at_same_revision(self):
        output = self.base / "site"
        build(self.root, output, REPOSITORY, REVISION)
        api = output / "api/v1"
        self.assertEqual((self.root / "explorer/index.html").read_bytes(),
                         (output / "index.html").read_bytes())
        index = json.loads((api / "index.json").read_bytes())
        self.assertEqual(1, index["version"])
        self.assertEqual(REVISION, index["revision"])
        for name in ("datasets", "sources", "questions", "reviews", "pilots"):
            self.assertEqual(index[name], json.loads((api / f"{name}.json").read_bytes()))
        self.assertEqual([], index["reviews"])  # synthetic schema examples stay outside the API
        native_paths = sorted((self.root / "catalog/datasets").glob("*.json"))
        self.assertEqual([p.stem for p in native_paths], [r["id"] for r in index["datasets"]])
        for path, row in zip(native_paths, index["datasets"]):
            self.assertEqual(f"https://github.com/{REPOSITORY}/blob/{REVISION}/catalog/datasets/{path.name}",
                             row.pop("record_url"))
            self.assertEqual(json.loads(path.read_bytes()), row)
        other = self.base / "site2"
        build(self.root, other, REPOSITORY, REVISION)
        for path in api.glob("*.json"):
            self.assertEqual(path.read_bytes(), (other / "api/v1" / path.name).read_bytes())

    def test_pilot_and_review_bind_to_exact_manifest(self):
        task = "egfr-luad-evidence-coverage"
        run = "pilot-20261110"
        result_dir = self.root / "results" / task / run
        digest = hashlib.sha256((result_dir / "manifest.json").read_bytes()).hexdigest()
        review_dir = self.root / "reviews" / task / run
        review_dir.mkdir(parents=True)
        review = json.loads((self.root / "examples/reviews/synthetic-methods-example.json").read_bytes())
        review.update(id="test-fixture", task_id=task, run_id=run,
                      synthetic=False, reviewed_manifest_sha256=digest)
        (review_dir / "test-fixture.json").write_text(json.dumps(review))

        index = make_index(self.root, REPOSITORY, REVISION)
        pilot = next(p for p in index["pilots"] if (p["task_id"], p["run_id"]) == (task, run))
        self.assertEqual(index["reviews"][0]["reviewed_manifest_sha256"], pilot["manifest_sha256"])
        self.assertEqual(digest, pilot["manifest_sha256"])
        self.assertEqual(json.loads((result_dir / "coverage.json").read_bytes()), pilot["coverage"])
        self.assertTrue(pilot["code_url"].endswith(f"/{run}/analyze.py"))
        for artifact in pilot["artifacts"]:
            self.assertTrue(artifact["url"].endswith(f"/{run}/{artifact['path']}"))

        # Reformatting the manifest changes the bytes that were reviewed.
        with (result_dir / "manifest.json").open("a") as stream:
            stream.write("\n")
        with self.assertRaisesRegex(ValueError, "stale or mismatched"):
            build(self.root, self.base / "site", REPOSITORY, REVISION)
        self.assertFalse((self.base / "site").exists())

    def test_only_declared_coverage_and_matching_checksums_can_publish(self):
        result_dir = self.root / "results/egfr-luad-evidence-coverage/pilot-20261110"
        path = result_dir / "coverage.json"
        path.write_text("{}")
        with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
            make_index(self.root, REPOSITORY, REVISION)
        manifest_path = result_dir / "manifest.json"
        manifest = json.loads(manifest_path.read_bytes())
        manifest["artifacts"] = [a for a in manifest["artifacts"] if a["path"] != "coverage.json"]
        manifest_path.write_text(json.dumps(manifest))
        index = make_index(self.root, REPOSITORY, REVISION)
        pilot = next(p for p in index["pilots"] if p["task_id"] == manifest["task_id"]
                     and p["run_id"] == manifest["run_id"])
        self.assertIsNone(pilot["coverage"])


if __name__ == "__main__":
    unittest.main()
