"""Offline tests for version-bound external review inputs and deterministic archives."""
from __future__ import annotations
from collections import Counter
import csv
import hashlib
import io
import json
from pathlib import Path
import shutil
import tempfile
import unittest
import zipfile

from scripts.build_review_packet import ROOT, RESULT, make_packet, review_sample, save_packet


class ReviewPacketTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.packet = make_packet(ROOT, "a" * 40)

    def test_exact_public_allowlist_and_manifest(self):
        packet = self.packet
        self.assertEqual(set(packet), {
            "manifest.json", "coverage.json", "civic_evidence_ids.csv",
            "cbio_protein_changes.csv", "analyze.py", "README.md",
            "SCIENTIFIC_CHECKLIST.md", "evidence_review_sample.csv",
            "packet-index.json", "SHA256SUMS",
        })
        index = json.loads(packet["packet-index.json"])
        self.assertEqual(index["kind"], "external-scientific-review-inputs-not-clinical-validation")
        self.assertEqual(index["original_manifest_sha256"], hashlib.sha256(
            (ROOT / RESULT / "manifest.json").read_bytes()).hexdigest())
        self.assertEqual(len(index["files"]), len(packet) - 2)
        self.assertEqual(json.loads(packet["manifest.json"])["findings"], [])
        self.assertEqual(json.loads(packet["coverage.json"])["status"], "unreviewed")
        for line in packet["SHA256SUMS"].decode().splitlines():
            checksum, filename = line.split("  ", 1)
            self.assertEqual(hashlib.sha256(packet[filename]).hexdigest(), checksum)

    def test_targeted_sample_is_twelve_unique_source_urls(self):
        sample = list(csv.DictReader(io.StringIO(
            self.packet["evidence_review_sample.csv"].decode())))
        self.assertEqual(len(sample), 12)
        self.assertEqual(len({row["record_id"] for row in sample}), 12)
        self.assertEqual(Counter(row["evidence_type"] for row in sample),
                         {"Predictive": 8, "Prognostic": 2, "Oncogenic": 2})
        self.assertTrue(all(row["record_url"].startswith("https://civicdb.org/links/evidence_items/")
                            for row in sample))

    def test_zip_determinism_and_file_allowlist(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            out1, out2 = base / "one", base / "two"
            zip1, zip2 = base / "one.zip", base / "two.zip"
            save_packet(out1, self.packet, zip1)
            save_packet(out2, self.packet, zip2)
            self.assertEqual(zip1.read_bytes(), zip2.read_bytes())
            with zipfile.ZipFile(zip1) as z:
                self.assertEqual(set(z.namelist()), set(self.packet))
                self.assertEqual(z.read("manifest.json"), self.packet["manifest.json"])
                self.assertTrue(all("/" not in name and ".." not in name
                                    for name in z.namelist()))
            with self.assertRaisesRegex(ValueError, "already exists"):
                save_packet(out1, self.packet, zip1)

    def test_corruption_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "repo"
            for dirname in ("schemas", "catalog", "tasks", "results",
                            "reviews", "examples", "review_packets"):
                path = ROOT / dirname
                if path.exists():
                    shutil.copytree(path, root / dirname)
            evidence = root / RESULT / "civic_evidence_ids.csv"
            evidence.write_text(evidence.read_text() + "unapproved-record\n")
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                make_packet(root, "a" * 40)
            evidence.unlink()
            evidence.symlink_to(ROOT / RESULT / "civic_evidence_ids.csv")
            with self.assertRaises(ValueError):
                make_packet(root, "a" * 40)

    def test_new_evidence_types_require_explicit_review(self):
        sample = [{"record_id": "42", "evidence_type": "Other",
                   "record_url": "https://civicdb.org/links/evidence_items/42/summary"}]
        with self.assertRaisesRegex(ValueError, "Evidence types changed"):
            review_sample(sample, "a" * 64)


if __name__ == "__main__":
    unittest.main()
