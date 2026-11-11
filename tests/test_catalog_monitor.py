"""Offline regression tests for scheduled public-metadata monitoring."""
import unittest
from urllib.error import HTTPError, URLError

from scripts.monitor_catalog import SOURCES, compare, managed, retry, scan


KEYS = {"gdc": "gdc-tcga-luad-project",
        "cbioportal": "cbioportal-luad-tcga-pancan-atlas-2018",
        "civic": "civic-nightly-clinical-evidence"}


def sample(source, key=None, release="v1", rights="verify-upstream"):
    return {"id": key or KEYS[source], "source_id": source,
            "cohort_definition": "test cohort",
            "release": {"label": release, "mutable": True}, "reuse": {"status": rights}}


class CatalogMonitorTests(unittest.TestCase):
    def test_detects_new_missing_and_changed_metadata(self):
        old = [sample("gdc"), sample("gdc", "gdc-tcga-brca-project")]
        new = [sample("gdc", release="v2", rights="restricted"),
               sample("gdc", "gdc-tcga-coad-project")]
        new[0]["cohort_definition"] = "different cohort"
        diff = compare("gdc", old, new)
        self.assertEqual(diff["added"], ["gdc-tcga-coad-project"])
        self.assertEqual(diff["removed"], ["gdc-tcga-brca-project"])
        self.assertEqual(diff["changed_fields"][KEYS["gdc"]],
                         ["cohort_definition", "release", "reuse"])
        self.assertFalse(diff["unexpected_shrinkage"])

    def test_empty_duplicates_and_wrong_identifiers_fail_closed(self):
        with self.assertRaisesRegex(ValueError, "zero upstream"):
            compare("gdc", [sample("gdc")], [])
        with self.assertRaisesRegex(ValueError, "duplicate"):
            compare("gdc", [], [sample("gdc"), sample("gdc")])
        with self.assertRaisesRegex(ValueError, "unexpected"):
            compare("gdc", [], [sample("civic", KEYS["gdc"])])

    def test_shrinkage_and_unrelated_pilot_exclusion(self):
        self.assertTrue(managed("gdc", KEYS["gdc"]))
        self.assertFalse(managed("gdc", "gdc-tcga-luad-egfr-ssm-occurrences"))
        difference = compare("gdc",
                             [sample("gdc"), sample("gdc", "gdc-tcga-brca-project")],
                             [sample("gdc")])
        self.assertTrue(difference["unexpected_shrinkage"])

    def test_source_outage_preserves_healthy_sources(self):
        current = {source: [sample(source)] for source in SOURCES}
        def discover(source):
            if source == "cbioportal":
                raise URLError("outage")
            return [sample(source)]
        report = scan(current=current, discover=discover, terms_urls={},
                      checked_at="2026-11-10T00:00:00+00:00", sleep=lambda _: None)
        self.assertTrue(report["sources"]["gdc"]["ok"])
        self.assertFalse(report["sources"]["cbioportal"]["ok"])
        self.assertIn("outage", report["sources"]["cbioportal"]["error"])
        self.assertTrue(report["sources"]["civic"]["ok"])
        self.assertEqual(report["sources"]["civic"]["last_checked"], report["checked_at"])

    def test_retry_transient_http_not_permanent_failures(self):
        attempts, sleeps = [], []
        def flaky():
            attempts.append(1)
            if len(attempts) < 3:
                raise HTTPError("https://example.org", 503, "busy", {}, None)
            return "ok"
        self.assertEqual("ok", retry(flaky, sleep=sleeps.append))
        self.assertEqual(len(attempts), 3)
        self.assertEqual(sleeps, [1, 2])
        with self.assertRaises(HTTPError):
            retry(lambda: (_ for _ in ()).throw(
                HTTPError("https://example.org", 404, "gone", {}, None)),
                sleep=lambda _: self.fail("404 must not retry"))

    def test_link_outage_and_drift_both_reported(self):
        current = {name: [sample(name)] for name in SOURCES}
        report = scan(current=current,
                      discover=lambda name: [sample(name, release="v2")],
                      terms_urls={"civic": "https://civicdb.org/terms"},
                      probe=lambda _: (_ for _ in ()).throw(URLError("terms offline")),
                      checked_at="2026-11-10T00:00:00+00:00", sleep=lambda _: None)
        self.assertTrue(report["sources"]["civic"]["ok"])
        self.assertIn("terms offline", report["sources"]["civic"]["terms_link_error"])
        self.assertEqual(report["sources"]["civic"]["drift"]["changed_fields"][KEYS["civic"]],
                         ["release"])


if __name__ == "__main__":
    unittest.main()
