"""Offline regression tests for scheduled public-metadata monitoring."""
import unittest
from urllib.error import HTTPError, URLError

from scripts.monitor_catalog import SOURCES, compare, managed, observation, retry, scan


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
        # Baseline comes from earlier connector observations, not manually
        # curated descriptions, which are expected to differ.
        previous = {row["id"]: observation(row) for row in old}
        diff = compare("gdc", old, new, previous_observations=previous)
        self.assertEqual(diff["added"], ["gdc-tcga-coad-project"])
        self.assertEqual(diff["removed"], ["gdc-tcga-brca-project"])
        self.assertEqual(diff["changed_fields"][KEYS["gdc"]], ["release"])
        self.assertTrue(diff["baseline_established"])
        self.assertFalse(diff["unexpected_shrinkage"])

    def test_initial_baseline_does_not_treat_curated_text_as_upstream_drift(self):
        curated = [sample("gdc")]
        curated[0]["title"] = "Editorial name intentionally differs"
        incoming = [sample("gdc")]
        incoming[0]["title"] = "Official upstream project name"
        first = compare("gdc", curated, incoming)
        self.assertEqual(first["added"], [])
        self.assertEqual(first["changed_fields"], {})
        self.assertFalse(first["baseline_established"])
        second = compare("gdc", curated, incoming,
                         previous_observations=first["observations"])
        self.assertEqual(second["changed_fields"], {})
        self.assertTrue(second["baseline_established"])
        altered = [dict(incoming[0], title="New upstream project name")]
        third = compare("gdc", curated, altered,
                        previous_observations=first["observations"])
        self.assertEqual(third["changed_fields"][KEYS["gdc"]], ["title"])

    def test_outage_preserves_previous_observations_after_recovery(self):
        current = {source: [sample(source)] for source in SOURCES}
        first = scan(current=current, discover=lambda source: [sample(source)],
                     checked_at="2026-11-10T00:00:00+00:00", sleep=lambda _: None)
        def discover(source):
            if source == "gdc":
                raise URLError("temporary outage")
            return [sample(source)]
        second = scan(current=current, discover=discover, previous_report=first,
                      checked_at="2026-11-11T00:00:00+00:00", sleep=lambda _: None)
        self.assertFalse(second["sources"]["gdc"]["ok"])
        # The scanner now preserves last-good observations in the Actions
        # report itself, without relying on a status-only Issue ledger.
        self.assertEqual(second["sources"]["gdc"]["observations"],
                         first["sources"]["gdc"]["observations"])
        self.assertEqual(second["sources"]["gdc"]["last_healthy"],
                         "2026-11-10T00:00:00+00:00")
        third = scan(current=current,
                     discover=lambda source: [sample(source, release="v2" if source == "gdc" else "v1")],
                     previous_report=second,
                     checked_at="2026-11-12T00:00:00+00:00", sleep=lambda _: None)
        self.assertEqual(third["sources"]["gdc"]["drift"]["changed_fields"][KEYS["gdc"]],
                         ["release"])

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

    def test_terms_review_warning_persists_until_human_review_date_advances(self):
        current = {source: [sample(source)] for source in SOURCES}
        terms = {"gdc": "https://example.org/terms"}
        reviewed = {"gdc": "2026-11-10"}

        def run(headers, prior=None, review_stamp=reviewed):
            return scan(
                current=current, discover=lambda source: [sample(source)],
                terms_urls=terms, probe=lambda url: {"url": url, **headers},
                previous_report=prior, rights_reviewed_on=review_stamp,
                checked_at="2026-11-12T00:00:00+00:00", sleep=lambda _: None,
            )

        first = run({"etag": "v1", "last_modified": "day-1"})
        self.assertEqual(first["sources"]["gdc"]["pending_terms_review"], [])
        changed = run({"etag": "v2", "last_modified": "day-2"}, prior=first)
        self.assertEqual(changed["sources"]["gdc"]["pending_terms_review"],
                         ["etag", "last_modified"])
        unchanged = run({"etag": "v2", "last_modified": "day-2"}, prior=changed)
        self.assertEqual(unchanged["sources"]["gdc"]["pending_terms_review"],
                         ["etag", "last_modified"])
        reviewed_now = run({"etag": "v2", "last_modified": "day-2"}, prior=unchanged,
                           review_stamp={"gdc": "2026-11-13"})
        self.assertEqual(reviewed_now["sources"]["gdc"]["pending_terms_review"], [])

    def test_previous_failed_run_terms_alerts_migrate_to_advisory_review(self):
        current = {source: [sample(source)] for source in SOURCES}
        prior = scan(
            current=current, discover=lambda source: [sample(source)],
            terms_urls={"gdc": "https://example.org/terms"},
            probe=lambda url: {"url": url, "etag": "v2"},
            rights_reviewed_on={"gdc": "2026-11-10"},
            checked_at="2026-11-11T00:00:00+00:00",
        )
        prior["alerts"] = ["gdc: source-terms etag changed; review reuse rights"]
        followup = scan(
            current=current, discover=lambda source: [sample(source)],
            terms_urls={"gdc": "https://example.org/terms"},
            probe=lambda url: {"url": url, "etag": "v2"},
            rights_reviewed_on={"gdc": "2026-11-10"},
            previous_report=prior, checked_at="2026-11-12T00:00:00+00:00",
        )
        self.assertEqual(followup["sources"]["gdc"]["pending_terms_review"], ["etag"])

    def test_changed_terms_url_is_reviewed_even_without_http_validators(self):
        current = {source: [sample(source)] for source in SOURCES}
        def run(url, previous=None):
            return scan(
                current=current, discover=lambda source: [sample(source)],
                terms_urls={"gdc": url},
                probe=lambda url: {"url": url, "http_status": 200},
                previous_report=previous,
                rights_reviewed_on={"gdc": "2026-11-10"},
                checked_at="2026-11-12T00:00:00+00:00",
            )
        first = run("https://example.org/terms-v1")
        second = run("https://example.org/terms-v2", first)
        self.assertEqual(second["sources"]["gdc"]["pending_terms_review"], ["url"])

    def test_link_outage_and_drift_both_reported(self):
        current = {name: [sample(name)] for name in SOURCES}
        report = scan(current=current,
                      discover=lambda name: [sample(name, release="v2")],
                      terms_urls={"civic": "https://civicdb.org/terms"},
                      probe=lambda _: (_ for _ in ()).throw(URLError("terms offline")),
                      checked_at="2026-11-10T00:00:00+00:00", sleep=lambda _: None)
        self.assertTrue(report["sources"]["civic"]["ok"])
        self.assertIn("terms offline", report["sources"]["civic"]["terms_link_error"])
        self.assertFalse(report["sources"]["civic"]["drift"]["baseline_established"])
        later = scan(current=current, discover=lambda name: [sample(name, release="v3")],
                     previous_report=report,
                     checked_at="2026-11-11T00:00:00+00:00", sleep=lambda _: None)
        self.assertEqual(later["sources"]["civic"]["drift"]["changed_fields"][KEYS["civic"]],
                         ["release"])


if __name__ == "__main__":
    unittest.main()
