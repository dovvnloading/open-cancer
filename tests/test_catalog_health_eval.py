"""Offline regression coverage for artifact-only catalog status monitoring."""
from __future__ import annotations

import unittest

from scripts.evaluate_catalog_health import evaluate, render
from scripts.monitor_catalog import scan


def healthy_report():
    return {
        "schema_version": 1,
        "checked_at": "2026-11-12T15:00:00Z",
        "sources": {
            key: {
                "ok": True,
                "last_checked": "2026-11-12T15:00:00Z",
                "last_healthy": "2026-11-12T15:00:00Z",
                "drift": {
                    "upstream_count": 2, "added": [], "removed": [],
                    "changed_fields": {}, "unexpected_shrinkage": False,
                },
            }
            for key in ("cbioportal", "civic", "gdc")
        },
    }


class ArtifactBackedCatalogHealthTests(unittest.TestCase):
    def test_clean_report_explicitly_marks_missing_comparison_baseline(self):
        report = evaluate(healthy_report())
        self.assertEqual(report["alerts"], [])
        self.assertFalse(report["previous_monitor_artifact_used"])
        self.assertIn("historical metadata drift is not assessed".lower(),
                      render(report).lower().replace("past drift not assessed",
                                                      "historical metadata drift is not assessed"))

    def test_drift_and_outages_fail_but_terms_headers_only_warn(self):
        report = healthy_report()
        report["sources"]["gdc"]["ok"] = False
        report["sources"]["gdc"]["error"] = "HTTP 503"
        report["sources"]["civic"]["drift"]["removed"] = ["civic-record"]
        report["sources"]["cbioportal"]["terms_link"] = {"etag": "new"}
        previous = healthy_report()
        previous["sources"]["cbioportal"]["last_verified_terms_link"] = {"etag": "old"}
        result = evaluate(report, previous)
        self.assertTrue(result["previous_monitor_artifact_used"])
        self.assertTrue(any("gdc: endpoint failure" in item for item in result["alerts"]))
        self.assertTrue(any("civic: metadata drift" in item for item in result["alerts"]))
        self.assertTrue(any("source-terms etag changed" in item for item in result["warnings"]))
        self.assertFalse(any("source-terms" in item for item in result["alerts"]))
        self.assertIn("ACTIONABLE DRIFT / FAILURE", render(result))

    def test_header_changes_alone_warn_without_failing_ci(self):
        previous = healthy_report()
        previous["sources"]["gdc"]["terms_link"] = {
            "url": "https://example.org/terms", "etag": "old", "last_modified": "day-1"
        }
        report = healthy_report()
        report["sources"]["gdc"]["terms_link"] = {
            "url": "https://example.org/terms", "etag": "new", "last_modified": "day-2"
        }
        result = evaluate(report, previous)
        self.assertEqual(result["alerts"], [])
        self.assertEqual(len(result["warnings"]), 2)
        self.assertTrue(all("gdc: source-terms" in w for w in result["warnings"]))
        self.assertIn("Advisory source-terms metadata reviews", render(result))

    def test_pending_review_does_not_hide_actionable_terms_link_outage(self):
        report = healthy_report()
        report["sources"]["cbioportal"]["pending_terms_review"] = ["etag"]
        report["sources"]["cbioportal"]["terms_link_error"] = "HTTP 404"
        result = evaluate(report)
        self.assertTrue(any("source-rights link issue" in a for a in result["alerts"]))
        self.assertEqual(len(result["warnings"]), 1)

    def test_invalid_pending_review_is_not_accepted(self):
        report = healthy_report()
        report["sources"]["gdc"]["pending_terms_review"] = ["unrecognized"]
        with self.assertRaisesRegex(ValueError, "Invalid terms review"):
            evaluate(report)

    def test_failed_source_keeps_last_non_sensitive_observation_for_recovery(self):
        prior = {
            "sources": {
                "gdc": {
                    "ok": True, "last_checked": "2026-11-11T00:00:00Z",
                    "observations": {"gdc-tcga-luad-project": {"accession": "TCGA-LUAD"}},
                    "terms_link": {"url": "https://example.org/terms", "etag": "v1"},
                }
            }
        }

        def down(_source):
            raise ValueError("Synthetic network failure")

        report = scan(
            current={"gdc": [], "cbioportal": [], "civic": []},
            discover=down,
            checked_at="2026-11-12T00:00:00Z",
            previous_report=prior,
        )
        gdc = report["sources"]["gdc"]
        self.assertFalse(gdc["ok"])
        self.assertEqual(gdc["last_healthy"], "2026-11-11T00:00:00Z")
        self.assertEqual(gdc["last_verified_terms_link"]["etag"], "v1")
        self.assertEqual(gdc["observations"], prior["sources"]["gdc"]["observations"])
        self.assertTrue(report["previous_monitor_artifact_used"])

    def test_invalid_source_record_and_schema_fail_closed(self):
        with self.assertRaisesRegex(ValueError, "version-1"):
            evaluate({"schema_version": 0, "sources": {}})
        report = healthy_report()
        report["sources"].pop("civic")
        with self.assertRaisesRegex(ValueError, "unexpected or missing"):
            evaluate(report)


if __name__ == "__main__":
    unittest.main()
