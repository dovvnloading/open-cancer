"""Offline regression checks for deployed Explorer release health."""
import io
import json
import unittest

from scripts.check_explorer_freshness import check_remote, validate_index


REV = "a" * 40
REPO = "dovvnloading/open-cancer"


def index():
    return {
        "version": 1, "repository": REPO, "revision": REV,
        "datasets": [{"id": "one"}, {"id": "two"}], "sources": [],
        "questions": [], "pilots": [], "reviews": [],
    }


class FreshnessTests(unittest.TestCase):
    def test_match(self):
        self.assertEqual([], validate_index(index(), revision=REV, repository=REPO,
                                            minimum_datasets=2))

    def test_fails_on_stale_or_mismatched_repo(self):
        errors = validate_index(index(), revision="b" * 40, repository="someone/other")
        self.assertTrue(any("STALE" in error for error in errors))
        self.assertTrue(any("repository" in error for error in errors))

    def test_fails_closed_on_missing_datasets_or_schema(self):
        data = index()
        data["datasets"] = None
        self.assertTrue(any("datasets" in error for error in
                            validate_index(data, revision=REV, repository=REPO)))
        data = index()
        data["datasets"] = []
        self.assertTrue(any("fewer" in error for error in
                            validate_index(data, revision=REV, repository=REPO,
                                           minimum_datasets=1)))
        self.assertTrue(validate_index([], revision=REV, repository=REPO))

    def test_https_and_bounded_response(self):
        request_url = "https://open-cancer-explorer.lagging.chatgpt.site/api/v1/index.json"
        self.assertTrue(check_remote("http://example.com/index.json", revision=REV,
                                     repository=REPO))
        self.assertTrue(check_remote("https://user:pass@example.com/index.json",
                                     revision=REV, repository=REPO))
        def good(request, timeout):
            self.assertEqual(timeout, 20)
            self.assertEqual(request.get_header("Accept"), "application/json")
            return io.BytesIO(json.dumps(index()).encode())
        self.assertEqual([], check_remote(request_url, revision=REV, repository=REPO,
                                          minimum_datasets=2, opener=good))
        self.assertTrue(check_remote(request_url, revision=REV, repository=REPO,
                                     opener=lambda request, timeout: io.BytesIO(b"<html>wrong</html>")))
        self.assertTrue(check_remote(request_url, revision=REV, repository=REPO,
                                     opener=lambda request, timeout: io.BytesIO(b"X" * (4*1024*1024 + 1))))


if __name__ == "__main__":
    unittest.main()
