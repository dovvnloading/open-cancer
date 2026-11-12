"""Offline regression checks for deployed Explorer release health."""
import io
import json
import unittest
from urllib.error import HTTPError

from scripts.release_polling import wait_for_healthy_release
from scripts.check_explorer_freshness import check_remote, validate_index


REV = "a" * 40
REPO = "dovvnloading/open-cancer"


def index():
    return {
        "version": 1, "repository": REPO, "revision": REV,
        "datasets": [{"id": "one"}, {"id": "two"}], "sources": [],
        "questions": [], "pilots": [], "reviews": [], "experimental_results": [],
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
        data = index()
        data.pop("experimental_results")
        self.assertTrue(any("experimental_results" in error for error in
                            validate_index(data, revision=REV, repository=REPO)))

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

    def test_short_lived_503_does_not_trigger_a_false_freshness_alarm(self):
        attempted = [0]
        clock = [0.0]

        def flaky(request, timeout):
            attempted[0] += 1
            if attempted[0] == 1:
                raise HTTPError(request.full_url, 503, "Service Unavailable", {}, None)
            return io.BytesIO(json.dumps(index()).encode())

        def advance(delay):
            clock[0] += delay

        errors = wait_for_healthy_release(
            lambda: check_remote("https://example.invalid/api/v1/index.json",
                                 revision=REV, repository=REPO, minimum_datasets=2,
                                 opener=flaky),
            wait_seconds=3, interval_seconds=1,
            clock=lambda: clock[0], sleep=advance)
        self.assertEqual([], errors)
        self.assertEqual(attempted[0], 2)


if __name__ == "__main__":
    unittest.main()
