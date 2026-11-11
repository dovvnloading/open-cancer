"""Catch accidental loss of the base-pinned publication review workflow.

These checks are not a substitute for enforced branch protection or GitHub
reviewer credentials. They prevent shipping an accidental CI omission.
"""
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]


class EditorialGateWorkflowTests(unittest.TestCase):
    def test_editorial_review_job_uses_trusted_base_evaluator(self):
        content = (ROOT / ".github/workflows/validate.yml").read_text(encoding="utf-8")
        self.assertIn("  publication-review:\n", content)
        self.assertIn("if: github.event_name == 'pull_request'", content)
        self.assertIn("pull-requests: read", content)
        self.assertIn("ref: ${{ github.event.pull_request.base.sha }}", content)
        self.assertIn("persist-credentials: false", content)
        self.assertIn("github.rest.pulls.listFiles", content)
        self.assertIn("github.rest.pulls.listReviews", content)
        self.assertIn("evaluateEditorialReview({", content)
        self.assertIn("author: request.user.login", content)
        self.assertIn("headSha: request.head.sha", content)
        self.assertIn("core.setFailed(verdict.reason", content)

    def test_review_policy_codeowners_routes_to_maintainer(self):
        owners = (ROOT / ".github/CODEOWNERS").read_text(encoding="utf-8")
        self.assertIn("/catalog/publication-decisions.json @dovvnloading", owners)
        self.assertIn("/scripts/publication_pr_gate.cjs @dovvnloading", owners)
        self.assertIn("/.github/workflows/validate.yml @dovvnloading", owners)


if __name__ == "__main__":
    unittest.main()
