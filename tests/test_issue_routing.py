"""Keep GitHub Issues strictly code/engineering scoped; no research recruitment tickets."""
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / '.github' / 'ISSUE_TEMPLATE'


class TechnicalIssueRoutingTests(unittest.TestCase):
    def test_only_engineering_issue_template_available(self):
        config = (TEMPLATE / 'config.yml').read_text(encoding='utf-8')
        self.assertIn('blank_issues_enabled: false', config)
        self.assertIn('docs/compute.md', config)
        self.assertIn('/actions', config)
        self.assertEqual(['engineering.yml'], sorted(p.name for p in TEMPLATE.glob('*.yml')
                                                      if p.name != 'config.yml'))
        self.assertNotIn('independent-reviewers', config)

    def test_engineering_issue_requires_reproduction_and_acceptance(self):
        engineering = (TEMPLATE / 'engineering.yml').read_text(encoding='utf-8')
        for part in ('id: component', 'id: problem', 'id: reproduction', 'id: acceptance'):
            self.assertIn(part, engineering)
        self.assertEqual(engineering.count('required: true'), 4)

    def test_monitoring_workflows_cannot_create_status_issues(self):
        for name in ('catalog-monitor.yml', 'explorer-freshness.yml', 'operations-report.yml'):
            with self.subTest(workflow=name):
                workflow = (ROOT / '.github/workflows' / name).read_text()
                self.assertNotIn('issues: write', workflow)
                for forbidden in ('github.rest.issues.create(', 'github.rest.issues.update(',
                                  'github.rest.issues.createComment(', 'issues: reopen'):
                    self.assertNotIn(forbidden, workflow)
                self.assertIn('actions/upload-artifact@v4', workflow)

    def test_routine_reports_use_actions_artifacts_not_issue_ledger(self):
        monitor = (ROOT / '.github/workflows/catalog-monitor.yml').read_text()
        public = (ROOT / '.github/workflows/explorer-freshness.yml').read_text()
        operations = (ROOT / '.github/workflows/operations-report.yml').read_text()
        self.assertIn('actions/download-artifact@v4', monitor)
        self.assertIn('catalog-health-report', monitor)
        self.assertIn('explorer-health-report', public)
        self.assertIn('actions/download-artifact@v4', operations)
        self.assertNotIn('Catalog metadata health ledger (automated)', monitor)
        self.assertNotIn('Public Explorer release health (automated)', public)
        self.assertNotIn('Open Cancer weekly operations report (automated)', operations)


if __name__ == '__main__':
    unittest.main()
