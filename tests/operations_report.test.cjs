"use strict";
const assert = require("node:assert/strict");
const test = require("node:test");
const {examine, render, latestRun, workflowState} = require("../scripts/operations_report.cjs");

const HEAD = "a".repeat(40);
const TIME = "2026-11-11T12:30:00Z";
const RUNS = [
  {id: 101, head_sha: HEAD, event: "push", status: "completed", conclusion: "success",
   created_at: TIME, name: "Validate catalog and research submissions"},
  {id: 102, head_sha: HEAD, event: "push", status: "completed", conclusion: "success",
   created_at: TIME, name: "Build and optionally deploy Open Cancer Explorer"},
  {id: 103, head_sha: HEAD, event: "push", status: "completed", conclusion: "failure",
   created_at: TIME, name: "Monitor public Explorer release freshness"},
  {id: 99, head_sha: "b".repeat(40), event: "push", status: "completed",
   conclusion: "success", created_at: "2026-11-11T13:00:00Z",
   name: "Validate catalog and research submissions"},
];

function fixture() {
  return {
    generated_at: TIME,
    revision: HEAD,
    branch_protected: false,
    workflow_runs: RUNS,
    deploy_jobs: [{name: "deploy", conclusion: "skipped"},
                  {name: "smoke-deployed", conclusion: "skipped"}],
    public_ledger: {checked_at: "2026-11-11T11:00:00Z",
                    passed: false, expected_revision: HEAD,
                    public_index_url: "https://example.org/api/v1/index.json",
                    output: "FAIL: stale"},
    catalog_ledger: {
      checked_at: "2026-11-11T10:00:00Z",
      alerts: [],
      sources: {
        gdc: {ok: true, last_checked: "2026-11-11T10:00:00Z",
              drift: {upstream_count: 33}},
        cbioportal: {ok: true, last_checked: "2026-11-11T10:00:00Z",
                     drift: {upstream_count: 32}},
        civic: {ok: true, last_checked: "2026-11-11T10:00:00Z",
                drift: {upstream_count: 3}},
      },
    },
    releases: [],
    issue_states: [
      {number: 14, state: "open", title: "Deployment"},
      {number: 18, state: "open", title: "Governance"},
      {number: 15, state: "closed", title: "Catalog monitoring"},
      {number: 11, state: "open", title: "Independent review"},
    ],
  };
}

test("passing build is never promoted to deployed or scientifically validated", () => {
  const summary = examine(fixture());
  assert.equal(summary.runs.validation.state, "passed");
  assert.equal(summary.runs.explorer.state, "passed");
  assert.equal(summary.deployment.deploy, "skipped");
  assert.equal(summary.deployment.live_smoke, "skipped");
  assert.equal(summary.public.current, false);
  assert.equal(summary.releases.review_kit_present, false);
  assert.equal(summary.catalog.healthy, true);
  assert.equal(summary.branch_protected, false);
  assert.deepEqual(summary.catalog.sources.map(x => x.count).sort((a, b) => a-b), [3,32,33]);
  const md = render(summary);
  assert.match(md, /Public Explorer API/);
  assert.match(md, /STALE \/ NOT VERIFIED/);
  assert.match(md, /Production deploy.*skipped/);
  assert.match(md, /not scientific or clinical validation/i);
  assert.match(md, /#14/);
  assert.doesNotMatch(md, /#15\]\(/);
});

test("older or unrelated successful runs cannot override current-head checks", () => {
  const r = fixture();
  r.workflow_runs = RUNS.filter(x => x.id !== 101);
  assert.equal(examine(r).runs.validation.state, "not verified at this revision");
  assert.equal(latestRun(RUNS, "Validate catalog and research submissions", HEAD).id, 101);
  assert.equal(workflowState({status: "in_progress"}), "pending");
  assert.equal(workflowState(null), "not verified at this revision");
});

test("green public checker only counts when current revision and checks are recent", () => {
  const r = fixture();
  r.public_ledger.passed = true;
  assert.equal(examine(r).public.current, true);
  r.public_ledger.expected_revision = "c".repeat(40);
  assert.equal(examine(r).public.current, false);
  r.public_ledger.expected_revision = HEAD;
  r.public_ledger.checked_at = "2026-11-01T00:00:00Z";
  assert.equal(examine(r).public.current, false);
});

test("stale or absent catalog ledger cannot be represented as healthy", () => {
  const r = fixture();
  r.catalog_ledger.checked_at = "2026-10-01T00:00:00Z";
  assert.equal(examine(r).catalog.healthy, false);
  r.catalog_ledger = null;
  assert.equal(examine(r).catalog.healthy, false);
  assert.match(render(examine(r)), /No durable source observations/);
});

test("release existence never certifies immutable inputs or scientific review", () => {
  const r = fixture();
  r.releases = [{tag_name: "reviewkit-egfr-luad-v1.0.0", draft: false, prerelease: false}];
  const summary = examine(r);
  assert.equal(summary.releases.review_kit_present, true);
  assert.ok(summary.blockers.some(x => /Independent scientific/.test(x)));
  assert.match(render(summary), /digest still requires independent verification/);
});

test("invalid revision and time fail fast", () => {
  assert.throws(() => examine({...fixture(), revision: "../../main"}), /full lowercase Git SHA/);
  assert.throws(() => examine({...fixture(), generated_at: "not a time"}), /ISO reporting time/);
});
