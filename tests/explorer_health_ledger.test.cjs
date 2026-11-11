"use strict";
const assert = require("node:assert/strict");
const test = require("node:test");
const {parsePrevious, renderHealth, shouldRecordTransition} =
  require("../scripts/explorer_health_ledger.cjs");

const status = {
  checked_at: "2026-11-11T10:20:00Z",
  expected_revision: "a".repeat(40),
  public_index_url: "https://example.org/api/v1/index.json",
  passed: false,
  output: "FAIL: STALE Explorer: deployed revision 'old', expected 'new'",
};

test("renders unambiguous production failure with Git revision and structured state", () => {
  const body = renderHealth(status);
  assert.match(body, /STALE OR UNAVAILABLE/);
  assert.match(body, /not scientific\/clinical validation/);
  assert.deepEqual(parsePrevious(body), {...status});
  assert.equal(shouldRecordTransition(null, status), true);
  assert.equal(shouldRecordTransition(status, status), false);
});

test("passing status and changing target produce transitions", () => {
  const good = {...status, passed: true, output: "OK: current revision"};
  assert.match(renderHealth(good), /CURRENT/);
  assert.equal(shouldRecordTransition(status, good), true);
  assert.equal(shouldRecordTransition(good, {...good, expected_revision: "b".repeat(40)}), true);
  assert.equal(shouldRecordTransition(good, {...good, public_index_url: "https://new.test/index.json"}), true);
});

test("malformed history and injection-like output fail safely", () => {
  assert.equal(parsePrevious("ordinary issue text"), null);
  assert.throws(() => parsePrevious("<!-- explorer-health-json\n{\"passed\":1}\n-->"), /Invalid/);
  assert.throws(() => renderHealth({...status, expected_revision: "../main"}), /Git SHA/);
  assert.throws(() => renderHealth({...status, public_index_url: "http://host"}), /HTTPS/);
  const body = renderHealth({...status, output: "~~~text\nfaked block"});
  assert.doesNotMatch(body, /~~~text\nfaked block/);
});
