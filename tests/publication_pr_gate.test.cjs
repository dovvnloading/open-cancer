"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");
const {REGISTRY, touchesRegistry, evaluateEditorialReview} =
  require("../scripts/publication_pr_gate.cjs");

const SHA = "a".repeat(40);
const AUTHOR = "dovvnloading";
const files = [{filename: REGISTRY, status: "modified"}];
function review(fields = {}) {
  return {
    id: 10,
    state: "APPROVED",
    commit_id: SHA,
    submitted_at: "2026-11-11T08:00:00Z",
    author_association: "COLLABORATOR",
    user: {login: "second-maintainer", type: "User"},
    ...fields,
  };
}
function decision(overrides = {}) {
  return evaluateEditorialReview({
    files, reviews: [review()], author: AUTHOR, headSha: SHA, ...overrides,
  });
}

test("no editorial registry change skips approval check, including unrelated research", () => {
  const result = decision({files: [{filename: "results/new-run/manifest.json"}], reviews: []});
  assert.deepEqual(result, {required: false, ok: true,
                           reason: "No change to the editorial publication decision registry"});
  assert.equal(touchesRegistry([{filename: "other", previous_filename: REGISTRY}]), true);
  assert.equal(touchesRegistry([{filename: REGISTRY}]), true);
  assert.throws(() => touchesRegistry(null), /changed-file list/);
});

test("current independent human owner/member/collaborator approval passes", () => {
  assert.equal(decision().ok, true);
  for (const association of ["OWNER", "MEMBER", "COLLABORATOR"]) {
    assert.equal(decision({reviews: [review({author_association: association})]}).ok, true);
  }
});

test("self-approval or bot/untrusted review cannot authorize publication", () => {
  assert.equal(decision({reviews: [review({user: {login: AUTHOR, type: "User"}})]}).ok, false);
  assert.equal(decision({reviews: [review({user: {login: "bot", type: "Bot"}})]}).ok, false);
  assert.equal(decision({reviews: [review({author_association: "NONE"})]}).ok, false);
  assert.equal(decision({reviews: [review({author_association: "FIRST_TIME_CONTRIBUTOR"})]}).ok, false);
  assert.equal(decision({reviews: []}).ok, false);
});

test("approval must target the exact current PR head revision", () => {
  assert.equal(decision({reviews: [review({commit_id: "b".repeat(40)})]}).ok, false);
  assert.equal(decision({headSha: "branch-name"}).ok, false);
  assert.equal(decision({author: null}).ok, false);
  assert.equal(decision({reviews: null}).ok, false);
});

test("a later dismissed review or changes request supersedes approval", () => {
  const dismissed = review({id: 20, state: "DISMISSED",
                            submitted_at: "2026-11-11T09:00:00Z"});
  const rejected = review({id: 21, state: "CHANGES_REQUESTED",
                           submitted_at: "2026-11-11T09:00:00Z"});
  assert.equal(decision({reviews: [review(), dismissed]}).ok, false);
  assert.equal(decision({reviews: [rejected, review()]}).ok, false);
});

test("comment-only review does not erase valid fresh approval", () => {
  const comment = review({id: 24, state: "COMMENTED",
                          submitted_at: "2026-11-11T09:00:00Z"});
  assert.equal(decision({reviews: [review(), comment]}).ok, true);
});
