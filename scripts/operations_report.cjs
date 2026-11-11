"use strict";

/* Pure, evidence-conservative status reporter. No networking or writes here. */
const REPO = "dovvnloading/open-cancer";
const RUN_NAMES = Object.freeze({
  validation: "Validate catalog and research submissions",
  explorer: "Build and optionally deploy Open Cancer Explorer",
  publicSite: "Monitor public Explorer release freshness",
  catalog: "Monitor catalog freshness and upstream source health",
});
const ACTIONS = "https://github.com/" + REPO + "/actions/runs/";
const ISSUE = "https://github.com/" + REPO + "/issues/";

function latestRun(runs, name, sha) {
  const matches = (Array.isArray(runs) ? runs : []).filter(run =>
    run && run.name === name && run.head_sha === sha && run.event === "push");
  return matches.sort((a, b) =>
    String(b.created_at || "").localeCompare(String(a.created_at || "")))[0] || null;
}

function workflowState(run) {
  if (!run) return "not verified at this revision";
  if (run.status !== "completed") return "pending";
  if (run.conclusion === "success") return "passed";
  return run.conclusion || "unknown";
}

function examine(input) {
  const revision = input?.revision;
  if (typeof revision !== "string" || !/^[a-f0-9]{40}$/.test(revision)) {
    throw new Error("Expected a full lowercase Git SHA");
  }
  const generated = input.generated_at;
  if (typeof generated !== "string" || !Number.isFinite(Date.parse(generated))) {
    throw new Error("Expected ISO reporting time");
  }
  const runs = {};
  for (const [key, name] of Object.entries(RUN_NAMES)) {
    const hit = latestRun(input.workflow_runs, name, revision);
    runs[key] = {
      state: workflowState(hit),
      id: hit?.id ?? null,
      url: hit?.id ? ACTIONS + hit.id : null,
    };
  }
  const deployJob = (input.deploy_jobs || []).find(job => job.name === "deploy");
  const smokeJob = (input.deploy_jobs || []).find(job => job.name === "smoke-deployed");
  const deployStatus = deployJob?.conclusion || (deployJob?.status === "in_progress" ? "pending" : "unverified");
  const smokeStatus = smokeJob?.conclusion || (smokeJob?.status === "in_progress" ? "pending" : "unverified");
  const publicLedger = input.public_ledger;
  const nowMs = Date.parse(generated);
  const siteAgeMs = nowMs - Date.parse(publicLedger?.checked_at || "");
  // A past green monitor becomes unknown again if daily checks stop.
  const publicCurrent = publicLedger?.passed === true &&
    publicLedger.expected_revision === revision &&
    Number.isFinite(siteAgeMs) && siteAgeMs >= 0 && siteAgeMs <= 2 * 86_400_000;
  const catalogLedger = input.catalog_ledger;
  const sources = catalogLedger?.sources;
  const sourceRows = sources && typeof sources === "object" ? Object.entries(sources) : [];
  const observedCounts = sourceRows.map(([name, row]) => ({
    source: name,
    count: row?.drift?.upstream_count ?? null,
    healthy: row?.ok === true && !row?.terms_link_error,
    checked_at: row?.last_checked || null,
  }));
  const catalogAgeMs = nowMs - Date.parse(catalogLedger?.checked_at || "");
  const catalogHealthy = sourceRows.length >= 3 &&
    observedCounts.every(source => source.healthy) &&
    Array.isArray(catalogLedger.alerts) && catalogLedger.alerts.length === 0 &&
    Number.isFinite(catalogAgeMs) && catalogAgeMs >= 0 && catalogAgeMs <= 10 * 86_400_000;
  const releases = Array.isArray(input.releases) ? input.releases : [];
  const kitRelease = releases.find(release =>
    /^reviewkit-egfr-luad-v1\./.test(release.tag_name || "") &&
    !release.draft && !release.prerelease) || null;
  const issueStates = Array.isArray(input.issue_states) ? input.issue_states : [];
  const openIssues = issueStates.filter(issue => issue.state === "open");
  const blockers = [];
  if (runs.validation.state !== "passed") blockers.push("Latest main validation not confirmed passing");
  if (runs.explorer.state !== "passed") blockers.push("Latest main Explorer build not confirmed passing");
  if (!publicCurrent) blockers.push("Public Explorer does not prove current revision");
  if (deployStatus !== "success" || smokeStatus !== "success") {
    blockers.push("Production deploy and live smoke not both confirmed successful");
  }
  if (input.branch_protected !== true) blockers.push("Main branch is not protected");
  if (!kitRelease) blockers.push("No verified durable reviewer-kit GitHub Release found");
  if (!catalogHealthy) blockers.push("Catalog health ledger is unavailable or has alerts");
  // Never infer clinical/scientific endorsement from GitHub issue state.
  blockers.push("Independent scientific interpretation and reviewer credentials require external human verification");
  return {
    generated_at: generated,
    revision,
    branch_protected: input.branch_protected === true,
    runs,
    deployment: {deploy: deployStatus, live_smoke: smokeStatus},
    public: {
      current: publicCurrent,
      checked_at: publicLedger?.checked_at || null,
      monitored_expected_revision: publicLedger?.expected_revision || null,
      url: publicLedger?.public_index_url || null,
      details: typeof publicLedger?.output === "string" ? publicLedger.output.slice(0, 500) : null,
    },
    catalog: {
      healthy: catalogHealthy,
      checked_at: catalogLedger?.checked_at || null,
      sources: observedCounts,
      alerts_count: catalogLedger?.alerts?.length ?? null,
    },
    releases: {review_kit_present: !!kitRelease, tag: kitRelease?.tag_name || null},
    open_issues: openIssues.map(issue => ({number: issue.number, title: issue.title})),
    blockers,
  };
}

function render(summary) {
  const sha = summary.revision;
  const toRun = (kind) => summary.runs[kind].url ?
    "[Actions run](" + summary.runs[kind].url + ")" : "No matching main run";
  const catalogCounts = summary.catalog.sources.length ?
    summary.catalog.sources.map(row => row.source + ": " +
      (row.count == null ? "unknown" : row.count) + " (" +
      (row.healthy ? "healthy check" : "needs attention") + ")").join("; ") :
    "No durable source observations";
  const siteStatus = summary.public.current ? "CURRENT — monitored match" :
    "STALE / NOT VERIFIED — investigate incident";
  const kitStatus = summary.releases.review_kit_present ?
    "Published tag " + summary.releases.tag + " (digest still requires independent verification)" :
    "No versioned reviewer-kit release found";
  const check = summary.generated_at;
  return [
    "# Open Cancer — Weekly Operations Readiness Report",
    "",
    "**Observed at (UTC):** " + check,
    "**Main Git revision:** [" + sha.slice(0, 12) +
      "](https://github.com/" + REPO + "/commit/" + sha + ")",
    "**Scope:** Engineering, public data metadata health, and release controls — **not** scientific or clinical validation.",
    "",
    "## Release status",
    "",
    "| Control | Evidence-based status | Link |",
    "| --- | --- | --- |",
    "| Source validation | **" + summary.runs.validation.state + "** | " + toRun("validation") + " |",
    "| Explorer artifact build | **" + summary.runs.explorer.state + "** | " + toRun("explorer") + " |",
    "| Production deploy | **" + summary.deployment.deploy + "** | " + toRun("explorer") + " |",
    "| Live post-deploy smoke | **" + summary.deployment.live_smoke + "** | " + toRun("explorer") + " |",
    "| Public Explorer API | **" + siteStatus + "** | [Release incident #29](" + ISSUE + "29) |",
    "| GitHub main protection | **" + (summary.branch_protected ? "protected" : "not protected") + "** | [Governance issue #18](" + ISSUE + "18) |",
    "| Reviewer-kit release | **" + kitStatus + "** | [Issue #17](" + ISSUE + "17) |",
    "",
    "## Upstream catalog health",
    "",
    "**Status:** " + (summary.catalog.healthy ? "Healthy bounded metadata check" :
      "Not verified healthy — investigate ledger") + ".",
    "**Last reported check:** " + (summary.catalog.checked_at || "not available") + ".",
    "**Source listings:** " + catalogCounts + ".",
    "Source health checks do **not** validate scientific claims, patient consent or redistribution rights.",
    "Full observation history: [catalog health ledger #27](" + ISSUE + "27).",
    "",
    "## Production API evidence",
    "",
    "Latest reported check: " + (summary.public.checked_at || "unavailable") +
      ". Monitored target: " + (summary.public.url || "not configured") + ".",
    "Reported expected revision: " +
      (summary.public.monitored_expected_revision || "not recorded") + ".",
    "When the checked revision differs from main, this report marks public readiness unverified.",
    "Public release and a successful **build** are separate states.",
    "",
    "## Blockers requiring action",
    "",
    ...summary.blockers.map(blocker => "- " + blocker),
    "",
    "## Issue triage",
    "",
    ...summary.open_issues
      .filter(issue => [3, 14, 17, 18, 29, 34, 35, 36].includes(issue.number))
      .sort((a, b) => a.number - b.number)
      .map(issue => "- [#" + issue.number + "](" + ISSUE + issue.number +
        "): " + issue.title),
    "",
    "### Acceptance gates",
    "",
    "1. **Production:** configure the actual hosting target and pass live API, revision and smoke checks before closing #14 or #29.",
    "2. **Governance:** activate main branch rules, required checks and author-independent CODEOWNERS approval before closing #18.",
    "3. **Reproducibility:** publish a rights-safe, tagged review archive and verify its independently pinned SHA before closing #17.",
    "4. **Scientific:** complete the CIViC QA (#34), GDC/cBioPortal cohort QA (#35), and an independently verified version-bound expert assessment (#36), tracked under milestone #3. The archived reviewer invitation #11 is not a scientific deliverable.",
    "",
    "### Auditability",
    "",
    "This is an automated snapshot of GitHub records and archived status ledgers, not a real-time medical or legal determination. It does not edit research results or catalog mappings. Each scheduled report adds a dated record; negative, disputed and pending findings remain visible.",
    "",
  ].join("\n");
}

module.exports = {examine, render, workflowState, latestRun};
