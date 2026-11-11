"use strict";
/* Pure formatting/decoding for public Explorer deployment health issue alerts. */
function parsePrevious(body) {
  if (typeof body !== "string") return null;
  const match = body.match(/<!-- explorer-health-json\n([\s\S]*?)\n-->/);
  if (!match) return null;
  const data = JSON.parse(match[1]);
  if (!data || typeof data !== "object" || typeof data.passed !== "boolean") {
    throw new Error("Invalid previous Explorer deployment health record");
  }
  return data;
}

function renderHealth(data) {
  const {checked_at, expected_revision, public_index_url, passed, output} = data;
  if (!/^[a-f0-9]{40}$/.test(expected_revision)) throw new Error("Invalid expected Git SHA");
  if (typeof public_index_url !== "string" || !public_index_url.startsWith("https://")) {
    throw new Error("Public API must have an HTTPS URL");
  }
  if (typeof passed !== "boolean") throw new Error("Invalid monitoring status");
  const sanitized = String(output || "(no stdout)").slice(0, 3000).replaceAll("~~~", "---");
  const compact = {
    checked_at, expected_revision, public_index_url, passed,
    output: sanitized.slice(0, 1500),
  };
  return [
    "# Open Cancer Explorer public release health",
    "",
    "Status: " + (passed ? "**CURRENT**" : "**STALE OR UNAVAILABLE — ACTION REQUIRED**"),
    "",
    "Expected revision: " + expected_revision,
    "Public API: " + public_index_url,
    "Checked (UTC): " + checked_at,
    "",
    "This checks deployed metadata, basic source catalog integrity and revision freshness.",
    "**It is not scientific/clinical validation, nor proof of upstream licensing.**",
    "",
    "### Latest check",
    "~~~text", sanitized, "~~~",
    "",
    "Builds may pass while public deployments are skipped. Inspect the Pages deployment",
    "and source revision before closing production release issues.",
    "",
    "<!-- explorer-health-json", JSON.stringify(compact), "-->",
  ].join("\n");
}

function shouldRecordTransition(previous, current) {
  return !previous || previous.passed !== current.passed ||
    previous.expected_revision !== current.expected_revision ||
    previous.public_index_url !== current.public_index_url;
}

module.exports = {parsePrevious, renderHealth, shouldRecordTransition};
