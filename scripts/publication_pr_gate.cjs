"use strict";

/*
 * GitHub PR editorial approval check for the separately curated featured
 * research registry. This authenticates GitHub review metadata only.
 * It cannot prove scientific interpretation, reviewer credentials, or a
 * trusted branch-protection setting. The protected base ref must supply
 * the implementation when used in Actions.
 */
const REGISTRY = "catalog/publication-decisions.json";
const AUTHORIZED_ASSOCIATIONS = new Set(["OWNER", "MEMBER", "COLLABORATOR"]);
const TERMINAL_REVIEW_STATES = new Set(["APPROVED", "CHANGES_REQUESTED", "DISMISSED"]);

function touchesRegistry(files) {
  if (!Array.isArray(files)) throw new TypeError("PR changed-file list was not available");
  return files.some(file => file && (file.filename === REGISTRY ||
                                    file.previous_filename === REGISTRY));
}

function evaluateEditorialReview({files, reviews, author, headSha}) {
  if (!touchesRegistry(files)) {
    return {required: false, ok: true,
            reason: "No change to the editorial publication decision registry"};
  }
  if (typeof author !== "string" || !author.trim() ||
      !/^[a-f0-9]{40}$/.test(headSha) || !Array.isArray(reviews)) {
    return {required: true, ok: false,
            reason: "Missing PR author, exact head SHA, or GitHub reviews"};
  }
  // A person's most recent substantive verdict supersedes earlier verdicts.
  // Comment-only submissions should not erase a valid earlier approval.
  const latest = new Map();
  for (const review of [...reviews].sort((a, b) => {
    const time = String(a?.submitted_at || "").localeCompare(String(b?.submitted_at || ""));
    return time || Number(a?.id || 0) - Number(b?.id || 0);
  })) {
    const login = review?.user?.login;
    if (typeof login !== "string" || !TERMINAL_REVIEW_STATES.has(review?.state)) continue;
    latest.set(login.toLowerCase(), review);
  }
  const authorLogin = author.toLowerCase();
  for (const [login, review] of latest) {
    if (login === authorLogin || review.state !== "APPROVED" ||
        review.commit_id !== headSha || review.user?.type !== "User" ||
        !AUTHORIZED_ASSOCIATIONS.has(review.author_association)) continue;
    return {required: true, ok: true,
            reason: "Independent trusted GitHub review of the current head was submitted",
            reviewer: review.user.login};
  }
  return {required: true, ok: false,
          reason: "Publication registry changes require an APPROVED review of the exact head commit " +
                  "from a separate human GitHub owner/member/collaborator. " +
                  "Self-approval, stale approvals, bots, or outside-only commenters do not qualify."};
}

module.exports = {REGISTRY, touchesRegistry, evaluateEditorialReview};
