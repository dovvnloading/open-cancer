# Publication governance and branch protection runbook

The public Open Cancer API contains externally citable research metadata. **JSON Schema validity, artifact checksums, a green workflow, and a GitHub merge do not prove scientific accuracy, redistribution rights, independent expert review, or clinical safety.**

## Current controls and known limitations

- The source tree records featured pilot decisions separately from submitted results in [catalog/publication-decisions.json](../catalog/publication-decisions.json). Each publication entry is pinned to the **exact original result-manifest SHA-256**. Unpinned submissions remain in the separate experimental collection.
- [schemas/publication.schema.json](../schemas/publication.schema.json) and [scripts/validate.py](../scripts/validate.py) check decision structure, uniqueness and digest freshness. A curated decision requires a self-declared curator name. Neither the schema nor the CI process authenticates that person's identity.
- [.github/CODEOWNERS](../.github/CODEOWNERS) routes registry, gate and CI changes to the owner, but **CODEOWNERS does not itself require approval**. As of the November 11, 2026 operations review, GitHub reported the main branch as **unprotected**. Human publication approval is not yet a fully enforceable repository requirement.
- The validation workflow also performs a **PR-only independent GitHub review gate** on changes to the publication-decision registry, using the evaluator checked out from the exact trusted **base SHA** (not the PR-submitted evaluator). It requires an `APPROVED` review for the **current head commit** by a separate human GitHub owner/member/collaborator. This authenticates GitHub account and review state, **not** scientific credentials, independent biomedical judgment or source rights. A malicious direct push or unauthorized workflow change can bypass it until branch ruleset checks are required.
- The original EGFR/LUAD demo is marked legacy-featured for continuity only. Its scientific interpretation has **not** been independently adjudicated, and it is not a clinical asset.

## Configure mandatory independent change review (repository owner)

1. Add a **second, trusted human maintainer** with appropriate repository permissions and expertise in the project's source and data policies. GitHub does not count the author approving their own pull request. A one-person CODEOWNERS roster can therefore block all owner-authored changes if reviews are made mandatory. Do not use a shared account or fictional reviewer to work around this.
2. In repository **Settings → Rules → Rulesets** (or the branch protection settings available in your GitHub plan), add an **active branch ruleset** targeting main. Require changes through pull requests, at least one approval from a reviewer other than the author, and **Code Owner review** for publication-decision changes. Require fresh approval after new commits and block force-pushes/deletion.
3. Require repository validation and Explorer build checks to pass. Select the actual status check names shown on successful Actions runs; check names can differ from workflow titles. Do not set a required check for the optional public deployment job while Pages publication is intentionally disabled.
4. Minimize/administer any bypass list. In particular, direct pushes and self-approval must not bypass changes to catalog/publication-decisions.json. Enable required review by a separate maintainer before the next curated-featured addition.
5. Test the policy safely with a harmless pull request to the publication registry: an unapproved change must be blocked, and an authorized second maintainer's approval must make it mergeable only after CI passes. Review current branch/ruleset status with GitHub settings or its protected-branch APIs and document the result on [issue #18](https://github.com/dovvnloading/open-cancer/issues/18).

## Editorial publication checklist

Before setting a new result to curated-featured, the maintainer should separately inspect: the original manifest SHA-256 and reproducibility method; source links and observed versions; source reuse/attribution restrictions and absence of patient-identifiable data; scientific limitations, negative evidence, and conflicts; and any existing independent expert reviews. Confirm any human credentials off-platform, rather than trusting the self-declared reviewer name. Only then record an auditable rationale, decision date and curator in a reviewed PR.

A curated publication is **discoverable**, not medically endorsed. An independent scientific assessment and clinical validation are different stages and are never implied by inclusion in the public API.


## When the publication-review job fails

If a PR modifies `catalog/publication-decisions.json` without fresh, independent authorized GitHub approval, the PR-only `publication-review` check fails **closed**. A separate trusted human collaborator must inspect the exact manifest digest, source-rights evidence and scientific caveats, approve the current PR head, and then rerun the failed Actions check (or use the host's approved re-evaluation workflow). Editing a self-declared `curator` JSON field, asking a bot, or approving your own PR does not satisfy the gate. An approver may reasonably submit `CHANGES_REQUESTED` or `DISMISSED`, which invalidates earlier approvals.

The gate is **defense in depth**. To make it enforceable, configure `publication-review` as a required PR status check together with the main `validate` job using GitHub rulesets/branch protection, enforce author-independent CODEOWNERS approval, and block direct pushes and bypasses. Without that owner-controlled configuration, a direct push to an unprotected branch can still avoid PR review altogether.

Do not close [issue #18](https://github.com/dovvnloading/open-cancer/issues/18) until the branch/ruleset reviewer gate has been enabled and tested. Publication policy must be enforced by the repository host in addition to application-level validation.
