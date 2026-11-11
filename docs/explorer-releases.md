# Explorer release engineering and freshness

The Explorer currently has **two distinct publication surfaces**:

1. The original public demonstration is hosted at [open-cancer-explorer.lagging.chatgpt.site](https://open-cancer-explorer.lagging.chatgpt.site/). This hosting provider is **not connected to the GitHub Actions deployment workflow**. A merged GitHub change does not update that public address. Check the index's revision against main before describing it as current. The daily [freshness monitor](../.github/workflows/explorer-freshness.yml) reports an old or unavailable revision as a release failure, not as a failure of the cancer data.
2. The [GitHub Pages release workflow](../.github/workflows/explorer-release.yml) builds a validated, version-pinned static site and JSON API on every main push and uploads an immutable artifact. GitHub Pages publication is **opt-in until the repository owner enables Pages**. It does not redirect or update the chatgpt.site deployment.

## Enable automatic GitHub Pages publishing (repository owner)

1. In [repository Settings → Pages](https://github.com/dovvnloading/open-cancer/settings/pages), set **Source: GitHub Actions**.
2. In **Settings → Secrets and variables → Actions → Variables**, create the repository variable `OPEN_CANCER_GITHUB_PAGES_ENABLED` with value `true`.
3. Run **Build and optionally deploy Open Cancer Explorer** in the GitHub Actions tab, or merge a reviewed change into main. Inspect the `github-pages` environment and the deployment URL. The usual GitHub Pages project URL is `https://dovvnloading.github.io/open-cancer/` unless custom settings override it.
4. To monitor the new GitHub Pages deployment as production, also set repository variable `OPEN_CANCER_EXPLORER_PUBLIC_INDEX_URL` to its exact HTTPS `api/v1/index.json` URL. If unset, the monitor continues checking the original chatgpt.site demonstration, which might remain stale.
5. Open the deployed `api/v1/index.json` and verify its `revision` matches the published Git commit. Inspect `api/v1/release.json` for the **actual artifact build timestamp** (not a scientific verification date). The footer shows both revision and build time when this metadata is available. After deployment the `smoke-deployed` job must succeed, checking the live index, separate collections, EGFR/LUAD crosswalk, pilot and version-pinned review links.
6. Keep an earlier, validated commit/artifact available for rollback. Never modify historical review-bound research snapshots when rolling back the site.

GitHub's documentation requires a one-time Pages source selection in Settings. This repository's connector cannot perform that account-level activation; the workflow intentionally skips deployment until explicit opt-in. Pages uses GitHub Actions short-lived ID tokens, not repository-stored deployment secrets.

## Check production release freshness

Run this command from the checked-out revision expected in production:

~~~bash
python scripts/check_explorer_freshness.py \
  --url https://open-cancer-explorer.lagging.chatgpt.site/api/v1/index.json \
  --expected-revision "$(git rev-parse HEAD)"
~~~

An API with a different revision is **STALE**, not silently treated as updated. After the next published site build, the footer will show the exact catalog revision. A green monitor proves only that the API revision and basic collections match, **not** that the scientific data or clinical interpretations are valid.

**Follow-up:** keep [QA issue #14](https://github.com/dovvnloading/open-cancer/issues/14) open until a production target is demonstrably current and its freshness monitor passes. Automated deployment to the current chatgpt.site host requires that host's deployment controls and is not implemented here.


## Post-deployment smoke checks and rollback

The build workflow writes `api/v1/release.json` with `built_at` in UTC, source revision and GitHub Actions run ID. This file is emitted **after** the deterministic source build, so `python scripts/build_explorer.py` remains bit-for-bit reproducible. The timestamp is when the artifact was built, not necessarily when the host served it. The original experimental evidence and reviewed research manifests are never rewritten.

Before upload, `scripts/smoke_explorer_release.py` checks the local output. After an opt-in Pages deployment succeeds, the same command fetches the **public** site and verifies the same SHA in the index and build metadata, cross-collection JSON roundtrips, discoverable GDC/cBioPortal/CIViC EGFR links, the unreviewed pilot, and manifest/review URLs bound to that Git revision. These checks validate public release consistency, **not** biology or clinical validity. They do not guarantee every upstream external website remains online.

For a manual smoke test:

~~~sh
python scripts/smoke_explorer_release.py \
  --url https://dovvnloading.github.io/open-cancer/ \
  --expected-revision "$(git rev-parse HEAD)" \
  --require-build-info
~~~

If the live Pages deploy fails or becomes stale, inspect its Actions `deploy` and `smoke-deployed` jobs and the published SHA. A maintainer can restore an earlier verified static revision by rerunning that successful tagged/commit-bound release workflow or by publishing the previously retained signed-off artifact through the approved Pages deployment process. Verify the reverted public API reports the **rollback revision**, not current `main`, and temporarily pin freshness expectations to the designated production revision if required. Never label a rollback as a new scientific result.

At present the Pages `deploy` job may be **skipped** by design until the repository owner enables Pages and `OPEN_CANCER_GITHUB_PAGES_ENABLED=true`. A green build does not imply a live production deployment. The separate chatgpt.site demonstration remains outside this workflow's publication control.


## Durable public-site freshness incidents

The [public Explorer freshness workflow](../.github/workflows/explorer-freshness.yml) checks the selected public JSON index on every push to `main`, daily at 11:37 UTC and when manually dispatched. It does **not** publish a site. If the live index is unreachable, malformed, points to the wrong Git repository or reports an unexpected revision/collection, the workflow fails and opens/reopens the durable GitHub issue **Public Explorer release health (automated)**. On a successful check, the ledger issue is closed with the checked SHA and time. Changes in expected production revision, endpoint URL or health state are recorded as issue comments, avoiding a flood of identical daily alerts.

The monitor defaults to comparing the original `chatgpt.site` demonstration with the latest `main` SHA. To switch production to a deployed GitHub Pages site, set `OPEN_CANCER_EXPLORER_PUBLIC_INDEX_URL` to the **actual** Pages index URL after the first deployment. For an intentional rollback, an operator may set `OPEN_CANCER_EXPLORER_DESIRED_SHA` to the known-good deployed full 40-character SHA, then restore it to match `main` once the replacement is deployed. Never set a desired SHA just to mask an unexplained stale site.

The issue records deployment **health**, not biomedical validation, scientific review, or legal permission. The release workflow separately distinguishes a successful build from a skipped Pages deployment. A green source-build workflow cannot close a stale public-site incident unless the live public index responds with the expected revision and required data collections.
