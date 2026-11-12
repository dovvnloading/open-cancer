# Explorer release engineering and freshness

The Explorer has **two distinct hosting surfaces**:

1. **Production:** [the GitHub Pages Explorer](https://dovvnloading.github.io/open-cancer/) with [JSON API](https://dovvnloading.github.io/open-cancer/api/v1/index.json) is the supported public release. The [Pages release workflow](../.github/workflows/explorer-release.yml) builds and deploys a version-pinned artifact on main after the owner explicitly opts in, and a live smoke test verifies the deployed API, metadata and research links. The [first verified successful deployment](https://github.com/dovvnloading/open-cancer/actions/runs/43041807429) and [independent freshness check](https://github.com/dovvnloading/open-cancer/actions/runs/43041948889) passed on November 11, 2026 at revision `71c58b3302fd0bb339570a5d160ce115f8367d11`.
2. **Historical demonstration:** [open-cancer-explorer.lagging.chatgpt.site](https://open-cancer-explorer.lagging.chatgpt.site/) is separately hosted, not updated by GitHub Actions, and was found serving an older snapshot. Do not describe this secondary URL as current production.

The Pages opt-in has been activated and independently exercised. Future changes must still pass build, deployment, and live smoke verification; a green source build alone is not a published release.

## Confirm or restore GitHub Pages publishing (repository owner)

1. In [repository Settings → Pages](https://github.com/dovvnloading/open-cancer/settings/pages), set **Source: GitHub Actions**.
2. In **Settings → Secrets and variables → Actions → Variables**, create the repository variable `OPEN_CANCER_GITHUB_PAGES_ENABLED` with value `true`.
3. Run **Build and optionally deploy Open Cancer Explorer** in the GitHub Actions tab, or merge a reviewed change into main. Inspect the `github-pages` environment and the deployment URL. The usual GitHub Pages project URL is `https://dovvnloading.github.io/open-cancer/` unless custom settings override it.
4. To monitor the new GitHub Pages deployment as production, also set repository variable `OPEN_CANCER_EXPLORER_PUBLIC_INDEX_URL` to its exact HTTPS `api/v1/index.json` URL. If unset, the monitor continues checking the original chatgpt.site demonstration, which might remain stale.
5. Open the deployed `api/v1/index.json` and verify its `revision` matches the published Git commit. Inspect `api/v1/release.json` for the **actual artifact build timestamp** (not a scientific verification date). The footer shows both revision and build time when this metadata is available. After deployment the `smoke-deployed` job must succeed, checking the live index, separate collections, EGFR/LUAD crosswalk, pilot and version-pinned review links.
6. Keep an earlier, validated commit/artifact available for rollback. Never modify historical review-bound research snapshots when rolling back the site.

GitHub's Pages source selection and opt-in were verified working by the successful November 11 release. These are owner-controlled settings that can later be changed or disabled. The deploy job still skips when Pages opt-in is absent. Pages uses GitHub Actions short-lived ID tokens, not stored deployment secrets.

## Check production release freshness

Run this command from the checked-out revision expected in production:

~~~bash
python scripts/check_explorer_freshness.py \
  --url https://dovvnloading.github.io/open-cancer/api/v1/index.json \
  --expected-revision "$(git rev-parse HEAD)"
~~~

An API with a different revision is **STALE**, not silently treated as updated. After the next published site build, the footer will show the exact catalog revision. A green monitor proves only that the API revision and basic collections match, **not** that the scientific data or clinical interpretations are valid.

**Acceptance verified:** [QA issue #14](https://github.com/dovvnloading/open-cancer/issues/14) has objective build, deployment, post-deploy smoke and public-freshness evidence from the November 11, 2026 GitHub Pages rollout. That validation covers the actual Pages production URL, **not** the separate historical chatgpt.site host. Continue monitoring future revisions and record real regressions as incidents.


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

The Pages `deploy` job **passed** in [the first production workflow](https://github.com/dovvnloading/open-cancer/actions/runs/43041807429), along with its `smoke-deployed` job; earlier runs had skipped both while hosting was disabled. If the owner removes `OPEN_CANCER_GITHUB_PAGES_ENABLED=true`, future builds may succeed without deployment and the public freshness monitor must not misrepresent them as published. The separate chatgpt.site demonstration remains outside this workflow's publication control.


### Transient Pages responses and release verification

GitHub Pages may briefly return a 503 or serve an earlier revision while a deployment
propagates, even after the deploy action reports success. The **live** `smoke-deployed`
job now retries the **entire** published API/UI verification for up to 180 seconds
(5-second intervals) and reports which asset failed. A valid, matching index, all
separate collections, the research crosswalk, the version-pinned review links and
the build metadata must pass **together** in one attempt. An HTTP 200 response alone
is not success. After the deadline, stale, malformed, missing or unavailable assets
still fail the job. Local builds remain immediately checked without network retries.

The independent production-health monitor similarly allows up to 90 seconds for
temporary Pages errors or propagation before opening an incident. It never changes
the expected production SHA, accepts mismatched JSON or suppresses a continuing
outage. These delays are intentional: the health ledger records durable release
problems rather than a single transient failed request.

## Durable public-site freshness incidents

The [public Explorer freshness workflow](../.github/workflows/explorer-freshness.yml) checks the selected public JSON index **after each completed Explorer release workflow** (for a non-PR `main` run), daily at 11:37 UTC, and when manually dispatched. Checking after deployment avoids false stale alarms caused by probing the previous Pages revision while a new release is still building. It does not silently accept a failed deployment; genuine live mismatches still fail and update the incident ledger. It does **not** publish a site. If the live index is unreachable, malformed, points to the wrong Git repository or reports an unexpected revision/collection, the workflow fails and opens/reopens the durable GitHub issue **Public Explorer release health (automated)**. On a successful check, the ledger issue is closed with the checked SHA and time. Changes in expected production revision, endpoint URL or health state are recorded as issue comments, avoiding a flood of identical daily alerts.

The monitor defaults to comparing the **GitHub Pages production API** to the latest `main` SHA: `https://dovvnloading.github.io/open-cancer/api/v1/index.json`. The repository variable `OPEN_CANCER_EXPLORER_PUBLIC_INDEX_URL` can override this for a verified replacement production address; do not point it to the known-stale historical demo unless intentionally auditing that separate host. For an intentional rollback, an operator may set `OPEN_CANCER_EXPLORER_DESIRED_SHA` to the known-good deployed full 40-character SHA, then restore it to match `main` once the replacement is deployed. Never set a desired SHA just to mask an unexplained stale site.

The issue records deployment **health**, not biomedical validation, scientific review, or legal permission. The release workflow separately distinguishes a successful build from a skipped Pages deployment. A green source-build workflow cannot close a stale public-site incident unless the live public index responds with the expected revision and required data collections.
