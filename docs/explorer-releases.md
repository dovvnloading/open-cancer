# Explorer release engineering and freshness

The Explorer currently has **two distinct publication surfaces**:

1. The original public demonstration is hosted at [open-cancer-explorer.lagging.chatgpt.site](https://open-cancer-explorer.lagging.chatgpt.site/). This hosting provider is **not connected to the GitHub Actions deployment workflow**. A merged GitHub change does not update that public address. Check the index's revision against main before describing it as current. The daily [freshness monitor](../.github/workflows/explorer-freshness.yml) reports an old or unavailable revision as a release failure, not as a failure of the cancer data.
2. The [GitHub Pages release workflow](../.github/workflows/explorer-release.yml) builds a validated, version-pinned static site and JSON API on every main push and uploads an immutable artifact. GitHub Pages publication is **opt-in until the repository owner enables Pages**. It does not redirect or update the chatgpt.site deployment.

## Enable automatic GitHub Pages publishing (repository owner)

1. In [repository Settings → Pages](https://github.com/dovvnloading/open-cancer/settings/pages), set **Source: GitHub Actions**.
2. In **Settings → Secrets and variables → Actions → Variables**, create the repository variable `OPEN_CANCER_GITHUB_PAGES_ENABLED` with value `true`.
3. Run **Build and optionally deploy Open Cancer Explorer** in the GitHub Actions tab, or merge a reviewed change into main. Inspect the `github-pages` environment and the deployment URL. The usual GitHub Pages project URL is `https://dovvnloading.github.io/open-cancer/` unless custom settings override it.
4. To monitor the new GitHub Pages deployment as production, also set repository variable `OPEN_CANCER_EXPLORER_PUBLIC_INDEX_URL` to its exact HTTPS `api/v1/index.json` URL. If unset, the monitor continues checking the original chatgpt.site demonstration, which might remain stale.
5. Open the deployed `api/v1/index.json` and verify its `revision` matches the published Git commit. Confirm catalog counts are nonempty, provenance links work and the pilot is correctly labeled unreviewed.
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
