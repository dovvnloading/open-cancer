# Durable EGFR/LUAD independent reviewer kit releases

The original [review packet workflow](../.github/workflows/reviewer-packet.yml) creates temporary GitHub Actions artifacts (retained up to 90 days). The separate [reviewer release workflow](../.github/workflows/reviewer-release.yml) publishes a versioned **GitHub Release** containing only the deterministic non-sensitive reviewer ZIP and its SHA-256 file. Releases have no 90-day Actions artifact expiration, but repository owners can still delete them, and public availability depends on GitHub hosting.

## Publication runbook (maintainer action required)

1. Review the exact source commit, its DATA_POLICY.md compliance and the review ZIP's allowlist. Do not include raw upstream mutation API responses, patient data, restricted sources, private tokens, or medical claims.
2. Protect tags matching reviewkit-egfr-luad-v1.* and, if desired, require a human approver in the reviewer-release GitHub Actions environment. Create and push a new immutable-by-policy annotated tag, e.g. reviewkit-egfr-luad-v1.0.0, **at the reviewed commit**. Never move a published tag or replace its release assets.
3. The tag-triggered workflow checks repository validation and packet tests, builds the ZIP from its exact Git SHA, emits egfr-luad-review-kit.zip.sha256, creates the release, downloads the public asset and checks its SHA. Read the workflow output before sharing the release link.
4. Copy the exact release tag and verified ZIP digest into repository Actions variables OPEN_CANCER_REVIEW_KIT_TAG and OPEN_CANCER_REVIEW_KIT_SHA256. A scheduled and manually dispatchable verifier checks that pinned release URL and digest each Monday. If either variable is absent, the scheduled/manual verifier **fails closed**, recording a configuration error instead of reporting a green run. A passing job requires a retrievable tag and a matching independently pinned digest.
5. Keep the public Release URL, checksum, source revision and original manifest digest together in external reviewer communications. Maintain a secondary permitted archive if institutional durability requires it; do not depend on mutable GitHub Actions run artifact links.

A sample workflow for a reviewed commit is:

~~~sh
git tag -a reviewkit-egfr-luad-v1.0.0 <REVIEWED_COMMIT_SHA> -m "Nonclinical EGFR/LUAD reviewer kit v1.0.0"
git push origin reviewkit-egfr-luad-v1.0.0
~~~

Do not present a candidate release URL as live until the release has actually been created and verified.

## What can and cannot be independently replayed

- **Method reproduction:** project collaborators reproduced the code and aggregate derived output on the original accessible upstream responses; external scientist interpretation review is still pending. Rebuilding the offline packet from the exact Git revision is deterministic and verifies the source manifest and artifact hashes.
- **Source-version drift:** the original manifest and packet-index.json retain source URLs, filters, accessed dates, response SHA-256, GDC project and gene identifiers, cBioPortal study ID and CIViC nightly export filename. A new upstream retrieval may differ; the reviewer should compare hashes and document any drift rather than silently treating changed source bytes as an exact rerun.
- **Exact-input replay:** the project does not archive the original raw upstream source responses, for privacy and data-use reasons. Their SHA fingerprints identify original content **only if the reviewer lawfully obtains matching bytes**. Source accessions identify studies/queries but are not guaranteed immutable releases; packet-index.json records immutable_release_url as null until an authoritative, rights-compatible fixed release can be verified.

A durable reviewer kit proves that *derived files and their linked research revision* are obtainable and integrity checked. It does not guarantee original byte-identical upstream replay, scientific accuracy, verified reviewer identity, clinical validity, or broad redistribution rights.
