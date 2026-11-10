# Data, privacy, and research integrity policy

Open Cancer is a public GitHub project. It **must not** be used as a place to deposit patient-level clinical data, identifying or potentially re-identifiable records, raw human genomes, health information, private access tokens, or any data under controlled-access conditions.

## Where data lives

The default is **metadata and links to the authoritative source**. A catalog record saying `access: "open"` describes discoverability, not a legal right to mirror its contents. Public access, licensing, consent, citation, and permitted redistribution are distinct questions.

Small derived data may be included only when (1) source terms permit redistribution, (2) privacy and consent conditions are met, (3) attribution and methodology are recorded, and (4) the result cannot reasonably expose person-level information. Prefer publishing large permitted data in a suitable research repository with persistent identifiers and linking from here.

The MIT license in this repository applies to the project's original code/documentation as appropriate; it **does not relicense third-party data**.

## Scientific claims

- Treat agent-generated outputs as **unreviewed hypotheses or evidence maps**, not established facts.
- Support claims with source-record URLs, source versions/access dates, transparent methods, limitations, and contradictory evidence where relevant.
- Distinguish observational associations, preclinical evidence, clinical trial results, and validated clinical guidance.
- Avoid patient-specific treatment, diagnosis, or drug recommendations.

## Agent and contributor safeguards

- Do not scrape sources in violation of terms or API rate limits; never circumvent authentication or access controls.
- Do not try to re-identify research participants or combine datasets to infer identities.
- Do not commit secrets or submit source materials without verified redistribution rights.
- Do not create misleading synthetic results presented as real experiments.
- When in doubt about confidentiality or licensing, **link to metadata instead of uploading data** and ask maintainers to review.

The catalog records are editorial snapshots; users must verify upstream conditions before use.
