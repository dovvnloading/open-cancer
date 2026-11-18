# Job 5 — Does SINE annotation add whole-order neoplasia transport beyond L1?

**Verdict: CONSISTENT under the frozen predictive protocol.** The primary nominal equal-order held-out binomial log-loss improvement from adding `log1p(SINEs_Counts)` to L1/longevity/assembly controls is **+0.010513743** (positive favors SINE), with a **complete-refit order-jackknife approximate 95% interval [+0.001361257, +0.048909796]**. The interval is centered on the bias-corrected jackknife center, **+0.025135526**, *not* the unadjusted point gain. **All 45 rounding-compatible neoplasia numerator assignments satisfy the frozen consistency gates**, including whole-order deletions and the prespecified assembly-quality restriction. See [REPORT.md](REPORT.md) for all 11 order scores, numerical checks and limitations.

**This is nonclinical, unreviewed research output — not medical or veterinary advice.** The source panel **excludes species with zero recorded neoplasia or malignancy**. Inherited reference-genome annotations are **not observations of ongoing somatic retrotransposition**. The apparent incremental portability is conditional on this selected 55-species panel, non-independent pathology sources, a fixed regression model and approximate order jackknife. It establishes **no causal mechanism, new biological association, validated clinical predictor, general-mammal risk claim, or independent outcome replication.** The original L1/SINE associations are prior published work.

## Immutable ordering and files

The [PROTOCOL.md](PROTOCOL.md) was separately committed at [`40521aa2`](https://github.com/dovvnloading/open-cancer/commit/40521aa2e96a801bca69fe5d2258aa90d1a9d786) **before any result computation**. This directory and the single new offline test are the only new files; no existing repository file was touched.

| File | Role |
| --- | --- |
| [PROTOCOL.md](PROTOCOL.md) | Frozen panel, exact single comparison, training transformations, denominator and rounding gates, full order-refit uncertainty, sensitivity and verdict thresholds |
| [data/species_aggregates.csv](data/species_aggregates.csv) | **Derived** minimal per-species covariates, denominator-linked nominal counts, all compatible numerator sets and ratio audit indicator (not mirrored source CSV) |
| [data/source_audit.json](data/source_audit.json) | Pinned source, upstream blob, licence decision, enumerated ambiguity and failed secondary gate |
| [analysis.py](analysis.py) | Pure-stdlib deterministic frozen Newton logistic regression, entire 45-assignment grid, complete jackknife and sensitivities, offline and without source fetching |
| [REPORT.md](REPORT.md) | Full result, order table, error analysis, sensitivity table, literature survey and hostile review |
| [../../tests/test_research_sine_l1_transport.py](../../tests/test_research_sine_l1_transport.py) | Ten offline integrity tests for counts, hashes, leakage, Newton conditions, 45 refits, and verdict boundaries; no numeric optimizer-trace snapshots |

## Reproduce

From the repository root with a stock Python 3.12 or 3.14 interpreter, without network or optional packages:

```sh
python -m unittest discover -s tests -p "test_research_sine_l1_transport.py" -v
python research/sine-l1-transport/analysis.py --output sine-l1-results.json
```

The script writes every held-out prediction, order score, order-deletion pseudovalue, all 45 assignment verdicts, diagnostics, sensitivity scores and **newline-normalized SHA256** keys with POSIX relative paths. The generated `sine-l1-results.json` is a local reproduction artifact, not a checked-in third-party dataset. Tests complete in approximately 16 seconds on the development environment; no network access is needed. Numerical snapshot tolerance is 3e-7 absolute/relative; optimizer gradients are checked against the **frozen analytic bounds**, not exact floating values or step counts.

## Attribution and redistribution

**Primary source:** V. Nikoonejad Fard and M. Tollis (2025), [*Recent Non-LTR Retrotransposon Activity Predicts Cancer Prevalence in Mammals*](https://pmc.ncbi.nlm.nih.gov/articles/PMC12424770/), preprint [DOI 10.1101/2025.08.30.673284](https://doi.org/10.1101/2025.08.30.673284). Data: [pinned supplementary CSV](https://github.com/VNF1981/Non-LTR-Retrotransposons-Cancer-Project/blob/792af204b4113238276b3ac43e511b62872dd368/Supplementary%20Data%20S1.csv), commit `792af204b4113238276b3ac43e511b62872dd368`, Git blob `257a6089179942e98e1cc454fa0e3f762c83218f`. Retrieved solely as text through the GitHub repository API; no browser download or download prompt was used. Article: **CC BY 4.0**. Pinned third-party repository: **no explicit LICENSE found; GitHub reports licence null**. Public availability cannot supply permission to mirror its content. Only original code and narrowly derived numerical facts are committed, subject to [DATA_POLICY.md](../../DATA_POLICY.md); the project's MIT licence does not relicense source materials.

**Related audit and overlapping outcome provenance:** [MHC-I mortality transport (#79)](../mhc-cancer-mortality-transport/) and [MHC-I pathology prevalence (#86)](../mhc-pathology-prevalence/). In particular, #86 identifies source dictionary disagreements between malignant counts and displayed prevalence. Here the secondary `MalignancyPrevalence` gate **fails** a stricter internal Ratio/count reconciliation and is not scored. Treat it as an unavailable secondary analysis, not a null result or an implied biological progression test.
