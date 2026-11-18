# Frozen job-5 result — SINE annotation proxy beyond L1

**Status (2026-11-18):** original, **nonclinical and unreviewed** statistical research. **Not medical or veterinary advice.** All numbers below are conditional on the source panel and [prefit frozen protocol](PROTOCOL.md) (commit [`40521aa2e96a`](https://github.com/dovvnloading/open-cancer/commit/40521aa2e96a801bca69fe5d2258aa90d1a9d786)). A favorable whole-order score is not a new biological finding, causal mechanism or clinical predictor.

## Literature survey, scientific novelty, and measurement boundaries

The prior Nikoonejad Fard and Tollis (2025) article/preprint [PMC12424770](https://pmc.ncbi.nlm.nih.gov/articles/PMC12424770/), [DOI 10.1101/2025.08.30.673284](https://doi.org/10.1101/2025.08.30.673284), reports associations of recent non-LTR retrotransposon annotations with comparative mammalian cancer prevalence. Its PGLS measures include individual L1 and combined L1/SINE summaries. Published abundance–malignancy associations were diminished after controlling for longevity. These are **prior associations** and do not independently demonstrate the **held-out incremental SINE-only information** beyond a full L1/longevity/assembly baseline. A targeted literature check including the publication and its cited preprint version did **not find** an already reported whole-order held-out conditional SINE vs L1 test with this complete baseline; this is a *limited survey*, not proof of historical priority or a new biological association. No pivot was required.

This work adds only an independently prespecified **cross-validation and denominator audit**. It mirrors the complete-refit order-jackknife and equal-order/species balancing in [MHC-I mortality transport (#79)](../mhc-cancer-mortality-transport/) and [MHC-I pathology prevalence (#86)](../mhc-pathology-prevalence/). Those studies and this source share pathology ascertainment/provenance families, so the data are **not an independent biological replication**. Pathology prevalence is neither lifetime incidence, functional driver activity, nor attributed cancer mortality. Species selection, necropsy sampling, phylogenetic relatedness and annotation/assembly quality remain unmeasured biases.

**Critical scope:** The published S1 panel **excluded species with zero recorded neoplasia or malignancy**. The 55 rows must not be generalized to the source-eligible mammalian universe. L1 and SINE counts are inherited *reference-genome annotations*, **not measurements of ongoing somatic retrotransposition**. The source `Ratio` does not identify a biological benign-to-malignant transition probability. No causal retrotransposition mechanism or new association is established.

## Sources, licence and data decisions

| Source | Version/licence decision | Use |
| --- | --- | --- |
| [Nikoonejad Fard & Tollis preprint/article](https://pmc.ncbi.nlm.nih.gov/articles/PMC12424770/) | PMC bioRxiv article marked **CC BY 4.0**; access 2026-11-18 | Citation and literature context; no article reproduction |
| [Supplementary Data S1.csv](https://github.com/VNF1981/Non-LTR-Retrotransposons-Cancer-Project/blob/792af204b4113238276b3ac43e511b62872dd368/Supplementary%20Data%20S1.csv) | Immutable commit `792af204b4113238276b3ac43e511b62872dd368`, source blob `257a6089179942e98e1cc454fa0e3f762c83218f`; repository has **no LICENSE** at that commit and GitHub API licence **null** | Accessed **only as GitHub API plain text**, not by browser download; no upstream CSV or scripts redistributed |
| [Compton pathology dictionary in PR #86](../mhc-pathology-prevalence/) | Dictionary/count arithmetic conflict, separate licences | Endpoint-provenance caution; no new Compton data copied |
| [Open Cancer data policy](../../DATA_POLICY.md) | Open Cancer MIT does not apply to third-party data | Only derived species-level numerical aggregates, original stdlib code, tests, documentation |

Because a licence to **redistribute the GitHub source data is not established**, `data/species_aggregates.csv` is deliberately *not* the raw S1 CSV or a reconstruction of all columns: it retains the minimum independently derived predictors and compatible integer-count sets needed to replay the fixed test. No individual necropsy, private human information, raw genome, reference annotation or article figure is included.

## Prefit data gate and frozen model

- **55 unique species, 11 source taxonomic orders**, **3,533 total necropsies**, nominal **691 neoplasia** and **405 malignancy** numerator equivalents; **9 L1-zero** species retained via `log1p(0)=0`.
- Exact decimal token-specific half-unit rounding, with frozen `2e-15` export allowance, admits **45 primary numerator assignments**: `Cercopithecus_mitis` (2–6 of 40; nominal 4), `Giraffa_camelopardalis` (2–4 of 30; nominal 3), `Panthera_pardus` (5–7 of 20; nominal 6). All other primary counts are unique. All 45 variants are fully refitted and jackknifed, **not** post-hoc sampled. No count was silently forced to its nominal fraction.
- The first prefit descriptive sum omitted the 3 rows with non-unique counts (**678/397**); hostile review caught and corrected the audit to **691/405** before final reporting. Source and panel count constraints were unchanged; the frozen protocol was **not rewritten**.
- **Secondary malignancy fails its prespecified data gate**: nonunique malignancy numerators in 2 rows, nonunique neoplasia numerators in 3, and **31 printed `Ratio` tokens incompatible at their literal displayed precision with integer malignancy/neoplasia ratios**. This likely reflects roundoff/provenance in ratios computed from displayed rates. It is not by itself a falsification of underlying clinical data. No numerator editing, borrowed `MalignancyKnown`, guessed provenance, exploratory second endpoint scoring or multiplicity-driven search was permitted. **Secondary: untestable / withheld.** The planned Bonferroni two-endpoint intervals do not apply without a usable second endpoint.
- Exactly one baseline regression: binomial logit ridge with `log1p(L1_Counts)`, maximum longevity (years), `Assembly_BUSCO_Scores`, `log1p(Contig_N50)`; candidate adds only `log1p(SINEs_Counts)`. Fixed slope penalty 1, unpenalized intercept, training-only unweighted z-scales, 11 leave-one-order-out held-out prediction refits; test score averages species within order and then orders equally. No held-out order contributes labels, transformations or tuning.
- Full one-order deletions repeat the *entire* ten-order prediction/refit procedure and generate (q_g=11\Delta-10\Delta_{(-g)}); 95% interval = bias-corrected jackknife mean +/- (t_{10,0.975}=2.2281388519649385) times its SE. This is an approximate related-order uncertainty interval, **not a confidence interval for independent animals**.

## Primary held-out neoplasia results

**Verdict: CONSISTENT** with incremental SINE information **under the exact frozen predictive rule**, and **only** on this selected source panel. The directional gate requires all 45 assignment gains and jackknife lower bounds positive, more than half the held-out orders improving, every full delete-order gain positive and all prespecified training/assembly sensitivities positive. All gates pass; no alternative annotation metric was searched after scoring.

| Scoring estimand | Baseline L1 + covariates | Added SINE | Gain (baseline − added) |
| --- | ---: | ---: | ---: |
| **Equal-order mean binomial log loss (primary)** | 0.453818394643 | 0.443304651628 | **+0.010513743016** |
| Necropsy-weighted log loss (descriptive) | 0.494446834943 | 0.489538157611 | **+0.004908677332** |

For nominal reconstructed numerators (not uniquely inferred underlying counts): full-refit **order jackknife bias-corrected center +0.025135526226**, SE **0.010670012464**, approximate **95% interval [+0.001361256904, +0.048909795547]**. The positive lower bound supports the *specific* conditional consistency verdict; it does **not** estimate a portable biological SINE effect. Of 11 held-out orders, **8** improve; no probabilities required numerical clipping.

### Whole-order out-of-sample gains and fully refitted order deletions

Positive gain means the SINE candidate has *lower* out-of-sample binomial log loss. Negative held-out-order gains are not hidden.

| Held-out order | Species | Nominal held-out gain | Refit gain when that order is deleted altogether |
| --- | ---: | ---: | ---: |
| Artiodactyla | 13 | +0.007483614 | +0.007164528 |
| Carnivora | 16 | **−0.006235896** | +0.011187543 |
| Chiroptera | 4 | +0.015994723 | +0.007496866 |
| Cingulata | 1 | +0.012304257 | +0.009768546 |
| Diprotodontia | 2 | **−0.019356375** | +0.014808514 |
| Hyracoidea | 1 | +0.000977538 | +0.011659477 |
| Lagomorpha | 1 | **−0.002411693** | +0.012490817 |
| Pilosa | 1 | +0.036128639 | +0.007170839 |
| Primates | 9 | +0.014834206 | +0.009334766 |
| Proboscidea | 1 | +0.044023862 | +0.006746552 |
| Rodentia | 6 | +0.011908300 | **+0.001738763** |

The full-refit deleted-order gain ranges from **+0.001738763** (Rodentia removed) to **+0.014808514** (Diprotodontia removed). Rodentia deletion is notably close to zero, and five orders have only one or two species, making apparent transport vulnerable to selection and phylogenetic structure despite the positive nominal interval.

### Frozen sensitivity checks (nominal counts)

| One-at-a-time sensitivity, with whole-order refits | Species / orders | Equal-order gain |
| --- | --- | ---: |
| Equal-species effective training counts | 55 / 11 | +0.006516192 |
| Ridge penalty 0.1 | 55 / 11 | +0.010541340 |
| Ridge penalty 10 | 55 / 11 | +0.010258681 |
| Training-fold N50 log-value winsorization at 10/90% | 55 / 11 | +0.010830157 |
| Assembly quality: BUSCO ≥ 90 and N50 ≥ 100,000 bases | **38 / 9** | +0.004578891 |

The restricted assembly-quality panel meets the fixed ≥8 order gate. It is **not** a tuned choice. Necropsy-weighted quality-restriction gain is actually **negative**: `0.5513320952−0.5525184959 = −0.0011855007`, an adverse descriptive observation that cannot be substituted for the frozen equal-order estimand. All listed robustness gains remain positive for **all 45 compatible numerator assignments**, including every complete-order deletion; the smallest quality-restricted gain across assignments is **+0.004188929**.

### Exhaustive primary rounding-compatibility uncertainty

| Quantity across 45 numerator assignments | Range |
| --- | --- |
| Full-panel primary gain | **[+0.009535222, +0.011529693]** |
| Approximate 95% order-jackknife **lower** endpoints | **[+0.000050994, +0.002355169]** |
| Approximate 95% order-jackknife **upper** endpoints | **[+0.045780597, +0.052566444]** |
| Improving orders | **8 of 11 in every assignment** |
| Full frozen verdict | **consistent in all 45**, no untestable fit |

The compatible-set envelope is **source display-rounding uncertainty**, not a new 95% confidence interval. The lower bound across assignments is only **+0.000050994**: numerical margins are close to the directional gate for one candidate count assignment. Approximate order-jackknife coverage and outcome-selected species ascertainment are additional limitations.

## Computation, hostile review and reproducibility

- **Prefit chronology:** protocol separately committed before extracting source endpoint counts or computing predictive scores; subsequent narrow derived panel and original scripts are new files on the same branch. No existing file was changed.
- **Data / transcription:** reread immutable S1 through GitHub **plain-text contents API** (no browser download, download prompt or binary repository mirror); verified source Git blob SHA; independently checked species order, 55/11 gate, 3,533 denominators, 9 zero-L1 species, 45 integer candidate assignments and 31 strict-ratio mismatches. Identified and corrected a three-ambiguous-row omission in an intermediate totals calculation. Stored only derived panel and explicit primary/secondary gates.
- **Leakage:** for each held-out order and each jackknife deletion, derive means, SDs, N50 winsor cutoffs and Newton fit **only from training orders**. Tests examine predictions under held-out label modification. Whole-order deletions are full refits, not cached-prediction subtraction.
- **Numerical fit:** original stdlib stable softplus/sigmoid and Cholesky Newton fit, analytic gradient ≤`max(1e-8,1e-10*effective_training_denominator)`, ≤100 accepted steps and ≤45 line-search steps; finite loss/probability/clipping and train/test membership audited. Diagnostics are checked **against protocol bounds**, not exact iteration counts/floating-point traces. An independent SciPy optimizer for the nominal comparison (review only, **not CI or production**) reproduced the equal-order gain within **3×10⁻¹⁰**.
- **Test environment:** ten stdlib-only, offline tests passed in approximately **16 seconds** in local Python 3.13; designed for Python 3.12/3.14 and Windows. The archived numerics are compared with absolute/relative `3e-7` tolerances. Every newline-normalized file SHA256 is obtained from path bytes after CRLF/CR→LF, with path keys from `Path.as_posix()`. Tests deny socket access; no dependency other than stdlib. Passing tests establish computational reproducibility, not source truth, permission or biological validity.
- **Adverse interpretation:** 3 of 11 orders worsen directly, quality-restricted necropsy weighting reverses the gain, rodent order deletion materially shrinks it, some single-species orders carry extreme gains, secondary cannot be scored and the panel excludes all zero-recorded cases. These are reported, not hidden. The study has no somatic measurement, causal intervention, individual-animal independence, phylogenetic random-effects model, prospective test, independent population or veterinary/clinical validation.
- **No post-score feature search:** this is the precommitted `SINEs_Counts` extension only; other SINE proximity and composite L1/SINE annotations in the original paper were intentionally not explored.

**Reproduction:** from repository root, `python -m unittest discover -s tests -p "test_research_sine_l1_transport.py" -v`; then `python research/sine-l1-transport/analysis.py --output sine-l1-results.json`. This offline JSON includes all held-out predictions, complete jackknife pseudovalues, convergence diagnostics, 45 count assignments and every sensitivity. Input-hash keys are `analysis.py` and `data/species_aggregates.csv` (POSIX paths). **Neither passing tests nor green CI establishes the mechanistic interpretation.**
