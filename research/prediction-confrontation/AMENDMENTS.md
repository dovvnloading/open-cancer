# Transparent clarifications and exploratory additions

**Nonclinical, unreviewed. No treatment advice.** The frozen protocol remains
unchanged in commit `29bedbd`. These notes record decisions made during the
evidence audit on 2026-11-17 (America/New_York), not before seeing all sources.

## Clarifications that do not change an empirical result

* R1's strict positive quadratic coefficient additionally requires finite
  second clone-size moment and positive probability of clone size at least
  two, as in PR57 Theorem 3. A singleton clone distribution is a failure
  stratum, not a counterexample to the conditional theorem.
* E4's explicit generation-two dose curves require **one-copy founders**,
  as stated in PR60 Theorem 6, not arbitrary monomorphic founders. The
  protocol's abbreviated row omitted this condition. It is enforced in the
  eligibility audit; no qualifying experiment was found and no numerical
  result changes.
* A box of ordinary marginal 95% intervals is not a simultaneous 95% box.
  Any future formal multi-input interval test must use simultaneous bounds
  or justified Bonferroni construction before study-level Holm adjustment.
  No qualifying R1 experiment was found, so this clarification changes no
  tested result.
* NoDerivatives alone does not prohibit unchanged small factual excerpts or
  a change of file format. The [Creative Commons FAQ](https://creativecommons.org/faq/)
  explicitly discusses unchanged excerpts, format changes, and factual data.
  We retain original licenses and attribution, transcribe only isolated
  numerical facts, and create our own analysis. No source figure, expressive
  table, article or underlying database is redistributed or adapted. Sources
  without an established affirmative reuse basis remain metadata-only.

## E3-X: exploratory rejection-frequency calculation

During the E3 literature search, Lundberg et al. (2008),
[DOI 10.1371/journal.pone.0003099](https://doi.org/10.1371/journal.pone.0003099),
was identified as reporting percentages of mitoses departing from a fair
binomial allocation test, with numbers of divisions and the test level.
Those summary outcomes were visible before this operationalization was
written. **The following calculation is exploratory, not preregistered or
part of the frozen E3 test.** It must not be presented as a prospective
confirmation, nor as a test of the dose laws or clustered extension.

Before computing its numerical result, specify: retain every reported
untreated cell-line panel, convert a reported rounded percentage to the full
set of compatible integer numerators using the explicitly reported
denominator, and use the largest conservative upper-tail probability among
those numerators. Under independently sampled mitoses and valid per-mitosis
tests of level a, the number of rejections is stochastically bounded by
Binomial(n,a), even if discreteness makes some rejection probabilities less
than a. Compute `Pr(Binomial(n,a)>=observed)` exactly with rational arithmetic;
report the observed rejection fraction, excess above a, and Wilson interval.
Apply Holm across all untreated panels, and also show the more conservative
combined family containing the immune association tests. An effect depends
on independence/calibration assumptions; the original per-mitosis p-values
cannot be reconstructed without forbidden individual observations. If those
assumptions are not established, the result is a conditional sensitivity
analysis and the empirical verdict remains untestable.

Do not use the irradiated panels without an explicit denominator, convert
reported p-value bounds to exact p-values, or use chromosomal amplification
controls as ecDNA observations. An excess supports clustering or another
departure from independent allocation; it cannot identify the beta/bias
kernel or its dose threshold.

After the conditional calculation, add a dependence sensitivity bound:
`Pr(K>=k) <= min(1,n*a/k)` by Markov's inequality. This needs valid marginal
per-mitosis tests but not independent mitoses. It is an exploratory upper
bound, not an estimated p-value or proof of compatibility. Wilson intervals
for rejection fractions also require independent-binomial sampling. These
assumptions remain unverified by the published summaries.

## Implementation and robustness decisions after extraction

* I2 is reported as a clonal-burden association proxy and a separate I2-H
  heterogeneity component. Separate clonal/subclonal coefficients do not test
  heterogeneity at fixed burden without their joint design/covariance. This
  prevents assigning the available association's verdict to a missing test.
* The formal immune family has all nine Yang locus/cohort contrasts and five
  planned clonal/HLA estimates: Litchfield clonal, Ravi clonal, Litchfield
  HLA, and both Lozac'hmeur histologies. Total/subclonal TMB are retained as
  context, not tested for significance or compared against each other.
* Yang intervals use central exact conditional OR inversion, Wilson group
  intervals, Newcombe risk-difference intervals and log-Wald risk-ratio
  intervals. The frozen protocol did not specify an OR interval construction.
  All source panels are retained; no interval method was selected to change
  a verdict. Central OR intervals need not invert probability-ordered Fisher
  p-values. All displayed intervals are marginal, not simultaneous.
* Ravi's pre-specified log-normal CI reconstruction is approximate. As a
  conservative sensitivity, an ordinary BH-adjusted p bounds its raw p
  before rounding. The printed value is not itself a guaranteed bound:
  use its upper nearest-last-digit rounding endpoint, explicitly conditional
  on that rounding convention. This is not a new raw p-value or an exact
  family-error guarantee. It changes no verdict.
* Shibata sensitivity calculations show both exclusion and inclusion of the
  internally disputed clone. Neither resolves the source's classification,
  sister dependence or fitted-null calibration. Lange supplementary tests
  remain outside formal inference because of the source-version discrepancy.

## Publication transport record

The command-line checkout could read the public repository but had no GitHub
write credentials. Publication therefore uses the connected GitHub app's
Git-data API. That transport creates new commit objects at publication time;
their timestamps must not be presented as external preregistration or as
the time the local analysis plan was frozen.

The original local protocol-only commit was
`29bedbd6ca3df0ecdd1fdb62cd52a4b3cad880c2`, created at
`2026-11-17T21:24:50-05:00`. Its complete repository tree is
`e27fc8dac0395aeaa4c205c5cd06d0ee9812975d` and its immutable protocol blob is
`1588e0f197ac8b568e2d47610e7e89045185f633`.
The [remote protocol-only commit](https://github.com/dovvnloading/open-cancer/commit/679e1f69da1546c1eb0f662471b1a5f11b1b8e38)
has exactly that same tree. The subsequent
[remote initial-amendment commit](https://github.com/dovvnloading/open-cancer/commit/ed4463556268b744ffe3f96f04a0b4c617c2d743)
likewise has the exact tree of local `a349bd2`. The final results are
published after these two commits, with their tree checked against the
locally tested files. Matching hashes establish content identity, not
independent evidence for the claimed analysis chronology or biological truth.
