# Transparent clarifications and exploratory additions

**Nonclinical, unreviewed. No treatment advice.** The frozen protocol remains
unchanged in commit `29bedbd`. These notes record decisions made during the
evidence audit on 2026-11-17 (America/New_York), not before seeing all sources.

## Clarifications that do not change an empirical result

* R1's strict positive quadratic coefficient additionally requires finite
  second clone-size moment and positive probability of clone size at least
  two, as in PR57 Theorem 3. A singleton clone distribution is a failure
  stratum, not a counterexample to the conditional theorem.
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
