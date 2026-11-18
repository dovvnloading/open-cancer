# Analysis protocol frozen before new evidence retrieval

**Nonclinical, unreviewed research. No treatment advice.**

Protocol date: 2026-11-17 (America/New_York). Base: main at
`864bd6aa09cbd2bd100a141750a99deb2688921c`. The first commit adding this
file precedes this work's literature retrieval and numerical extraction.
This is a prospective analysis plan for a retrospective, targeted evidence
audit, **not a registered or blinded study**. The merged notes already expose
some published results (including McGranahan, Ravi, Reinert and TRACERx).
We read those notes first as requested; we cannot claim ignorance of those
results. Future changes to a test must be identified as amendments, with
their timing relative to extraction. Do not silently replace this protocol.

## Eligibility and extraction

Use primary peer-reviewed papers or explicitly labelled preprints, their
permitted manuscripts and aggregate supplements. Extract small published
aggregate numbers only, with source URL, DOI, exact figure/table/text locator,
access date, endpoint/population, units, uncertainty type and reuse basis on
every CSV row. No patient records, cell/culture-level records, controlled
access, digitized individual points, reconstructed survival records, or
inferred denominators from rounded percentages. Missing is not zero.
An abstract can establish a reported aggregate but cannot supply missing
methods. Restrict redistribution to sources with an established reuse basis;
otherwise retain a bibliographic exclusion without copying numerical data.

Search the user-named sources and follow their relevant primary references,
including contrary/null evidence. Keep an inclusion/exclusion search ledger.
Do not select a cancer type, assay, endpoint, dose transformation or cutoff
because its effect agrees with a model. Do not pool overlapping cohorts or
different estimands. Published fitted models and simulations are not new
experimental observations. A statistic already used to construct a theorem
is not an independent replication of that theorem.

## Predictions and tests, specified before extraction

| ID | Prediction and assumptions | Required published aggregates and planned test |
| --- | --- | --- |
| I1 | PR59 Theorem 6: single-founder `k*p_escape -> sigma*h/kappa`; homogeneous founders `Q -> exp(-sigma*h*N/(kappa*k))`, heterogeneous founders depend on `sum(1/k_i)`; the half-k/2, half-3k/2 comparison has exponent 4/3. All rates fixed, source suppressed, minimum functional load increasing, no pre-existing null cells. | Matched groups with functional load distribution, effective founder count, independently measured rates and escape/extinction counts. For a response bridge, independently fixed `f,s` and justified finite-load errors, at least three disjoint finite positive dose intervals and response counts. Use the note's binomial acceptance-interval/intersection test for a common A at family alpha=.05; never substitute total TMB, an open-ended bin midpoint, an unreported per-SD scale, or survival HR for these inputs. Otherwise untestable. |
| I2 | Directional coverage proxy: higher clonal burden / lower heterogeneity predicts greater checkpoint benefit under matched background conditions. This is weaker than I1 and is already prior biological literature. | Extract all eligible reported clonal-response ORs/CIs and published contingency counts from McGranahan 2016, Litchfield 2021 and Ravi 2023. Recompute crude OR, risk difference, risk ratio, Wilson group intervals and Fisher's exact test when counts exist; retain adjusted published ORs separately. Compare to OR=1, never treat association as a test of the reciprocal coefficient. |
| I3 | Deliberately strong observation bridge: any HLA LOH predicts uniformly worse CPI outcome. PR59 itself warns against this bridge; this tests its tempting extrapolation, not its theorem. | Montesion and other eligible primary HLA-LOH studies: reported adjusted/unadjusted OR/HR, CI, endpoint, histology, TMB adjustment and population. Orient effects explicitly. An opposite-direction interval can contradict the blanket bridge; a null result does not prove absence. No LOH prevalence can estimate transition hazard h. |
| C1 | PR58 Theorem 6: on an exact d-spaced schedule with persistent detectability and fixed recurrence dates, parity-averaged alternate-test removal gives `Lthin=(L+max(0,L-d))/2` in expectation. Sharp mean/support bounds follow Eq25. | Published same-cohort original mean/support, d and both-parity-averaged, zero-imputed thinned mean with unchanged endpoints. Propagate reported decimal rounding outward; reject only if observed interval is outside the prediction envelope. Without thinned aggregate, calculate a labelled forecast and mark the prediction untestable, never confirmed. |
| C2 | PR58: fixed residual burden, perfect threshold, unimodal log-growth density and stated horizon conditions imply optimal geometric times with ratio `q=log(C/B*)/log(D/B*)`. | Calibrated threshold/residual quantities and aggregate capture counts under comparable geometric/equal-calendar schedules with the same number of tests and horizon. Compare capture risk difference and interval. Observational contrasts with changed assay/imaging or a different cancer do not test the optimum. No claim about an adaptive policy. |
| C3 | Scope stress test only: common detection window W plus exact periodic perfect testing and continuous recurrence observation implies lead-range width at most d (or d+I with bounded imaging lag I). | Published lead range, cadence and imaging cadence. Compare rounding-expanded range width to support bound only if all conditioning assumptions are documented. Heterogeneous rates and intermittent imaging are candidate minimal repairs. This conjunction is not PR58's general model. |
| E1 | PR55 generation-two Fano for monomorphic founders: `F2=(3+s)/(2*(1+s))`, independent of founder copy k and u, with `1<F2<1.5` for `0<s<1`. All cells including zeros, pooled by expected cell number, exact genealogical generation. | Published eligible means/variances or Fano with replicate uncertainty and independently measured s. Difference from predicted F2, propagate s interval. No substitution of cross-sectional calendar-time positive-only distributions. Without eligible moments/uncertainty, untestable. |
| E2 | PR55 threshold `lambda=2*u*(1-sqrt(1-s))=1`; at u=1, s*=3/4. PR60 linear-hazard binomial threshold `c*=0.378474994554...`, c=gamma/b. | Independent aggregate division/death-by-copy measurements, segregation law, founder state and finite-generation surviving-well counts. Test finite-generation survival against exact model probabilities with binomial intervals. A finite observation of persistence cannot falsify eventual extinction; drug concentration is not c. Missing calibration means untestable. |
| E3 | PR60 conditional segregation: unbiased binomial gives daughter mean k, variance k/2 and Fano 1/2. Its symmetric clustering/bias family gives `F(k)=(1-R)/2+R*k`, `0<=R<=1`. | Published conditional aggregate daughter moments by maternal k, or published fitted allocation p and intraclass correlation rho with uncertainty. Effective sister-pooled `R=(2p-1)^2+4*p*(1-p)*rho`. Test baseline p=1/2,rho=0 using reported intervals/tests; separately assess whether the broader mean/Fano constraint is testable. Rejecting baseline is not rejecting the clustered extension. |
| E4 | PR60 generation-two linear-dose laws: A=2+(3-R)c, D=2+4c+(1+R)c^2, B=3+R+4c; mean=A/D and Fano=(1+R)B/A+(1-R)/2-A/D. At fixed R, mean falls and Fano rises with c. Equal second-moment curves need not identify threshold. | Published generation-two aggregate moments with calibrated c,R and uncertainty. Compare residuals/direction. Distinguishing beta versus bias at R=1/4 requires conditional fourth moment (1/8 vs1/16), not Fano alone. Missing generation/rate/moment inputs means untestable. |
| R1 | PR57 independent switching, shared positive clone-size law, no holiday growth/death/new acquisition, unchanged readout: `mu_j=B+C*r^j`, C>0,0<r<1; `W_j=var_j-mu_j=A0+A1*r^j+A2*r^(2j)`, A1>=0,A2>0. | Independent-culture n, mean, variance and uncertainty at holidays 0,H,2H,3H. Estimate r=(mu1-mu2)/(mu0-mu1); predict held-out mu3=mu2-(mu1-mu2)^2/(mu0-mu1). Propagate published CI boxes conservatively; reject only if every compatible vector violates required constraints or held-out CI misses forecast. Check three-arm mean/Fano inverse admissibility, not only a fitted exponential. Without matched holiday arms and uncertainty, untestable. |
| R2 | Assumption stress test: independent cellular retention conditional on initial state predicts zero residual within-lineage retention covariance. | Published conditioned lineage-correlation coefficients/intervals or paired aggregate counts, with initial-state/exposure controls. Unconditioned sister similarity can reflect shared latent state and is not a falsification of PR57. Treat it as contextual evidence when conditioning is absent. |

## Error control, effect sizes and verdicts

The unit of evidence is the study/experimental panel, not every reported
p-value. No pooled meta-analysis is planned because cohorts may overlap and
estimands differ. Report every eligible contrast, including opposite and
imprecise findings. For eligible formal tests in this audit, collect all
two-sided p-values into one Holm family at .05. Keep source-reported p-values
distinct from recomputed ones. A p-value reconstructed from a rounded CI
under a log-normal approximation must be labelled approximate and cannot
upgrade an otherwise ineligible experiment. Exact binomial intersection
tests already protect their within-study bins; include any resulting study
p-value in the family if a calibrated one is available.

Report effect sizes and the source interval or reproducible interval method.
Ranges, IQRs, bootstrap CIs, SDs and SEs must never be interchanged. No interval
is invented where an aggregate omits required uncertainty. Rounding envelopes
are numerical uncertainty, not sampling confidence intervals.

* **Consistent:** an eligible observation agrees with the specified
  direction/envelope within available precision; this is compatibility only,
  not validation, causal identification or proof of a universal claim.
* **Inconsistent (falsified):** an eligible contrast violates the specified
  prediction after its stated uncertainty/error control. State exactly which
  model-plus-observation assumptions are rejected and the minimal repair.
* **Untestable with current public aggregates:** required measurements,
  design controls, bridge calibration or uncertainty are absent. A forecast
  generated from an input statistic is not an observed test result.

Non-significance is not equivalence. An assumption challenge must be labelled
separately from the conditional theorem. CI, unit tests, reproducibility,
and agreement among AI agents do not establish source truth or biology.

## Implementation and audit commitments

Use Python standard library only, no network at analysis/test time. One new
test file must recompute every reported statistic from the CSV inputs in
under 60 seconds, exercise missing-data/endpoint guards, and compare generated
result tables with README text. Source verification is a manual provenance
audit, not something unit tests can prove. Freeze an explicit search ledger,
check transcriptions against primary sources, check overlap and multiplicity,
and preserve contrary evidence before pushing. Only this directory and
`tests/test_research_prediction_confrontation.py` may be added.
