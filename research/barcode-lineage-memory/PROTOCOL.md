# Reconstructed analysis specification

Read [RECONSTRUCTION.md](RECONSTRUCTION.md) first. This documents a known-result
reproduction and preserves the primary-versus-exploratory distinction of the
earlier record; it is not a fresh preregistration. All results are nonclinical
and unreviewed. No biological p-value or read-level standard error is used.

## Data and observation unit

Use Whiting et al., Nature Communications 2025, DOI
[10.1038/s41467-025-59479-7](https://doi.org/10.1038/s41467-025-59479-7),
publisher Supplementary Data 2, sheets `HCTbc_Barcode_Data` and
`SW6bc_Barcode_Data`. Expected shapes: 13,962 and 688,439 barcode rows, each
with eight count columns ordered DT1_P1, DT1_P2, ..., DT4_P1, DT4_P2. Barcodes
are artificial lineage identifiers, not human genomes. Counts are amplified
sequencing reads, not independent cells. Four treated cultures per cell line
share one barcoded parent and are the reporting units.

P1 is post-treatment (source reports day 67 HCT and day 37 SW), not a parental
pretreatment baseline. The inspected processed release does not supply full
POT/untreated CO distributions. Their effects are not adjusted away. The
source has already filtered barcodes; without its complete filtering history,
P1-only support selection does not certify raw-data outcome independence.

## Primary support and probability model

At each abundance floor, set entries with `count / original_sample_total <
floor` to zero. If any required sample becomes empty, mark the scenario
unavailable. Keep every barcode present in any of the four thresholded P1
samples, and collapse all remaining P2 counts to one `OTHER_UNSEEN` category.
Retain that category even if its mass is zero. Normalize sample frequencies.

Smooth each early vector as `(1-delta)*p + delta/K`, primary delta `1e-4`,
sensitivity `1e-6` and `1e-3`. Targets are empirical late probability vectors.
For target r, own early vector x, and the arithmetic mean z of eligible
normalized/smoothed donor early vectors:

`q0_i(a) = x_i^a / sum(x^a)`;
`q1_i(a,b) = x_i^a*z_i^b / sum(x^a*z^b)`.

Tune `a` over `[0.5,1,1.5,2,3,4]`, `b` over `[0,0.25,0.5,1,2]`. Baseline and
augmented models tune separately. For each outer target, use the remaining
three cultures as inner validation targets; each inner donor set contains
the other two training cultures, never the outer target. Maximize the equally
weighted mean inner log score. Scores within `1e-12` of the maximum are tied;
choose smaller b, then a closest to one, then smaller a. Predict the outer P2
using its own P1 and the other three P1 donors. Never tune using outer P2.

Report `sum(p2*log(q1/q0))`, nats per read of the evaluated empirical
distribution, and equal-culture mean. Exponentiating gives a geometric
probability ratio, not accuracy or explained resistance. Report all four
folds, both selected exponents, baseline/full log scores, unsharpened own-only
and donor-only diagnostics, read totals, retained masses, and support size.

## Primary decision and declared sensitivities

The descriptive conjunction requires: mean gain >0.01; all four primary gains
positive; `1e-3` abundance-floor complete refit available with positive mean
and no negative fold; fixed-predictor 1% adverse target-deletion lower gain
positive in all folds. Unavailable or zero required-floor gain fails. This
is a descriptive robustness rule, not a significance or biological test.

Also compute complete-refit floors `1e-4` and `1e-3`, smoothing sensitivities,
and a highest-P2-category removal diagnostic. For the latter, exclude the
category, renormalize both predictions on the surviving support, and score
the renormalized target without refitting parameters.

For fixed per-category contrast d=`log(q1/q0)` and target probabilities p,
the sharp lower retained score after deleting at most epsilon mass removes
the highest d mass first, with fractional deletion at the boundary. Divide
the retained weighted sum by `1-epsilon`. Evaluate epsilon 0, 0.001, 0.01,
0.05 and 0.1; report first zero-crossing tipping mass if attainable below one.
This concerns outcomes only; predictors and parameter fits stay fixed. It
does not estimate an assay error rate or bound predictor contamination.

## Exploratory A2: isolate intervention location

After the primary result, the earlier investigation declared cutoffs `1e-5`,
`1e-4`, `1e-3` in four modes: early-only, late-only, both, donor-feature-only.
Keep the original primary P1 dictionary and its smoothing cardinality fixed
in every mode. Apply thresholds before renormalization, using each original
sample frequency. Refit the entire nested procedure. Empty required inputs
make the scenario unavailable, without dropping cultures.

Donor-feature-only changes only the donor view of early data. It preserves
own-P1, every P2 observation, the baseline fit/predictions, and the dictionary.
The donor view is separately thresholded, renormalized and smoothed and is
used in both inner and outer donor means. Record donor retained mass per
source culture. Removed read mass is a sensitivity witness, not measured
contamination, and selecting b=0 does not establish universal absence of
biological information or performance of other model classes.

## Exploratory A2: temporal profile diagnostic

For each home culture select original barcodes with home P1 frequency >=1e-3
and all other P1 frequencies <1e-4. Split them deterministically by parity of
the final byte of SHA256(barcode ASCII). Require at least ten selected
barcodes in each half; this does not guarantee occupied receiver support.
For each other culture's minor P1 profile, compare normalized home P1 with
`(1-w)*home_P1 + w*home_P2` profiles. Normalize home profiles separately
within each half and smooth with delta=1e-4 within that half. Tune w over
`0,0.05,...,1` on the opposite half, choosing smaller w on ties at 1e-12.
Evaluate held-out-half log-score gain relative to home P1. If any required
profile has no observations, mark that comparison unavailable.

Report every home/receiver/half, selected and occupied barcode counts,
receiver read support, fitted shape weights, gains, per-home and equal-home
means. Scores are nats per receiver read CONDITIONAL on the selected half;
they are not comparable directly with the primary full-distribution score.
These comparisons deliberately use later information and are not prospective
forecasts. Weights are profile-shape parameters, not contamination fractions.
Shared observations, sparse support and post hoc design preclude interpreting
the 24 comparisons as independent biological replicates or significance.

## Evidence and publication limits

Authenticate count and barcode bytes against pinned canonical hashes before
either implementation runs. Independently implement the probability-space
forecast and count-weighted score; verify important numerical and exact
mathematical results. Package only original aggregate calculations and
metadata, never licensed source tables or raw reads. Preserve contrary
evidence and all failed/unavailable scenarios. A first-header flowcell match
can support sequencing opportunity for spillover, not causation or its rate.
The published source already discusses residual index hopping and robust
diversity metrics; no rejection of its full mechanistic model follows.
