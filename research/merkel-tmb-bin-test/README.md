# Merkel TMB bins: an outcome-matched test of frozen response curves

**Nonclinical, unreviewed computational research. No treatment, diagnosis,
surveillance or patient-selection advice.** These are published aggregates.
Software checks, hashes and AI critique do not establish source truth,
external scientific review, biological validity or clinical usefulness.

## Verdict

**The primary frozen observation bridge survives with stated latitude under
the declared printed-interval compatibility rule.** It fails the identity
scale and the inherited narrower assumed band, while surviving in the upper
part of the wider assumed band. All required curve sensitivities have the
same qualitative domain verdict. The complete sharp sets, including their
open boundaries, are generated below.

This is a constrained survival result. It allows a common uncalibrated assay
multiplier, a different arbitrary TMB distribution in every cohort/bin, and
any mean within each broad printed marginal interval. The point-ORR check
is adverse: both reciprocal curves have an empty common point-compatible
region because of their finite-load response ceilings; both log curves need
point-compatible scales outside the assumed bands. Interval survival is
neither exact point agreement nor confirmation of the bridge.

The matched bin endpoints remove the particular median substitution and
broad-ORR/TMB-subset denominator mismatch in [PR #72](../merkel-viral-antigen-test/).
They do not establish the within-bin interpretation of an ecological curve,
correct sequencing-subset selection, calibrate the assay multiplier, or make
the treatment populations exchangeable. **This does not isolate viral
antigens or test [PR #59's functional-load theorem](../immune-escape-neoantigen-clonality/).**
The same JAVELIN cohorts are reused evidence. The earlier
[PR #67](../immune-escape-pancancer-test/) pan-cancer comparison is unchanged.

## Protocol chronology and scope

This implements [ranked job 3](../README.md) on main at
[`37c37e2`](https://github.com/dovvnloading/open-cancer/commit/37c37e2b4cfc711ce938201123ef3085b8868f2f).
The first local branch commit, `6d6d4088c11099c73fd23b79c24ffe0c5c781b79`,
contains **only [PROTOCOL.md](PROTOCOL.md)**. It precedes every new fit replay,
curve evaluation, scale calculation and score, and the protocol remains
unchanged. The original results and source outcomes were already visible;
this is a retrospective pre-computation freeze, not blinding or an
independent preregistration.

CLI push later failed because the workspace lacked GitHub credentials.
The authenticated GitHub connector published the **identical protocol-only
tree** as
[`892df395`](https://github.com/dovvnloading/open-cancer/commit/892df395d5b3fa0c293e76c03d8312f79832dcf3),
with new server commit metadata. That publication timestamp is after the
calculations; it is not presented as a second pre-computation freeze. The
original local commit remains preserved. For durable reconstruction, its
complete commit payload is also retained here, as UTF-8 with a trailing LF;
the tree is available in the published protocol commit. The test checks that
the Git commit-object hash of this payload is the original local commit ID.

<!-- BEGIN ORIGINAL PROTOCOL COMMIT -->
```text
tree 23572b1081642a8063d23a0261c0b150df5e2cb6
parent 37c37e2b4cfc711ce938201123ef3085b8868f2f
author Codex <codex@openai.com> 1794991072 -0500
committer Codex <codex@openai.com> 1794991072 -0500

Freeze outcome-matched Merkel TMB-bin protocol before computation
```
<!-- END ORIGINAL PROTOCOL COMMIT -->

The primary curve, all required sensitivities, scale domains, point and
rounding checks, endpoint rules and verdict rule were fixed in that commit.
The code reads the existing full-precision frozen curve CSV directly.
Historical fit replay is an identity check only. **No parameter is fitted to
either new bin outcome, including Part B.** Feasible scale regions are
reported in full rather than selecting a best scale. Mathematical witnesses
in JSON are existence certificates, not measured TMB values or fitted
population distributions.

Only this directory and
`tests/test_research_merkel_tmb_bin_test.py` are added. No prior file,
research index, workflow, dependency, source table or other research directory
is modified. This follows the user-authorized scope and
[DATA_POLICY.md](../../DATA_POLICY.md).

## Sources, denominators and HTML-only audit

The small [aggregate CSV](data/tmb_bins.csv) contains one row per cohort/bin.
Every row carries source, DOI, section, paragraph, denominator locator,
access date, actual inspection route, reuse terms and limitations. Access
date is **2026-11-18, America/New_York**. Printed decimal precision is
preserved. Transformations are categorical labels, factual numerical
transcription and percentage-to-probability conversion in code.

| Evidence | Precise article-text locator | Reuse basis |
| --- | --- | --- |
| [D'Angelo et al. 2020, Part A](https://pmc.ncbi.nlm.nih.gov/articles/PMC7239697/), DOI [10.1136/jitc-2020-000674](https://doi.org/10.1136/jitc-2020-000674) | Results, "Exploratory biomarker analyses", first prose paragraph: TMB-evaluable sentence and high-versus-low TMB ORR sentence | Article HTML copyright notice explicitly states [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Attributed factual aggregate transcription. |
| [D'Angelo et al. 2021, Part B](https://pmc.ncbi.nlm.nih.gov/articles/PMC8311489/), DOI [10.1136/jitc-2021-002646](https://doi.org/10.1136/jitc-2021-002646) | Results, "Exploratory genomic and transcriptomic analyses" > "Tumor mutational burden", first/sole paragraph: opening evaluability sentence and high-versus-low TMB ORR sentence | CC BY-NC/no commercial reuse notice independently corroborated in [PubMed HTML](https://pubmed.ncbi.nlm.nih.gov/34301810/). The [version-4.0 URI](https://creativecommons.org/licenses/by-nc/4.0/) is inherited from PR #72's audit; this HTML-only audit did not independently recover the versioned URI. |

Source licences remain attached to these noncommercial factual aggregate
transcriptions. Repository MIT does not relicense them. Source articles,
figures, supplements and participant data are not redistributed.

The numeric passages were inspected as **indexed primary PMC article HTML
text**, independently cross-checked, and linked to their canonical article
pages. An initial direct Part A HTML request also displayed its article and
copyright notice. Subsequent direct PMC requests returned browser-check
text; publisher HTML routes returned access errors. Ordinary alternate HTML
routes did not resolve the Part B versioned-licence limitation. Search
results also included unsolicited PDF links; those were not opened or used.
**No browser was used, no file-download prompt occurred, and no file download,
PDF reading, supplementary-file reading or figure digitization was performed.**
This is a text audit, not a visual audit of the original printed pages.

The selected paragraphs supply the total TMB-evaluable cohort denominator,
but do not print a separate denominator for either bin. Those fields and
responder counts remain blank. The total is repeated as cohort context,
never used as either bin denominator or summed across the two rows. Nothing
is inferred from rounded percentages or interval endpoints. The analysis
uses unweighted set feasibility and thus does not require bin-specific
counts. It neither pools cohorts nor constructs a binomial likelihood.

Part A concerns previously chemotherapy-treated metastatic disease; Part B
concerns first-line metastatic therapy. Both are from the same trial
programme and were already used in PR #72. The source authors' exploratory
cutoff is inherited, not newly validated. Part A's biomarker paragraph has
its own follow-up/cutoff; the broad response paragraph from a later follow-up
is not substituted. The ORRs here correspond to the reported TMB bins within
the measured subset. Missingness and generalization beyond that subset remain
unidentified. No viral-status-by-bin cross-tabulation is inferred.

## Observation bridge and why the bounds are sharp

Let T be finite and nonnegative in NSSV/Mb, and let c be one finite positive
assay multiplier shared across all bins and cohorts. The added bridge is

\[
\operatorname{ORR}_{j,b}=\mathbb E[f(cT)\mid j,b],\qquad
b\in\{[0,2),[2,\infty)\}.
\]

An ecological regression on cancer-type summaries does not justify this
within-bin response interpretation. It is a testable additional observation
assumption. The code does not substitute a median or assume that a nonlinear
curve commutes with averaging. Each cohort/bin can have its own distribution;
no observed within-bin distribution, mean, median or maximum is supplied as
an extra restriction in this frozen question.

Write `s=f(2c)`. For a reciprocal curve with positive A and B,

\[
f(x)=A\exp(-B/x),\quad f(0)=0,
\qquad
\mathcal M_{\mathrm{low}}=[0,s),\quad
\mathcal M_{\mathrm{high}}=[s,A).
\]

The low-bin upper endpoint is excluded because every T is strictly below
the threshold and f is strictly increasing. The high-bin ceiling is excluded
because f is strictly below A at every finite T. This remains true for an
unbounded distribution on finite T: the nonnegative gap from the ceiling is
positive almost surely, so its expectation is positive. Taking a limit at
infinite T does not produce an admissible distribution with mean A.

For clipped log curves, set `x0=exp(-a/b)` and `x1=exp((1-a)/b)`.
The high-bin range is `[s,1]`. The exact low-bin range is:

| Threshold position | Attainable low-bin means |
| --- | --- |
| `2c <= x0` | `{0}` |
| `x0 < 2c <= x1` | `[0,s)` |
| `2c > x1` | `[0,1]` |

At the upper clipping threshold itself, low-bin mean one remains unattained.
Above it, a finite T strictly within the low bin can already lie on the
plateau. The high bin always permits response one at a sufficiently large
finite T. These plateaus distinguish the log curves from a reciprocal
ceiling.

The bounds are **sharp**, rather than only necessary. Every permitted mean
can be realized by a single finite T in the correct bin, using the continuous
curve's inverse or a clipping plateau. Hence arbitrary mixtures cannot add
anything outside these ranges, while finite point masses attain everything
inside them. These are existence arguments, not reconstructed patients.

For a target interval `[L,U]` on the probability scale, the exact inverse
rules used by the code are:

| Curve/bin | Case | Feasible c |
| --- | --- | --- |
| Reciprocal low | `L=0` | `(0,infinity)` |
| Reciprocal low | `0<L<A` | `(B/(2 ln(A/L)),infinity)` |
| Reciprocal low | `L>=A` | empty |
| Reciprocal high | `L>=A` or `U=0` | empty |
| Reciprocal high | `L<A` and `U>=A` | `(0,infinity)` |
| Reciprocal high | `L<A` and `0<U<A` | `(0,B/(2 ln(A/U))]` |
| Clipped log low | `L=0` | `(0,infinity)` |
| Clipped log low | `0<L<=1` | `(exp((L-a)/b)/2,infinity)` |
| Clipped log high | `U=1` | `(0,infinity)` |
| Clipped log high | `0<=U<1` | `(0,exp((U-a)/b)/2]` |

Intersect these sets within cohort and across cohorts, retaining endpoint
flags. An open lower bound touching a closed upper bound gives an empty
intersection. Neither zero scale nor infinite scale belongs to the domain.
The analytic inversion covers the full domain; the fixed grid is only a
descriptive display and does not determine the result.

## Scoring and uncertainty

The primary rule asks whether every bin's attainable mean range intersects
its **unchanged printed marginal interval** at a common scale. The fixed-scale
score counts compatible bins; the domain score asks whether the common set
intersects the stated domain. The two inherited bands are assumptions, not
measured calibration intervals. In particular, their numerical direction
and width do not follow from a synonymous-count adjustment.

This intersection is **not a simultaneous confidence test, joint confidence
region, p-value or calibrated significance result**. We do not assign a
coverage probability to the common set or assume independence of the source
intervals. An empty set falsifies the specified bridge only under this
descriptive rule and its particular domain. Nonemptiness means not excluded
under that rule, not that the bridge is true.

The point-ORR sensitivity uses exact printed points separately. The rounding
sensitivity widens interval endpoints by half the last printed unit, without
changing the bin rule. Point-distance scores are infima over the permitted
mean ranges; a zero infimum at an excluded boundary is still incompatible.
No surrogate response or within-bin distribution is fitted. Selected feasible
means used for JSON witnesses certify existence only and cannot change the
computed regions.

<!-- BEGIN RECOMPUTED RESULTS -->

### Transcribed outcome-matched aggregates

All percentages and intervals below are printed source values. The n column is the **total TMB-evaluable cohort**, repeated for context; it is not either bin denominator. Bin denominators and responder counts are not printed in the selected prose and remain blank in the CSV.

| Part | TMB bin (NSSV/Mb) | ORR % | Printed marginal 95% interval % | TMB-evaluable cohort n |
| --- | --- | ---: | --- | ---: |
| A | <2 | 28.0 | [12.1, 49.4] | 36 |
| A | >=2 | 45.5 | [16.7, 76.6] | 36 |
| B | <2 | 41.2 | [24.6, 59.3] | 52 |
| B | >=2 | 50.0 | [26.0, 74.0] | 52 |

### Frozen curve identity and full response range

Coefficients are the unchanged full-precision PR #72 CSV values. Its original fit replay is checked with absolute coefficient tolerance 2e-6; optimizer diagnostics are bounded, not compared with exact saved outputs. These curves are never fitted to the new bins.

| Curve | Coefficients | Training rows / parents | Response range for finite nonnegative input |
| --- | --- | --- | --- |
| Reciprocal holdout (primary) | A=0.48718323459062135, B=5.5353933555249 | 20 / 19 | [0, 48.718323%) |
| Reciprocal full panel | A=0.44912947779826107, B=4.199233915543866 | 21 / 20 | [0, 44.912948%) |
| Clipped log holdout | a=-0.03576193945209477, b=0.12639822421480662 | 20 / 19 | [0%, 100%] |
| Clipped log full panel | a=-0.014484116854756118, b=0.1252672387866171 | 21 / 20 | [0%, 100%] |

### Primary printed-interval scale regions

All sets concern finite positive c. Parentheses are scientifically meaningful: a low-bin mean cannot attain the threshold response while the curve is increasing. Values are displayed to six decimals; endpoint flags apply to the unrounded analytic values. These sets do not have simultaneous confidence coverage.

| Curve | Part | Low-bin scales | High-bin scales | Both bins |
| --- | --- | --- | --- | --- |
| Reciprocal holdout (primary) | A | (1.987075, infinity) | (0.000000, infinity) | (1.987075, infinity) |
| Reciprocal holdout (primary) | B | (4.050433, infinity) | (0.000000, infinity) | (4.050433, infinity) |
| Reciprocal full panel | A | (1.600903, infinity) | (0.000000, infinity) | (1.600903, infinity) |
| Reciprocal full panel | B | (3.487854, infinity) | (0.000000, infinity) | (3.487854, infinity) |
| Clipped log holdout | A | (1.728191, infinity) | (0.000000, 284.290116] | (1.728191, 284.290116] |
| Clipped log holdout | B | (4.646031, infinity) | (0.000000, 231.434330] | (4.646031, 231.434330] |
| Clipped log full panel | A | (1.474639, infinity) | (0.000000, 254.018066] | (1.474639, 254.018066] |
| Clipped log full panel | B | (3.999943, infinity) | (0.000000, 206.406831] | (3.999943, 206.406831] |

| Curve | Both cohorts, all c>0 | In assumed [1,2] | In assumed [0.5,5] | Primary rule |
| --- | --- | --- | --- | --- |
| Reciprocal holdout (primary) | (4.050433, infinity) | empty | (4.050433, 5.000000] | survives with stated latitude |
| Reciprocal full panel | (3.487854, infinity) | empty | (3.487854, 5.000000] | survives with stated latitude |
| Clipped log holdout | (4.646031, 231.434330] | empty | (4.646031, 5.000000] | survives with stated latitude |
| Clipped log full panel | (3.999943, 206.406831] | empty | (3.999943, 5.000000] | survives with stated latitude |

### Exact printed point-ORR sensitivity

Points are treated as exact means only in this secondary diagnostic. A point failure is not a population rejection when the printed-interval rule survives.

| Curve | Part A | Part B | Both cohorts | In [1,2] | In [0.5,5] |
| --- | --- | --- | --- | --- | --- |
| Reciprocal holdout (primary) | (4.997189, 40.497217] | empty | empty | empty | empty |
| Reciprocal full panel | empty | empty | empty | empty | empty |
| Clipped log holdout | (6.079999, 24.276715] | (17.276082, 34.658233] | (17.276082, 24.276715] | empty | empty |
| Clipped log full panel | (5.247228, 21.215099] | (15.051039, 30.384881] | (15.051039, 21.215099] | empty | empty |

Per-bin sets (including ceiling failures):

| Curve | A low | A high | B low | B high |
| --- | --- | --- | --- | --- |
| Reciprocal holdout (primary) | (4.997189, infinity) | (0.000000, 40.497217] | (16.512033, infinity) | empty |
| Reciprocal full panel | (4.443431, infinity) | empty | (24.332702, infinity) | empty |
| Clipped log holdout | (6.079999, infinity) | (0.000000, 24.276715] | (17.276082, infinity) | (0.000000, 34.658233] |
| Clipped log full panel | (5.247228, infinity) | (0.000000, 21.215099] | (15.051039, infinity) | (0.000000, 30.384881] |

### Source-rounding sensitivity

Only interval endpoints are widened by half their last printed unit (0.05 percentage points here). Thresholds and coefficients are unchanged. This is not added sampling uncertainty.

| Curve | Part A | Part B | Both cohorts | In [1,2] | In [0.5,5] |
| --- | --- | --- | --- | --- | --- |
| Reciprocal holdout (primary) | (1.981185, infinity) | (4.038409, infinity) | (4.038409, infinity) | empty | (4.038409, 5.000000] |
| Reciprocal full panel | (1.595864, infinity) | (3.476105, infinity) | (3.476105, infinity) | empty | (3.476105, 5.000000] |
| Clipped log holdout | (1.721369, 285.416924] | (4.627689, 232.351639] | (4.627689, 232.351639] | empty | (4.627689, 5.000000] |
| Clipped log full panel | (1.468765, 255.033997] | (3.984009, 207.232343] | (3.984009, 207.232343] | empty | (3.984009, 5.000000] |

Per-bin sets (including ceiling failures):

| Curve | A low | A high | B low | B high |
| --- | --- | --- | --- | --- |
| Reciprocal holdout (primary) | (1.981185, infinity) | (0.000000, infinity) | (4.038409, infinity) | (0.000000, infinity) |
| Reciprocal full panel | (1.595864, infinity) | (0.000000, infinity) | (3.476105, infinity) | (0.000000, infinity) |
| Clipped log holdout | (1.721369, infinity) | (0.000000, 285.416924] | (4.627689, infinity) | (0.000000, 232.351639] |
| Clipped log full panel | (1.468765, infinity) | (0.000000, 255.033997] | (3.984009, infinity) | (0.000000, 207.232343] |

### Scale-independent reciprocal ceiling deficits

Positive deficits are infima at infinite load: the finite-T ceiling is excluded. Zero deficits here do not decide low-bin feasibility.

| Reciprocal curve | A high point above ceiling (pp) | B high point above ceiling (pp) |
| --- | ---: | ---: |
| Reciprocal holdout (primary) | 0.000000 | 1.281677 |
| Reciprocal full panel | 0.587052 | 5.087052 |

### Every inherited fixed scale: attainable means and feasibility

Ranges are percentages, identical across cohorts at a given scale because the cutoff is shared. Counts are numbers of compatible bins out of the four published bins. Each bin may use a different admissible distribution.

| Curve | c | Low mean range % | High mean range % | Interval-compatible bins | Point-compatible bins |
| --- | ---: | --- | --- | ---: | ---: |
| Reciprocal holdout (primary) | 0.5 | [0.000000, 0.192177) | [0.192177, 48.718323) | 2 | 1 |
| Reciprocal holdout (primary) | 1 | [0.000000, 3.059827) | [3.059827, 48.718323) | 2 | 1 |
| Reciprocal holdout (primary) | 4/3 | [0.000000, 6.112175) | [6.112175, 48.718323) | 2 | 1 |
| Reciprocal holdout (primary) | 2 | [0.000000, 12.209409) | [12.209409, 48.718323) | 3 | 1 |
| Reciprocal holdout (primary) | 5 | [0.000000, 28.008720) | [28.008720, 48.718323) | 4 | 2 |
| Reciprocal full panel | 0.5 | [0.000000, 0.674012) | [0.674012, 44.912948) | 2 | 0 |
| Reciprocal full panel | 1 | [0.000000, 5.501986) | [5.501986, 44.912948) | 2 | 0 |
| Reciprocal full panel | 4/3 | [0.000000, 9.299991) | [9.299991, 44.912948) | 2 | 0 |
| Reciprocal full panel | 2 | [0.000000, 15.719746) | [15.719746, 44.912948) | 3 | 0 |
| Reciprocal full panel | 5 | [0.000000, 29.512170) | [29.512170, 44.912948) | 4 | 1 |
| Clipped log holdout | 0.5 | [0.000000, 0.000000] | [0.000000, 100.000000] | 2 | 2 |
| Clipped log holdout | 1 | [0.000000, 5.185063) | [5.185063, 100.000000] | 2 | 2 |
| Clipped log holdout | 4/3 | [0.000000, 8.821314) | [8.821314, 100.000000] | 2 | 2 |
| Clipped log holdout | 2 | [0.000000, 13.946321) | [13.946321, 100.000000] | 3 | 2 |
| Clipped log holdout | 5 | [0.000000, 25.528073) | [25.528073, 100.000000] | 4 | 2 |
| Clipped log full panel | 0.5 | [0.000000, 0.000000] | [0.000000, 100.000000] | 2 | 2 |
| Clipped log full panel | 1 | [0.000000, 7.234452) | [7.234452, 100.000000] | 2 | 2 |
| Clipped log full panel | 4/3 | [0.000000, 10.838166) | [10.838166, 100.000000] | 2 | 2 |
| Clipped log full panel | 2 | [0.000000, 15.917315) | [15.917315, 100.000000] | 3 | 2 |
| Clipped log full panel | 5 | [0.000000, 27.395436) | [27.395436, 100.000000] | 4 | 2 |

### Point-distance infima at every fixed scale

Distances in percentage points grant arbitrary within-bin distributions. RMS is the root mean of the four squared distances. These are lower bounds on aggregate discrepancy, **not fitted prediction errors**. An excluded endpoint can have zero distance but fail exact compatibility, which is evaluated separately above.

| Curve | c | A low gap | A high gap | B low gap | B high gap | RMS gap |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Reciprocal holdout (primary) | 0.5 | 27.807823 | 0.000000 | 41.007823 | 1.281677 | 24.781844 |
| Reciprocal holdout (primary) | 1 | 24.940173 | 0.000000 | 38.140173 | 1.281677 | 22.794340 |
| Reciprocal holdout (primary) | 4/3 | 21.887825 | 0.000000 | 35.087825 | 1.281677 | 20.687406 |
| Reciprocal holdout (primary) | 2 | 15.790591 | 0.000000 | 28.990591 | 1.281677 | 16.518473 |
| Reciprocal holdout (primary) | 5 | 0.000000 | 0.000000 | 13.191280 | 1.281677 | 6.626699 |
| Reciprocal full panel | 0.5 | 27.325988 | 0.587052 | 40.525988 | 5.087052 | 24.572790 |
| Reciprocal full panel | 1 | 22.498014 | 0.587052 | 35.698014 | 5.087052 | 21.252832 |
| Reciprocal full panel | 4/3 | 18.700009 | 0.587052 | 31.900009 | 5.087052 | 18.664965 |
| Reciprocal full panel | 2 | 12.280254 | 0.587052 | 25.480254 | 5.087052 | 14.372462 |
| Reciprocal full panel | 5 | 0.000000 | 0.587052 | 11.687830 | 5.087052 | 6.380206 |
| Clipped log holdout | 0.5 | 28.000000 | 0.000000 | 41.200000 | 0.000000 | 24.907027 |
| Clipped log holdout | 1 | 22.814937 | 0.000000 | 36.014937 | 0.000000 | 21.316643 |
| Clipped log holdout | 4/3 | 19.178686 | 0.000000 | 32.378686 | 0.000000 | 18.816225 |
| Clipped log holdout | 2 | 14.053679 | 0.000000 | 27.253679 | 0.000000 | 15.331903 |
| Clipped log holdout | 5 | 2.471927 | 0.000000 | 15.671927 | 0.000000 | 7.932839 |
| Clipped log full panel | 0.5 | 28.000000 | 0.000000 | 41.200000 | 0.000000 | 24.907027 |
| Clipped log full panel | 1 | 20.765548 | 0.000000 | 33.965548 | 0.000000 | 19.905191 |
| Clipped log full panel | 4/3 | 17.161834 | 0.000000 | 30.361834 | 0.000000 | 17.438245 |
| Clipped log full panel | 2 | 12.082685 | 0.000000 | 25.282685 | 0.000000 | 14.010759 |
| Clipped log full panel | 5 | 0.604564 | 0.000000 | 13.804564 | 0.000000 | 6.908898 |

### Canonical input and implementation fingerprints

SHA-256 of UTF-8 after CRLF/CR normalization to LF. Every path key is built with `Path.as_posix()`. README is excluded from its own fingerprint. Original source files are read without modification; the protocol hash is pinned in code.

| Path | SHA-256 |
| --- | --- |
| `research/immune-escape-pancancer-test/analysis.py` | `06a7e15a24fa23b52f4c74b3659c678dd2778a88fdb429b0712e170e4b564b1c` |
| `research/immune-escape-pancancer-test/data/tmb_orr.csv` | `b34c1dd8d8d2c403f7251cf02611283d5326cbdc2caf2c72a6afadcf5ddc540b` |
| `research/merkel-tmb-bin-test/PROTOCOL.md` | `4789e3562692bf563d0eed4754e519b5f5ae45424092516e074117fd7a1ae1d6` |
| `research/merkel-tmb-bin-test/analysis.py` | `859cd4c9e536feafecbcd5d7e30fcac22237c400d99d8a359eb5412741aeafbc` |
| `research/merkel-tmb-bin-test/data/tmb_bins.csv` | `71c8fa0ecee9cbf4ea3ab30a53dd3307caad0fbf201e1d7e041a7c2deeb0adfb` |
| `research/merkel-viral-antigen-test/PROTOCOL.md` | `22923506eff606d81284a476a56813ec0773d73d238a5dcf633067344e711ea1` |
| `research/merkel-viral-antigen-test/analysis.py` | `39c6d7be8b81a02c8b6ee5c9ffacff37df30e9d705a6745039424af4cebcde31` |
| `research/merkel-viral-antigen-test/data/frozen_curves.csv` | `1acc4eff8810b2afd59c6e5cac250705f2e4e4f50b1c5f079a0c18308bb893f9` |
| `tests/test_research_merkel_tmb_bin_test.py` | `f991d40e41c41762f97e41ec8dc0ecbc5d4c9530ab778aedee14783425547964` |

<!-- END RECOMPUTED RESULTS -->

## What the result changes, and what it leaves unresolved

For the primary curve, the resulting common-scale set is narrower than
PR #72's reported median-envelope set. These use different aggregate
stratifications, so this compares the reported compatible sets rather than
establishing that one model is nested inside the other. Here the low-bin
intervals impose a positive lower scale boundary even after arbitrary
within-bin distributions are allowed. Passing the wider assumed band does
not mean passing every scale in it. That is a constraint on this observation
bridge, not an assay calibration.

The reciprocal high-bin intervals are broad enough to straddle each curve's
ceiling, so they place no additional scale restriction under the primary
interval rule. Their point estimates expose the ceiling restriction instead.
The log curves retain a finite upper scale boundary from the high-bin
intervals. All of these effects are retained in the tables rather than
choosing the most favorable curve or uncertainty convention.

There is no measured functional antigen load, presentation quality, founder
burden, immune killing rate, viral-antigen contrast or eradication endpoint
in the new CSV. Consequently, rejecting a scale scenario cannot identify the
missing biological input, and survival cannot identify a viral mechanism.
PR #59's conditional mathematical theorem remains outside the tested
observation model. Reusing the same trial cohorts also prevents claiming an
independent empirical replication of PR #72.

## Falsification attempts and hostile review

| Attempt | Outcome and boundary of the conclusion |
| --- | --- |
| Reuse the actual primary held-out curve rather than a favorable new fit | Frozen source hashes and an independent derivative-based training replay check the original coefficients. No JAVELIN outcome changes them. |
| Require one scale for both bins and both cohorts | A nonempty exact common set remains for every frozen curve. Part B supplies the binding lower boundary; for log curves it also supplies the binding upper boundary. |
| Retain the identity reference and both inherited bands | Identity and the narrower band fail. Only an upper portion of the wider band survives; no surviving scale is claimed to be calibrated. |
| Enforce the full response range, not only the threshold inequality | Both reciprocal curves fail exact common point matching through their ceilings. Their printed-interval survival remains a separate result. |
| Challenge strict low-bin and finite-T endpoints | Low-bin boundary contact is excluded; reciprocal ceilings cannot be reached with finite T. Log plateau attainment is checked separately on each side of both clipping thresholds. |
| Use point ORRs and source rounding instead of only broad intervals | Point sensitivity is adverse; widening only source rounding does not change the domain verdicts. Both remain displayed. |
| Replace one proxy assumption with a hidden distribution assumption | Rejected. Arbitrary within-bin distributions are explicit nuisance latitude. Sufficiency witnesses are mathematical examples, not recovered TMB distributions. |
| Reconstruct missing counts, pool totals or substitute viral-status ORRs | Rejected. Missing bin denominators/responders stay blank. The matched TMB-bin ORRs are the only new outcome data. |
| Treat the marginal-interval intersection as a joint significance test | Rejected. No joint coverage, independence, likelihood or p-value is asserted. |
| Turn transport or reused cohorts into causal or independent validation | Rejected. Treatment context, subset selection, assay uncertainty and reused evidence remain limitations. |
| Infer viral attribution or a clinical decision rule | Untestable here. No qualifying antigen measurement or causal comparison is supplied, and no treatment advice follows. |

Before computation, a separate mathematical review derived the open-endpoint
rules symbolically. After the protocol commit, an independent high-precision
decimal calculation checked all primary, point and rounding regions and the
fixed-grid diagnostics. A separate source audit checked every transcribed
number against the specified HTML passages. No numeric transcription
correction was needed. The licence-version limitation is retained explicitly
rather than promoted to a new direct verification. Implementation tests use
independent bisection, training fits, finite-distribution witnesses and
synthetic boundary counterexamples. These checks are internal AI/software
review, not external scientific validation.

No post-result empirical sensitivity or protocol amendment was added to
rescue the result. Numerical and wording corrections made during review are
confined to implementing the frozen rule. The protocol stays unchanged.

Hostile implementation review caught and fixed a floating-point witness
artifact at an excluded low-bin scale boundary: a forward log evaluation
could round above the target and suggest T exactly on the forbidden bin
edge. Witness generation now first requires membership in the analytic
scale set and then checks finite T and bin membership. A regression test
preserves the failure case. The scale sets, compatibility scores and source
values were already correct and did not change. Comparisons with PR #72's
median-envelope set were also qualified as different aggregate
stratifications, rather than a proof of model nesting.

## Honest novelty assessment

The trial authors already published the TMB-bin response analyses. Earlier
ecological TMB-response regressions and their limitations are credited in
[PR #67's source audit](../immune-escape-pancancer-test/TMB_ORR_SOURCES.md);
the curve and assay scenarios here are inherited from PR #72. The question
itself is the already-published repository roadmap job.

The contribution is a **repository-specific, auditable sharp-feasibility
calculation for a frozen observation bridge**, removing a specific median and
denominator substitution while preserving failures and uncertainty. The proof
uses elementary monotonicity, interval intersection and inverse functions.
No new biological mechanism, biomarker, general probability theorem, clinical
efficacy estimate or bibliographic priority is claimed. This targeted audit
is not a systematic literature review. A passing calculation does not reverse
PR #67's earlier predictive comparison.

## Reproduce and inspect

From the repository root, offline with Python's standard library:

```bash
python -B research/merkel-tmb-bin-test/analysis.py --check
python -B research/merkel-tmb-bin-test/analysis.py --json
python -B -m unittest discover -s tests -p 'test_research_merkel_tmb_bin_test.py' -v
```

The JSON output includes every per-bin set, cohort/common intersection,
assumed-band intersection, fixed-grid range, distance and finite mathematical
witness. All reported numbers are recomputed from the new aggregate CSV and
the unchanged curve/training CSVs. `--check` verifies the marked README block
without writing. `--write` updates only that block after an authorized change.
The analysis performs no numerical optimization; the independent historical
fit check compares coefficients within the inherited tolerance and checks
gradient/objective diagnostics against declared bounds rather than saved
optimizer values: absolute coefficient tolerance `2e-6`, absolute gradient
bound `1e-7`, and projected objective gap at most `1e-12`.

The single new test runs without network access, uses newline-normalized
hashes and POSIX path keys, and is designed for Python 3.12/3.14 on Windows
and Linux. It checks all source rows, missing counts, exact endpoint flags,
ceiling failures, log plateaus, printed/point/rounding results, fixed-grid
scores, witnesses, frozen inputs and the generated report. It does not need
git history in a shallow CI checkout. Path handling uses `Path` throughout;
network operations and nonstandard packages are absent from the analysis.

Repository integration commands are:

```bash
python -B scripts/validate.py
python -B -m unittest discover -s tests -v
node --test tests/*.test.*
```

The existing pull-request workflow discovers this test. CI establishes
software/bookkeeping checks, without scientific or clinical endorsement.
