# Does the older-age breast-incidence pattern have a common period trend within each stage?

**Nonclinical, unreviewed computational research. No screening or medical
advice. No patient-level data, clinical prediction, or biological validation.**

**Verdict: the conjunction of four separate within-stage common-average
trend restrictions fails the frozen descriptive compatibility rule.
Localized and distant disease fail; regional and unknown stage pass.** Localized failure survives the declared
rounding stress and every single-age omission. Distant failure is a narrow
boundary result: the wider rounding stress makes its intervals touch, and
some age omissions restore compatibility. We therefore do not describe the
distant result as robust evidence of population-trend heterogeneity.

Passing means that a common average log-period change is compatible with
these marginal intervals under this rule. It does not establish equal true
slopes, equal annual trajectories, or full age-period factorisation. Failure
means an empty intersection of the selected published-summary intervals;
**interval intersection is a descriptive compatibility rule, not a
simultaneous significance test**. There is no claimed simultaneous confidence
level or calibrated rejection of a biological model.

This executes the highest-ranked question proposed in
[PR #84](https://github.com/dovvnloading/open-cancer/pull/84), which was still
open when this work began; the user's brief supplies the specification.
It follows [PR #80's cross-sectional stage analysis][3] and
[PR #73's incidence-strata audit][4]. It keeps passing stages and the unknown
category visible and does not modify the central research index.

## Frozen design and chronology

The branch starts from `main` at
[`2a6d357224a286a05dbae03d6d3ad88ed622d552`](https://github.com/dovvnloading/open-cancer/commit/2a6d357224a286a05dbae03d6d3ad88ed622d552).
The local checkout and GitHub branch metadata agreed on that base after main
advanced during initial reading. [PROTOCOL.md](PROTOCOL.md) was committed
alone on the GitHub branch as
[`04ad47e180bb6ee1444784f2e0f96e0ef449aafb`](https://github.com/dovvnloading/open-cancer/commit/04ad47e180bb6ee1444784f2e0f96e0ef449aafb)
**before any transformed-trend, interval-intersection, allowance, or
sensitivity computation**. The source licence, table cells, and method design
were read beforehand. The committed protocol was then fetched into the local
working branch before the CSV and analysis were created.

There is **no protocol amendment or analytical deviation**. The wider
rounding check, point-rounding check, every age omission, and all passing
stages were fixed before computation. Published patterns were already known,
so this is a committed computational specification, not blinding, independent
preregistration, independent data replication, or an independent timestamp
attestation. The preserved Git history documents the work sequence.

## Source, licence, and what was actually read

The only numerical input is
[data/lee_argov_table2.csv](data/lee_argov_table2.csv), transcribed from
**Lee Argov EJ, Lui ML, Karr AG, Tehranifar P, Kehm RD (2025),
*Breast Cancer Incidence Trends in Older US Women by Race, Ethnicity,
Geography, and Stage*, JAMA Network Open,
[doi:10.1001/jamanetworkopen.2025.16947](https://doi.org/10.1001/jamanetworkopen.2025.16947)** [1].
The selected source is the **main article's HTML Table 2, Stage at diagnosis**,
with Localized, Regional, Distant, and Unknown rows and the age columns
65-74, 75-84, and 85+. The overlapping All ages column and the In situ row
are outside the requested panel. They were not omitted because of a result.

**Licence checked before transcription:** the article's Article Information /
Open Access notice explicitly states **CC-BY**. Its linked
[JAMA CC-BY permissions page][2] permits reuse with author and journal credit.
Neither readable page specifies a licence version; the record therefore says
**CC-BY, version not stated**. Attribution is retained on every CSV row.
The CSV reformats selected facts, normalizes minus signs and age labels, and
adds provenance. It does not reproduce article prose, figures, or layout.
The repository's MIT licence does not relicense third-party material.
This follows [DATA_POLICY.md](../../DATA_POLICY.md).

**Access: 2026-11-18, America/New_York.** The main article and linked
permissions page were read as HTML text. No browser file download was
initiated, no download prompt appeared, and no PDF, supplementary file, or
raw registry data were acquired. The text reader exposed the same selected
Table 2 cells in two places; both representations agreed. These are duplicate
representations of one source, not independent datasets. A separate source
audit checked every printed numeric string, sign, trailing zero, row, and
age header against the HTML before publication.

Every CSV row carries the source citation, DOI, URL, table/section/row/column,
human-readable locator, period, measure, units, marginal-interval label,
source standardization and database-vintage labels, access date/timezone,
and reuse basis/URL. There are no derived numerical columns. The test pins
the independently audited transcription; it cannot verify the registry
behind the article.

### Definitions and data limits

Table 1 provides period-pooled age-adjusted rate and rate-ratio context,
not annual observations. Its rates are not used to manufacture a time series.
The article Methods describes female USCS incidence, standardized to the
2000 US standard population using the stated Census age-group scheme, and
analyses using the 2023 NPCR/SEER submission [1, Methods]. This standardization
label is preserved, without claiming to reproduce the unpublished SEER*Stat
session.

Recorded stage uses a combined summary-stage variable, with all forms of
regional disease combined. Unknown includes unknown, unstaged, unspecified,
and death-certificate-only cases. The HTML text reader omitted the table
footnote prose; these definitions were checked in the main Methods, not
claimed as an observed footnote transcription. The paper names different
staging systems across calendar periods but does not provide the complete
reclassification histories needed to harmonize diagnoses ourselves [1,
Methods]. Unknown is retained as an observation category, not interpreted
as a separate biological stage or redistributed into known stages.

The selected data are published full-period joinpoint AAPCs with marginal
intervals. Annual counts, annual denominators, fitted segment covariance,
cross-age/stage covariance, and original model-selection sessions are not
available in the permitted inputs. We do not reconstruct them or infer them
from interval widths. No annual counts, joint covariance, 85+ midpoints,
or stage-reclassification histories may be invented.

## Restriction and allowance

For a fixed recorded stage, the factorised model is

```text
I_s(age, period) = A_s(age) * P_s(period).
```

It allows a different period factor in every stage. With fixed nonnegative
age weights within a bin, aggregation preserves the period factor, so the
bin rate also has the form `J_sg(t) = B_sg * P_s(t)`. All age groups then
share the same average log-period change over the same calendar interval.
We use the specified operational transform

```text
beta_sg = log(1 + AAPC_sg / 100).
```

The transform is applied to the published estimates, which are not exact
population slopes. The source describes joinpoint fits on logarithmic rates;
this audit does not refit them or equate a summary with a constant slope in
every year [1, Statistical Analysis]. The bin-weight and observation bridge
is an explicit condition, not an independently verified demographic model.

Equality of full-period averages is **necessary, not sufficient** for the
factorised surface. For example, an age-specific log-rate perturbation that
is zero at both period endpoints but varies between them preserves the
average endpoint change while breaking common annual trajectories. This is
a mathematical counterexample, not invented annual cancer data. A compatible
average therefore cannot verify transport of the full age pattern over time.

For each printed decimal, the primary rounding allowance is half its last
displayed unit. Widen each printed marginal interval endpoint outward **on
the AAPC percentage-point scale**, then apply the increasing logarithm.
Also retain the analogous rounding interval around each printed estimate.
The last-digit convention is an assumption about printing, not a sampling
distribution, validated publisher rounding algorithm, or added standard error.

Let `C_sg = [l_sg, u_sg]` be the transformed widened interval. The intersection
has limiting bounds

```text
L_s = max_g l_sg,  U_s = min_g u_sg.
```

Compatibility passes exactly when `L_s <= U_s`; equality passes. The code
makes this decision using the exact decimal AAPC endpoints before applying
the monotone transform, so a rounded display cannot change the verdict.
Natural-log values are evaluated in floating point and independently checked
with high-precision decimal logarithms.

The reported effect size is the smallest **symmetric allowance per age**
added to the log-trend intervals:

```text
epsilon_s = inf_b max_g distance(b, C_sg)
          = max(0, (L_s - U_s) / 2).
```

After widening by epsilon, compatibility requires and is equivalent to
`L_s - epsilon <= U_s + epsilon`. This proves sharpness and attainment.
The midpoint of the limiting log bounds supplies an attaining reference;
for a pass it is just one compatible value. It is not an estimate of a
biological common slope. The allowance is **half the positive gap**, since
both opposing endpoints can move. A common translation of every interval
cannot remove separation. These are established interval-geometry facts,
not new theorems.

The units are natural-log change per calendar year. The allowance is not an
AAPC percentage-point difference, an annual incidence-rate difference, a
confidence bound on true heterogeneity, or an estimate of screening effects.
The overall rule is the conjunction of the stage-specific rules, with
separate reference slopes by stage. Its common allowance is the largest
stage-specific allowance.

<!-- BEGIN RECOMPUTED RESULTS -->

## Recomputed results

Panel: **12 stage-age cells**, with all selected stages and ages retained.
**Overall descriptive verdict: FAIL.** Failing stages: **localized, distant**. Passing stages: **regional, unknown**.

The smallest allowance shared across stages (with a separate reference trend per stage) is **0.005638219609 natural-log units per calendar year**.

### Printed inputs and primary transformed intervals

AAPCs and printed marginal intervals are percentages per year from Table 2 [1]. The log estimate and widened intervals are derived. Rounding is applied to percent endpoints before the logarithm. The point-rounding intervals, which omit sampling uncertainty, are also retained in the JSON output.

| Stage | Age | Printed AAPC | Printed marginal interval | Log estimate | Widened log interval |
|---|---|---:|---|---:|---|
| localized | 65-74 | 0.95 | [0.78, 1.18] | 0.009455158771 | [0.007720123015, 0.011780338536] |
| localized | 75-84 | 0.38 | [0.20, 0.74] | 0.003792798239 | [0.001948101218, 0.007422385816] |
| localized | 85+ | -0.66 | [-0.94, -0.36] | -0.006621876309 | [-0.009494934562, -0.003556316203] |
| regional | 65-74 | -0.80 | [-1.07, -0.48] | -0.008032171697 | [-0.010808198717, -0.004761317102] |
| regional | 75-84 | -0.92 | [-1.23, -0.60] | -0.009242581367 | [-0.012426895008, -0.005967771780] |
| regional | 85+ | -0.89 | [-1.23, -0.53] | -0.008939841569 | [-0.012426895008, -0.005263829675] |
| distant | 65-74 | 1.13 | [0.92, 1.34] | 0.011236631926 | [0.009108392364, 0.013360351702] |
| distant | 75-84 | 1.64 | [1.33, 1.96] | 0.016266972464 | [0.013162986526, 0.019459431156] |
| distant | 85+ | 1.64 | [1.36, 1.95] | 0.016266972464 | [0.013459019684, 0.019361353479] |
| unknown | 65-74 | -3.99 | [-4.61, -3.34] | -0.040717833279 | [-0.047248852601, -0.033918793204] |
| unknown | 75-84 | -3.64 | [-4.37, -3.07] | -0.037079008241 | [-0.044735893841, -0.031129535197] |
| unknown | 85+ | -3.80 | [-4.32, -3.11] | -0.038740828316 | [-0.044213154676, -0.031542268010] |

### Primary compatibility and minimum allowance

L is the greatest lower log bound; U is the least upper log bound. For PASS, [L,U] is the common-average compatibility set. For FAIL, L and U are separating bounds, not a reversed confidence interval. All log quantities have units per calendar year. Decisions use exact decimal percent endpoints, not the rounded display.

| Stage | Verdict | L (attaining ages) | U (attaining ages) | Full positive gap | Symmetric allowance | Attaining reference |
|---|---|---|---|---:|---:|---:|
| localized | **FAIL** | 0.007720123015 (65-74) | -0.003556316203 (85+) | 0.011276439218 | 0.005638219609 | 0.002081903406 |
| regional | **PASS** | -0.010808198717 (65-74) | -0.005967771780 (75-84) | 0.000000000000 | 0.000000000000 | -0.008387985248 |
| distant | **FAIL** | 0.013459019684 (85+) | 0.013360351702 (65-74) | 0.000098667982 | 0.000049333991 | 0.013409685693 |
| unknown | **PASS** | -0.044213154676 (85+) | -0.033918793204 (65-74) | 0.000000000000 | 0.000000000000 | -0.039065973940 |

The exact percent-scale endpoints below determine the same verdict because the logarithm is strictly increasing. They are the limiting widened source endpoints, not new uncertainty estimates.

| Stage | Greatest lower endpoint (%) | Least upper endpoint (%) |
|---|---:|---:|
| localized | 0.775 | -0.355 |
| regional | -1.075 | -0.595 |
| distant | 1.355 | 1.345 |
| unknown | -4.325 | -3.335 |

### Frozen rounding and uncertainty challenges

Each cell shows verdict / symmetric log allowance. Printed CI uses no rounding widening; full-unit CI uses twice the primary half-unit; point rounding omits the published sampling uncertainty. These sensitivities do not replace the primary result.

| Stage | Printed CI | Primary half-unit CI | Full-unit CI stress | Point rounding only |
|---|---|---|---|---|
| localized | FAIL / 0.005688116429 | FAIL / 0.005638219609 | FAIL / 0.005588322817 | FAIL / 0.007988586729 |
| regional | PASS / 0.000000000000 | PASS / 0.000000000000 | PASS / 0.000000000000 | FAIL / 0.000554771088 |
| distant | FAIL / 0.000098667983 | FAIL / 0.000049333991 | PASS / 0.000000000000 | FAIL / 0.002465853003 |
| unknown | PASS / 0.000000000000 | PASS / 0.000000000000 | PASS / 0.000000000000 | FAIL / 0.001767429194 |

### Every leave-one-age-out diagnostic

The omitted age is removed only in this diagnostic. Each row uses the primary rounding rule on the remaining pair. Passing rows remain visible.

| Stage | Omitted age | Remaining ages | Verdict | L | U | Symmetric allowance |
|---|---|---|---|---:|---:|---:|
| localized | 65-74 | 75-84, 85+ | FAIL | 0.001948101218 | -0.003556316203 | 0.002752208710 |
| localized | 75-84 | 65-74, 85+ | FAIL | 0.007720123015 | -0.003556316203 | 0.005638219609 |
| localized | 85+ | 65-74, 75-84 | FAIL | 0.007720123015 | 0.007422385816 | 0.000148868600 |
| regional | 65-74 | 75-84, 85+ | PASS | -0.012426895008 | -0.005967771780 | 0.000000000000 |
| regional | 75-84 | 65-74, 85+ | PASS | -0.010808198717 | -0.005263829675 | 0.000000000000 |
| regional | 85+ | 65-74, 75-84 | PASS | -0.010808198717 | -0.005967771780 | 0.000000000000 |
| distant | 65-74 | 75-84, 85+ | PASS | 0.013459019684 | 0.019361353479 | 0.000000000000 |
| distant | 75-84 | 65-74, 85+ | FAIL | 0.013459019684 | 0.013360351702 | 0.000049333991 |
| distant | 85+ | 65-74, 75-84 | PASS | 0.013162986526 | 0.013360351702 | 0.000000000000 |
| unknown | 65-74 | 75-84, 85+ | PASS | -0.044213154676 | -0.031542268010 | 0.000000000000 |
| unknown | 75-84 | 65-74, 85+ | PASS | -0.044213154676 | -0.033918793204 | 0.000000000000 |
| unknown | 85+ | 65-74, 75-84 | PASS | -0.044735893841 | -0.033918793204 | 0.000000000000 |

### Reproduction hashes

SHA-256 of UTF-8 content after CRLF/CR-to-LF normalization, with repository-relative POSIX path keys. Hashes document exact inputs/code; they do not prove source truth or scientific validity.

| Path | SHA-256 |
|---|---|
| `research/breast-age-period-stage/PROTOCOL.md` | `d8ad921768f0e387687f17b9a82faf0ba22e0789b9adf833e4584019ddd58d40` |
| `research/breast-age-period-stage/data/lee_argov_table2.csv` | `856f7b1ab0c7bab4a47c005247efb70abc75dacf430380eba0868ce6c58a5f47` |
| `research/breast-age-period-stage/analysis.py` | `90341a6407ff49944b8e2e442cb9bc866260e6a0b8b4c9cfa42deef844b4ebda` |

<!-- END RECOMPUTED RESULTS -->

## Explicit stage verdicts and interpretation

**Localized: FAIL.** The youngest selected age group's lower trend bound
exceeds the oldest group's upper bound. The failure survives both rounding
sensitivities and every single-age omission. Thus the observed localized
age pattern cannot be described by one common average period multiplier
under this published-interval rule. This does not identify why the
stage-by-age summaries differ or establish a statistically calibrated
interaction.

**Regional: PASS.** All selected age groups have a nonempty common-average
interval, with zero additional allowance. The primary pass survives all
declared CI-endpoint rounding and age-omission checks. Point-rounding intervals alone fail,
so equality of the printed point estimates is not being claimed. Discarding
the published uncertainty would produce a different and unjustified primary
conclusion.

**Distant: FAIL under the primary rule, with substantial qualification.** The
oldest group's lower endpoint is just beyond the youngest group's upper
endpoint after the frozen half-unit widening. The full-unit stress restores
compatibility at a touching boundary. Omitting either member of that
constraining age pair also restores compatibility; omitting the middle group
does not. This is a fragile incompatibility of the selected printed
summaries, not robust evidence for a general population-age interaction.
The article's already published increase in distant disease across age
groups is a different claim from equality of its cross-age average trends.

**Unknown: PASS.** The interval intersection is nonempty and the allowance
is zero. The primary pass survives all declared CI-endpoint rounding and
age-omission checks,
while the point-rounding-only restriction fails. This retained pass concerns
recorded unknown-stage incidence. It does not establish stable classification
or independence from the trends in named stages.

**Overall: FAIL the frozen descriptive conjunction.** Localized alone is
sufficient for this result and its failure persists in the wider-rounding
stress. The overall outcome therefore does not depend on promoting the
fragile distant result. Regional and unknown provide explicit compatible
counterexamples to a statement that every stage needs a nonzero allowance.
No stage is removed from the primary report because it passes.

## Falsification attempts and hostile review

| Challenge | Outcome and limit |
|---|---|
| Misread a sign, age column, or pooled column | Every selected printed estimate and endpoint was checked against the main HTML by a separate source audit and again in the transcription invariants. No discrepancy was found. This does not audit registry records. |
| Use inaccessible or unlicensed material | The article-specific CC-BY notice and its linked reuse page were read before transcription. Only HTML facts are used; no file download or prompt occurred. The unspecified licence version remains unspecified. |
| A rounding decision creates the primary failure | Printed-CI and half-unit results agree on the set of failing stages. Full-unit widening retains localized failure but restores distant compatibility at equality. The adverse distant sensitivity remains prominent. |
| The oldest open-ended group is solely responsible | Every age omission is retained. Localized still fails after omitting the oldest group. Distant becomes compatible when either constraining age is omitted. No age midpoint is used. These are diagnostics, not revised primary populations. |
| Passing intervals conceal unequal printed estimates | Point-rounding-only checks fail for all stages. This illustrates how much the result depends on including source uncertainty; it does not override the primary regional and unknown passes. |
| The effect size mistakenly reports the full gap | The independent arithmetic checks both the half-gap minimum and an attaining reference, then demonstrates incompatibility below each positive minimum. Common-shift and equality fixtures check the convention. |
| Compatible averages prove full factorisation | The endpoint-preserving log-rate perturbation above breaks annual factorisation while preserving the full-period average. Passing is deliberately weaker than validating the surface. |
| The result is a calibrated cross-age significance test | The source intervals are marginal and the joint covariance is unavailable. No p-values, joint confidence region, multiple-comparison correction, or claimed simultaneous error rate are manufactured. |
| Unknown stage can be reassigned to explain known-stage trends | No joint classification history or annual crosswalk is provided. The recorded unknown trend is retained, but neither cancellation nor a reclassification mechanism is identified. |
| The analysis independently replicates the age decline or identifies screening | The source populations overlap and the estimands/age bins differ. Screening history, cohorts, and stage-reclassification histories are unavailable. Those claims are explicitly withheld. |

The source audit, independent arithmetic implementation, and hostile methods
review were separate AI investigations using the same evidence. Their
agreement is not independent biological replication, external human peer
review, or clinical validation. Frozen hashes and software tests document
reproduction, not truth of the source measurements.

## Relationship to the earlier age pattern and honest novelty

PR #80 used a single calendar year's USCS incidence counts and population
denominators in finite, narrower age bins, with the open-ended oldest group
descriptive. Its **cross-sectional age slopes** compare different ages at
the same period and its stage decomposition locates which recorded stages
carry a net decline [3]. Here **period slopes** compare calendar changes
within each broader age group. A stage can decline across ages in a
cross-section while increasing over calendar time at each age. Neither
direction logically determines the other.

The present use of the reported 85+ period slope requires neither a midpoint
nor a finite width for that age bin. It does not retroactively supply an age
slope for the open-ended bin in PR #80. The conditional count model and
simultaneous inference of that earlier note are not inherited here. PR #73's
rounded-mean shape restriction is also a different estimand; its contribution
to this design is transparent restriction checking and retention of passing
strata, not a transferable statistical calibration [4].

**The published trends are prior findings.** Lee Argov et al. already describe
age-varying localized trends, regional and unknown declines, and distant
increases [1, Abstract and Results]. The present contribution is an
incremental, reproducible compatibility audit with a committed rounding rule,
sharp symmetric allowance, exact witnesses, and disclosed fragility. It is
not discovery of those epidemiological directions or a new mathematical
method. The novelty comparison is against the identified source and prior
repository notes; no systematic literature review or bibliographic priority
claim was attempted.

Both this source and PR #80 draw on USCS populations, with overlapping
calendar coverage and different age bins; earlier SEER-based notes also
overlap the USCS system. Broader time coverage does not make this independent
replication [1; 3; 4]. Different age bins and standardization further prevent
literal transport of PR #80's finite-bin contrast. The new summaries cannot
separate age, calendar-period, and birth-cohort effects, nor supply a screening
intervention comparison. The original article itself states that its
screening hypothesis could not be tested with the available USCS data and
notes changing unknown-classification proportions [1, Discussion and
Limitations].

Thus no **age, cohort, or screening attribution** is made. Detection,
selection, staging, competing risks, and biology remain possible explanations
requiring other observations; none is identified by the allowance. Stage at
diagnosis is not detection mode, and recorded incidence is not a directly
observed latent-onset hazard. No screening, diagnosis, surveillance,
treatment, or other medical advice follows.

## Reproduce and verify

From the repository root, using only Python's standard library:

```bash
python -B research/breast-age-period-stage/analysis.py --check
python -B research/breast-age-period-stage/analysis.py --format json
python -B -m unittest discover -s tests -p 'test_research_breast_age_period_stage.py' -v
```

[analysis.py](analysis.py) derives the entire numerical block from the CSV.
`--write-readme` deliberately regenerates that block. JSON exposes every
intermediate log interval and every sensitivity's limiting ages and minimum
allowance, including compatible cases. No network or downloaded source file
is needed. An incomplete or invalid panel stops with an explicit error;
it is never imputed or silently reduced.

The [dedicated test](../../tests/test_research_breast_age_period_stage.py)
independently parses the CSV and recomputes the transformed bounds and every
primary, sensitivity, omission, and overall result. It checks the complete
generated report, frozen protocol/data hashes, source/panel invariants,
attainment/minimality, equality/ties, domain/malformed inputs, and path/newline
portability. All hash path keys are constructed with `Path.as_posix()`;
content hashes normalize CRLF and CR to LF. Compatibility targets are Python
3.12 and 3.14 on Windows and Linux. Simulated Windows path/newline checks
are not a claim that a Windows runner was executed.

The existing PR workflow discovers the new test and runs repository
validation. No workflow or existing file is changed. The PR's live checks
are the authority for CI status; passing CI establishes neither independent
scientific replication nor clinical validity.

## Sources

1. Lee Argov et al. (2025), [main article HTML][1]. Table 2 is the sole
   numerical source; Table 1 and Methods supply rate/definition context;
   Results, Discussion, and Limitations constrain novelty and interpretation.
2. JAMA Network, [CC-BY License Permissions][2]. Article-linked reuse basis;
   version unspecified in the readable notice and permissions page.
3. Open Cancer, [PR #80 note][3], *Which older-age breast-incidence declines
   survive uncertainty and stage accounting?* Its source, Turbow et al.,
   supplies the earlier USCS cross-sectional analysis, not input data here.
4. Open Cancer, [PR #73 note][4], *Incidence shape within published demographic
   strata*. Background restriction-checking discipline; no data imported.

[1]: https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2835519
[2]: https://jamanetwork.com/pages/cc-by-license-permissions
[3]: ../breast-incidence-decline-stage/README.md
[4]: ../incidence-strata-test/README.md
