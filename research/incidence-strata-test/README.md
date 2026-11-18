# Incidence shape within published demographic strata

**Nonclinical, unreviewed computational research. No medical advice.**
This is an aggregate-data shape audit, not biological validation, a clinical
prediction tool, or an estimate of cancer drivers. It executes ranked job 4
in [the research programme](../README.md), extending
[PR #64's incidence analysis](../multistage-incidence-identifiability/).

The result below is conditional on treating published rounding intervals as
sets of possible mean rates. Sampling fluctuation is not represented by those
intervals. A failure of that deterministic condition is not a statistical
rejection of the true population hazard or of multistage carcinogenesis.

<!-- BEGIN RECOMPUTED RESULTS -->

## Recomputed results

**Verdict: shape failure persists within strata.**

The CSVs contain **330 published rates** in **30 prespecified curves**; **30** are testable. **6** fail the rounded-interval monotone-mean condition.

Breast curves failing the condition: **All Races, Whites, Blacks**.

The failure persists within **both** printed race-specific breast strata. Pooling across those race categories is therefore insufficient, by itself, to explain away this descriptive failure.

### Primary and contrast monotonicity frontiers

`rho` is the minimum for any nondecreasing positive mean vector. The exact rational `rho^2` is included so a rounded display is never used as a cutoff. Ages in a witness are earlier -> later. No witness is assigned to a passing curve.

| Curve | Rounded rho | Exact rho^2 | Printed-point rho | Limiting age pair(s) | First failing prefix ends |
|---|---:|---|---:|---|---|
| female breast / All Races / Females | 1.066238 | 9403/8271 | 1.066359 | 70-74 -> 80-84 | 75-79 |
| female breast / Whites / Females | 1.062759 | 3263/2889 | 1.062875 | 70-74 -> 80-84 | 75-79 |
| female breast / Blacks / Females | 1.056962 | 8791/7869 | 1.057089 | 75-79 -> 80-84 | 80-84 |
| colon and rectum / All Races / Both Sexes | 1.000000 | 1 | 1.000000 | none | none |
| colon and rectum / All Races / Males | 1.000000 | 1 | 1.000000 | none | none |
| colon and rectum / All Races / Females | 1.000000 | 1 | 1.000000 | none | none |
| colon and rectum / Whites / Both Sexes | 1.000000 | 1 | 1.000000 | none | none |
| colon and rectum / Whites / Males | 1.000000 | 1 | 1.000000 | none | none |
| colon and rectum / Whites / Females | 1.000000 | 1 | 1.000000 | none | none |
| colon and rectum / Blacks / Both Sexes | 1.000000 | 1 | 1.000000 | none | none |
| colon and rectum / Blacks / Males | 1.000000 | 1 | 1.000000 | none | none |
| colon and rectum / Blacks / Females | 1.000000 | 1 | 1.000000 | none | none |
| lung and bronchus / All Races / Both Sexes | 1.000000 | 1 | 1.000000 | none | none |
| lung and bronchus / All Races / Males | 1.000000 | 1 | 1.000000 | none | none |
| lung and bronchus / All Races / Females | 1.000000 | 1 | 1.000000 | none | none |
| lung and bronchus / Whites / Both Sexes | 1.000000 | 1 | 1.000000 | none | none |
| lung and bronchus / Whites / Males | 1.000000 | 1 | 1.000000 | none | none |
| lung and bronchus / Whites / Females | 1.000000 | 1 | 1.000000 | none | none |
| lung and bronchus / Blacks / Both Sexes | 1.023930 | 7187/6855 | 1.024076 | 75-79 -> 80-84 | 80-84 |
| lung and bronchus / Blacks / Males | 1.004756 | 9529/9439 | 1.004862 | 75-79 -> 80-84 | 80-84 |
| lung and bronchus / Blacks / Females | 1.024663 | 5677/5407 | 1.024849 | 75-79 -> 80-84 | 80-84 |
| pancreas / All Races / Both Sexes | 1.000000 | 1 | 1.000000 | none | none |
| pancreas / All Races / Males | 1.000000 | 1 | 1.000000 | none | none |
| pancreas / All Races / Females | 1.000000 | 1 | 1.000000 | none | none |
| pancreas / Whites / Both Sexes | 1.000000 | 1 | 1.000000 | none | none |
| pancreas / Whites / Males | 1.000000 | 1 | 1.000000 | none | none |
| pancreas / Whites / Females | 1.000000 | 1 | 1.000000 | none | none |
| pancreas / Blacks / Both Sexes | 1.000000 | 1 | 1.000000 | none | none |
| pancreas / Blacks / Males | 1.000000 | 1 | 1.000000 | none | none |
| pancreas / Blacks / Females | 1.000000 | 1 | 1.000000 | none | none |

### Every rounding-robust decrease

Rates and minimum declines are per 100,000. Every violating earlier/later pair is shown, including nonadjacent ages; absence from this table means no violating pair.

| Curve | Earlier age | Printed rate | Later age | Printed rate | Minimum decline |
|---|---|---:|---|---:|---:|
| female breast / All Races / Females | 65-69 | 426.5 | 80-84 | 413.5 | 12.9 |
| female breast / All Races / Females | 70-74 | 470.2 | 75-79 | 458.9 | 11.2 |
| female breast / All Races / Females | 70-74 | 470.2 | 80-84 | 413.5 | 56.6 |
| female breast / All Races / Females | 75-79 | 458.9 | 80-84 | 413.5 | 45.3 |
| female breast / Whites / Females | 65-69 | 441.3 | 80-84 | 433.3 | 7.9 |
| female breast / Whites / Females | 70-74 | 489.5 | 75-79 | 480.8 | 8.6 |
| female breast / Whites / Females | 70-74 | 489.5 | 80-84 | 433.3 | 56.1 |
| female breast / Whites / Females | 75-79 | 480.8 | 80-84 | 433.3 | 47.4 |
| female breast / Blacks / Females | 65-69 | 402.1 | 80-84 | 393.4 | 8.6 |
| female breast / Blacks / Females | 70-74 | 428.1 | 80-84 | 393.4 | 34.6 |
| female breast / Blacks / Females | 75-79 | 439.6 | 80-84 | 393.4 | 46.1 |
| lung and bronchus / Blacks / Both Sexes | 75-79 | 359.4 | 80-84 | 342.7 | 16.6 |
| lung and bronchus / Blacks / Males | 75-79 | 476.5 | 80-84 | 471.9 | 4.5 |
| lung and bronchus / Blacks / Females | 75-79 | 283.9 | 80-84 | 270.3 | 13.5 |

### Age-removal and truncation attempts

The full panel remains primary. These are diagnostic deletions, not a revised fitting panel. All tied smallest deletion sets are listed; `{}` denotes the empty set. Deleting the oldest bin equals the prefix ending at the preceding bin.

| Curve | Through 65-69 rho | Through 70-74 rho | Through 75-79 rho | Minimum deletion sets restoring compatibility |
|---|---:|---:|---:|---|
| female breast / All Races / Females | 1.000000 | 1.000000 | 1.012128 | {70-74, 80-84}; {75-79, 80-84} |
| female breast / Whites / Females | 1.000000 | 1.000000 | 1.008903 | {70-74, 80-84}; {75-79, 80-84} |
| female breast / Blacks / Females | 1.000000 | 1.000000 | 1.000000 | {80-84} |
| colon and rectum / All Races / Both Sexes | 1.000000 | 1.000000 | 1.000000 | {} |
| colon and rectum / All Races / Males | 1.000000 | 1.000000 | 1.000000 | {} |
| colon and rectum / All Races / Females | 1.000000 | 1.000000 | 1.000000 | {} |
| colon and rectum / Whites / Both Sexes | 1.000000 | 1.000000 | 1.000000 | {} |
| colon and rectum / Whites / Males | 1.000000 | 1.000000 | 1.000000 | {} |
| colon and rectum / Whites / Females | 1.000000 | 1.000000 | 1.000000 | {} |
| colon and rectum / Blacks / Both Sexes | 1.000000 | 1.000000 | 1.000000 | {} |
| colon and rectum / Blacks / Males | 1.000000 | 1.000000 | 1.000000 | {} |
| colon and rectum / Blacks / Females | 1.000000 | 1.000000 | 1.000000 | {} |
| lung and bronchus / All Races / Both Sexes | 1.000000 | 1.000000 | 1.000000 | {} |
| lung and bronchus / All Races / Males | 1.000000 | 1.000000 | 1.000000 | {} |
| lung and bronchus / All Races / Females | 1.000000 | 1.000000 | 1.000000 | {} |
| lung and bronchus / Whites / Both Sexes | 1.000000 | 1.000000 | 1.000000 | {} |
| lung and bronchus / Whites / Males | 1.000000 | 1.000000 | 1.000000 | {} |
| lung and bronchus / Whites / Females | 1.000000 | 1.000000 | 1.000000 | {} |
| lung and bronchus / Blacks / Both Sexes | 1.000000 | 1.000000 | 1.000000 | {75-79}; {80-84} |
| lung and bronchus / Blacks / Males | 1.000000 | 1.000000 | 1.000000 | {75-79}; {80-84} |
| lung and bronchus / Blacks / Females | 1.000000 | 1.000000 | 1.000000 | {75-79}; {80-84} |
| pancreas / All Races / Both Sexes | 1.000000 | 1.000000 | 1.000000 | {} |
| pancreas / All Races / Males | 1.000000 | 1.000000 | 1.000000 | {} |
| pancreas / All Races / Females | 1.000000 | 1.000000 | 1.000000 | {} |
| pancreas / Whites / Both Sexes | 1.000000 | 1.000000 | 1.000000 | {} |
| pancreas / Whites / Males | 1.000000 | 1.000000 | 1.000000 | {} |
| pancreas / Whites / Females | 1.000000 | 1.000000 | 1.000000 | {} |
| pancreas / Blacks / Both Sexes | 1.000000 | 1.000000 | 1.000000 | {} |
| pancreas / Blacks / Males | 1.000000 | 1.000000 | 1.000000 | {} |
| pancreas / Blacks / Females | 1.000000 | 1.000000 | 1.000000 | {} |

<details>
<summary>Every prefix and leave-one-age-out frontier</summary>

Vectors follow ages 30-34, 35-39, 40-44, 45-49, 50-54, 55-59, 60-64, 65-69, 70-74, 75-79, 80-84. Prefix entries end at the indicated age; deletion entries remove only that age bin. All entries are rho.

| Curve | Prefix vector | Leave-one-age-out vector |
|---|---|---|
| female breast / All Races / Females | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.012128, 1.066238 | 1.066238, 1.066238, 1.066238, 1.066238, 1.066238, 1.066238, 1.066238, 1.066238, 1.053347, 1.066238, 1.012128 |
| female breast / Whites / Females | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.008903, 1.062759 | 1.062759, 1.062759, 1.062759, 1.062759, 1.062759, 1.062759, 1.062759, 1.062759, 1.053271, 1.062759, 1.008903 |
| female breast / Blacks / Females | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.056962 | 1.056962, 1.056962, 1.056962, 1.056962, 1.056962, 1.056962, 1.056962, 1.056962, 1.056962, 1.043044, 1.000000 |
| colon and rectum / All Races / Both Sexes | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| colon and rectum / All Races / Males | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| colon and rectum / All Races / Females | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| colon and rectum / Whites / Both Sexes | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| colon and rectum / Whites / Males | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| colon and rectum / Whites / Females | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| colon and rectum / Blacks / Both Sexes | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| colon and rectum / Blacks / Males | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| colon and rectum / Blacks / Females | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| lung and bronchus / All Races / Both Sexes | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| lung and bronchus / All Races / Males | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| lung and bronchus / All Races / Females | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| lung and bronchus / Whites / Both Sexes | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| lung and bronchus / Whites / Males | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| lung and bronchus / Whites / Females | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| lung and bronchus / Blacks / Both Sexes | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.023930 | 1.023930, 1.023930, 1.023930, 1.023930, 1.023930, 1.023930, 1.023930, 1.023930, 1.023930, 1.000000, 1.000000 |
| lung and bronchus / Blacks / Males | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.004756 | 1.004756, 1.004756, 1.004756, 1.004756, 1.004756, 1.004756, 1.004756, 1.004756, 1.004756, 1.000000, 1.000000 |
| lung and bronchus / Blacks / Females | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.024663 | 1.024663, 1.024663, 1.024663, 1.024663, 1.024663, 1.024663, 1.024663, 1.024663, 1.024663, 1.000000, 1.000000 |
| pancreas / All Races / Both Sexes | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| pancreas / All Races / Males | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| pancreas / All Races / Females | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| pancreas / Whites / Both Sexes | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| pancreas / Whites / Males | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| pancreas / Whites / Females | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| pancreas / Blacks / Both Sexes | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| pancreas / Blacks / Males | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |
| pancreas / Blacks / Females | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 | 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000 |

</details>

### Specified-shape calibration frontiers

These are finite-grid minima under uniform within-bin averaging. `L / U ages` name the lower and upper amplitude constraints, in that order, not necessarily chronological order. The grid index k is not a biological driver estimate. Negative g is net loss.

| Curve | Pure-power k | Pure rho | Pure L / U ages | One-growth k | g/year | Grid rho | Grid L / U ages | Unique grid minimum? |
|---|---:|---:|---|---:|---:|---:|---|---|
| female breast / All Races / Females | 4 | 1.578568 | 45-49 / 80-84 | 4 | -0.015 | 1.504293 | 45-49 / 30-34 | yes |
| female breast / Whites / Females | 4 | 1.549890 | 45-49 / 80-84 | 4 | -0.010 | 1.512726 | 45-49 / 30-34 | yes |
| female breast / Blacks / Females | 4 | 1.554288 | 45-49 / 80-84 | 4 | -0.030 | 1.408107 | 45-49 / 80-84 | yes |
| colon and rectum / All Races / Both Sexes | 5 | 1.265066 | 50-54 / 80-84 | 5 | -0.010 | 1.241860 | 50-54 / 30-34 | yes |
| colon and rectum / All Races / Males | 5 | 1.300289 | 50-54 / 30-34 | 5 | 0.005 | 1.287041 | 50-54 / 30-34 | yes |
| colon and rectum / All Races / Females | 5 | 1.267344 | 50-54 / 70-74 | 5 | -0.030 | 1.212991 | 50-54 / 35-39 | yes |
| colon and rectum / Whites / Both Sexes | 5 | 1.221759 | 50-54 / 80-84 | 5 | -0.005 | 1.204393 | 50-54 / 80-84 | yes |
| colon and rectum / Whites / Males | 5 | 1.249364 | 50-54 / 30-34 | 5 | 0.005 | 1.241608 | 50-54 / 80-84 | yes |
| colon and rectum / Whites / Females | 5 | 1.240077 | 50-54 / 70-74 | 5 | -0.045 | 1.190490 | 50-54 / 35-39 | yes |
| colon and rectum / Blacks / Both Sexes | 5 | 1.377748 | 50-54 / 80-84 | 5 | -0.005 | 1.361139 | 50-54 / 30-34 | yes |
| colon and rectum / Blacks / Males | 5 | 1.477013 | 50-54 / 30-34 | 5 | 0.020 | 1.412458 | 50-54 / 30-34 | yes |
| colon and rectum / Blacks / Females | 5 | 1.425479 | 50-54 / 80-84 | 5 | -0.025 | 1.342072 | 50-54 / 80-84 | yes |
| lung and bronchus / All Races / Both Sexes | 7 | 1.576755 | 55-59 / 35-39 | 7 | 0.040 | 1.482398 | 55-59 / 80-84 | yes |
| lung and bronchus / All Races / Males | 8 | 1.502307 | 55-59 / 80-84 | 8 | -0.060 | 1.413102 | 55-59 / 35-39 | yes |
| lung and bronchus / All Races / Females | 7 | 1.509827 | 55-59 / 35-39 | 7 | 0.020 | 1.475344 | 55-59 / 80-84 | yes |
| lung and bronchus / Whites / Both Sexes | 8 | 1.576140 | 55-59 / 80-84 | 8 | -0.060 | 1.477024 | 55-59 / 80-84 | yes |
| lung and bronchus / Whites / Males | 8 | 1.467828 | 55-59 / 80-84 | 8 | -0.015 | 1.437103 | 55-59 / 80-84 | yes |
| lung and bronchus / Whites / Females | 7 | 1.565526 | 55-59 / 35-39 | 7 | 0.030 | 1.494009 | 55-59 / 80-84 | yes |
| lung and bronchus / Blacks / Both Sexes | 7 | 1.664506 | 55-59 / 40-44 | 7 | 0.000 | 1.664506 | 55-59 / 40-44 | yes |
| lung and bronchus / Blacks / Males | 7 | 1.710494 | 60-64 / 40-44 | 7 | 0.030 | 1.644055 | 55-59 / 40-44 | yes |
| lung and bronchus / Blacks / Females | 7 | 1.730492 | 55-59 / 80-84 | 7 | -0.015 | 1.689551 | 55-59 / 80-84 | yes |
| pancreas / All Races / Both Sexes | 6 | 1.333798 | 65-69 / 35-39 | 6 | 0.030 | 1.245370 | 60-64 / 80-84 | yes |
| pancreas / All Races / Males | 7 | 1.389806 | 55-59 / 80-84 | 6 | 0.035 | 1.297958 | 60-64 / 80-84 | yes |
| pancreas / All Races / Females | 6 | 1.239332 | 65-69 / 40-44 | 6 | 0.025 | 1.174503 | 65-69 / 80-84 | yes |
| pancreas / Whites / Both Sexes | 6 | 1.343438 | 65-69 / 35-39 | 6 | 0.030 | 1.235369 | 65-69 / 80-84 | yes |
| pancreas / Whites / Males | 7 | 1.366955 | 55-59 / 80-84 | 7 | -0.060 | 1.275768 | 60-64 / 80-84 | yes |
| pancreas / Whites / Females | 6 | 1.231760 | 65-69 / 35-39 | 6 | 0.020 | 1.172564 | 65-69 / 80-84 | yes |
| pancreas / Blacks / Both Sexes | 6 | 1.393231 | 60-64 / 40-44 | 6 | 0.015 | 1.358507 | 60-64 / 80-84 | yes |
| pancreas / Blacks / Males | 6 | 1.462279 | 60-64 / 40-44 | 6 | 0.000 | 1.462279 | 60-64 / 40-44 | yes |
| pancreas / Blacks / Females | 6 | 1.357306 | 60-64 / 40-44 | 6 | 0.030 | 1.284923 | 60-64 / 80-84 | yes |

Certified unique joint grid minima: **30/30**. Per-k selected-grid rows with certified single L and U ages: **270/270**. Every reported grid threshold has a rational squared-threshold enclosure narrower than `1e-24`; both bounds round to the displayed rho. Full rational enclosures are emitted by the CLI.

<details>
<summary>Complete per-k frontiers, bin-weight sensitivity and constraint ages</summary>

`Weight rho` is the pure-power infimum over freely varying within-bin age weights. It is not a confidence bound or a guaranteed attained fit. A list in an age cell retains tied/unresolved constraints. g contenders are retained at every k.

| Curve | k | Pure rho | Pure L / U ages | Weight rho | Weight L / U ages | Grid rho | g/year | g contenders | Grid L / U ages |
|---|---:|---:|---|---:|---|---:|---:|---|---|
| female breast / All Races / Females | 2 | 2.675151 | 70-74 / 30-34 | 2.527002 | 70-74 / 30-34 | 1.722975 | 0.045 | 0.045 | 45-49 / 30-34 |
| female breast / All Races / Females | 3 | 1.833579 | 65-69 / 30-34 | 1.630852 | 65-69 / 30-34 | 1.601425 | 0.035 | 0.035 | 45-49 / 80-84 |
| female breast / All Races / Females | 4 | 1.578568 | 45-49 / 80-84 | 1.397025 | 45-49 / 80-84 | 1.504293 | -0.015 | -0.015 | 45-49 / 30-34 |
| female breast / All Races / Females | 5 | 2.112619 | 40-44 / 80-84 | 1.776420 | 40-44 / 80-84 | 1.822619 | -0.060 | -0.060 | 45-49 / 80-84 |
| female breast / All Races / Females | 6 | 2.938490 | 40-44 / 80-84 | 2.368560 | 40-44 / 80-84 | 2.570781 | -0.060 | -0.060 | 40-44 / 80-84 |
| female breast / All Races / Females | 7 | 4.304212 | 30-34 / 80-84 | 3.181309 | 30-34 / 80-84 | 3.693617 | -0.060 | -0.060 | 30-34 / 80-84 |
| female breast / All Races / Females | 8 | 6.824417 | 30-34 / 80-84 | 4.809687 | 30-34 / 80-84 | 5.945817 | -0.060 | -0.060 | 30-34 / 80-84 |
| female breast / All Races / Females | 9 | 10.811953 | 30-34 / 80-84 | 7.271564 | 30-34 / 80-84 | 9.536809 | -0.060 | -0.060 | 30-34 / 80-84 |
| female breast / All Races / Females | 10 | 17.116569 | 30-34 / 80-84 | 10.993571 | 30-34 / 80-84 | 15.253224 | -0.060 | -0.060 | 30-34 / 80-84 |
| female breast / Whites / Females | 2 | 2.752977 | 70-74 / 30-34 | 2.600518 | 70-74 / 30-34 | 1.696127 | 0.050 | 0.050 | 45-49 / 30-34 |
| female breast / Whites / Females | 3 | 1.881162 | 65-69 / 30-34 | 1.673175 | 65-69 / 30-34 | 1.609144 | 0.035 | 0.035 | 45-49 / 30-34 |
| female breast / Whites / Females | 4 | 1.549890 | 45-49 / 80-84 | 1.371645 | 45-49 / 80-84 | 1.512726 | -0.010 | -0.010 | 45-49 / 30-34 |
| female breast / Whites / Females | 5 | 2.066159 | 40-44 / 80-84 | 1.737354 | 40-44 / 80-84 | 1.789507 | -0.060 | -0.060 | 45-49 / 80-84 |
| female breast / Whites / Females | 6 | 2.873869 | 40-44 / 80-84 | 2.316472 | 40-44 / 80-84 | 2.514246 | -0.060 | -0.060 | 40-44 / 80-84 |
| female breast / Whites / Females | 7 | 4.168762 | 30-34 / 80-84 | 3.088629 | 40-44 / 80-84 | 3.577382 | -0.060 | -0.060 | 30-34 / 80-84 |
| female breast / Whites / Females | 8 | 6.609659 | 30-34 / 80-84 | 4.658331 | 30-34 / 80-84 | 5.758707 | -0.060 | -0.060 | 30-34 / 80-84 |
| female breast / Whites / Females | 9 | 10.471710 | 30-34 / 80-84 | 7.042734 | 30-34 / 80-84 | 9.236694 | -0.060 | -0.060 | 30-34 / 80-84 |
| female breast / Whites / Females | 10 | 16.577925 | 30-34 / 80-84 | 10.647613 | 30-34 / 80-84 | 14.773218 | -0.060 | -0.060 | 30-34 / 80-84 |
| female breast / Blacks / Females | 2 | 2.438677 | 65-69 / 30-34 | 2.300786 | 65-69 / 30-34 | 1.573965 | 0.045 | 0.045 | 45-49 / 30-34 |
| female breast / Blacks / Females | 3 | 1.693451 | 65-69 / 30-34 | 1.506218 | 65-69 / 30-34 | 1.496317 | 0.030 | 0.030 | 45-49 / 80-84 |
| female breast / Blacks / Females | 4 | 1.554288 | 45-49 / 80-84 | 1.375537 | 45-49 / 80-84 | 1.408107 | -0.030 | -0.030 | 45-49 / 80-84 |
| female breast / Blacks / Females | 5 | 2.127439 | 40-44 / 80-84 | 1.788882 | 40-44 / 80-84 | 1.826229 | -0.060 | -0.060 | 40-44 / 80-84 |
| female breast / Blacks / Females | 6 | 3.030779 | 35-39 / 80-84 | 2.402185 | 35-39 / 80-84 | 2.601688 | -0.060 | -0.060 | 35-39 / 80-84 |
| female breast / Blacks / Females | 7 | 4.639981 | 30-34 / 80-84 | 3.429481 | 30-34 / 80-84 | 3.981753 | -0.060 | -0.060 | 30-34 / 80-84 |
| female breast / Blacks / Females | 8 | 7.356786 | 30-34 / 80-84 | 5.184888 | 30-34 / 80-84 | 6.409646 | -0.060 | -0.060 | 30-34 / 80-84 |
| female breast / Blacks / Females | 9 | 11.655387 | 30-34 / 80-84 | 7.838814 | 30-34 / 80-84 | 10.280770 | -0.060 | -0.060 | 30-34 / 80-84 |
| female breast / Blacks / Females | 10 | 18.451822 | 30-34 / 80-84 | 11.851172 | 30-34 / 80-84 | 16.443118 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / All Races / Both Sexes | 2 | 3.876099 | 80-84 / 30-34 | 3.668862 | 80-84 / 30-34 | 1.646597 | 0.060 | 0.060 | 50-54 / 30-34 |
| colon and rectum / All Races / Both Sexes | 3 | 2.434843 | 80-84 / 30-34 | 2.179629 | 80-84 / 30-34 | 1.462377 | 0.060 | 0.060 | 50-54 / 30-34 |
| colon and rectum / All Races / Both Sexes | 4 | 1.546167 | 50-54 / 30-34 | 1.294892 | 80-84 / 30-34 | 1.331451 | 0.045 | 0.045 | 50-54 / 30-34 |
| colon and rectum / All Races / Both Sexes | 5 | 1.265066 | 50-54 / 80-84 | 1.085332 | 50-54 / 80-84 | 1.241860 | -0.010 | -0.010 | 50-54 / 30-34 |
| colon and rectum / All Races / Both Sexes | 6 | 1.634967 | 30-34 / 80-84 | 1.308960 | 50-54 / 80-84 | 1.437426 | -0.060 | -0.060 | 50-54 / 80-84 |
| colon and rectum / All Races / Both Sexes | 7 | 2.594305 | 30-34 / 80-84 | 1.917491 | 30-34 / 80-84 | 2.226277 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / All Races / Both Sexes | 8 | 4.113324 | 30-34 / 80-84 | 2.898973 | 30-34 / 80-84 | 3.583760 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / All Races / Both Sexes | 9 | 6.516757 | 30-34 / 80-84 | 4.382836 | 30-34 / 80-84 | 5.748181 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / All Races / Both Sexes | 10 | 10.316778 | 30-34 / 80-84 | 6.626225 | 30-34 / 80-84 | 9.193672 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / All Races / Males | 2 | 4.157463 | 80-84 / 30-34 | 3.935182 | 80-84 / 30-34 | 1.756810 | 0.060 | 0.060 | 50-54 / 30-34 |
| colon and rectum / All Races / Males | 3 | 2.611587 | 80-84 / 30-34 | 2.337847 | 80-84 / 30-34 | 1.560259 | 0.060 | 0.060 | 50-54 / 30-34 |
| colon and rectum / All Races / Males | 4 | 1.649657 | 50-54 / 30-34 | 1.388888 | 80-84 / 30-34 | 1.389835 | 0.050 | 0.050 | 50-54 / 30-34 |
| colon and rectum / All Races / Males | 5 | 1.300289 | 50-54 / 30-34 | 1.079636 | 50-54 / 80-84 | 1.287041 | 0.005 | 0.005 | 50-54 / 30-34 |
| colon and rectum / All Races / Males | 6 | 1.576112 | 50-54 / 80-84 | 1.302091 | 50-54 / 80-84 | 1.429882 | -0.060 | -0.060 | 50-54 / 80-84 |
| colon and rectum / All Races / Males | 7 | 2.418460 | 30-34 / 80-84 | 1.787521 | 30-34 / 80-84 | 2.075377 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / All Races / Males | 8 | 3.834519 | 30-34 / 80-84 | 2.702478 | 30-34 / 80-84 | 3.340849 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / All Races / Males | 9 | 6.075044 | 30-34 / 80-84 | 4.085762 | 30-34 / 80-84 | 5.358564 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / All Races / Males | 10 | 9.617496 | 30-34 / 80-84 | 6.177092 | 30-34 / 80-84 | 8.570515 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / All Races / Females | 2 | 3.694179 | 80-84 / 30-34 | 3.496668 | 80-84 / 30-34 | 1.547627 | 0.060 | 0.060 | 50-54 / 30-34 |
| colon and rectum / All Races / Females | 3 | 2.320567 | 80-84 / 30-34 | 2.077330 | 80-84 / 30-34 | 1.374479 | 0.060 | 0.060 | 50-54 / 30-34 |
| colon and rectum / All Races / Females | 4 | 1.458912 | 80-84 / 30-34 | 1.234118 | 80-84 / 30-34 | 1.277637 | 0.040 | 0.040 | 50-54 / 30-34 |
| colon and rectum / All Races / Females | 5 | 1.267344 | 50-54 / 70-74 | 1.077642 | 50-54 / 70-74 | 1.212991 | -0.030 | -0.030 | 50-54 / 35-39 |
| colon and rectum / All Races / Females | 6 | 1.715444 | 30-34 / 80-84 | 1.330730 | 30-34 / 80-84 | 1.444083 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / All Races / Females | 7 | 2.722003 | 30-34 / 80-84 | 2.011874 | 30-34 / 80-84 | 2.335860 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / All Races / Females | 8 | 4.315793 | 30-34 / 80-84 | 3.041668 | 30-34 / 80-84 | 3.760161 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / All Races / Females | 9 | 6.837528 | 30-34 / 80-84 | 4.598570 | 30-34 / 80-84 | 6.031121 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / All Races / Females | 10 | 10.824596 | 30-34 / 80-84 | 6.952384 | 30-34 / 80-84 | 9.646208 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Whites / Both Sexes | 2 | 3.849377 | 80-84 / 30-34 | 3.643568 | 80-84 / 30-34 | 1.579260 | 0.060 | 0.060 | 50-54 / 30-34 |
| colon and rectum / Whites / Both Sexes | 3 | 2.418057 | 80-84 / 30-34 | 2.164602 | 80-84 / 30-34 | 1.402573 | 0.060 | 0.060 | 50-54 / 30-34 |
| colon and rectum / Whites / Both Sexes | 4 | 1.520203 | 80-84 / 30-34 | 1.285965 | 80-84 / 30-34 | 1.277001 | 0.045 | 0.045 | 50-54 / 30-34 |
| colon and rectum / Whites / Both Sexes | 5 | 1.221759 | 50-54 / 80-84 | 1.048178 | 50-54 / 80-84 | 1.204393 | -0.005 | -0.005 | 50-54 / 80-84 |
| colon and rectum / Whites / Both Sexes | 6 | 1.646751 | 30-34 / 80-84 | 1.277442 | 30-34 / 80-84 | 1.388218 | -0.060 | -0.060 | 50-54 / 80-84 |
| colon and rectum / Whites / Both Sexes | 7 | 2.613004 | 30-34 / 80-84 | 1.931311 | 30-34 / 80-84 | 2.242323 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Whites / Both Sexes | 8 | 4.142972 | 30-34 / 80-84 | 2.919868 | 30-34 / 80-84 | 3.609591 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Whites / Both Sexes | 9 | 6.563728 | 30-34 / 80-84 | 4.414426 | 30-34 / 80-84 | 5.789613 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Whites / Both Sexes | 10 | 10.391139 | 30-34 / 80-84 | 6.673985 | 30-34 / 80-84 | 9.259938 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Whites / Males | 2 | 4.112685 | 80-84 / 30-34 | 3.892799 | 80-84 / 30-34 | 1.688005 | 0.060 | 0.060 | 50-54 / 30-34 |
| colon and rectum / Whites / Males | 3 | 2.583459 | 80-84 / 30-34 | 2.312667 | 80-84 / 30-34 | 1.499152 | 0.060 | 0.060 | 50-54 / 30-34 |
| colon and rectum / Whites / Males | 4 | 1.624190 | 80-84 / 30-34 | 1.373929 | 80-84 / 30-34 | 1.335402 | 0.050 | 0.050 | 50-54 / 30-34 |
| colon and rectum / Whites / Males | 5 | 1.249364 | 50-54 / 30-34 | 1.048649 | 50-54 / 80-84 | 1.241608 | 0.005 | 0.005 | 50-54 / 80-84 |
| colon and rectum / Whites / Males | 6 | 1.541159 | 30-34 / 80-84 | 1.264719 | 50-54 / 80-84 | 1.388842 | -0.060 | -0.060 | 50-54 / 80-84 |
| colon and rectum / Whites / Males | 7 | 2.445454 | 30-34 / 80-84 | 1.807472 | 30-34 / 80-84 | 2.098542 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Whites / Males | 8 | 3.877318 | 30-34 / 80-84 | 2.732641 | 30-34 / 80-84 | 3.378137 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Whites / Males | 9 | 6.142851 | 30-34 / 80-84 | 4.131366 | 30-34 / 80-84 | 5.418373 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Whites / Males | 10 | 9.724841 | 30-34 / 80-84 | 6.246038 | 30-34 / 80-84 | 8.666175 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Whites / Females | 2 | 3.676263 | 80-84 / 30-34 | 3.479710 | 80-84 / 30-34 | 1.480587 | 0.060 | 0.060 | 50-54 / 30-34 |
| colon and rectum / Whites / Females | 3 | 2.309313 | 80-84 / 30-34 | 2.067256 | 80-84 / 30-34 | 1.314939 | 0.060 | 0.060 | 50-54 / 30-34 |
| colon and rectum / Whites / Females | 4 | 1.451837 | 80-84 / 30-34 | 1.228133 | 80-84 / 30-34 | 1.222293 | 0.040 | 0.040 | 50-54 / 30-34 |
| colon and rectum / Whites / Females | 5 | 1.240077 | 50-54 / 70-74 | 1.054457 | 50-54 / 70-74 | 1.190490 | -0.045 | -0.045 | 50-54 / 35-39 |
| colon and rectum / Whites / Females | 6 | 1.724261 | 30-34 / 80-84 | 1.337570 | 30-34 / 80-84 | 1.451506 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Whites / Females | 7 | 2.735994 | 30-34 / 80-84 | 2.022215 | 30-34 / 80-84 | 2.347866 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Whites / Females | 8 | 4.337976 | 30-34 / 80-84 | 3.057302 | 30-34 / 80-84 | 3.779489 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Whites / Females | 9 | 6.872673 | 30-34 / 80-84 | 4.622206 | 30-34 / 80-84 | 6.062121 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Whites / Females | 10 | 10.880234 | 30-34 / 80-84 | 6.988119 | 30-34 / 80-84 | 9.695789 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Blacks / Both Sexes | 2 | 3.936769 | 80-84 / 30-34 | 3.726288 | 80-84 / 30-34 | 1.821319 | 0.060 | 0.060 | 50-54 / 30-34 |
| colon and rectum / Blacks / Both Sexes | 3 | 2.472954 | 80-84 / 30-34 | 2.213745 | 80-84 / 30-34 | 1.617551 | 0.060 | 0.060 | 50-54 / 30-34 |
| colon and rectum / Blacks / Both Sexes | 4 | 1.710232 | 50-54 / 30-34 | 1.411944 | 50-54 / 30-34 | 1.472733 | 0.045 | 0.045 | 50-54 / 30-34 |
| colon and rectum / Blacks / Both Sexes | 5 | 1.377748 | 50-54 / 80-84 | 1.182005 | 50-54 / 80-84 | 1.361139 | -0.005 | -0.005 | 50-54 / 30-34 |
| colon and rectum / Blacks / Both Sexes | 6 | 1.725555 | 50-54 / 80-84 | 1.425552 | 50-54 / 80-84 | 1.565460 | -0.060 | -0.060 | 50-54 / 80-84 |
| colon and rectum / Blacks / Both Sexes | 7 | 2.554340 | 30-34 / 80-84 | 1.887952 | 30-34 / 80-84 | 2.191982 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Blacks / Both Sexes | 8 | 4.049960 | 30-34 / 80-84 | 2.854315 | 30-34 / 80-84 | 3.528553 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Blacks / Both Sexes | 9 | 6.416368 | 30-34 / 80-84 | 4.315319 | 30-34 / 80-84 | 5.659632 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Blacks / Both Sexes | 10 | 10.157850 | 30-34 / 80-84 | 6.524149 | 30-34 / 80-84 | 9.052046 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Blacks / Males | 2 | 4.609543 | 80-84 / 30-34 | 4.363091 | 80-84 / 30-34 | 1.995580 | 0.060 | 0.060 | 50-54 / 30-34 |
| colon and rectum / Blacks / Males | 3 | 2.895569 | 80-84 / 30-34 | 2.592062 | 80-84 / 30-34 | 1.772316 | 0.060 | 0.060 | 50-54 / 30-34 |
| colon and rectum / Blacks / Males | 4 | 1.933007 | 65-69 / 30-34 | 1.619629 | 65-69 / 30-34 | 1.535964 | 0.060 | 0.060 | 50-54 / 80-84 |
| colon and rectum / Blacks / Males | 5 | 1.477013 | 50-54 / 30-34 | 1.142563 | 50-54 / 30-34 | 1.412458 | 0.020 | 0.020 | 50-54 / 30-34 |
| colon and rectum / Blacks / Males | 6 | 1.614773 | 50-54 / 80-84 | 1.334029 | 50-54 / 80-84 | 1.464956 | -0.060 | -0.060 | 50-54 / 80-84 |
| colon and rectum / Blacks / Males | 7 | 2.179665 | 30-34 / 80-84 | 1.611024 | 30-34 / 80-84 | 1.870458 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Blacks / Males | 8 | 3.455904 | 30-34 / 80-84 | 2.435639 | 30-34 / 80-84 | 3.010978 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Blacks / Males | 9 | 5.475203 | 30-34 / 80-84 | 3.682340 | 30-34 / 80-84 | 4.829466 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Blacks / Males | 10 | 8.667878 | 30-34 / 80-84 | 5.567175 | 30-34 / 80-84 | 7.724274 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Blacks / Females | 2 | 3.499591 | 80-84 / 30-34 | 3.312484 | 80-84 / 30-34 | 1.675206 | 0.060 | 0.060 | 50-54 / 30-34 |
| colon and rectum / Blacks / Females | 3 | 2.198333 | 80-84 / 30-34 | 1.967909 | 80-84 / 30-34 | 1.534092 | 0.060 | 0.060 | 50-54 / 80-84 |
| colon and rectum / Blacks / Females | 4 | 1.573031 | 50-54 / 30-34 | 1.298673 | 50-54 / 30-34 | 1.410358 | 0.035 | 0.035 | 50-54 / 30-34 |
| colon and rectum / Blacks / Females | 5 | 1.425479 | 50-54 / 80-84 | 1.222955 | 50-54 / 80-84 | 1.342072 | -0.025 | -0.025 | 50-54 / 80-84 |
| colon and rectum / Blacks / Females | 6 | 1.812151 | 30-34 / 80-84 | 1.474939 | 50-54 / 80-84 | 1.619695 | -0.060 | -0.060 | 50-54 / 80-84 |
| colon and rectum / Blacks / Females | 7 | 2.875454 | 30-34 / 80-84 | 2.125293 | 30-34 / 80-84 | 2.467543 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Blacks / Females | 8 | 4.559093 | 30-34 / 80-84 | 3.213140 | 30-34 / 80-84 | 3.972138 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Blacks / Females | 9 | 7.222990 | 30-34 / 80-84 | 4.857811 | 30-34 / 80-84 | 6.371122 | -0.060 | -0.060 | 30-34 / 80-84 |
| colon and rectum / Blacks / Females | 10 | 11.434826 | 30-34 / 80-84 | 7.344321 | 30-34 / 80-84 | 10.190007 | -0.060 | -0.060 | 30-34 / 80-84 |
| lung and bronchus / All Races / Both Sexes | 2 | 11.445072 | 75-79 / 30-34 | 10.822892 | 75-79 / 30-34 | 4.484726 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / All Races / Both Sexes | 3 | 7.417575 | 75-79 / 30-34 | 6.627641 | 75-79 / 30-34 | 3.801923 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / All Races / Both Sexes | 4 | 4.811223 | 75-79 / 30-34 | 4.058585 | 75-79 / 30-34 | 2.985218 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / All Races / Both Sexes | 5 | 3.224045 | 70-74 / 30-34 | 2.554995 | 70-74 / 30-34 | 2.226626 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / All Races / Both Sexes | 6 | 2.165324 | 70-74 / 30-34 | 1.631620 | 70-74 / 35-39 | 1.689266 | 0.060 | 0.060 | 55-59 / 35-39 |
| lung and bronchus / All Races / Both Sexes | 7 | 1.576755 | 55-59 / 35-39 | 1.127094 | 65-69 / 35-39 | 1.482398 | 0.040 | 0.040 | 55-59 / 80-84 |
| lung and bronchus / All Races / Both Sexes | 8 | 1.607992 | 55-59 / 80-84 | 1.248204 | 55-59 / 80-84 | 1.506873 | -0.060 | -0.060 | 55-59 / 80-84 |
| lung and bronchus / All Races / Both Sexes | 9 | 2.158520 | 30-34 / 80-84 | 1.451709 | 30-34 / 80-84 | 1.903948 | -0.060 | -0.060 | 30-34 / 80-84 |
| lung and bronchus / All Races / Both Sexes | 10 | 3.417186 | 30-34 / 80-84 | 2.194778 | 30-34 / 80-84 | 3.045184 | -0.060 | -0.060 | 30-34 / 80-84 |
| lung and bronchus / All Races / Males | 2 | 12.947192 | 80-84 / 30-34 | 12.254964 | 80-84 / 30-34 | 4.936100 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / All Races / Males | 3 | 8.311849 | 75-79 / 30-34 | 7.426680 | 75-79 / 30-34 | 4.184575 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / All Races / Males | 4 | 5.391271 | 75-79 / 30-34 | 4.547894 | 75-79 / 30-34 | 3.285671 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / All Races / Males | 5 | 3.548536 | 70-74 / 30-34 | 2.812148 | 70-74 / 30-34 | 2.461468 | 0.060 | 0.060 | 65-69 / 30-34 |
| lung and bronchus / All Races / Males | 6 | 2.383258 | 70-74 / 30-34 | 1.778558 | 70-74 / 30-34 | 1.822476 | 0.060 | 0.060 | 65-69 / 30-34 |
| lung and bronchus / All Races / Males | 7 | 1.684008 | 65-69 / 35-39 | 1.218349 | 65-69 / 35-39 | 1.482836 | 0.055 | 0.055 | 55-59 / 35-39 |
| lung and bronchus / All Races / Males | 8 | 1.502307 | 55-59 / 80-84 | 1.166165 | 55-59 / 80-84 | 1.413102 | -0.060 | -0.060 | 55-59 / 35-39 |
| lung and bronchus / All Races / Males | 9 | 1.879769 | 30-34 / 80-84 | 1.346572 | 55-59 / 80-84 | 1.693728 | -0.060 | -0.060 | 55-59 / 80-84 |
| lung and bronchus / All Races / Males | 10 | 2.975892 | 30-34 / 80-84 | 1.911346 | 30-34 / 80-84 | 2.651930 | -0.060 | -0.060 | 30-34 / 80-84 |
| lung and bronchus / All Races / Females | 2 | 10.729459 | 75-79 / 30-34 | 10.146182 | 75-79 / 30-34 | 4.264029 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / All Races / Females | 3 | 6.953784 | 75-79 / 30-34 | 6.213242 | 75-79 / 30-34 | 3.614828 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / All Races / Females | 4 | 4.567677 | 70-74 / 30-34 | 3.841000 | 70-74 / 30-34 | 2.838313 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / All Races / Females | 5 | 3.065387 | 70-74 / 30-34 | 2.429262 | 70-74 / 30-34 | 2.117052 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / All Races / Females | 6 | 2.058767 | 70-74 / 30-34 | 1.536400 | 70-74 / 30-34 | 1.628330 | 0.060 | 0.060 | 55-59 / 30-34 |
| lung and bronchus / All Races / Females | 7 | 1.509827 | 55-59 / 35-39 | 1.138471 | 55-59 / 80-84 | 1.475344 | 0.020 | 0.020 | 55-59 / 80-84 |
| lung and bronchus / All Races / Females | 8 | 1.703059 | 50-54 / 80-84 | 1.314594 | 55-59 / 80-84 | 1.587021 | -0.060 | -0.060 | 55-59 / 80-84 |
| lung and bronchus / All Races / Females | 9 | 2.330537 | 30-34 / 80-84 | 1.574499 | 50-54 / 80-84 | 2.055677 | -0.060 | -0.060 | 30-34 / 80-84 |
| lung and bronchus / All Races / Females | 10 | 3.689508 | 30-34 / 80-84 | 2.369685 | 30-34 / 80-84 | 3.287861 | -0.060 | -0.060 | 30-34 / 80-84 |
| lung and bronchus / Whites / Both Sexes | 2 | 12.238571 | 75-79 / 30-34 | 11.573254 | 75-79 / 30-34 | 4.788700 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / Whites / Both Sexes | 3 | 7.931843 | 75-79 / 30-34 | 7.087142 | 75-79 / 30-34 | 4.059617 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / Whites / Both Sexes | 4 | 5.144790 | 75-79 / 30-34 | 4.339970 | 75-79 / 30-34 | 3.187555 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / Whites / Both Sexes | 5 | 3.442570 | 70-74 / 30-34 | 2.728172 | 70-74 / 30-34 | 2.377546 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / Whites / Both Sexes | 6 | 2.312089 | 70-74 / 30-34 | 1.737939 | 70-74 / 35-39 | 1.763750 | 0.060 | 0.060 | 65-69 / 35-39 |
| lung and bronchus / Whites / Both Sexes | 7 | 1.647064 | 65-69 / 35-39 | 1.191621 | 65-69 / 35-39 | 1.505494 | 0.050 | 0.050 | 55-59 / 80-84 |
| lung and bronchus / Whites / Both Sexes | 8 | 1.576140 | 55-59 / 80-84 | 1.223479 | 55-59 / 80-84 | 1.477024 | -0.060 | -0.060 | 55-59 / 80-84 |
| lung and bronchus / Whites / Both Sexes | 9 | 2.015730 | 30-34 / 80-84 | 1.420229 | 50-54 / 80-84 | 1.788964 | -0.060 | -0.060 | 50-54 / 80-84 |
| lung and bronchus / Whites / Both Sexes | 10 | 3.191134 | 30-34 / 80-84 | 2.049590 | 30-34 / 80-84 | 2.843740 | -0.060 | -0.060 | 30-34 / 80-84 |
| lung and bronchus / Whites / Males | 2 | 13.718937 | 80-84 / 30-34 | 12.985448 | 80-84 / 30-34 | 5.202034 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / Whites / Males | 3 | 8.813604 | 75-79 / 30-34 | 7.875000 | 75-79 / 30-34 | 4.410021 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / Whites / Males | 4 | 5.716722 | 75-79 / 30-34 | 4.822433 | 75-79 / 30-34 | 3.462687 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / Whites / Males | 5 | 3.739714 | 70-74 / 30-34 | 2.963653 | 70-74 / 30-34 | 2.582763 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / Whites / Males | 6 | 2.511656 | 70-74 / 30-34 | 1.882468 | 70-74 / 35-39 | 1.922146 | 0.060 | 0.060 | 65-69 / 35-39 |
| lung and bronchus / Whites / Males | 7 | 1.794982 | 65-69 / 35-39 | 1.298637 | 65-69 / 35-39 | 1.541545 | 0.060 | 0.060 | 55-59 / 35-39 |
| lung and bronchus / Whites / Males | 8 | 1.467828 | 55-59 / 80-84 | 1.139401 | 55-59 / 80-84 | 1.437103 | -0.015 | -0.015 | 55-59 / 80-84 |
| lung and bronchus / Whites / Males | 9 | 1.765966 | 30-34 / 80-84 | 1.315667 | 55-59 / 80-84 | 1.654856 | -0.060 | -0.060 | 55-59 / 80-84 |
| lung and bronchus / Whites / Males | 10 | 2.795727 | 30-34 / 80-84 | 1.795630 | 30-34 / 80-84 | 2.491379 | -0.060 | -0.060 | 30-34 / 80-84 |
| lung and bronchus / Whites / Females | 2 | 11.125501 | 75-79 / 30-34 | 10.520694 | 75-79 / 30-34 | 4.438054 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / Whites / Females | 3 | 7.210460 | 75-79 / 30-34 | 6.442583 | 75-79 / 30-34 | 3.762356 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / Whites / Females | 4 | 4.754094 | 70-74 / 30-34 | 3.997759 | 70-74 / 30-34 | 2.954150 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / Whites / Females | 5 | 3.190492 | 70-74 / 30-34 | 2.528405 | 70-74 / 30-34 | 2.203454 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / Whites / Females | 6 | 2.142790 | 70-74 / 30-34 | 1.614640 | 70-74 / 35-39 | 1.677235 | 0.060 | 0.060 | 55-59 / 35-39 |
| lung and bronchus / Whites / Females | 7 | 1.565526 | 55-59 / 35-39 | 1.123145 | 55-59 / 80-84 | 1.494009 | 0.030 | 0.030 | 55-59 / 80-84 |
| lung and bronchus / Whites / Females | 8 | 1.682846 | 50-54 / 80-84 | 1.296896 | 55-59 / 80-84 | 1.565656 | -0.060 | -0.060 | 55-59 / 80-84 |
| lung and bronchus / Whites / Females | 9 | 2.258810 | 30-34 / 80-84 | 1.555812 | 50-54 / 80-84 | 1.992410 | -0.060 | -0.060 | 30-34 / 80-84 |
| lung and bronchus / Whites / Females | 10 | 3.575958 | 30-34 / 80-84 | 2.296754 | 30-34 / 80-84 | 3.186671 | -0.060 | -0.060 | 30-34 / 80-84 |
| lung and bronchus / Blacks / Both Sexes | 2 | 9.860156 | 75-79 / 30-34 | 9.324136 | 75-79 / 30-34 | 3.990313 | 0.060 | 0.060 | 65-69 / 30-34 |
| lung and bronchus / Blacks / Both Sexes | 3 | 6.390388 | 75-79 / 30-34 | 5.709844 | 75-79 / 30-34 | 3.409775 | 0.060 | 0.060 | 65-69 / 30-34 |
| lung and bronchus / Blacks / Both Sexes | 4 | 4.256289 | 70-74 / 30-34 | 3.579151 | 70-74 / 30-34 | 2.742802 | 0.060 | 0.060 | 60-64 / 30-34 |
| lung and bronchus / Blacks / Both Sexes | 5 | 2.856414 | 70-74 / 30-34 | 2.263654 | 70-74 / 30-34 | 2.162544 | 0.060 | 0.060 | 55-59 / 30-34 |
| lung and bronchus / Blacks / Both Sexes | 6 | 2.026753 | 55-59 / 35-39 | 1.533988 | 60-64 / 35-39 | 1.757209 | 0.060 | 0.060 | 55-59 / 35-39 |
| lung and bronchus / Blacks / Both Sexes | 7 | 1.664506 | 55-59 / 40-44 | 1.335187 | 55-59 / 80-84 | 1.664506 | 0.000 | 0.000 | 55-59 / 40-44 |
| lung and bronchus / Blacks / Both Sexes | 8 | 1.986140 | 55-59 / 80-84 | 1.541741 | 55-59 / 80-84 | 1.861241 | -0.060 | -0.060 | 55-59 / 80-84 |
| lung and bronchus / Blacks / Both Sexes | 9 | 2.639729 | 30-34 / 80-84 | 1.780249 | 55-59 / 80-84 | 2.328404 | -0.060 | -0.060 | 30-34 / 80-84 |
| lung and bronchus / Blacks / Both Sexes | 10 | 4.178996 | 30-34 / 80-84 | 2.684071 | 30-34 / 80-84 | 3.724062 | -0.060 | -0.060 | 30-34 / 80-84 |
| lung and bronchus / Blacks / Males | 2 | 10.685146 | 75-79 / 30-34 | 10.104278 | 75-79 / 30-34 | 4.338635 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / Blacks / Males | 3 | 6.930657 | 70-74 / 30-34 | 6.187581 | 75-79 / 30-34 | 3.678074 | 0.060 | 0.060 | 70-74 / 30-34 |
| lung and bronchus / Blacks / Males | 4 | 4.647595 | 70-74 / 30-34 | 3.908204 | 70-74 / 30-34 | 2.923921 | 0.060 | 0.060 | 65-69 / 30-34 |
| lung and bronchus / Blacks / Males | 5 | 3.119021 | 70-74 / 30-34 | 2.471765 | 70-74 / 30-34 | 2.277979 | 0.060 | 0.060 | 60-64 / 30-34 |
| lung and bronchus / Blacks / Males | 6 | 2.149463 | 60-64 / 30-34 | 1.608648 | 60-64 / 40-44 | 1.785383 | 0.060 | 0.060 | 60-64 / 40-44 |
| lung and bronchus / Blacks / Males | 7 | 1.710494 | 60-64 / 40-44 | 1.261927 | 60-64 / 40-44 | 1.644055 | 0.030 | 0.030 | 55-59 / 40-44 |
| lung and bronchus / Blacks / Males | 8 | 1.828021 | 55-59 / 80-84 | 1.419001 | 55-59 / 80-84 | 1.713066 | -0.060 | -0.060 | 55-59 / 80-84 |
| lung and bronchus / Blacks / Males | 9 | 2.399706 | 30-34 / 80-84 | 1.638522 | 55-59 / 80-84 | 2.116689 | -0.060 | -0.060 | 30-34 / 80-84 |
| lung and bronchus / Blacks / Males | 10 | 3.799012 | 30-34 / 80-84 | 2.440016 | 30-34 / 80-84 | 3.385444 | -0.060 | -0.060 | 30-34 / 80-84 |
| lung and bronchus / Blacks / Females | 2 | 9.390061 | 75-79 / 30-34 | 8.879596 | 75-79 / 30-34 | 3.743153 | 0.060 | 0.060 | 65-69 / 30-34 |
| lung and bronchus / Blacks / Females | 3 | 6.085718 | 75-79 / 30-34 | 5.437620 | 75-79 / 30-34 | 3.230606 | 0.060 | 0.060 | 55-59 / 30-34 |
| lung and bronchus / Blacks / Females | 4 | 3.963754 | 70-74 / 30-34 | 3.333156 | 70-74 / 30-34 | 2.676508 | 0.060 | 0.060 | 55-59 / 30-34 |
| lung and bronchus / Blacks / Females | 5 | 2.665785 | 55-59 / 30-34 | 2.108073 | 70-74 / 30-34 | 2.145723 | 0.060 | 0.060 | 55-59 / 30-34 |
| lung and bronchus / Blacks / Females | 6 | 2.046176 | 55-59 / 35-39 | 1.541677 | 55-59 / 35-39 | 1.788651 | 0.060 | 0.060 | 55-59 / 80-84 |
| lung and bronchus / Blacks / Females | 7 | 1.730492 | 55-59 / 80-84 | 1.392122 | 55-59 / 80-84 | 1.689551 | -0.015 | -0.015 | 55-59 / 80-84 |
| lung and bronchus / Blacks / Females | 8 | 2.070834 | 55-59 / 80-84 | 1.607485 | 55-59 / 80-84 | 1.940609 | -0.060 | -0.060 | 55-59 / 80-84 |
| lung and bronchus / Blacks / Females | 9 | 2.759662 | 30-34 / 80-84 | 1.881190 | 50-54 / 80-84 | 2.434192 | -0.060 | -0.060 | 30-34 / 80-84 |
| lung and bronchus / Blacks / Females | 10 | 4.368864 | 30-34 / 80-84 | 2.806019 | 30-34 / 80-84 | 3.893260 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / All Races / Both Sexes | 2 | 6.654579 | 80-84 / 30-34 | 6.298789 | 80-84 / 30-34 | 2.464002 | 0.060 | 0.060 | 65-69 / 30-34 |
| pancreas / All Races / Both Sexes | 3 | 4.180197 | 80-84 / 30-34 | 3.742038 | 80-84 / 30-34 | 2.105523 | 0.060 | 0.060 | 65-69 / 30-34 |
| pancreas / All Races / Both Sexes | 4 | 2.679477 | 75-79 / 30-34 | 2.260316 | 75-79 / 30-34 | 1.677924 | 0.060 | 0.060 | 65-69 / 30-34 |
| pancreas / All Races / Both Sexes | 5 | 1.786676 | 70-74 / 35-39 | 1.449660 | 70-74 / 35-39 | 1.364172 | 0.060 | 0.060 | 60-64 / 35-39 |
| pancreas / All Races / Both Sexes | 6 | 1.333798 | 65-69 / 35-39 | 1.019741 | 65-69 / 35-39 | 1.245370 | 0.030 | 0.030 | 60-64 / 80-84 |
| pancreas / All Races / Both Sexes | 7 | 1.430677 | 30-34 / 80-84 | 1.098701 | 55-59 / 80-84 | 1.271340 | -0.060 | -0.060 | 55-59 / 80-84 |
| pancreas / All Races / Both Sexes | 8 | 2.268369 | 30-34 / 80-84 | 1.598693 | 30-34 / 80-84 | 1.976331 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / All Races / Both Sexes | 9 | 3.593786 | 30-34 / 80-84 | 2.416996 | 30-34 / 80-84 | 3.169941 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / All Races / Both Sexes | 10 | 5.689378 | 30-34 / 80-84 | 3.654154 | 30-34 / 80-84 | 5.070021 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / All Races / Males | 2 | 7.516043 | 80-84 / 30-34 | 7.114195 | 80-84 / 30-34 | 2.844816 | 0.060 | 0.060 | 65-69 / 30-34 |
| pancreas / All Races / Males | 3 | 4.721341 | 80-84 / 30-34 | 4.226461 | 80-84 / 30-34 | 2.430933 | 0.060 | 0.060 | 65-69 / 30-34 |
| pancreas / All Races / Males | 4 | 3.004356 | 75-79 / 30-34 | 2.534373 | 75-79 / 30-34 | 1.937249 | 0.060 | 0.060 | 65-69 / 30-34 |
| pancreas / All Races / Males | 5 | 2.025689 | 65-69 / 30-34 | 1.597713 | 65-69 / 30-34 | 1.494051 | 0.060 | 0.060 | 60-64 / 30-34 |
| pancreas / All Races / Males | 6 | 1.409805 | 65-69 / 30-34 | 1.069651 | 65-69 / 35-39 | 1.297958 | 0.035 | 0.035 | 60-64 / 80-84 |
| pancreas / All Races / Males | 7 | 1.389806 | 55-59 / 80-84 | 1.118052 | 55-59 / 80-84 | 1.300104 | -0.060 | -0.060 | 60-64 / 80-84 |
| pancreas / All Races / Males | 8 | 1.990560 | 30-34 / 80-84 | 1.402899 | 30-34 / 80-84 | 1.734288 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / All Races / Males | 9 | 3.153652 | 30-34 / 80-84 | 2.120984 | 30-34 / 80-84 | 2.781716 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / All Races / Males | 10 | 4.992595 | 30-34 / 80-84 | 3.206627 | 30-34 / 80-84 | 4.449091 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / All Races / Females | 2 | 6.355342 | 80-84 / 30-34 | 6.015551 | 80-84 / 30-34 | 2.272422 | 0.060 | 0.060 | 70-74 / 30-34 |
| pancreas / All Races / Females | 3 | 3.992226 | 80-84 / 30-34 | 3.573770 | 80-84 / 30-34 | 1.937346 | 0.060 | 0.060 | 65-69 / 30-34 |
| pancreas / All Races / Females | 4 | 2.561365 | 75-79 / 30-34 | 2.160681 | 75-79 / 30-34 | 1.543902 | 0.060 | 0.060 | 65-69 / 30-34 |
| pancreas / All Races / Females | 5 | 1.692386 | 75-79 / 35-39 | 1.378876 | 75-79 / 35-39 | 1.245317 | 0.060 | 0.060 | 65-69 / 35-39 |
| pancreas / All Races / Females | 6 | 1.239332 | 65-69 / 40-44 | 1.000000 | 65-69 / 40-44 | 1.174503 | 0.025 | 0.025 | 65-69 / 80-84 |
| pancreas / All Races / Females | 7 | 1.497964 | 30-34 / 80-84 | 1.107169 | 30-34 / 80-84 | 1.285463 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / All Races / Females | 8 | 2.375054 | 30-34 / 80-84 | 1.673882 | 30-34 / 80-84 | 2.069281 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / All Races / Females | 9 | 3.762808 | 30-34 / 80-84 | 2.530671 | 30-34 / 80-84 | 3.319028 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / All Races / Females | 10 | 5.956959 | 30-34 / 80-84 | 3.826015 | 30-34 / 80-84 | 5.308472 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / Whites / Both Sexes | 2 | 7.113936 | 80-84 / 30-34 | 6.733586 | 80-84 / 30-34 | 2.642090 | 0.060 | 0.060 | 65-69 / 30-34 |
| pancreas / Whites / Both Sexes | 3 | 4.468750 | 80-84 / 30-34 | 4.000346 | 80-84 / 30-34 | 2.257701 | 0.060 | 0.060 | 65-69 / 30-34 |
| pancreas / Whites / Both Sexes | 4 | 2.861168 | 75-79 / 30-34 | 2.413585 | 75-79 / 30-34 | 1.799197 | 0.060 | 0.060 | 65-69 / 30-34 |
| pancreas / Whites / Both Sexes | 5 | 1.881335 | 65-69 / 30-34 | 1.487799 | 70-74 / 30-34 | 1.370115 | 0.060 | 0.060 | 65-69 / 30-34 |
| pancreas / Whites / Both Sexes | 6 | 1.343438 | 65-69 / 35-39 | 1.027111 | 65-69 / 35-39 | 1.235369 | 0.030 | 0.030 | 65-69 / 80-84 |
| pancreas / Whites / Both Sexes | 7 | 1.353442 | 55-59 / 80-84 | 1.088798 | 55-59 / 80-84 | 1.259880 | -0.060 | -0.060 | 55-59 / 80-84 |
| pancreas / Whites / Both Sexes | 8 | 2.102960 | 30-34 / 80-84 | 1.482117 | 30-34 / 80-84 | 1.832218 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / Whites / Both Sexes | 9 | 3.331729 | 30-34 / 80-84 | 2.240750 | 30-34 / 80-84 | 2.938790 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / Whites / Both Sexes | 10 | 5.274511 | 30-34 / 80-84 | 3.387695 | 30-34 / 80-84 | 4.700317 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / Whites / Males | 2 | 8.174239 | 80-84 / 30-34 | 7.737199 | 80-84 / 30-34 | 3.074611 | 0.060 | 0.060 | 65-69 / 30-34 |
| pancreas / Whites / Males | 3 | 5.134799 | 80-84 / 30-34 | 4.596581 | 80-84 / 30-34 | 2.627296 | 0.060 | 0.060 | 65-69 / 30-34 |
| pancreas / Whites / Males | 4 | 3.260585 | 75-79 / 30-34 | 2.750519 | 75-79 / 30-34 | 2.093733 | 0.060 | 0.060 | 65-69 / 30-34 |
| pancreas / Whites / Males | 5 | 2.189317 | 65-69 / 30-34 | 1.728387 | 70-74 / 30-34 | 1.594453 | 0.060 | 0.060 | 60-64 / 30-34 |
| pancreas / Whites / Males | 6 | 1.523684 | 65-69 / 30-34 | 1.130437 | 65-69 / 30-34 | 1.316419 | 0.045 | 0.045 | 60-64 / 80-84 |
| pancreas / Whites / Males | 7 | 1.366955 | 55-59 / 80-84 | 1.099669 | 55-59 / 80-84 | 1.275768 | -0.060 | -0.060 | 60-64 / 80-84 |
| pancreas / Whites / Males | 8 | 1.808509 | 30-34 / 80-84 | 1.305346 | 50-54 / 80-84 | 1.575675 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / Whites / Males | 9 | 2.865229 | 30-34 / 80-84 | 1.927005 | 30-34 / 80-84 | 2.527308 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / Whites / Males | 10 | 4.535987 | 30-34 / 80-84 | 2.913358 | 30-34 / 80-84 | 4.042190 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / Whites / Females | 2 | 6.738552 | 80-84 / 30-34 | 6.378272 | 80-84 / 30-34 | 2.422458 | 0.060 | 0.060 | 65-69 / 30-34 |
| pancreas / Whites / Females | 3 | 4.232946 | 80-84 / 30-34 | 3.789258 | 80-84 / 30-34 | 2.070022 | 0.060 | 0.060 | 65-69 / 30-34 |
| pancreas / Whites / Females | 4 | 2.714069 | 75-79 / 30-34 | 2.289497 | 75-79 / 30-34 | 1.649633 | 0.060 | 0.060 | 65-69 / 30-34 |
| pancreas / Whites / Females | 5 | 1.761819 | 75-79 / 30-34 | 1.402025 | 75-79 / 30-34 | 1.256219 | 0.060 | 0.060 | 65-69 / 30-34 |
| pancreas / Whites / Females | 6 | 1.231760 | 65-69 / 35-39 | 1.000000 | 65-69 / 80-84 | 1.172564 | 0.020 | 0.020 | 65-69 / 80-84 |
| pancreas / Whites / Females | 7 | 1.400157 | 30-34 / 80-84 | 1.052071 | 55-59 / 80-84 | 1.217382 | -0.060 | -0.060 | 55-59 / 80-84 |
| pancreas / Whites / Females | 8 | 2.219978 | 30-34 / 80-84 | 1.564588 | 30-34 / 80-84 | 1.934170 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / Whites / Females | 9 | 3.517121 | 30-34 / 80-84 | 2.365434 | 30-34 / 80-84 | 3.102317 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / Whites / Females | 10 | 5.568008 | 30-34 / 80-84 | 3.576201 | 30-34 / 80-84 | 4.961863 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / Blacks / Both Sexes | 2 | 5.906719 | 80-84 / 30-34 | 5.590913 | 80-84 / 30-34 | 2.369622 | 0.060 | 0.060 | 60-64 / 30-34 |
| pancreas / Blacks / Both Sexes | 3 | 3.765074 | 75-79 / 30-34 | 3.364113 | 75-79 / 30-34 | 2.045276 | 0.060 | 0.060 | 60-64 / 30-34 |
| pancreas / Blacks / Both Sexes | 4 | 2.442120 | 75-79 / 30-34 | 2.060090 | 75-79 / 30-34 | 1.660771 | 0.060 | 0.060 | 60-64 / 35-39 |
| pancreas / Blacks / Both Sexes | 5 | 1.723653 | 60-64 / 35-39 | 1.384293 | 60-64 / 35-39 | 1.427577 | 0.060 | 0.060 | 60-64 / 80-84 |
| pancreas / Blacks / Both Sexes | 6 | 1.393231 | 60-64 / 40-44 | 1.109644 | 60-64 / 80-84 | 1.358507 | 0.015 | 0.015 | 60-64 / 80-84 |
| pancreas / Blacks / Both Sexes | 7 | 1.639663 | 30-34 / 80-84 | 1.231040 | 60-64 / 80-84 | 1.432692 | -0.060 | -0.060 | 60-64 / 80-84 |
| pancreas / Blacks / Both Sexes | 8 | 2.599720 | 30-34 / 80-84 | 1.832221 | 30-34 / 80-84 | 2.265022 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / Blacks / Both Sexes | 9 | 4.118747 | 30-34 / 80-84 | 2.770058 | 30-34 / 80-84 | 3.632989 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / Blacks / Both Sexes | 10 | 6.520452 | 30-34 / 80-84 | 4.187933 | 30-34 / 80-84 | 5.810622 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / Blacks / Males | 2 | 6.151567 | 75-79 / 30-34 | 5.819411 | 80-84 / 30-34 | 2.716844 | 0.060 | 0.060 | 60-64 / 30-34 |
| pancreas / Blacks / Males | 3 | 4.009074 | 70-74 / 30-34 | 3.574513 | 70-74 / 30-34 | 2.344972 | 0.060 | 0.060 | 60-64 / 30-34 |
| pancreas / Blacks / Males | 4 | 2.703615 | 65-69 / 30-34 | 2.265306 | 65-69 / 30-34 | 1.902139 | 0.060 | 0.060 | 60-64 / 30-34 |
| pancreas / Blacks / Males | 5 | 1.935716 | 60-64 / 30-34 | 1.518408 | 60-64 / 30-34 | 1.530447 | 0.055 | 0.055 | 60-64 / 80-84 |
| pancreas / Blacks / Males | 6 | 1.462279 | 60-64 / 40-44 | 1.222280 | 60-64 / 80-84 | 1.462279 | 0.000 | 0.000 | 60-64 / 40-44 |
| pancreas / Blacks / Males | 7 | 1.669992 | 60-64 / 80-84 | 1.355998 | 60-64 / 80-84 | 1.578119 | -0.060 | -0.060 | 60-64 / 80-84 |
| pancreas / Blacks / Males | 8 | 2.486276 | 30-34 / 80-84 | 1.752268 | 30-34 / 80-84 | 2.166184 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / Blacks / Males | 9 | 3.939018 | 30-34 / 80-84 | 2.649181 | 30-34 / 80-84 | 3.474457 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / Blacks / Males | 10 | 6.235920 | 30-34 / 80-84 | 4.005185 | 30-34 / 80-84 | 5.557065 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / Blacks / Females | 2 | 5.682189 | 80-84 / 30-34 | 5.378388 | 80-84 / 30-34 | 2.077861 | 0.060 | 0.060 | 60-64 / 30-34 |
| pancreas / Blacks / Females | 3 | 3.581947 | 75-79 / 30-34 | 3.200488 | 75-79 / 30-34 | 1.822418 | 0.060 | 0.060 | 60-64 / 35-39 |
| pancreas / Blacks / Females | 4 | 2.365296 | 75-79 / 35-39 | 2.030118 | 75-79 / 35-39 | 1.570376 | 0.060 | 0.060 | 60-64 / 40-44 |
| pancreas / Blacks / Females | 5 | 1.674750 | 75-79 / 40-44 | 1.388894 | 75-79 / 40-44 | 1.366744 | 0.060 | 0.060 | 60-64 / 40-44 |
| pancreas / Blacks / Females | 6 | 1.357306 | 60-64 / 40-44 | 1.054229 | 60-64 / 40-44 | 1.284923 | 0.030 | 0.030 | 60-64 / 80-84 |
| pancreas / Blacks / Females | 7 | 1.710940 | 30-34 / 80-84 | 1.264582 | 30-34 / 80-84 | 1.468226 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / Blacks / Females | 8 | 2.712730 | 30-34 / 80-84 | 1.911868 | 30-34 / 80-84 | 2.363483 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / Blacks / Females | 9 | 4.297790 | 30-34 / 80-84 | 2.890473 | 30-34 / 80-84 | 3.790916 | -0.060 | -0.060 | 30-34 / 80-84 |
| pancreas / Blacks / Females | 10 | 6.803897 | 30-34 / 80-84 | 4.369984 | 30-34 / 80-84 | 6.063211 | -0.060 | -0.060 | 30-34 / 80-84 |

Joint-grid contender sets and exact-power ties:

| Curve | Pure-power k contenders | Joint (k, g/year) contenders |
|---|---|---|
| female breast / All Races / Females | 4 | (4, -0.015) |
| female breast / Whites / Females | 4 | (4, -0.010) |
| female breast / Blacks / Females | 4 | (4, -0.030) |
| colon and rectum / All Races / Both Sexes | 5 | (5, -0.010) |
| colon and rectum / All Races / Males | 5 | (5, 0.005) |
| colon and rectum / All Races / Females | 5 | (5, -0.030) |
| colon and rectum / Whites / Both Sexes | 5 | (5, -0.005) |
| colon and rectum / Whites / Males | 5 | (5, 0.005) |
| colon and rectum / Whites / Females | 5 | (5, -0.045) |
| colon and rectum / Blacks / Both Sexes | 5 | (5, -0.005) |
| colon and rectum / Blacks / Males | 5 | (5, 0.020) |
| colon and rectum / Blacks / Females | 5 | (5, -0.025) |
| lung and bronchus / All Races / Both Sexes | 7 | (7, 0.040) |
| lung and bronchus / All Races / Males | 8 | (8, -0.060) |
| lung and bronchus / All Races / Females | 7 | (7, 0.020) |
| lung and bronchus / Whites / Both Sexes | 8 | (8, -0.060) |
| lung and bronchus / Whites / Males | 8 | (8, -0.015) |
| lung and bronchus / Whites / Females | 7 | (7, 0.030) |
| lung and bronchus / Blacks / Both Sexes | 7 | (7, 0.000) |
| lung and bronchus / Blacks / Males | 7 | (7, 0.030) |
| lung and bronchus / Blacks / Females | 7 | (7, -0.015) |
| pancreas / All Races / Both Sexes | 6 | (6, 0.030) |
| pancreas / All Races / Males | 7 | (6, 0.035) |
| pancreas / All Races / Females | 6 | (6, 0.025) |
| pancreas / Whites / Both Sexes | 6 | (6, 0.030) |
| pancreas / Whites / Males | 7 | (7, -0.060) |
| pancreas / Whites / Females | 6 | (6, 0.020) |
| pancreas / Blacks / Both Sexes | 6 | (6, 0.015) |
| pancreas / Blacks / Males | 6 | (6, 0.000) |
| pancreas / Blacks / Females | 6 | (6, 0.030) |

</details>

### Content integrity

SHA-256 after UTF-8 universal-newline normalization; these establish local content integrity, not upstream source truth. Relative path keys use `Path.as_posix()`.

| Input | Normalized SHA-256 |
|---|---|
| `PROTOCOL.md` | `daa9c8b3d838d79043c253cc9c16d679ec3dabfb39938ae16d3ce3cc786ef822` |
| `data/colon_and_rectum.csv` | `7180847c206fff5a97734fc6999817efc9ed5d8f3150edc77d72079511a35358` |
| `data/female_breast.csv` | `dd23cb15fd17063a02451901e63a68937328492c8987bce30c48429be162361d` |
| `data/lung_and_bronchus.csv` | `c94ba12c47034b633f6b7b2a9658390599caeedf7b9451b2126cb4b744f6ebba` |
| `data/pancreas.csv` | `aa2d0ae30e43dcd6e83bfda57f4e1b6a04d47b07f49f96ce57eb19a87d446cef` |

Unmodified inherited primitive file `research/multistage-incidence-identifiability/verifier.py`: `670ab2a88124f42633105b063a0d945910eb7f63cc92fbc85b1935fe0c11d66e`.

<!-- END RECOMPUTED RESULTS -->

## What the result establishes

The female-breast failure remains within both source-defined race-specific
columns. Pooling across those race categories is insufficient by itself to
explain away the observed shape failure. The age audit is more informative
than a single threshold: the White and All Races breast curves already fail
before the oldest retained bin is added, whereas the Black breast curve's
failure first appears when that bin enters. Removing the oldest bin therefore
does not rescue all breast curves. The generated minimum-deletion table
retains every tied smallest set and makes this difference explicit.

The contrast findings are not uniformly negative. All prespecified colon and
rectum and pancreas curves pass the necessary monotonicity condition. The
All Races and White lung/bronchus curves pass it, while the Black lung/bronchus
curves fail for each printed sex column. Thus a passing pooled curve need not
imply passing race-specific curves. There is no weighting information here
with which to attribute either pattern to demographic mixing.

The shape-family frontiers answer a different question. After division by a
specified model shape, an earlier or middle-age bin can impose the largest
required amplitude while another bin imposes the smallest permitted amplitude.
Those competing constraints, reported as `L / U ages`, are not necessarily the
older-age decline that determines the abstract monotonicity frontier. A model
family may require substantial calibration even when the necessary monotone
condition passes. The fitted grid index and net-growth parameter are shape
descriptors under the frozen limiting family, not biological estimates.

## Frozen protocol and chronology

[PROTOCOL.md](PROTOCOL.md) was committed and published alone in
[`f06822e9f2427cbec8fef16c48589290561b40c3`](https://github.com/dovvnloading/open-cancer/commit/f06822e9f2427cbec8fef16c48589290561b40c3)
before the first new monotonicity, calibration, age-deletion or grid calculation.
Its parent is the then-current `main`,
[`88c7e07668887f6e941a3fe7592d53002199669a`](https://github.com/dovvnloading/open-cancer/commit/88c7e07668887f6e941a3fe7592d53002199669a),
which includes PR #70. The protocol content has not changed. There were no
post-computation changes to the panel, rounding convention, grid, sensitivity
checks, or verdict rule. The separate protocol commit is an auditable sequence,
not independent registration or a blinded study: the published breast result
and public source tables were known beforehand.

The frozen primary rule requires complete All Races, Whites and Blacks breast
curves. Given a reproduced All Races failure, failure in either race-specific
curve means persistence, with the breadth stated explicitly. Passing both
means disappearance of this particular shape failure; missing required
inputs or an unreproduced All Races benchmark means the stated question is
untestable. Site/sex contrasts do not determine or override that primary rule.

## Source, extraction and reuse audit

The sole numeric source is Howlader N and colleagues, editors, **SEER Cancer
Statistics Review, 1975-2017**, National Cancer Institute, Bethesda, Maryland;
November 2019 SEER submission, released April 15, 2020. We use the
**Age-Specific SEER Incidence Rates, 2013-2017** subsection for **SEER 21
areas**, invasive disease, ages **30-34 through 80-84**. The period and release
are deliberately historical; this is not a statement of current incidence.

| Selected CSV | Source table | Source race/sex columns | Official CSV | PDF, page 1 | Public HTML |
|---|---|---|---|---|---|
| [female_breast.csv](data/female_breast.csv) | 4.11 | All Races / Whites / Blacks; Females | [CSV](https://seer.cancer.gov/archive/csr/1975_2017/download_csr_datafile.php/sect_04_table.11.csv) | [PDF](https://seer.cancer.gov/archive/csr/1975_2017/results_single/sect_04_table.11.pdf) | [HTML](https://seer.cancer.gov/archive/csr/1975_2017/browse_csr.php?pageSEL=sect_04_table.11&sectionSEL=4) |
| [colon_and_rectum.csv](data/colon_and_rectum.csv) | 6.10 | All Races / Whites / Blacks; Both Sexes / Males / Females | [CSV](https://seer.cancer.gov/archive/csr/1975_2017/download_csr_datafile.php/sect_06_table.10.csv) | [PDF](https://seer.cancer.gov/archive/csr/1975_2017/results_single/sect_06_table.10.pdf) | [HTML](https://seer.cancer.gov/archive/csr/1975_2017/browse_csr.php?pageSEL=sect_06_table.10&sectionSEL=6) |
| [lung_and_bronchus.csv](data/lung_and_bronchus.csv) | 15.10 | All Races / Whites / Blacks; Both Sexes / Males / Females | [CSV](https://seer.cancer.gov/archive/csr/1975_2017/download_csr_datafile.php/sect_15_table.10.csv) | [PDF](https://seer.cancer.gov/archive/csr/1975_2017/results_single/sect_15_table.10.pdf) | [HTML](https://seer.cancer.gov/archive/csr/1975_2017/browse_csr.php?pageSEL=sect_15_table.10&sectionSEL=15) |
| [pancreas.csv](data/pancreas.csv) | 22.7 | All Races / Whites / Blacks; Both Sexes / Males / Females | [CSV](https://seer.cancer.gov/archive/csr/1975_2017/download_csr_datafile.php/sect_22_table.07.csv) | [PDF](https://seer.cancer.gov/archive/csr/1975_2017/results_single/sect_22_table.07.pdf) | [HTML](https://seer.cancer.gov/archive/csr/1975_2017/browse_csr.php?pageSEL=sect_22_table.07&sectionSEL=22) |

Every data row records the printed rate, source/version, table, subsection,
column, PDF page, URLs, access date **2026-11-18**, units, display precision,
and reuse basis. `case_count` and `person_years` are deliberately blank.
These selected tables do not print the age-by-site-by-sex-by-race counts or
exposures required for a sampling analysis. A blank is not zero, and a rounded
rate cannot uniquely reconstruct a denominator. No data were joined from
another release, graph, registry record, or patient source.

The [archive's Suggested Citation](https://seer.cancer.gov/archive/csr/1975_2017/)
expressly places the report in the public domain and permits reproduction,
requesting attribution. This exact-source permission is recorded on each row.
The repository's MIT license does not relicense the NCI source material.
Only the small selected aggregates are committed; there are no person-level
records, raw genomes, or controlled-access inputs. This follows
[DATA_POLICY.md](../../DATA_POLICY.md).

**Acquisition limitation, resolved without a browser download:** the public
web text reader rejected the official CSV response's
`text/comma-separated-values` content type. Ordinary public HTTP requests
read the official CSVs successfully, and those CSVs supplied the rate strings.
No browser file download was used and no download prompt appeared. Public
HTML and the PDFs were readable. It would be incorrect to describe the CSVs
as unavailable or the primary data as manually reconstructed from a plot.

[extract_sources.py](extract_sources.py) is an offline selector that checks
source title, table, subsection, column order, ages, and decimal strings before
writing the selected long-format CSVs. It also supports an independent PDF-text
comparison. Every selected rate string was checked against its PDF cell with
no mismatch; every one-page PDF was visually inspected to confirm headings,
row/column alignment and incidence-versus-mortality boundaries. The PDFs label
the combined-sex column `Total`; the CSVs label it `Both Sexes`. We preserve
the CSV labels. The table footnote suppresses statistics based on fewer than
16 cases; none of the selected cells is suppressed. This does not reveal the
actual case counts. CSV/PDF agreement checks transcription between two
representations of one source, not the truth of the registry estimates or
independent epidemiological replication.

## Methods and exact inherited statistics

These are direct applications of PR #64, Theorems 7 and 8. Under constant
Poisson initiation of independent identically distributed clones, the direct
first-event hazard is `h(t)=lambda*F(t)` for a clone completion CDF. An
independent **clone-level** diagnosis lag changes the CDF but preserves its
monotonicity. A nondecreasing hazard has ordered means under any nonnegative
normalized age weights supported on ordered nonoverlapping bins. This is the
homogeneous observation model being challenged, not a universal law about
published cancer incidence.

For a printed rate `y_j`, set `l_j=y_j-0.05` and `u_j=y_j+0.05` using exact
rational arithmetic. Closed endpoints are conservative for the assumed
round-to-nearest display convention. The archive's exact rounding algorithm
has not been established. The intervals represent rounding sensitivity only.

The complete ordered-pair check and sharp abstract relaxation are

\[
  l_i>u_j\quad(i<j),\qquad
  \rho_{\mathrm{mono}}^2=
  \max\left\{1,\max_{i\le j}\frac{l_i}{u_j}\right\}.
\]

At the threshold, `z_j=max_{i<=j} l_i/rho` is a feasible nondecreasing vector
inside `[l_j/rho,rho*u_j]`. Below a nontrivial threshold, a maximizing ordered
pair proves impossibility. Sufficiency applies only to arbitrary positive
nondecreasing vectors; such a vector is not necessarily realizable by a
particular branching process. Passing is therefore only a necessary-condition
pass. No empirical rate is assigned a significance level.

For a specified positive shape vector `f_j`, free common amplitude `A`, and
age-specific calibration multipliers `a_j in [1/rho,rho]`, define

\[
  L_f=\max_j l_j/f_j,\quad U_f=\min_j u_j/f_j,\qquad
  \rho_f^2=\max\{1,L_f/U_f\}.
\]

The amplitude interval is `[L_f/rho,rho*U_f]`. This is a sharp compatibility
frontier under the stated allowance, not estimated registry incompleteness.
The constraining ages are all maximizers of `l/f` and minimizers of `u/f`.
The symmetric factor is multiplicative, after fitting amplitude; it is not an
additive error percentage or a demographic risk ratio.

The frozen shape grid is exactly #64's: `k=2,...,10`, pure powers
`t^(k-1)/(k-1)!`, and one-net-growth-pole limiting shapes with
`g=-0.060,...,+0.060` per year at spacing `0.005`. Uniform means on each
half-open bin are used for this comparison. The inherited rational Taylor
enclosures use indices zero through 70 and an explicit absolute tail bound;
the implementation compares bounds before formatting any decimal output.
The code retains tied or unresolved contenders. Its thresholds are finite-grid
minima, not continuous optima, priors, full-MSCE optima, or estimates of
biological stage counts. Several nonzero growth rates and off-grid parameters
are outside this declared comparison. Negative net growth is loss, and the
limiting family is distinct from a fitted finite-rate branching process.

For pure powers, the prespecified unknown-weight sensitivity permits each bin
shape over its endpoint closure. It optimizes a nuisance allowance, not
numerical error. The result is a sharp **infimum** over that closure; the
excluded upper endpoint of a half-open bin need not be attained by actual age
weights. Changing weights cannot remove an ordered-bin decline. The abstract
monotonicity test already allows arbitrary normalized nonnegative weights.

## Falsification ledger

| Attempt to defeat or narrow the result | Outcome and its limit |
|---|---|
| Reproduce the original All Races breast result from official source CSVs | The original monotonicity threshold and declared-grid minima reproduce. This is a computational replication using the same historical source, not independent epidemiological evidence. |
| Allow half-last-unit rounding instead of using printed point values exactly | The breast failures remain; both point and rounded thresholds are retained. Rounding intervals cannot stand in for sampling uncertainty. |
| Separate the source-defined White and Black female populations | Both retain the shape failure. Broad within-stratum heterogeneity remains, so this does not establish a race-specific mechanism. |
| Apply all frozen race-by-sex contrasts at the other sites | Colon/rectum and pancreas pass the necessary condition; Black lung/bronchus curves fail while the other printed lung strata pass. A universal all-site failure claim is contradicted. |
| Truncate later ages, remove each single age, and enumerate every smallest rescue set | Failure depends on retained ages, with different rescue requirements across breast strata. All results are reported; deletion is a diagnostic and does not replace the primary panel. |
| Permit arbitrary within-bin age weights | This can reduce a specific pure-power threshold but cannot cure a genuine ordered-bin decrease. Endpoint results are closure infima. |
| Add the frozen one-growth grid | The grid includes zero growth and can improve specific-shape compatibility. Its curves remain nondecreasing; it cannot explain a rounding-robust decline without calibration. This does not exhaust biological alternatives. |
| Allow unrestricted age-specific multipliers | Every positive finite shape can fit after amplitude/multiplier adjustment. This removes identifiability by construction; it is not a successful causal explanation or evidence for a selected k. |
| Treat all source strata as independent replications or reconstruct pooling weights | Rejected: All Races overlaps its race strata, Both Sexes overlaps sex strata, and the table supplies no mixing denominators. No independent-study count or mixture decomposition is reported. |
| Demand a calibrated statistical rejection or an identified cause | Untestable with these selected aggregates: age-cell counts, person-years and a defensible sampling/observation model are absent. No artificial sample size, p-value, confidence interval or causal fraction is supplied. |

## What this cannot say about cause

The direct conclusion is a failure of a conjunction: homogeneous constant
initiation/clone laws, the stated observation bridge, and interpretation of the
rounding intervals as possible mean rates. The exercise does not identify
which part fails. Sampling fluctuation can produce a nonmonotone estimated
rate even when the true hazard is monotone. Period rates also combine people
born in different years; changes across age can reflect cohort composition,
exposure histories, detection or screening, subtype composition, selection,
competing mortality, and other within-stratum heterogeneity. Those are
possible explanations, not effects estimated by this analysis. In particular,
screening cessation, cohort effects, and diagnosis delays are **not identified
causes** here: these aggregate rates do not distinguish them, establish that
they occurred, or estimate their contributions.

The broad source race and sex labels do not establish biological homogeneity.
Differences between their calibration factors are descriptive differences in
curve shape, not statistically established demographic effects. Nothing here
supplies missing denominators, separates biological from observation effects,
estimates mutation/driver counts, measures reporting bias, or identifies a
preventable fraction. The result neither recommends screening changes nor
supports diagnosis, treatment, surveillance, or patient selection.

## Novelty assessment

**A reproducible extension and constraint-age audit; no new theorem or
biological discovery.** The monotonicity argument, interval-intersection
frontier, shape families, calibration interpretation and driver-count warning
are inherited from [#64](../multistage-incidence-identifiability/). The
elementary optimization formulas and their prior mathematical foundations are
not new here. [#70](../stem-cell-divisions-risk/) is a methodological predecessor
for preserving compatibility/causation distinctions and portable reproduction,
not a new source of incidence data.

The contribution to this repository is the frozen extension to the specified
published race/sex/site curves, complete limiting-age attribution, all tied
minimum rescue sets, retained passing contrasts, official-CSV/PDF audit, and
offline recomputation. The printed patterns predate this analysis. We do not
claim the first observation of an older-age decline, an independent population
replication, global bibliographic priority, or an explanation of cancer's
causes. No new literature-wide priority search was needed to support that
deliberately limited contribution claim.

## Reproduction and hostile review

From the repository root, with standard-library Python only:

```bash
python -B research/incidence-strata-test/analysis.py --check
python -B -m unittest discover -s tests -p 'test_research_incidence_strata_test.py' -v
```

Running `analysis.py` without arguments emits the complete result object,
including exact rational squared thresholds and grid enclosures, as JSON on
stdout. `--write-readme` replaces only the generated results block. The tests
recompute the entire block from CSVs and independently check the arithmetic;
no cached numerical results or network requests are required. Input hashes
normalize CRLF and CR-only newlines to LF, preserve all other text, and use
portable POSIX path keys. The frozen protocol is also hash-checked.

To repeat source extraction with separately acquired official source CSVs and
`pdftotext -layout` outputs named by the table stems:

```bash
python -B research/incidence-strata-test/extract_sources.py \
  --source-dir /path/to/nci-source-files \
  --output-dir research/incidence-strata-test/data \
  --pdf-text-dir /path/to/nci-source-files --check
```

The offline source check requires those original documents to be acquired
separately. The ordinary test suite proves consistency with the checked-in
aggregates; it cannot authenticate an upstream source. The original documents
are intentionally not redistributed in the repository.

The hostile review checked every extracted rate against the source, complete
age and column selection, source reuse, nonadjacent/equality edge cases,
monotonicity feasibility at and below the exact threshold, model extrema,
grid nesting and tail bounds, all age-deletion outcomes, missing-cell
abstention, newline/path portability, allowed file scope and interpretation.
Independent arithmetic checks and passed CI assess the implementation, not
clinical validity, external scientific review, source truth or priority.

Repository gates are unchanged and the existing PR workflow discovers the
new test. No existing research note, research index, policy, catalog, or
workflow file is modified by this contribution.
