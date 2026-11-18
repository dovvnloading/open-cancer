# Persister-colony resistance-origin reanalysis

**Nonclinical computational research. Unreviewed. No treatment advice.**

## Completion status and pre-fit deviation record

**The requested raw Mendeley V2 reanalysis is blocked and incomplete.** The
public dataset's licence and advertised archive checksum were inspected, but
the archive could not be downloaded in this environment. Its bytes have not
been verified, its internal raw-data locators have not been inspected, and no
raw population trajectory has been recalibrated. A successful offline test or
CI run must not be interpreted as completing those missing scientific steps.

The protocol was frozen in commit `cb76f03` before any empirical fitting. Its
normalized SHA-256 is
`a2a2b44d24977626c5cb395e8832e45378d0db2f9e898c9db23a2cda5b6b35e2`.
The protocol itself remains unchanged.

**Deviation D1, recorded before fitting:** perform a separately labeled,
partial **published-count analysis** using the exact annotated Fig. 3e values
for DiFi B6 and WiDr B7. These are real experimental aggregate observations,
not synthetic data, but they are **not a fetched raw-deposit dataset**. They
are literal printed integers, not pixel-height estimates. The fetched
supplementary workbook supplies context and published model-dependent
sensitivity anchors only. Do not substitute its fitted parameters for raw
population measurements or treat Table 5's inferred rates as count outcomes.
The additional clones cannot enter a raw-data fit without their source data.

The planned whole-biological-replicate scoring, competing observation laws,
priors, sensitivity choices, and verdict rule remain those in
[PROTOCOL.md](PROTOCOL.md). This deviation is a restricted deliverable while
the archive is inaccessible; it does not assert that the user's requested
raw-data job is finished. The eventual PR is to remain draft for that reason.

The acquisition failures were an HTTP 403 / Cloudflare 1010 response at the
documented Mendeley public download and API endpoints, and timeouts without a
saved-file result from supported browser download operations. No access
control was bypassed. [fetch_sources.py](fetch_sources.py) accepts an official
locally downloaded archive and checks its entire advertised digest before
listing candidate data members. A source-access success must precede the
remaining extraction, raw-data cross-check, and population calibration work.

## Pre-fit interpretation constraints

The allowed families can have exactly the same observable occupancy laws:
late visibility may arise from an already competent but delayed lineage, or
from competence acquired after exposure. A fit of that common law cannot
identify the latent allocation. The relevant exact obstruction is the
vector-readout extension of PR #57, Theorem 5.3, **“Induced clonal
amplification”**, with the singleton special case in Theorem 5.1,
**“Different unknown clone laws”**. This is a structural statement conditional
on the broad kernel family, not empirical proof of either mechanism.

The independent-holiday positive identification theorem in PR #57 is not
applicable to repeated readouts during continuous drug exposure. Drug-caused
excess acquisition is a further causal question; Theorem 5.6, **“Drug causation
needs an exposure control”**, remains relevant even if event timing becomes
identifiable.

Numerical results, source-audit details, falsification results, follow-up
design, and novelty assessment will be added after this deviation record is
committed and the partial published-count calculations are run.
