# Prospective source-correction addendum

**Nonclinical computational research; no medical advice.**

Date: **2026-11-17, America/New_York**. This addendum was written and published
after the protocol-only commit and before the first fit or inverse-rate
calculation on the empirical data. No empirical results motivated it.

The continuing literature/source audit identified a specific later correction
that the original protocol did not name: Nunney and Thai (2020), Methods,
the paragraph correcting the total-stem-cell-division formula,
[doi:10.1098/rspb.2020.2291](https://doi.org/10.1098/rspb.2020.2291), replace
the TV developmental-count convention `N*(2+L)-2` by `N*(1+L)-1`. Their reason
is that binary growth from one cell to N cells uses N-1 division events.
This is a published mathematical correction, not a formal TV table erratum,
and it does not make the underlying tissue estimates observed truths.

Add a separately labeled small derived CSV using the original TV2015 risk, N,
annual-rate and lifetime-L entries, and recompute D by `N*(1+L)-1` using exact
decimal arithmetic. Attribute both the original inputs' permitted CC0
republication and the corrective formula. Do not replace the original CSV.

Repeat the primary-panel descriptive fits, all fixed-k no-expansion fits,
both prediction schemes, and the exact point/stress-box rate frontiers using
this corrected D and the effective clock D/N. Do not run an expansion grid
or a further developmental-clock analysis for this correction: N and L are
unchanged, so the latter would duplicate the original alternative-clock
analysis. No primary selection, model, metric, parameter grid, outlier rule or
verdict criterion is changed. Report the correction's results even if they
weaken the proposed interpretation.

The later-source contingency in PROTOCOL.md still applies to TV2017 coverage:
missing fields or unavailable reuse terms must be shown, not silently filled
from a different publication or endpoint.
