# Breast receptor-by-age proportionality: frozen aggregate audit

**Nonclinical, unreviewed research. Not medical advice or a screening recommendation.**

**Verdict: INCONCLUSIVE.** See [the complete generated analysis](REPORT.md).
HER2 is **UNTESTABLE** in the source Table 1: its two recorded classes
are ER-positive and/or PR-positive versus ER-negative **and** PR-negative.

The qualitative receptor-dependent age pattern has already been described
in the literature. This is a narrower, conditional reproducibility audit
of the simultaneous two-class proportionality restriction in **2016 USCS**
data, not a claim of epidemiological discovery or independent replication.
The strict failure gate asks whether an interval excludes Q=1 after
allowing **every** allocation of unknown and residual receptor cases.
No pair passed that stronger failure gate; a minimum restoration completion
was therefore not triggered or computed.

The [PROTOCOL.md](PROTOCOL.md) was committed by itself before numerical
work (GitHub commit `ec33191bec49a945149ab1fc998ed588a41715d3`).
It freezes data, exact binomial uncertainty, 10 contrasts, missingness
bounds, rounding stress tests, verdicts and reporting. It was not
blinded or registered outside this repository.

Source: [Turbow et al. 2021 Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC10846538/#T1), 2016 USCS.
This is an author manuscript with no explicit CC licence visible;
the derived CSV contains attributed factual counts/rates only.

## Reproduce, offline (stdlib Python only)

```sh
python research/breast-receptor-age/analysis.py --check
python -m unittest discover -s tests -p test_research_breast_receptor_age.py -v
```

Use `--write` only if regenerating README.md and REPORT.md from the
committed CSV/protocol. Hashes normalize all newline conventions to LF,
and path keys use `Path.as_posix()` for cross-platform reproduction.

Results concern one-year population-level group contrasts; no stage
cross-tabulation or HER2 split is available. They identify no screening,
birth-cohort, biological, competing-risk, or clinical cause.
