# Contributing to Open Cancer

Open Cancer is a compute-first AI software project. Contribute code, tests, source connectors, model-worker adapters, profiling, reproducible computational runs or bug fixes. It does **not** require lab access, volunteer recruitment, professors, medical reviewer teams or formal research programs.

## GitHub Issues: engineering only

Create an issue only for an **actionable code or technical defect/task** with an affected file/endpoint, evidence or reproduction, a deliverable, and testable acceptance criteria. Do **not** create Issues for invitations, governance, research committees, broad hypotheses, status, or speculative projects. The issue chooser provides a single engineering form. Routine monitor outputs live in Actions artifacts, not issues.

## Working on AI computation

1. Start with a specific task and public-source metadata, then use [the compute protocol](docs/compute.md) to plan deterministic jobs.
2. Supply adapters for an investigator and a separate critic model/process. Run on local hardware, cloud GPU inference or your chosen cluster; the repository ships no default paid provider.
3. Keep worker output in `.compute/` and review the unverified candidate, critique, failures, source URLs and reproducibility before selecting anything to check in.
4. For a real checked-in computational output, submit only licensed non-sensitive derived artifacts and a schema-valid, checksum-bound `results/<task-id>/<run-id>/manifest.json`; label each claim `unreviewed` with limitations and contradictory evidence when known.
5. No AI or code validator is a clinical expert. Do not infer medical effectiveness, identify patients or republish restricted upstream data.

## PR checklist

- [ ] Concrete executable behavior or defect fixed, with tests and failure cases.
- [ ] All inputs/sources, commands, model versions and outcomes are traceable.
- [ ] No patient information, credentials, controlled data or unlicensed source dumps.
- [ ] Findings remain explicitly unreviewed and nonclinical; negative results are retained.
- [ ] `python scripts/validate.py` and `python -m unittest discover -s tests -v` pass.

Read [AGENTS.md](AGENTS.md) and [DATA_POLICY.md](DATA_POLICY.md). Historical manual review formats remain available but are **optional** and are not project milestones or requirements for AI compute work.
