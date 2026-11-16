"""Render a local, nonclinical technical report from the measured benchmark."""
from __future__ import annotations

import json


def literal(value):
    """Render source/model strings as literal Markdown, including on multiline responses."""
    return "\n".join("    " + line for line in str(value).replace("\r", "\n").split("\n"))


def amount(value):
    return "unavailable" if value is None else str(value)


def render(report, bodies):
    metrics = report["metrics"]
    cfg = report["configuration"]
    civic = report["sources"]["civic"]
    execution = metrics["execution"]
    jobs = metrics["jobs"]
    checks = metrics["deterministic_checks"]
    observations = metrics["observations"]
    fixture = report["synthetic"]
    lines = [
        "# EGFR/LUAD technical research report", "",
        "**Unreviewed · unpublished · nonclinical**", "",
        ("Synthetic fixture experiment. The records and both role responses were generated "
         "for software validation. No live inference was performed." if fixture else
         "Real-mode experiment using a pinned public CIViC projection and operator-configured "
         "inference endpoints. The checkpoints identify which requests succeeded."),
        "", f"Run `{report['run_id']}` · {report['reported_at']}", "",
        f"{jobs['completed']} of {jobs['total']} jobs completed both roles. "
        f"{jobs['failed']} failed. {jobs['with_candidate']} investigator and "
        f"{jobs['with_critique']} critic responses are available.", "",
        "## Evidence available to the models", "", report["selection"], "",
        "| Source | Status | Jobs with records | Unique IDs supplied | Unique IDs cited by investigator |",
        "| --- | --- | ---: | ---: | ---: |"]
    for source, coverage in metrics["source_coverage"].items():
        lines.append(f"| {source} | {report['sources'][source]['status']} "
                     f"| {coverage['jobs_with_records']} | {coverage['unique_records_supplied']} "
                     f"| {coverage['unique_records_cited_by_investigator']} |")
    lines.extend(["", "GDC and cBioPortal have no record inputs in this pilot. Their cohort "
                  "membership, mutation calls and profiled denominators were unexamined.", "",
                  f"CIViC input · {civic['attribution']} · {civic['license']} · {civic['release']}",
                  "", f"Input `{civic['url']}` · {civic['connector']} · {civic['retrieval_status']} · "
                  f"retrieved {civic['retrieved_at'] or 'offline fixture'}",
                  "", f"Source fingerprint (native GraphQL bundle or raw fallback TSV SHA-256) `{civic['sha256']}`",
                  "", f"Projected snapshot SHA-256 `{cfg['snapshot_sha256']}`",
                  "", "## Benchmark", "",
                  f"- Citation occurrences resolving to the job packet: {metrics['citations']['valid']}. "
                  f"Unsupported occurrences: {metrics['citations']['unsupported']}.",
                  f"- Deterministic field and count checks: {checks['passed']} passed, "
                  f"{checks['failed']} failed, {checks['unassessed']} unassessed.",
                  f"- Observations flagged by the evaluator: {observations['flagged_unsupported']} "
                  f"of {observations['total']}. Observations without checkable assertions: "
                  f"{observations['without_checkable_assertions']}.",
                  f"- Investigator counterevidence entries: {observations['contradictory_evidence_entries']}. "
                  f"Observations with a cohort quantity for review: {observations['with_cohort_quantity_to_review']}.",
                  "", "| Critic issue category | Issues raised |", "| --- | ---: |"])
    for category, n in metrics["critic_issues_by_category"].items():
        lines.append(f"| {category} | {n} |")
    native = civic["agent1"]
    lines += ["", "### CIViC request lineage", "",
              f"Agent 1 source status: {native['status']}. Failure code: {native['failure'] or 'none'}.",
              "", "| Page | Nodes | Response-byte SHA-256 |", "| ---: | ---: | --- |"]
    for i, page in enumerate(native["pages"], start=1):
        lines.append(f"| {i} | {page['node_count']} | `{page['response_sha256']}` |")
    lines += ["", "The saved evidence snapshot includes the original GraphQL query, pagination "
              "variables, annotation fingerprints and source capture time. All match counts are "
              "limited to the pages or export actually retrieved.",
              "", "## Execution and reproducibility", "",
              f"Plan SHA-256 `{cfg['plan_sha256']}` · configuration SHA-256 `{report['config_sha256']}`",
              "", f"Git HEAD `{cfg['git_head']}`. The file fingerprints in run-config.json identify "
              f"the actual runtime inputs, including local edits. Python {cfg['python']}, jsonschema {cfg['jsonschema']}.",
              "", "Requested models and options (endpoint and prompt identified by SHA-256)", "",
              literal(json.dumps(cfg["models"], indent=2, sort_keys=True)), "",
              f"{execution['attempts']} role attempts, {execution['failed_attempts']} failed. "
              f"Confirmed endpoint attempts: {execution['endpoint_attempts']}. "
              f"Uncertain dispatches: {execution['attempts_of_unknown_dispatch']}.",
              "", f"This invocation took {report['wall_seconds_this_invocation']:.3f} seconds. "
              f"Accumulated measured role time in seconds: {amount(execution['elapsed_attempt_seconds'])}.",
              "", "| Role | Prompt tokens | Completion tokens | Estimated USD |",
              "| --- | ---: | ---: | ---: |"]
    for stage, row in execution["by_stage"].items():
        lines.append(f"| {stage} | {amount(row['usage']['prompt_tokens'])} "
                     f"| {amount(row['usage']['completion_tokens'])} | {amount(row['estimated_cost_usd'])} |")
    lines += ["", f"Total estimated token cost in USD: {amount(execution['estimated_cost_usd'])}. "
              "Per-attempt usage, returned model identifiers, available provider fingerprints and partial "
              "subtotals are retained in the checkpoints and benchmark.json.",
              "", "Single-run output consistency is unmeasured. Use the compare command on two independent "
              "launches with the same pinned configuration for paired output and cited-record agreement.",
              "", "## Job results", "",
              "| Alteration | Question | Status | Critic verdict | Checkpoint |", "| --- | --- | --- | --- | --- |"]
    for b in bodies:
        dims = b["job"]["dimensions"]
        job_id = b["job"]["id"]
        verdict = (b["critique"] or {}).get("verdict", "unavailable")
        lines.append(f"| {dims['alteration']} | {dims['check']} | {b['status']} | {verdict} "
                     f"| [JSON](records/{job_id}.json) |")
    lines += ["", "## Investigator observations and critic challenges", "",
              "The following entries reproduce generated text. Each checkpoint retains the exact "
              "citation triples, field checks, assumptions, limitations and attempt history."]
    for b in bodies:
        dims = b["job"]["dimensions"]
        lines += ["", f"### {dims['alteration']} · {dims['check']}", ""]
        for attempt in b["attempts"]:
            if attempt["error"]:
                lines += [literal(f"{attempt['stage']} attempt {attempt['number']}: {attempt['error']}"), ""]
        candidate = b["candidate"]
        if candidate:
            for section in ("observations", "contradictory_evidence"):
                for o in candidate[section]:
                    ids = ", ".join(f"{c['source_id']}/{c['record_id']}" for c in o["citations"])
                    lines += [literal(f"{section} [{o['id']}]: {o['text']} ({ids})"), ""]
            if not candidate["observations"]:
                lines += ["No investigator observations supplied.", ""]
        for issue in (b["critique"] or {}).get("issues", []):
            lines += [literal(f"Critic, {issue['category']}: {issue['detail']}"), ""]
    lines += ["## Limitations", "",
              "The projection contains structured CIViC annotations. Original papers, evidence statements, "
              "compound profiles and aliases outside the exact filter were outside scope. Selection favors "
              "diversity of direction/significance labels and cannot estimate population frequencies. "
              "No cross-cohort denominator is available.",
              "", "A matching citation identifies a supplied row. Literal field and packet-count checks "
              "evaluate only those assertions. Free-text meaning, source quality and whether differing "
              "annotations constitute a biological contradiction require separate assessment. Matching "
              "model outputs cannot establish efficacy or scientific validity.", ""]
    return "\n".join(lines)
