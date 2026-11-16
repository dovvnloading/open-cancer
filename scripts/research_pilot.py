#!/usr/bin/env python3
"""Execute, replay and compare the fixed nine-job EGFR/LUAD research experiment."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
from hashlib import sha256
from http.client import HTTPException
from importlib.metadata import version
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import time
import uuid

if __package__:
    from . import compute_fanout as compute, pilot_evaluation as evaluator, pilot_worker
    from .pilot_report import render
    from .pilot_evidence import (EvidenceUnavailable, FIXTURE_PATH, evidence_for_job, from_agent1_bundle,
                                 prepare_snapshot, validate_snapshot)
else:
    import compute_fanout as compute
    import pilot_evaluation as evaluator
    import pilot_worker
    from pilot_report import render
    from pilot_evidence import (EvidenceUnavailable, FIXTURE_PATH, evidence_for_job, from_agent1_bundle,
                                prepare_snapshot, validate_snapshot)

ROOT = compute.ROOT
WORKLOAD = "workloads/egfr-luad-compute.json"
STAGES = ("investigate", "critique")
MAX_STAGE_ATTEMPTS = 2  # second attempt requires an explicit resume
STATUS = {"review_status": "unreviewed", "publication_status": "not-published",
          "clinical_validation_status": "not-assessed", "nonclinical": True}
INPUT_FILES = (
    WORKLOAD, "tasks/egfr-luad-evidence-coverage.json", "DATA_POLICY.md",
    "catalog/sources/civic.json", "catalog/sources/gdc.json", "catalog/sources/cbioportal.json",
    "catalog/datasets/civic-nightly-clinical-evidence.json", "scripts/compute_fanout.py",
    "scripts/egfr_luad_pilot.py", "scripts/egfr_live_source_audit.py",
    "scripts/openai_compatible_worker.py", "scripts/pilot_worker.py",
    "scripts/pilot_evidence.py", "scripts/pilot_evaluation.py", "scripts/research_pilot.py",
    "scripts/pilot_report.py", "scripts/validate.py", "scripts/civic_evidence.py",
    "scripts/grounded_claims.py",
)


def now():
    return datetime.now(timezone.utc).isoformat()


def local_output(path):
    path = Path(path).resolve()
    if path.is_relative_to(ROOT):
        compute.require(path.is_relative_to(ROOT / ".compute") and path != ROOT / ".compute",
                        "Use a run directory under .compute/ or outside the repository")
    return path


@contextmanager
def exclusive(out):
    """A crashed process releases this lock. Concurrent resumes must not duplicate requests."""
    out.mkdir(parents=True, exist_ok=True)
    fd = os.open(out, os.O_RDONLY)
    try:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError("Another pilot process owns this output directory") from exc
        yield
    finally:
        os.close(fd)


def atomic_text(path, content):
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent,
                                     prefix=".tmp-", delete=False) as f:
        tmp = Path(f.name)
        try:
            f.write(content)
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise
    tmp.replace(path)


def pricing(mode, env):
    prices = {}
    for stage, role in zip(STAGES, ("INVESTIGATOR", "CRITIC")):
        prices[stage] = {}
        for side in ("input", "output"):
            name = f"OPEN_CANCER_{role}_{side.upper()}_USD_PER_MILLION"
            price = env.get(name) if mode == "real" else None
            if price is not None:
                try:
                    price = float(price)
                except ValueError:
                    raise ValueError(f"Invalid nonnegative price in {name}") from None
                compute.require(math.isfinite(price) and price >= 0, f"Invalid nonnegative price in {name}")
            prices[stage][side + "_per_million_usd"] = price
    return prices


def config_for(mode, jobs, snapshot, env):
    files = (*INPUT_FILES, FIXTURE_PATH) if mode == "fixture" else INPUT_FILES
    try:
        commit = subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                                         stderr=subprocess.DEVNULL, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = None
    return {"protocol": "open-cancer-egfr-pilot-v1", "mode": mode, "synthetic": mode == "fixture",
            "models": pilot_worker.configuration(mode, env), "prices_usd": pricing(mode, env),
            "python": platform.python_version(), "jsonschema": version("jsonschema"),
            "git_head": commit,
            "file_sha256": {p: sha256((ROOT / p).read_bytes()).hexdigest() for p in files},
            "workload": WORKLOAD, "plan_sha256": compute.digest(jobs),
            "snapshot_sha256": snapshot["snapshot_sha256"],
            "max_attempts_per_stage": MAX_STAGE_ATTEMPTS, "workers": 1, **STATUS}


def new_body(job, identity, packet):
    return {"job": job, "job_sha256": compute.digest(job), **identity, **STATUS,
            "packet_sha256": compute.digest(packet), "status": "pending", "attempts": [],
            "candidate": None, "critique": None}


def record_path(out, job):
    return out / "records" / (job["id"] + ".json")


def load_body(out, job, identity, packet, required=False):
    expected = new_body(job, identity, packet)
    path = record_path(out, job)
    if not required and not path.exists():
        return expected
    body = compute.read_envelope(path)
    keys = ("job", "job_sha256", *identity, *STATUS, "packet_sha256")
    compute.require(all(body[k] == expected[k] for k in keys),
                    "Checkpoint belongs to different job, evidence, mode or run")
    return body


def run_job(out, body, packet, mode, env):
    path = record_path(out, body["job"])
    for stage, field in zip(STAGES, ("candidate", "critique")):
        if body[field] is not None:
            continue
        stage_calls = sum(a["stage"] == stage for a in body["attempts"])
        if stage_calls >= MAX_STAGE_ATTEMPTS:
            break
        # Reserve and charge the attempt budget before dispatch. A process killed during
        # HTTP can never establish whether the provider handled or billed the call.
        attempt = {"stage": stage, "number": stage_calls + 1, "started_at": now(),
                   "synthetic": mode == "fixture", "endpoint_attempted": None if mode == "real" else False,
                   "error": "interrupted-attempt-outcome-unknown", "usage": None,
                   "elapsed_seconds": None, "returned_model": None}
        body["attempts"].append(attempt)
        body["status"] = "in-progress"
        compute.atomic_json(path, compute.envelope(body))
        result = pilot_worker.invoke(stage, body["job"], packet, body["candidate"], mode, env)
        body[field] = result.pop("response")
        attempt.update(result)
        body["status"] = "failed" if attempt["error"] else "in-progress"
        compute.atomic_json(path, compute.envelope(body))
        if attempt["error"]:
            break
    body["status"] = "completed" if body["critique"] is not None else "failed"
    body["evaluation"] = evaluator.evaluate_job(body["job"], packet, body["candidate"], body["critique"])
    compute.atomic_json(path, compute.envelope(body))
    return body


def snapshot_at(mode, supplied, pinned):
    imported = compute.read_json(supplied) if supplied is not None else None
    if isinstance(imported, dict) and imported.get("source_id") == "civic":
        imported = from_agent1_bundle(imported, mode)
    if pinned.is_file():
        snapshot = compute.read_json(pinned)
        if supplied is not None:
            compute.require(snapshot == imported, "Resume evidence differs")
    elif supplied is not None:
        snapshot = imported
    else:
        snapshot = prepare_snapshot(mode)
    validate_snapshot(snapshot, mode)
    return snapshot


def preflight_block(out, mode, jobs, reason, snapshot, elapsed, unavailable=None):
    sources = snapshot["sources"] if snapshot else unavailable
    report = {"status": "blocked-before-inference", "mode": mode, "synthetic": mode == "fixture",
              **STATUS, "reason": reason, "planned_jobs": len(jobs), "jobs_attempted": 0,
              "jobs_completed": 0, "jobs_failed": 0, "jobs_blocked": len(jobs),
              "sources": sources, "snapshot_sha256": snapshot["snapshot_sha256"] if snapshot else None,
              "plan_sha256": compute.digest(jobs), "confirmed_endpoint_requests": 0,
              "wall_seconds": elapsed, "timestamp": now()}
    compute.atomic_json(out / "benchmark.json", report)
    source_note = ("CIViC retrieval did not produce a usable projection." if snapshot is None else
                   f"CIViC source fingerprint `{sources['civic']['sha256']}`. Projected snapshot "
                   f"SHA-256 `{snapshot['snapshot_sha256']}`. See evidence.json for the "
                   "record selection, license and capture time.")
    atomic_text(out / "report.md", "# EGFR/LUAD research pilot\n\n"
                f"**Blocked before inference.** {reason}\n\n"
                f"Mode `{mode}`. {len(jobs)} jobs planned, zero attempted. "
                "No model calls, findings or token measurements. Unreviewed, unpublished, nonclinical.\n\n"
                f"{source_note} GDC and cBioPortal records were unexamined.\n")
    return report


def run(mode, output, supplied=None, resume=False, env=None):
    env = dict(os.environ if env is None else env)
    started = time.monotonic()
    out = local_output(output)
    jobs = compute.plan(compute.read_json(ROOT / WORKLOAD), max_jobs=9)
    compute.require(len(jobs) == 9, "Expected the nine-job EGFR workload")
    with exclusive(out):
        files = set(p.name for p in out.iterdir() if not p.name.startswith(".tmp-"))
        meta_path = out / "run-config.json"
        if files:
            compute.require(resume, "Run directory exists. Use --resume or choose a new directory")
            # Allow recovery across initial writes, including source/config preflight failure.
            compute.require(meta_path.is_file() or files.issubset(
                {"plan.json", "evidence.json", "benchmark.json", "report.md"}), "Unexpected output files")
        plan_path = out / "plan.json"
        if plan_path.exists():
            compute.require(compute.read_json(plan_path) == jobs, "Plan changed")
        else:
            compute.atomic_json(plan_path, jobs)
        pinned = out / "evidence.json"
        if meta_path.exists():
            compute.require(pinned.is_file(), "Pinned evidence is missing")
        try:
            snapshot = snapshot_at(mode, supplied, pinned)
        except EvidenceUnavailable as exc:
            if meta_path.exists():
                raise
            return preflight_block(out, mode, jobs, str(exc), None, time.monotonic() - started,
                                   unavailable=exc.sources)
        except (OSError, HTTPException, ValueError, RuntimeError, TypeError, KeyError):
            if meta_path.exists() or pinned.exists() or supplied is not None:
                raise
            return preflight_block(out, mode, jobs, "CIViC retrieval or projection failed.",
                                   None, time.monotonic() - started)
        if not pinned.exists():
            compute.atomic_json(pinned, snapshot)
        try:
            config = config_for(mode, jobs, snapshot, env)
        except ValueError as exc:
            if meta_path.exists():
                raise
            # Adapter diagnostics may refer to environment variable names, never values.
            return preflight_block(out, mode, jobs, str(exc), snapshot, time.monotonic() - started)
        config_sha = compute.digest(config)
        if meta_path.exists():
            metadata = compute.read_envelope(meta_path)
            compute.require(metadata["configuration"] == config, "Resume configuration differs")
        else:
            metadata = {"run_id": str(uuid.uuid4()), "created_at": now(),
                        "checkpoints_initialized": False, "configuration": config}
            compute.atomic_json(meta_path, compute.envelope(metadata))
        identity = {"run_id": metadata["run_id"], "config_sha256": config_sha,
                    "mode": mode, "synthetic": config["synthetic"]}
        packets = [evidence_for_job(snapshot, job) for job in jobs]
        # Initialize all nine files before allowing a single dispatch. A missing
        # checkpoint on a later resume must never imply that no request was made.
        initialized = metadata["checkpoints_initialized"]
        bodies = [load_body(out, job, identity, p, required=initialized)
                  for job, p in zip(jobs, packets)]
        if not initialized:
            for body in bodies:
                compute.atomic_json(record_path(out, body["job"]), compute.envelope(body))
            metadata["checkpoints_initialized"] = True
            compute.atomic_json(meta_path, compute.envelope(metadata))
        for body, packet in zip(bodies, packets):
            if body["status"] != "completed":
                run_job(out, body, packet, mode, env)
        metrics = evaluator.aggregate(bodies, config["prices_usd"], mode=mode)
        report = {**identity, **STATUS, "status": "finished-with-failures" if metrics["jobs"]["failed"] else "finished",
                  "created_at": metadata["created_at"], "reported_at": now(),
                  "configuration": config, "sources": snapshot["sources"], "selection": snapshot["selection"],
                  "metrics": metrics, "consistency": {"status": "requires-independent-repeat"},
                  "wall_seconds_this_invocation": time.monotonic() - started,
                  "records": [{"job_id": b["job"]["id"], "status": b["status"],
                               "sha256": compute.digest(b)} for b in bodies]}
        compute.atomic_json(out / "benchmark.json", report)
        atomic_text(out / "report.md", render(report, bodies))
        return report


def read_run(path, *, allow_missing=False):
    meta = compute.read_envelope(path / "run-config.json")
    config = meta["configuration"]
    snapshot = compute.read_json(path / "evidence.json")
    validate_snapshot(snapshot, config["mode"])
    jobs = compute.read_json(path / "plan.json")
    compute.require(config["snapshot_sha256"] == snapshot["snapshot_sha256"] and
                    config["plan_sha256"] == compute.digest(jobs), "Plan or evidence changed")
    identity = {"run_id": meta["run_id"], "mode": config["mode"], "synthetic": config["synthetic"],
                "config_sha256": compute.digest(config)}
    bodies = []
    for job in jobs:
        compute.validate_job(job)
        if allow_missing and not record_path(path, job).is_file():
            bodies.append({"job": job, "status": "missing"})
        else:
            bodies.append(load_body(path, job, identity, evidence_for_job(snapshot, job), required=True))
    return identity, bodies


def compare(first, second):
    a, left = read_run(first, allow_missing=True)
    b, right = read_run(second, allow_missing=True)
    compute.require(a["run_id"] != b["run_id"], "Compare independently launched runs")
    compute.require(a["config_sha256"] == b["config_sha256"], "Run configurations or evidence differ")
    return {**STATUS, "mode": a["mode"], "synthetic": a["synthetic"],
            "run_ids": [a["run_id"], b["run_id"]], "config_sha256": a["config_sha256"],
            "consistency": evaluator.compare_outputs(left, right)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="command", required=True)
    r = subs.add_parser("run", help="Execute the fixed plan using pinned evidence")
    p = subs.add_parser("prepare", help="Retrieve and fingerprint evidence only, without inference")
    for sub in (r, p):
        sub.add_argument("--mode", choices=["fixture", "real"], required=True)
        sub.add_argument("--output", type=Path, required=True)
    r.add_argument("--snapshot", type=Path, help="Use a previously captured evidence.json")
    r.add_argument("--resume", action="store_true", help="Retry unfinished roles once within the attempt bound")
    c = subs.add_parser("compare", help="Compare independently launched runs on identical pinned inputs")
    c.add_argument("first", type=Path)
    c.add_argument("second", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.command == "compare":
            print(json.dumps(compare(args.first, args.second), sort_keys=True, indent=2))
        elif args.command == "prepare":
            path = local_output(args.output)
            path.parent.mkdir(parents=True, exist_ok=True)
            with exclusive(path.parent):
                compute.require(not path.exists(), "Evidence output already exists")
                snap = prepare_snapshot(args.mode)
                compute.atomic_json(path, snap)
            print(json.dumps({"mode": args.mode, "snapshot_sha256": snap["snapshot_sha256"],
                              "inference_calls": 0}, sort_keys=True))
        else:
            report = run(args.mode, args.output, args.snapshot, args.resume)
            print(json.dumps({"status": report["status"], "mode": args.mode,
                              "jobs": report.get("metrics", {}).get("jobs"), "output": str(args.output),
                              "reason": report.get("reason")}, sort_keys=True))
            return 0 if report["status"] == "finished" else 1
        return 0
    except (OSError, HTTPException, ValueError, RuntimeError, KeyError, TypeError) as exc:
        print(f"research-pilot: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
