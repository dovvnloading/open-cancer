#!/usr/bin/env python3
"""Provider-neutral, resumable AI job fan-out. No claims are scientifically verified."""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile

if __package__:
    from . import civic_evidence, grounded_claims
else:
    import civic_evidence
    import grounded_claims

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = "open-cancer-compute-v1"
SLUG = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")
SHA = re.compile(r"^[0-9a-f]{64}$")
MAX_RESPONSE_BYTES = 262_144
MAX_QUESTION_LENGTH = 4096


def encoded(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":"),
                      allow_nan=False).encode("utf-8")


def digest(value: object) -> str:
    return hashlib.sha256(encoded(value)).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def read_json(path: Path) -> object:
    return json.loads(path.read_text(encoding="utf-8"))


def stable_id(job: dict) -> str:
    payload = {k: v for k, v in job.items() if k != "id"}
    return f"{job['task_id']}-{digest(payload)[:16]}"


def validate_job(job: object) -> dict:
    require(isinstance(job, dict), "Job must be a JSON object")
    require(set(job) == {"protocol", "id", "task_id", "question", "source_ids", "dimensions"},
            "Job has missing or unsupported fields")
    require(job["protocol"] == PROTOCOL, "Unsupported job protocol")
    require(isinstance(job["task_id"], str) and bool(SLUG.fullmatch(job["task_id"])),
            "Invalid task ID")
    require(isinstance(job["question"], str) and 30 <= len(job["question"]) <= MAX_QUESTION_LENGTH,
            "Question must be 30-4096 characters")
    ids = job["source_ids"]
    require(isinstance(ids, list) and bool(ids) and
            all(isinstance(s, str) and SLUG.fullmatch(s) for s in ids) and
            ids == sorted(set(ids)), "source_ids must be unique, sorted source slugs")
    dims = job["dimensions"]
    require(isinstance(dims, dict) and bool(dims) and
            all(isinstance(k, str) and SLUG.fullmatch(k) and
                isinstance(v, str) and 1 <= len(v) <= 160 for k, v in dims.items()),
            "Dimensions must be a nonempty mapping of small strings")
    require(job["id"] == stable_id(job), "Job ID does not match content digest")
    return job


def plan(spec: object, max_jobs: int, root: Path = ROOT) -> list[dict]:
    require(isinstance(spec, dict), "Workload must be a JSON object")
    require(set(spec) == {"schema_version", "task_id", "question", "source_ids", "axes"},
            "Workload fields must be schema_version, task_id, question, source_ids and axes")
    require(spec["schema_version"] == 1, "Workload schema_version must be 1")
    task_id = spec["task_id"]
    require(isinstance(task_id, str) and bool(SLUG.fullmatch(task_id)), "Invalid task ID")
    task_path = root / "tasks" / (task_id + ".json")
    require(task_path.is_file(), f"Unknown checked-in task: {task_id}")
    task = read_json(task_path)
    require(isinstance(task, dict) and task.get("id") == task_id, "Task ID does not match file")
    ids = spec["source_ids"]
    require(isinstance(ids, list) and bool(ids) and all(isinstance(s, str) for s in ids),
            "source_ids must be a nonempty list")
    require(set(ids).issubset(task.get("input_source_ids", [])),
            "Workload source_ids must belong to the checked-in task")
    for source in ids:
        require(bool(SLUG.fullmatch(source)) and
                (root / "catalog" / "sources" / (source + ".json")).is_file(),
                f"Unknown catalog source: {source}")
    axes = spec["axes"]
    require(isinstance(axes, dict) and bool(axes) and len(axes) <= 8,
            "axes must contain 1-8 dimensions")
    size = 1
    for key, values in axes.items():
        require(isinstance(key, str) and bool(SLUG.fullmatch(key)), "Invalid axis slug")
        require(isinstance(values, list) and bool(values) and len(values) <= 1000 and
                all(isinstance(v, str) and 1 <= len(v) <= 160 for v in values) and
                len(set(values)) == len(values), f"Axis {key} has invalid/duplicate values")
        size *= len(values)
        require(size <= max_jobs, f"Workload has {size}+ jobs, exceeds --max-jobs {max_jobs}")
    keys = sorted(axes)
    jobs = []
    for values in itertools.product(*(sorted(axes[k]) for k in keys)):
        job = {"protocol": PROTOCOL, "task_id": task_id, "question": spec["question"],
               "source_ids": sorted(set(ids)), "dimensions": dict(zip(keys, values))}
        job["id"] = stable_id(job)
        jobs.append(validate_job(job))
    return sorted(jobs, key=lambda x: x["id"])


def load_jobs(path: Path) -> list[dict]:
    jobs = []
    ids = set()
    for n, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        try:
            job = validate_job(json.loads(line))
        except (ValueError, TypeError, KeyError) as exc:
            raise ValueError(f"Invalid plan line {n}: {exc}") from exc
        require(job["id"] not in ids, f"Duplicate job ID on line {n}")
        ids.add(job["id"])
        jobs.append(job)
    require(bool(jobs), "Plan has no jobs")
    return jobs


def check_response(stage: str, payload: object, job: dict) -> dict:
    require(isinstance(payload, dict), f"{stage} response must be a JSON object")
    if stage == "investigate":
        require(set(payload) == {"summary", "evidence", "limitations", "synthetic"},
                "Investigator response fields are invalid")
        require(isinstance(payload["summary"], str) and
                20 <= len(payload["summary"]) <= 10000, "Invalid summary")
        require(type(payload["synthetic"]) is bool, "synthetic must be a boolean")
        evidence = payload["evidence"]
        require(isinstance(evidence, list) and bool(evidence) and len(evidence) <= 100,
                "Investigator must cite 1-100 evidence/context URLs")
        for item in evidence:
            require(isinstance(item, dict) and
                    set(item) == {"source_id", "url", "relation"} and
                    item["source_id"] in job["source_ids"] and
                    isinstance(item["url"], str) and
                    item["url"].startswith("https://") and len(item["url"]) <= 2048 and
                    item["relation"] in {"supports", "contradicts", "context"},
                    "Evidence needs a task source, HTTPS URL and relation")
        limits = payload["limitations"]
        require(isinstance(limits, list) and bool(limits) and len(limits) <= 40 and
                all(isinstance(x, str) and 10 <= len(x) <= 1000 for x in limits),
                "Investigator must describe limitations")
    else:
        require(set(payload) == {"verdict", "issues"}, "Critic response fields are invalid")
        require(payload["verdict"] in {"no-internal-issues-found", "inconsistent", "insufficient-evidence"},
                "Invalid critique verdict")
        issues = payload["issues"]
        require(isinstance(issues, list) and len(issues) <= 40 and
                all(isinstance(x, str) and 10 <= len(x) <= 1000 for x in issues),
                "Invalid critique issues")
        require(payload["verdict"] == "no-internal-issues-found" or bool(issues),
                "Nonpassing critique must explain issues")
    return payload


def check_directory_response(stage: str, payload: object, job: dict, root: Path = ROOT) -> dict:
    answer = check_response(stage, payload, job)
    if stage == "investigate":
        for item in answer["evidence"]:
            source = read_json(root / "catalog/sources" / (item["source_id"] + ".json"))
            if item["relation"] != "context" or item["url"] not in {
                    source.get("homepage"), source.get("documentation_url")}:
                raise ValueError("Source directory only: cite a contextual catalog URL")
    elif answer["verdict"] == "no-internal-issues-found":
        return {"verdict": "insufficient-evidence", "issues": [
            "Source directory metadata cannot validate record-level or biological claims.", *answer["issues"][:39]]}
    return answer


def invoke(stage: str, command: list[str], job: dict, candidate: dict | None,
           timeout: int, retries: int, bundle: dict | None = None) -> tuple[dict, int]:
    request = {"protocol": PROTOCOL, "stage": stage, "job": job}
    if bundle is not None:
        request.update(protocol=grounded_claims.GROUNDED_PROTOCOL, evidence_bundle=bundle)
    if candidate is not None:
        request["candidate"] = candidate
    last_error = "unknown-error"
    for attempt in range(1, retries + 2):
        try:
            process = subprocess.run(command, input=encoded(request), capture_output=True,
                                     timeout=timeout, check=False)
            if process.returncode != 0:
                last_error = f"worker-exit-{process.returncode}"
                continue
            if len(process.stdout) > MAX_RESPONSE_BYTES:
                last_error = "response-too-large"
                continue
            try:
                payload = json.loads(process.stdout)
                if bundle is None:
                    response = check_directory_response(stage, payload, job)
                else:
                    response = grounded_claims.check_grounded_response(stage, payload, bundle,
                                                                       job, candidate)
            except (ValueError, TypeError, KeyError, OSError, RecursionError):
                last_error = "invalid-worker-response"
                continue
            return response, attempt
        except subprocess.TimeoutExpired:
            last_error = "worker-timeout"
        except OSError:
            last_error = "worker-execution-error"
    raise RuntimeError(f"{stage}: {last_error} after {retries + 1} attempts")


def atomic_json(path: Path, value: object, *, exclusive: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=path.parent, prefix=".tmp-", suffix=".json",
                                     encoding="utf-8", delete=False) as output:
        temp_path = Path(output.name)
        try:
            json.dump(value, output, sort_keys=True, indent=2, ensure_ascii=False,
                      allow_nan=False)
            output.write("\n")
        except BaseException:
            temp_path.unlink(missing_ok=True)
            raise
    try:
        if exclusive:
            os.link(temp_path, path)  # Publish only a complete file, without replacing another run.
        else:
            temp_path.replace(path)
    finally:
        temp_path.unlink(missing_ok=True)


def envelope(body: dict) -> dict:
    return {"body": body, "sha256": digest(body)}


def read_envelope(path: Path) -> dict:
    saved = read_json(path)
    require(isinstance(saved, dict) and set(saved) == {"body", "sha256"} and
            isinstance(saved["body"], dict) and saved["sha256"] == digest(saved["body"]),
            f"Corrupt or modified checkpoint: {path}")
    return saved["body"]


def investigate(job: dict, commands: dict, opts: dict, bundle: dict | None) -> dict:
    if bundle is not None and not civic_evidence.records_for_job(bundle, job):
        candidate = grounded_claims.insufficient_response("investigate", bundle)
        critique = grounded_claims.insufficient_response("critique", bundle)
        attempts_i = attempts_c = 0
    else:
        candidate, attempts_i = invoke("investigate", commands["investigate"], job, None,
                                       opts["timeout"], opts["retries"], bundle)
        critique, attempts_c = invoke("critique", commands["critique"], job, candidate,
                                      opts["timeout"], opts["retries"], bundle)
    if bundle is not None and critique["verdict"] == "insufficient-evidence":
        triage = "insufficient-evidence"
    elif candidate["synthetic"]:
        triage = "synthetic"
    elif critique["verdict"] == "no-internal-issues-found":
        triage = "candidate-unverified"
    else:
        triage = "needs-investigation"
    return {"status": "completed", "triage": triage, "candidate": candidate,
            "critique": critique, "attempts": {"investigate": attempts_i, "critique": attempts_c}}


def run_one(job: dict, out: Path, commands: dict, opts: dict, config_sha: str,
            bundle: dict | None = None) -> dict:
    try:
        result = investigate(job, commands, opts, bundle)
    except RuntimeError as exc:
        result = {"status": "failed", "triage": "failed", "error": str(exc)}
    body = {"protocol": grounded_claims.GROUNDED_PROTOCOL if bundle is not None else PROTOCOL,
            "job": job, "job_sha256": digest(job), "config_sha256": config_sha,
            "review_status": "unreviewed", "clinical_validation_status": "not-assessed",
            "publication_status": "not-published", **result}
    if bundle is not None:
        body["evidence_sha256"] = bundle["sha256"]
    atomic_json(out / "records" / (job["id"] + ".json"), envelope(body))
    return body


def civic_catalog() -> dict:
    return read_json(ROOT / "catalog/sources/civic.json")


def private_output(path: Path) -> None:
    resolved = path.resolve()
    require(not resolved.is_relative_to(ROOT) or resolved.is_relative_to(ROOT / ".compute"),
            "Keep evidence snapshots and grounded runs under .compute/ or outside the repository")


def load_evidence(path: Path) -> dict:
    with path.open("rb") as source:
        raw = source.read(2_097_153)
    require(len(raw) <= 2_097_152, "Evidence snapshot exceeds the input limit")
    try:
        bundle = json.loads(raw)
    except (ValueError, RecursionError) as exc:
        raise ValueError("Malformed evidence snapshot") from exc
    return civic_evidence.validate_bundle(bundle, civic_catalog())


def run_batch(jobs: list[dict], out: Path, investigator: str, critic: str,
              workers: int, timeout: int, retries: int, max_jobs: int,
              shard_index: int, shard_count: int, resume: bool, evidence: Path | None = None) -> dict:
    require(workers >= 1 and timeout >= 1 and retries >= 0 and max_jobs >= 1,
            "workers, timeout and max-jobs must be positive; retries >= 0")
    require(shard_count >= 1 and 0 <= shard_index < shard_count,
            "Invalid shard index/count")
    commands = {"investigate": shlex.split(investigator), "critique": shlex.split(critic)}
    require(all(commands.values()), "Both worker commands are required")
    selected = [job for job in jobs if int(job["id"].rsplit("-", 1)[1], 16) % shard_count
                == shard_index]
    require(len(selected) <= max_jobs, "Selected jobs exceed --max-jobs")
    configuration = {"protocol": PROTOCOL, "plan_sha256": digest(jobs),
                     "worker_fingerprint": digest(commands), "timeout": timeout,
                     "retries": retries, "shard_index": shard_index, "shard_count": shard_count}
    bundle = None
    if evidence is not None:
        private_output(out)
        for job in selected:
            require(civic_evidence.supports_job(job), "Job outside the CIViC EGFR/LUAD scope")
        bundle = load_evidence(evidence)
        configuration.update(protocol=grounded_claims.GROUNDED_PROTOCOL,
                             evidence_bundle=bundle, evidence_sha256=bundle["sha256"])
        # Distinguish model settings across resumes without storing endpoint URLs or keys.
        settings = {k: os.environ.get("OPEN_CANCER_" + k) for k in (
            "LLM_URL", "LLM_MODEL", "INVESTIGATOR_MODEL", "CRITIC_MODEL",
            "LLM_TIMEOUT", "MAX_TOKENS")}
        configuration["inference_config_sha256"] = digest(settings)
    config_sha = digest(configuration)
    meta = out / "run-config.json"
    if meta.exists():
        require(resume, "Output exists; pass --resume to reuse matching checkpoints")
        require(read_json(meta) == configuration, "Resume configuration differs from checkpoint")
    else:
        require(not out.exists() or not any(out.iterdir()), "Output directory must be empty")
        atomic_json(meta, configuration, exclusive=True)
    bodies = {}
    pending = []
    for job in selected:
        checkpoint = out / "records" / (job["id"] + ".json")
        if resume and checkpoint.is_file():
            saved = read_envelope(checkpoint)
            require(saved.get("job") == job and saved.get("config_sha256") == config_sha,
                    "Checkpoint does not match job or worker configuration")
            if saved.get("status") == "completed":
                if bundle is not None:
                    require(saved.get("evidence_sha256") == configuration["evidence_sha256"],
                            "Checkpoint evidence differs from the pinned snapshot")
                    candidate = grounded_claims.check_grounded_response(
                        "investigate", saved.get("candidate"), bundle, job)
                    grounded_claims.check_grounded_response(
                        "critique", saved.get("critique"), bundle, job, candidate)
                else:
                    check_directory_response("investigate", saved.get("candidate"), job)
                    require(check_directory_response("critique", saved.get("critique"), job)
                            == saved.get("critique"), "Legacy checkpoint lacks record evidence")
                bodies[job["id"]] = saved
                continue
        pending.append(job)
    opts = {"timeout": timeout, "retries": retries}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(run_one, job, out, commands, opts, config_sha, bundle): job["id"]
                   for job in pending}
        for future in as_completed(futures):
            bodies[futures[future]] = future.result()
    counts = Counter(body["status"] for body in bodies.values())
    triage_counts = Counter(body["triage"] for body in bodies.values())
    summary = {"protocol": configuration["protocol"], "plan_sha256": configuration["plan_sha256"],
               "config_sha256": config_sha, "shard_index": shard_index,
               "shard_count": shard_count, "selected_jobs": len(selected),
               "completed": counts["completed"], "failed": counts["failed"],
               "triage": dict(sorted(triage_counts.items())),
               "status": "unreviewed-computational-output-not-clinical-evidence",
               "records": [{"id": job_id, "status": bodies[job_id]["status"],
                            "triage": bodies[job_id]["triage"]} for job_id in sorted(bodies)]}
    atomic_json(out / "summary.json", summary)
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("plan", help="Build a deterministic JSONL worklist")
    p.add_argument("--spec", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--max-jobs", type=int, default=256)
    fetch = sub.add_parser("retrieve", help="Snapshot permitted CIViC EGFR/LUAD evidence")
    fetch.add_argument("--source", choices=["civic"], required=True)
    fetch.add_argument("--output", type=Path, required=True)
    fetch.add_argument("--fixture", action="store_true", help="Replay the labeled offline fixture")
    fetch.add_argument("--max-pages", type=int, help="Live page budget (default 4, maximum 8)")
    fetch.add_argument("--page-size", type=int, default=20, help="Items per page (maximum 25)")
    r = sub.add_parser("run", help="Fan out to subprocess worker adapters")
    r.add_argument("--jobs", type=Path, required=True)
    r.add_argument("--output", type=Path, required=True)
    r.add_argument("--evidence", type=Path, help="Snapshot created by the retrieve command")
    r.add_argument("--investigator", required=True, help="Subprocess command, reads JSON stdin")
    r.add_argument("--critic", required=True, help="Separate subprocess command, reads JSON stdin")
    r.add_argument("--workers", type=int, default=4)
    r.add_argument("--timeout", type=int, default=120)
    r.add_argument("--retries", type=int, default=1)
    r.add_argument("--max-jobs", type=int, default=100)
    r.add_argument("--shard-index", type=int, default=0)
    r.add_argument("--shard-count", type=int, default=1)
    r.add_argument("--resume", action="store_true")
    args = parser.parse_args(argv)
    try:
        if args.command == "plan":
            require(args.max_jobs >= 1, "max-jobs must be positive")
            require(not args.output.exists(), "Refusing to overwrite existing job plan")
            jobs = plan(read_json(args.spec), args.max_jobs)
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text("".join(json.dumps(j, sort_keys=True) + "\n" for j in jobs),
                                   encoding="utf-8")
            print(f"Planned {len(jobs)} content-addressed jobs: {args.output}")
        elif args.command == "retrieve":
            private_output(args.output)
            require(not args.output.exists(), "Refusing to overwrite existing evidence snapshot")
            catalog = civic_catalog()
            if args.fixture:
                require(args.max_pages in (None, 1) and args.page_size == 20,
                        "Fixture replay uses its recorded page budget")
                bundle = civic_evidence.retrieve_fixture(catalog)
            else:
                bundle = civic_evidence.retrieve(catalog, max_pages=4 if args.max_pages is None else args.max_pages,
                                                 page_size=args.page_size)
            civic_evidence.validate_bundle(bundle, catalog)
            atomic_json(args.output, bundle, exclusive=True)
            print(json.dumps({"status": bundle["status"], "synthetic": bundle["synthetic"],
                              "records": len(bundle["records"]), "sha256": bundle["sha256"]},
                             sort_keys=True))
            return 1 if bundle["status"] == "unavailable" else 0
        else:
            summary = run_batch(load_jobs(args.jobs), args.output, args.investigator,
                                args.critic, args.workers, args.timeout, args.retries,
                                args.max_jobs, args.shard_index, args.shard_count, args.resume,
                                args.evidence)
            print(json.dumps({k: summary[k] for k in ("selected_jobs", "completed", "failed", "triage")},
                             sort_keys=True))
            return 1 if summary["failed"] else 0
        return 0
    except (ValueError, OSError, json.JSONDecodeError) as exc:
        print(f"compute-fanout: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
