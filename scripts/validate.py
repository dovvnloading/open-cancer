#!/usr/bin/env python3
"""Offline structural checks for catalog, tasks, and reproducible submissions.

This does NOT verify medical or biological claims or fetch network resources.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath
import sys

from jsonschema import Draft202012Validator, FormatChecker

ROOT = Path(__file__).resolve().parents[1]
MAX_ARTIFACT_BYTES = 10_000_000


def load_json(path: Path, errors: list[str]):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        errors.append(f"{path}: cannot read JSON: {exc}")
        return None


def validate_repository(root: Path = ROOT) -> list[str]:
    """Return all structural and local integrity errors (empty list means pass)."""
    errors: list[str] = []
    schemas = {}
    for kind in ("source", "dataset", "task", "result", "review", "publication"):
        schema_path = root / "schemas" / f"{kind}.schema.json"
        schema = load_json(schema_path, errors)
        if isinstance(schema, dict):
            try:
                Draft202012Validator.check_schema(schema)
                schemas[kind] = Draft202012Validator(schema, format_checker=FormatChecker())
            except Exception as exc:
                errors.append(f"{schema_path}: invalid JSON Schema: {exc}")

    def validate(path: Path, kind: str):
        data = load_json(path, errors)
        if not isinstance(data, dict):
            if data is not None:
                errors.append(f"{path}: expected JSON object")
            return None
        validator = schemas.get(kind)
        if validator:
            for err in sorted(validator.iter_errors(data), key=lambda e: (list(map(str, e.path)), e.message)):
                where = ".".join(map(str, err.path)) or "<root>"
                errors.append(f"{path}:{where}: {err.message}")
        return data

    source_paths = sorted((root / "catalog" / "sources").glob("*.json"))
    dataset_paths = sorted((root / "catalog" / "datasets").glob("*.json"))
    task_paths = sorted((root / "tasks").glob("*.json"))
    if not source_paths:
        errors.append("catalog/sources/: no source records")
    if not dataset_paths:
        errors.append("catalog/datasets/: no dataset records")
    if not task_paths:
        errors.append("tasks/: no research tasks")

    sources = {}
    for path in source_paths:
        obj = validate(path, "source")
        if obj:
            key = obj.get("id")
            if key != path.stem:
                errors.append(f"{path}: id must equal file stem")
            if key in sources:
                errors.append(f"{path}: duplicate source id {key}")
            if isinstance(key, str):
                sources[key] = obj


    dataset_ids = set()
    datasets_with_context = []
    for path in dataset_paths:
        obj = validate(path, "dataset")
        if not obj:
            continue
        key = obj.get("id")
        if key != path.stem:
            errors.append(f"{path}: dataset id must equal file stem")
        if key in dataset_ids:
            errors.append(f"{path}: duplicate dataset id {key}")
        if isinstance(key, str):
            dataset_ids.add(key)
        if obj.get("source_id") not in sources:
            errors.append(f"{path}: dataset references unknown source {obj.get('source_id')}")
        if obj.get("research_context"):
            datasets_with_context.append((path, obj))

    tasks = {}
    for path in task_paths:
        obj = validate(path, "task")
        if obj:
            key = obj.get("id")
            if key != path.stem:
                errors.append(f"{path}: id must equal file stem")
            if key in tasks:
                errors.append(f"{path}: duplicate task id {key}")
            if isinstance(key, str):
                tasks[key] = obj
            for source_id in obj.get("input_source_ids", []):
                if source_id not in sources:
                    errors.append(f"{path}: unknown source id {source_id}")

    # Research context is an explicit relationship to a particular submitted run,
    # not an inferred claim that every item in the upstream dataset is gene-specific.
    for dataset_path, record in datasets_with_context:
        for relation in record.get("research_context", []):
            task_id, run_id = relation.get("task_id"), relation.get("run_id")
            if task_id not in tasks:
                errors.append(f"{dataset_path}: research context references unknown task {task_id}")
                continue
            if not (isinstance(task_id, str) and isinstance(run_id, str)
                    and task_id.replace("-", "").isalnum() and run_id.replace("-", "").isalnum()):
                errors.append(f"{dataset_path}: research context has invalid task or run identifier")
                continue
            result_path = root / "results" / task_id / run_id / "manifest.json"
            if not result_path.is_file():
                errors.append(f"{dataset_path}: research context references missing result {task_id}/{run_id}")
                continue
            result_obj = load_json(result_path, errors)
            if not isinstance(result_obj, dict) or (result_obj.get("task_id"), result_obj.get("run_id")) != (task_id, run_id):
                errors.append(f"{dataset_path}: research context points to mismatched result")
                continue
            if record.get("source_id") not in [item.get("source_id") for item in result_obj.get("inputs", [])
                                                 if isinstance(item, dict)]:
                errors.append(f"{dataset_path}: linked research result did not use this source")

    for path in sorted((root / "results").glob("*/*/manifest.json")):
        obj = validate(path, "result")
        if not obj:
            continue
        if obj.get("task_id") != path.parent.parent.name:
            errors.append(f"{path}: task_id must match parent task directory")
        if obj.get("run_id") != path.parent.name:
            errors.append(f"{path}: run_id must match run directory")
        if obj.get("task_id") not in tasks:
            errors.append(f"{path}: unknown task_id {obj.get('task_id')}")

        input_ids = set()
        for record in obj.get("inputs", []):
            if not isinstance(record, dict):
                continue
            source_id = record.get("source_id")
            input_ids.add(source_id)
            if source_id not in sources:
                errors.append(f"{path}: input references unknown source {source_id}")
        for finding in obj.get("findings", []):
            if not isinstance(finding, dict):
                continue
            for citation in finding.get("evidence", []):
                if isinstance(citation, dict) and citation.get("source_id") not in input_ids:
                    errors.append(f"{path}: evidence source must also appear in inputs: {citation.get('source_id')}")

        def local_file(relative: str):
            if not isinstance(relative, str) or not relative or "\\" in relative:
                errors.append(f"{path}: invalid local path {relative!r}")
                return None
            rel = PurePosixPath(relative)
            if rel.is_absolute() or any(x in (".", "..") for x in rel.parts):
                errors.append(f"{path}: unsafe local path {relative!r}")
                return None
            resolved = (path.parent / relative).resolve()
            base = path.parent.resolve()
            if not resolved.is_relative_to(base) or not resolved.is_file():
                errors.append(f"{path}: missing or out-of-directory file {relative!r}")
                return None
            return resolved

        execution = obj.get("execution", {})
        if isinstance(execution, dict) and isinstance(execution.get("code_path"), str):
            local_file(execution["code_path"])
        for artifact in obj.get("artifacts", []):
            if not isinstance(artifact, dict):
                continue
            actual_path = local_file(artifact.get("path"))
            if actual_path is None:
                continue
            if actual_path.stat().st_size > MAX_ARTIFACT_BYTES:
                errors.append(f"{path}: artifact exceeds {MAX_ARTIFACT_BYTES} bytes: {actual_path.name}")
                continue
            expected = artifact.get("sha256")
            actual = hashlib.sha256(actual_path.read_bytes()).hexdigest()
            if actual != expected:
                errors.append(f"{path}: SHA-256 mismatch for {actual_path.name}")

    # Only an auditable, SHA-pinned editorial decision can feature a submitted result.
    publication_path = root / "catalog" / "publication-decisions.json"
    decisions = load_json(publication_path, errors)
    if not isinstance(decisions, list):
        errors.append(f"{publication_path}: expected an array of publication decisions")
    else:
        seen = set()
        for index, decision in enumerate(decisions):
            if not isinstance(decision, dict):
                errors.append(f"{publication_path}[{index}]: expected a publication decision object")
                continue
            validator = schemas.get("publication")
            if validator:
                for err in validator.iter_errors(decision):
                    errors.append(f"{publication_path}[{index}]: {err.message}")
            task_id, run_id = decision.get("task_id"), decision.get("run_id")
            # Schema errors must be reported, not turned into an unhashable-value
            # exception by malformed JSON. Never use unvalidated IDs in paths.
            if not isinstance(task_id, str) or not isinstance(run_id, str):
                continue
            pair = (task_id, run_id)
            if pair in seen:
                errors.append(f"{publication_path}: duplicate publication decision {pair}")
            seen.add(pair)
            if (decision.get("publication_status") == "legacy-featured"
                    and pair != ("egfr-luad-evidence-coverage", "pilot-20261110")):
                errors.append(f"{publication_path}: legacy status is reserved for the original EGFR/LUAD pilot")
            manifest = root / "results" / task_id / run_id / "manifest.json"
            if not manifest.is_file():
                errors.append(f"{publication_path}: missing manifest for publication decision {pair}")
                continue
            if hashlib.sha256(manifest.read_bytes()).hexdigest() != decision.get("manifest_sha256"):
                errors.append(f"{publication_path}: stale publication pin for {pair}")

    # Synthetic examples demonstrate the schema, not a real reviewer attestation.
    for path in sorted((root / "examples" / "reviews").glob("*.json")):
        obj = validate(path, "review")
        if obj and obj.get("synthetic") is not True:
            errors.append(f"{path}: example review must have synthetic=true")

    # Bind every real review to exactly the result manifest bytes it evaluates.
    for path in sorted((root / "reviews").glob("*/*/*.json")):
        obj = validate(path, "review")
        if not obj:
            continue
        if obj.get("synthetic") is not False:
            errors.append(f"{path}: live reviews may not be synthetic")
        if obj.get("id") != path.stem:
            errors.append(f"{path}: review id must equal file stem")
        if obj.get("task_id") != path.parent.parent.name or obj.get("run_id") != path.parent.name:
            errors.append(f"{path}: review target does not match directory")
        result_path = root / "results" / path.parent.parent.name / path.parent.name / "manifest.json"
        if not result_path.is_file():
            errors.append(f"{path}: no matching submitted result manifest")
            continue
        actual_sha = hashlib.sha256(result_path.read_bytes()).hexdigest()
        if obj.get("reviewed_manifest_sha256") != actual_sha:
            errors.append(f"{path}: stale or mismatched reviewed manifest SHA-256")
        result_data = load_json(result_path, errors)
        if not isinstance(result_data, dict):
            continue
        known_findings = {finding.get("id") for finding in result_data.get("findings", [])
                          if isinstance(finding, dict)}
        for finding_id in obj.get("reviewed_finding_ids", []):
            if finding_id not in known_findings:
                errors.append(f"{path}: unknown reviewed finding id {finding_id}")

    return sorted(errors)


def main() -> int:
    errors = validate_repository()
    if errors:
        print("Validation failed:")
        for err in errors:
            print(f" - {err}")
        return 1
    print("OK: source catalog, tasks, and submission integrity checks passed.")
    print("Note: structural validation does not prove source truth, biomedical accuracy, or clinical suitability.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
