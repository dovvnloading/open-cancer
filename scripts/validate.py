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
    for kind in ("source", "task", "result"):
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
    task_paths = sorted((root / "tasks").glob("*.json"))
    if not source_paths:
        errors.append("catalog/sources/: no source records")
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

    return sorted(errors)


def main() -> int:
    errors = validate_repository()
    if errors:
        print("Validation failed:")
        for err in errors:
            print(f" - {err}")
        return 1
    print("OK: source catalog, tasks, and submission integrity checks passed.")
    print("Note: scientific claims and upstream rights still require human review.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
