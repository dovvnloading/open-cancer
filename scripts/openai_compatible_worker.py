#!/usr/bin/env python3
"""Optional worker adapter for an operator-controlled OpenAI-compatible inference endpoint.

This worker receives only checked-in source DIRECTORY metadata, never source
record bodies. It can propose hypotheses but cannot verify cancer data or claims.
Secrets are read exclusively from environment variables, not command arguments.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import sys
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

if __package__:
    from .compute_fanout import PROTOCOL, ROOT, check_response, validate_job
else:
    from compute_fanout import PROTOCOL, ROOT, check_response, validate_job

MAX_HTTP_BYTES = 1_048_576


def source_context(job: dict, root: Path = ROOT) -> list[dict]:
    """Only trusted source directory links and limitations, not record data."""
    sources = []
    for source_id in job["source_ids"]:
        path = root / "catalog" / "sources" / (source_id + ".json")
        record = json.loads(path.read_text(encoding="utf-8"))
        if record.get("id") != source_id:
            raise ValueError(f"Catalog source mismatch: {source_id}")
        sources.append({key: record[key] for key in ("id", "name", "summary",
                           "homepage", "documentation_url") if key in record})
    return sources


def prompt(request: dict, sources: list[dict]) -> list[dict]:
    stage = request["stage"]
    if stage == "investigate":
        instruction = (
            "You are an AI computational research ideation worker. You have source DIRECTORY "
            "metadata ONLY, not retrieved evidence records. Propose an explicitly unverified "
            "hypothesis or data-quality check; never state you fetched or confirmed source "
            "records. No treatments, clinical advice, expert review or biological proof. "
            "Return ONLY a JSON object with keys: summary (string >=20 chars), evidence "
            "(nonempty list of {source_id,url,relation} using ONLY supplied source "
            "homepage/documentation_url and relation='context'), limitations (nonempty list "
            "of detailed strings), synthetic (false). No markdown."
        )
    else:
        instruction = (
            "You are the adversarial AI critique worker. Identify unsupported assumptions, "
            "missing records, denominators and conflicts. Source directory metadata alone "
            "does not verify a hypothesis: the verdict MUST be insufficient-evidence. "
            "Return ONLY a JSON object with verdict='insufficient-evidence' and issues "
            "as a nonempty list of detailed strings. No clinical guidance or markdown."
        )
    return [{"role": "system", "content": instruction},
            {"role": "user", "content": json.dumps({"question": request["job"]["question"],
                "dimensions": request["job"]["dimensions"], "source_directory": sources,
                "candidate": request.get("candidate") if stage == "critique" else None},
                sort_keys=True)}]


def validate_endpoint(endpoint: str) -> None:
    parsed = urlsplit(endpoint)
    if parsed.scheme not in {"http", "https"} or not parsed.hostname or parsed.username or \
            parsed.password or parsed.query or parsed.fragment:
        raise ValueError("OPEN_CANCER_LLM_URL must be a clean HTTP(S) inference endpoint")
    if parsed.scheme == "http" and parsed.hostname not in {"localhost", "127.0.0.1", "::1"}:
        raise ValueError("Non-local inference endpoints must use HTTPS")
    if not parsed.path.endswith("/chat/completions"):
        raise ValueError("Endpoint must end with /chat/completions")


def extract_response(payload: object) -> dict:
    if not isinstance(payload, dict):
        raise ValueError("Inference response must be JSON")
    try:
        content = payload["choices"][0]["message"]["content"]
        if not isinstance(content, str):
            raise ValueError("Model response content must be a JSON string")
        answer = json.loads(content)
    except (KeyError, IndexError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError("Model response must contain JSON message content") from exc
    if not isinstance(answer, dict):
        raise ValueError("Model response must contain a JSON object")
    return answer


def work(request: object, env: dict[str, str] | None = None, root: Path = ROOT,
         opener=urlopen) -> dict:
    env = os.environ if env is None else env
    if not isinstance(request, dict) or request.get("protocol") != PROTOCOL or \
            request.get("stage") not in {"investigate", "critique"}:
        raise ValueError("Unsupported worker request")
    job = validate_job(request.get("job"))
    stage = request["stage"]
    if stage == "critique" and not isinstance(request.get("candidate"), dict):
        raise ValueError("Critique requires an investigator candidate")
    endpoint = env.get("OPEN_CANCER_LLM_URL", "")
    validate_endpoint(endpoint)
    model = (env.get("OPEN_CANCER_INVESTIGATOR_MODEL") if stage == "investigate" else
             env.get("OPEN_CANCER_CRITIC_MODEL")) or env.get("OPEN_CANCER_LLM_MODEL", "")
    if not model:
        raise ValueError("Set OPEN_CANCER_LLM_MODEL (or stage-specific model names)")
    timeout = int(env.get("OPEN_CANCER_LLM_TIMEOUT", "120"))
    max_tokens = int(env.get("OPEN_CANCER_MAX_TOKENS", "1200"))
    if not 1 <= timeout <= 3600 or not 128 <= max_tokens <= 32768:
        raise ValueError("Invalid timeout/token budget")
    sources = source_context(job, root)
    body = json.dumps({"model": model, "messages": prompt(request, sources),
                       "temperature": 0, "max_tokens": max_tokens}).encode("utf-8")
    headers = {"Content-Type": "application/json"}
    if env.get("OPEN_CANCER_LLM_KEY"):
        headers["Authorization"] = "Bearer " + env["OPEN_CANCER_LLM_KEY"]
    http_request = Request(endpoint, data=body, headers=headers, method="POST")
    try:
        with opener(http_request, timeout=timeout) as response:
            raw = response.read(MAX_HTTP_BYTES + 1)
    except (HTTPError, URLError, TimeoutError, OSError) as exc:
        # Never echo URLs, response bodies, credentials, or source payloads to logs.
        raise RuntimeError("Inference endpoint request failed") from exc
    if len(raw) > MAX_HTTP_BYTES:
        raise ValueError("Inference response exceeds size limit")
    candidate = extract_response(json.loads(raw))
    answer = check_response(stage, candidate, job)
    if stage == "investigate":
        allowed = {(entry["id"], entry[key]) for entry in sources
                   for key in ("homepage", "documentation_url") if key in entry}
        if answer["synthetic"] or any(x["relation"] != "context" or
                                     (x["source_id"], x["url"]) not in allowed
                                     for x in answer["evidence"]):
            raise ValueError("Source directory only: cite contextual catalog URLs, no claims of checked evidence")
    else:
        # The adapter has no record-level inputs. Never treat AI agreement as verification.
        if answer["verdict"] != "insufficient-evidence":
            answer = {"verdict": "insufficient-evidence", "issues": [
                "Source directory metadata cannot validate record-level or biological claims."]}
    return answer


def main() -> int:
    try:
        response = work(json.load(sys.stdin))
        print(json.dumps(response, sort_keys=True))
        return 0
    except (ValueError, TypeError, KeyError, OSError, RuntimeError, json.JSONDecodeError):
        print("AI worker failed; check endpoint config, response schema and source context.", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
