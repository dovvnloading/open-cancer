#!/usr/bin/env python3
"""Optional worker adapter for an operator-controlled OpenAI-compatible inference endpoint.

The v1 protocol supplies source directory metadata. The v2 protocol supplies
validated, permission-checked CIViC record projections pinned by the coordinator.
Credentials are read exclusively from the environment.
"""
from __future__ import annotations

from http.client import HTTPException
import json
import os
from pathlib import Path
import sys
from urllib.error import HTTPError
from urllib.parse import urlsplit
from urllib.request import HTTPRedirectHandler, Request, build_opener

if __package__:
    from . import civic_evidence, grounded_claims
    from .compute_fanout import PROTOCOL, ROOT, check_directory_response, validate_job
else:
    import civic_evidence
    import grounded_claims
    from compute_fanout import PROTOCOL, ROOT, check_directory_response, validate_job

MAX_HTTP_BYTES = 1_048_576


class NoRedirect(HTTPRedirectHandler):
    """Keep the inference input and bearer credential on the configured endpoint."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


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


def infer(stage: str, messages: list[dict], env: dict[str, str], opener) -> dict:
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
    body = json.dumps({"model": model, "messages": messages,
                       "temperature": 0, "max_tokens": max_tokens}).encode("utf-8")
    headers = {"Content-Type": "application/json"}
    if env.get("OPEN_CANCER_LLM_KEY"):
        headers["Authorization"] = "Bearer " + env["OPEN_CANCER_LLM_KEY"]
    http_request = Request(endpoint, data=body, headers=headers, method="POST")
    try:
        open_request = build_opener(NoRedirect).open if opener is None else opener
        with open_request(http_request, timeout=timeout) as response:
            raw = response.read(MAX_HTTP_BYTES + 1)
    except HTTPError as exc:
        exc.close()
        raise RuntimeError("Inference endpoint request failed") from None
    except (HTTPException, OSError):
        # Server status lines can include untrusted text.
        raise RuntimeError("Inference endpoint request failed") from None
    if len(raw) > MAX_HTTP_BYTES:
        raise ValueError("Inference response exceeds size limit")
    return extract_response(json.loads(raw))


def grounded_work(request: dict, job: dict, env: dict[str, str], root: Path, opener) -> dict:
    stage = request["stage"]
    fields = {"protocol", "stage", "job", "evidence_bundle"}
    if stage == "critique":
        fields.add("candidate")
    if set(request) != fields or not civic_evidence.supports_job(job):
        raise ValueError("Unsupported grounded worker request")
    catalog = json.loads((root / "catalog/sources/civic.json").read_text(encoding="utf-8"))
    bundle = request["evidence_bundle"]
    civic_evidence.validate_bundle(bundle, catalog)
    records = civic_evidence.records_for_job(bundle, job)
    candidate = request.get("candidate")
    if stage == "critique":
        grounded_claims.check_grounded_response("investigate", candidate, bundle, job)
    if not records or (stage == "critique" and candidate["outcome"] == "insufficient-evidence"):
        return grounded_claims.insufficient_response(stage, bundle)

    context = {"question": job["question"], "dimensions": job["dimensions"],
               "retrieval_status": bundle["status"], "retrieved_at": bundle["retrieved_at"],
               "synthetic": bundle["synthetic"], "records": records}
    if stage == "critique":
        context["candidate"] = candidate
    messages = [{"role": "system", "content": grounded_claims.stage_instruction(stage)},
                {"role": "user", "content": json.dumps(context, sort_keys=True)}]
    answer = infer(stage, messages, env, opener)
    return grounded_claims.check_grounded_response(stage, answer, bundle, job, candidate)


def work(request: object, env: dict[str, str] | None = None, root: Path = ROOT,
         opener=None) -> dict:
    env = os.environ if env is None else env
    if not isinstance(request, dict) or request.get("stage") not in {"investigate", "critique"}:
        raise ValueError("Unsupported worker request")
    job = validate_job(request.get("job"))
    if request.get("protocol") == grounded_claims.GROUNDED_PROTOCOL:
        return grounded_work(request, job, env, root, opener)
    if request.get("protocol") != PROTOCOL or "evidence_bundle" in request:
        raise ValueError("Unsupported worker request")
    stage = request["stage"]
    if stage == "critique" and not isinstance(request.get("candidate"), dict):
        raise ValueError("Critique requires an investigator candidate")
    sources = source_context(job, root)
    answer = check_directory_response(stage, infer(stage, prompt(request, sources), env, opener),
                                      job, root)
    if stage == "investigate" and answer["synthetic"]:
        raise ValueError("Source directory model output must set synthetic=false")
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
