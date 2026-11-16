"""Score the job-specific evidence packet and compare independent pilot executions.

A valid citation resolves within the supplied packet. A passing check establishes
only the asserted field or packet count. Free-text interpretations are unassessed.
"""
from __future__ import annotations

from collections import Counter
from decimal import Decimal
import re

if __package__:
    from .compute_fanout import digest
    from .pilot_evidence import FIELDS as RECORD_FIELDS
else:
    from compute_fanout import digest
    from pilot_evidence import FIELDS as RECORD_FIELDS

CATEGORIES = ("unsupported-claim", "terminology", "source-interpretation",
              "missing-denominator", "contradiction")
STAGES = ("investigate", "critique")
TOKEN_FIELDS = ("prompt_tokens", "completion_tokens", "total_tokens")


# Deliberately constrained to affirmative, exact packet counts. Numbers in variant
# names, evidence IDs, exon numbers and citations are not quantitative assertions.
PACKET_COUNT = re.compile(
    r"\b(?:this|the) (?:supplied )?packet (?:contains|includes|has) "
    r"(?:exactly )?(\d+) (?:CIViC )?evidence items?\b", re.IGNORECASE)
SOURCE_EXAMINATION = re.compile(
    r"\b(?:we|I|this investigation) (?:directly )?(?:examined|retrieved|analy[sz]ed|inspected) "
    r"(?:the )?(?:public )?(GDC|cBioPortal) "
    r"(?:(?:mutation.call|mutation|occurrence|cohort) )?(?:records|rows|calls)\b", re.IGNORECASE)
COHORT_QUANTITY = re.compile(
    r"\b\d+(?:\.\d+)?(?:\s*%|\s+(?:patients?|cases?|samples?)(?:\b))", re.IGNORECASE)


def citation_key(record: dict) -> tuple:
    return record["source_id"], record["record_id"], record["url"]


def _check(assertion, cited_records, evidence):
    """Return status, recomputed value (if known), and an explicit decision rule."""
    if assertion["kind"] == "record-field":
        if assertion["field"] not in RECORD_FIELDS:
            return "unassessed", None, "field outside projected biomedical fields"
        matches = [r for r in cited_records if r["record_id"] == assertion["record_id"]]
        if len(matches) != 1:
            return "unassessed", None, "requires one resolvable citation for this record ID"
        actual = matches[0].get(assertion["field"])
        if actual is None or actual == "":
            return "unassessed", None, "source field missing"
        passed = type(actual) is type(assertion["value"]) and actual == assertion["value"]
        return "passed" if passed else "failed", actual, "exact projected field equality"

    if assertion["kind"] == "packet-count":
        source = assertion["source_id"]
        if evidence["sources"].get(source, {}).get("status") not in {"retrieved", "fixture"}:
            return "unassessed", None, "source records unavailable"
        if source != "civic" or assertion.get("unit") != "evidence-items":
            return "failed", None, "CIViC packet measures evidence items only"
        actual = sum(r["source_id"] == source for r in evidence["records"])
        passed = type(assertion["value"]) is int and assertion["value"] == actual
        return "passed" if passed else "failed", actual, "count of supplied rows for source"
    return "unassessed", None, "unsupported deterministic check"


def evaluate_job(job: dict, evidence: dict | None, candidate: dict | None = None,
                 critique: dict | None = None) -> dict:
    """Score structurally validated role outputs, retaining valid partial stages on failure."""
    packet = evidence if evidence is not None else {"sources": {}, "records": []}
    packet_records = {citation_key(record): record for record in packet["records"]}
    citations, checks, observations = [], [], []

    def resolve(refs, role, item):
        resolved = {}
        for ref in refs:
            key = citation_key(ref)
            valid = key in packet_records and ref["source_id"] in job["source_ids"]
            citations.append({**ref, "role": role, "item": item, "valid": valid})
            if valid:
                resolved[key] = packet_records[key]
        return list(resolved.values())

    if candidate is not None:
        for section in ("observations", "contradictory_evidence"):
            for observation in candidate[section]:
                item = observation["id"]
                citation_start = len(citations)
                cited = resolve(observation["citations"], "investigate", item)
                reasons = []
                check_start = len(checks)
                if not cited:
                    reasons.append("no resolvable record citation")
                if any(not c["valid"] for c in citations[citation_start:]):
                    reasons.append("contains unsupported citation")
                assertions = [("model", check) for check in observation["checks"]]
                for match in PACKET_COUNT.finditer(observation["text"]):
                    assertion = {"kind": "packet-count", "source_id": "civic",
                                 "unit": "evidence-items", "value": int(match[1])}
                    if assertion not in observation["checks"]:
                        assertions.append(("prose", assertion))
                for origin, assertion in assertions:
                    status, actual, rule = _check(assertion, cited, packet)
                    checks.append({"item": item, "origin": origin, "assertion": assertion,
                                   "status": status,
                                   "recomputed": actual, "rule": rule})
                    if status != "passed":
                        reasons.append(f"machine check {status}")
                for match in SOURCE_EXAMINATION.finditer(observation["text"]):
                    source = match[1].lower()
                    status = packet["sources"].get(source, {}).get("status")
                    if status == "not-retrieved":
                        checks.append({"item": item, "origin": "prose",
                                       "assertion": {"kind": "source-examined", "source_id": source},
                                       "status": "failed", "recomputed": status,
                                       "rule": "source not retrieved by this investigation"})
                        reasons.append("claims examination of unavailable source")
                observations.append({"item": item, "section": section,
                                     "has_valid_citation": bool(cited), "checks": len(checks) - check_start,
                                     "cohort_quantity_requires_review": bool(
                                         COHORT_QUANTITY.search(observation["text"])),
                                     "unsupported_reasons": sorted(set(reasons))})

    issues = Counter()
    if critique is not None:
        for number, issue in enumerate(critique["issues"], 1):
            issues[issue["category"]] += 1
            resolve(issue["citations"], "critique", f"issue-{number}")
    sources = {}
    for source in job["source_ids"]:
        metadata = packet["sources"].get(source, {})
        sources[source] = {
            "status": metadata.get("status", "unknown"),
            "retrieval_status": metadata.get("retrieval_status"),
            "matching_records": metadata.get("matching_records"),
            "supplied_ids": sorted(k[1] for k in packet_records if k[0] == source),
            "cited_ids": sorted({r["record_id"] for r in citations if
                                 r["source_id"] == source and r["valid"] and
                                 r["role"] == "investigate"}),
        }
    return {"citations": citations, "checks": checks, "observations": observations,
            "critic_issues_by_category": {key: issues[key] for key in CATEGORIES},
            "source_coverage": sources,
            "output_sha256": {"investigate": digest(candidate) if candidate is not None else None,
                              "critique": digest(critique) if critique is not None else None}}


def _coverage(records):
    source_ids = sorted({source for r in records for source in r["evaluation"]["source_coverage"]})
    coverage = {}
    for source in source_ids:
        by_job = [r["evaluation"]["source_coverage"][source] for r in records
                  if source in r["evaluation"]["source_coverage"]]
        supplied = {rid for row in by_job for rid in row["supplied_ids"]}
        cited = {rid for row in by_job for rid in row["cited_ids"]}
        coverage[source] = {
            "jobs_requested": len(by_job),
            "jobs_with_records": sum(bool(row["supplied_ids"]) for row in by_job),
            "jobs_live_retrieved": sum(row["status"] == "retrieved" for row in by_job),
            "jobs_with_fixture": sum(row["status"] == "fixture" for row in by_job),
            "jobs_without_retrieval": sum(row["status"] not in {"retrieved", "fixture"}
                                          for row in by_job),
            "jobs_without_matches": sum(row["matching_records"] == 0 and
                                        row["retrieval_status"] == "complete" for row in by_job),
            "jobs_with_truncated_retrieval": sum(row["retrieval_status"] == "truncated"
                                                 for row in by_job),
            "jobs_without_matches_before_truncation": sum(row["matching_records"] == 0 and
                                                          row["retrieval_status"] == "truncated"
                                                          for row in by_job),
            "unique_records_supplied": len(supplied),
            "unique_records_cited_by_investigator": len(cited),
            "fraction_supplied_cited": len(cited) / len(supplied) if supplied else None,
        }
    return coverage


def validate_prices(prices):
    """Reject mistaken prices rather than emitting a misleading cost calculation."""
    if not isinstance(prices, dict) or not prices.keys() <= set(STAGES):
        raise ValueError("prices must be keyed by stage")
    for rates in prices.values():
        if not isinstance(rates, dict) or not rates.keys() <= {
            "input_per_million_usd", "output_per_million_usd"
        }:
            raise ValueError("expected per-million input/output USD rates")
        for price in rates.values():
            if price is not None and (type(price) not in {int, float} or
                                      not Decimal(str(price)).is_finite() or price < 0):
                raise ValueError("token prices must be finite, nonnegative numbers")


def _tokens(usage):
    counts = {}
    for field in TOKEN_FIELDS:
        value = (usage or {}).get(field)
        counts[field] = value if type(value) is int and value >= 0 else None
    if counts["total_tokens"] is None and None not in (counts["prompt_tokens"],
                                                       counts["completion_tokens"]):
        counts["total_tokens"] = counts["prompt_tokens"] + counts["completion_tokens"]
    return counts


def _resources(attempts, prices, mode):
    # A failed request may still be billable. Price all endpoint attempts, including failures.
    possible_calls = [a for a in attempts if mode == "real" and a["endpoint_attempted"] is not False]
    dispatched = sum(a["endpoint_attempted"] is True for a in possible_calls)
    totals = Counter()
    known = Counter()
    subtotal = Decimal(0)
    priced_components = 0
    for attempt in possible_calls:
        usage = _tokens(attempt.get("usage"))
        for field, value in usage.items():
            if value is not None:
                totals[field] += value
                known[field] += 1
        rates = prices.get(attempt["stage"], {})
        for field, price_key in (("prompt_tokens", "input_per_million_usd"),
                                 ("completion_tokens", "output_per_million_usd")):
            if usage[field] is not None and rates.get(price_key) is not None:
                subtotal += Decimal(str(rates[price_key])) * usage[field] / 1_000_000
                priced_components += 1
    elapsed = [a["elapsed_seconds"] for a in attempts if a["elapsed_seconds"] is not None]
    complete_cost = possible_calls and priced_components == 2 * len(possible_calls)
    return {
        "attempts": len(attempts),
        "failed_attempts": sum(a.get("error") is not None for a in attempts),
        "endpoint_attempts": dispatched,
        "attempts_of_unknown_dispatch": len(possible_calls) - dispatched,
        "elapsed_attempt_seconds": sum(elapsed) if len(elapsed) == len(attempts) else None,
        "known_elapsed_attempt_seconds": sum(elapsed) if elapsed or not attempts else None,
        "attempts_with_elapsed_seconds": len(elapsed),
        "usage": {k: totals[k] if possible_calls and known[k] == len(possible_calls)
                  else None for k in TOKEN_FIELDS},
        "known_token_subtotals": {k: totals[k] if known[k] else None for k in TOKEN_FIELDS},
        "attempts_with_token_counts": {k: known[k] for k in TOKEN_FIELDS},
        "estimated_cost_usd": float(subtotal) if complete_cost else None,
        "known_estimated_cost_subtotal_usd": float(subtotal) if priced_components else None,
        "priced_input_or_output_components": priced_components,
    }


def aggregate(records: list[dict], prices: dict | None = None, mode: str = "fixture") -> dict:
    """Count both roles' citation occurrences, distinct supplied records, and all attempts."""
    if mode not in {"fixture", "real"}:
        raise ValueError("unknown pilot mode")
    prices = {} if prices is None else prices
    validate_prices(prices)
    status = Counter(r["status"] for r in records)
    if status.keys() - {"failed", "completed"}:
        raise ValueError("unexpected pilot job status")
    attempts = [a for r in records for a in r["attempts"]]
    if any(a["stage"] not in STAGES or a["synthetic"] is not (mode == "fixture") for a in attempts):
        raise ValueError("attempt stage or synthetic flag contradicts run configuration")
    evaluations = [r["evaluation"] for r in records]
    observations = [o for e in evaluations for o in e["observations"]]
    citations = [c for e in evaluations for c in e["citations"]]
    checks = Counter(x["status"] for e in evaluations for x in e["checks"])
    execution = _resources(attempts, prices, mode)
    execution["by_stage"] = {s: _resources([a for a in attempts if a["stage"] == s], prices, mode)
                             for s in STAGES}
    return {
        "jobs": {"total": len(records), "completed": status["completed"], "failed": status["failed"],
                 "with_candidate": sum(e["output_sha256"]["investigate"] is not None
                                       for e in evaluations),
                 "with_critique": sum(e["output_sha256"]["critique"] is not None
                                      for e in evaluations)},
        "citations": {"valid": sum(c["valid"] for c in citations),
                      "unsupported": sum(not c["valid"] for c in citations)},
        "deterministic_checks": {k: checks[k] for k in ("passed", "failed", "unassessed")},
        "observations": {"total": len(observations),
                         "without_checkable_assertions": sum(not o["checks"] for o in observations),
                         "with_cohort_quantity_to_review": sum(
                             o["cohort_quantity_requires_review"] for o in observations),
                         "flagged_unsupported": sum(bool(o["unsupported_reasons"])
                                                    for o in observations),
                         "with_valid_citation": sum(o["has_valid_citation"] for o in observations),
                         "contradictory_evidence_entries":
                             sum(o["section"] == "contradictory_evidence" for o in observations)},
        "critic_issues_by_category": {k: sum(e["critic_issues_by_category"][k] for e in evaluations)
                                      for k in CATEGORIES},
        "source_coverage": _coverage(records),
        "execution": execution,
    }


def _provider_comparison(first, second):
    """Compare reported identifiers on the two successful attempts that produced each output."""
    by_stage = {}
    for stage in STAGES:
        served = []
        for body in (first, second):
            attempt = next((a for a in reversed(body["attempts"])
                            if a["stage"] == stage and a["error"] is None and
                            a["synthetic"] is False and a["endpoint_attempted"] is True), {})
            served.append(attempt)
        by_stage[stage] = {}
        for field in ("returned_model", "system_fingerprint"):
            left, right = (a.get(field) for a in served)
            by_stage[stage][field] = {"first": left, "second": right,
                                      "differ": left != right if left and right else None}
    return by_stage


def compare_outputs(first: list[dict], second: list[dict]) -> dict:
    """Compare outputs only after the caller verifies independent, configuration-matched runs."""
    def index(records):
        by_id = {r["job"]["id"]: r for r in records}
        if len(by_id) != len(records):
            raise ValueError("duplicate job in comparison")
        return by_id

    left, right = index(first), index(second)
    pairs = []
    overlaps = []
    for job_id in sorted(left.keys() | right.keys()):
        a, b = left.get(job_id), right.get(job_id)
        pair = {"id": job_id, "first_status": a["status"] if a else "missing",
                "second_status": b["status"] if b else "missing"}
        if pair["first_status"] == pair["second_status"] == "completed":
            for role in ("candidate", "critique"):
                pair[f"{role}_exact"] = digest(a[role]) == digest(b[role])
            pair["provider_identifiers"] = _provider_comparison(a, b)
            pair["verdict_match"] = a["critique"]["verdict"] == b["critique"]["verdict"]
            sets = [{citation_key(c) for c in r["evaluation"]["citations"]
                     if c["valid"] and c["role"] == "investigate"} for r in (a, b)]
            union = sets[0] | sets[1]
            jaccard = len(sets[0] & sets[1]) / len(union) if union else None
            pair["cited_record_jaccard"] = jaccard
            if jaccard is not None:
                overlaps.append(jaccard)
        pairs.append(pair)
    compared = [p for p in pairs if "candidate_exact" in p]
    exact = sum(p["candidate_exact"] and p["critique_exact"] for p in compared)
    return {"job_ids_total": len(pairs), "paired_completed": len(compared),
            "missing_or_failed_pairs": len(pairs) - len(compared), "both_outputs_exact": exact,
            "exact_fraction": exact / len(compared) if compared else None,
            "verdict_matches": sum(p["verdict_match"] for p in compared),
            "pairs_with_cited_records": len(overlaps),
            "mean_cited_record_jaccard": sum(overlaps) / len(overlaps) if overlaps else None,
            "pairs": pairs}
