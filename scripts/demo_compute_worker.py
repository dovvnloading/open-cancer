#!/usr/bin/env python3
"""Synthetic demo of the compute worker JSON protocol. NEVER performs AI research."""
import json
import sys

if __package__:
    from . import grounded_claims
else:
    import grounded_claims

SOURCE_URLS = {
    "cbioportal": "https://www.cbioportal.org/",
    "civic": "https://civicdb.org/",
    "gdc": "https://gdc.cancer.gov/",
}


def main() -> None:
    request = json.load(sys.stdin)
    job = request["job"]
    if request["protocol"] == grounded_claims.GROUNDED_PROTOCOL:
        response = grounded_claims.fixture_response(request["stage"], request["evidence_bundle"],
                                                  job, request.get("candidate"))
    elif request["stage"] == "investigate":
        source = job["source_ids"][0]
        response = {
            "summary": "Synthetic execution fixture only; no source records were fetched or analyzed.",
            "evidence": [{"source_id": source, "url": SOURCE_URLS[source], "relation": "context"}],
            "limitations": ["This demonstration made no external requests and produced no biological findings."],
            "synthetic": True,
        }
    elif request["stage"] == "critique":
        response = {
            "verdict": "insufficient-evidence",
            "issues": ["Synthetic demo contains no source-level observations to verify or reproduce."],
        }
    else:
        raise SystemExit("Unknown stage")
    print(json.dumps(response))


if __name__ == "__main__":
    main()
