#!/usr/bin/env python3
"""Summarize public catalog metadata monitoring in Actions, never an issue ledger.

The source observation history is a previous *workflow artifact*. Missing or
expired history is explicitly reported, not treated as evidence of no drift.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

SOURCES = ("cbioportal", "civic", "gdc")


def evaluate(report: dict, previous: dict | None = None) -> dict:
    if not isinstance(report, dict) or report.get("schema_version") != 1:
        raise ValueError("Expected version-1 catalog health report")
    if not isinstance(report.get("sources"), dict) or set(report["sources"]) != set(SOURCES):
        raise ValueError("Catalog health report has unexpected or missing sources")
    previous_sources = previous.get("sources", {}) if isinstance(previous, dict) else {}
    alerts = []
    warnings = []
    for name in SOURCES:
        row = report["sources"][name]
        if not isinstance(row, dict):
            raise ValueError(f"Invalid health row: {name}")
        old = previous_sources.get(name, {}) or {}
        if not isinstance(old, dict):
            raise ValueError(f"Invalid previous health row: {name}")
        if not row.get("ok"):
            alerts.append(f"{name}: endpoint failure: {row.get('error', 'unknown error')}")
        if row.get("terms_link_error"):
            alerts.append(f"{name}: source-rights link issue: {row['terms_link_error']}")
        drift = row.get("drift") or {}
        if not isinstance(drift, dict):
            raise ValueError(f"Invalid drift row: {name}")
        if drift.get("unexpected_shrinkage"):
            alerts.append(f"{name}: unexpected source-count shrinkage")
        if drift.get("added") or drift.get("removed") or drift.get("changed_fields"):
            alerts.append(
                f"{name}: metadata drift; additions={drift.get('added', [])}, "
                f"removals={drift.get('removed', [])}, "
                f"changed_fields={drift.get('changed_fields', {})}"
            )
        now_terms = row.get("terms_link") or {}
        old_terms = old.get("last_verified_terms_link") or old.get("terms_link") or {}
        pending = row.get("pending_terms_review")
        if pending is None:
            # Compatibility with reports that predate the durable review flag.
            changed = [
                field for field in ("url", "etag", "last_modified")
                if now_terms.get(field) and old_terms.get(field)
                and now_terms[field] != old_terms[field]
            ]
        else:
            if (not isinstance(pending, list)
                    or any(field not in ("url", "etag", "last_modified") for field in pending)):
                raise ValueError(f"Invalid terms review fields: {name}")
            changed = sorted(set(pending))
        for field in changed:
            # HTTP validators can change without a legal terms change.
            # Keep these visible for manual review, but do not fail CI.
            warnings.append(f"{name}: source-terms {field} changed; review reuse rights")
    report["alerts"] = alerts
    report["warnings"] = warnings
    report["previous_monitor_artifact_used"] = previous is not None
    return report


def render(report: dict) -> str:
    lines = [
        "# Public catalog metadata monitoring", "",
        f"Checked: {report['checked_at']}",
        "**Scope:** public metadata only; no patient records, scientific review or reuse-rights adjudication.",
        "",
        f"**Outcome:** {'ACTIONABLE DRIFT / FAILURE' if report['alerts'] else 'no actionable technical failures'}",
        f"**Historical comparison:** {'previous Actions report' if report.get('previous_monitor_artifact_used') else 'unavailable — initial/new baseline; past drift not assessed'}",
        "",
        "## Source checks", "",
    ]
    for name in SOURCES:
        row = report["sources"][name]
        drift = row.get("drift") or {}
        count = drift.get("upstream_count", "unavailable")
        lines.append(
            f"- {name}: {'reachable' if row.get('ok') else 'FAILED'}; "
            f"upstream public records {count}; "
            f"last healthy {row.get('last_healthy') or 'unknown'}"
        )
    lines += ["", "## Alerts", ""]
    lines.extend("- " + line.replace("\n", " ")[:600] for line in report["alerts"])
    if not report["alerts"]:
        lines.append("- None from these bounded checks.")
    lines += ["", "## Advisory source-terms metadata reviews", ""]
    lines.extend("- " + line.replace("\n", " ")[:600] for line in report.get("warnings", []))
    if not report.get("warnings"):
        lines.append("- None.")
    lines += [
        "", "HTTP ETag and Last-Modified changes are not proof of a legal terms change. "
        "Review the source's actual reuse conditions and update its `reviewed_on` "
        "date in `catalog/sources/` once reviewed; these notices do not fail CI.",
        "", "See the attached `catalog-health-report.json` artifact for exact "
        "public-source observations, query drift and terms-URL metadata. "
        "An artifact expires after 90 days; absence of an older artifact "
        "means historical changes cannot be checked for that interval.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--previous-report", type=Path)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--check", action="store_true", help="Exit nonzero if saved report has actionable alerts")
    args = parser.parse_args()

    report = json.loads(args.report.read_text(encoding="utf-8"))
    if args.check:
        if not isinstance(report.get("alerts"), list):
            parser.error("Report must be evaluated before checking alert status")
        if report["alerts"]:
            print("\n".join(report["alerts"]))
            return 1
        print("No actionable catalog health alerts in the current Actions artifact.")
        return 0

    previous = (
        json.loads(args.previous_report.read_text(encoding="utf-8"))
        if args.previous_report and args.previous_report.is_file() else None
    )
    evaluate(report, previous)
    args.report.write_text(json.dumps(report, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    if args.summary:
        args.summary.write_text(render(report), encoding="utf-8")
    print(f"Catalog monitoring summary: {len(report['alerts'])} alerts, "
          f"{len(report['warnings'])} advisory reviews; "
          f"prior artifact {'available' if previous is not None else 'missing'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
