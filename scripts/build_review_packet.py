#!/usr/bin/env python3
"""Create a version-bound, deterministic, public-data-only EGFR reviewer kit.

No network requests are made. Clinical interpretations are outside the scope.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

try:
    from .validate import ROOT, validate_repository, MAX_ARTIFACT_BYTES
except ImportError:
    from validate import ROOT, validate_repository, MAX_ARTIFACT_BYTES

TASK = "egfr-luad-evidence-coverage"
RUN = "pilot-20261110"
RESULT = Path("results") / TASK / RUN
DOCS = Path("review_packets") / "egfr-luad-v1"
RAW_FILES = ("manifest.json", "coverage.json", "civic_evidence_ids.csv",
             "cbio_protein_changes.csv", "analyze.py")
DOC_FILES = ("README.md", "SCIENTIFIC_CHECKLIST.md")
SAMPLE_FIELDS = ("record_id", "evidence_type", "molecular_profile",
                 "disease", "doid", "record_url")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def checked_public_file(path: Path) -> bytes:
    if path.is_symlink() or not path.is_file() or path.stat().st_size > MAX_ARTIFACT_BYTES:
        raise ValueError(f"Invalid, symlinked or oversized review input: {path}")
    return path.read_bytes()


def csv_records(data: bytes, fields: list[str]) -> list[dict]:
    reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig")))
    if reader.fieldnames != fields:
        raise ValueError(f"Unexpected CSV header: {reader.fieldnames}")
    rows = list(reader)
    if not rows or any(None in row or any(v is None for v in row.values()) for row in rows):
        raise ValueError("CSV empty, incomplete or unexpectedly wide")
    return rows


def review_sample(rows: list[dict], manifest_sha: str) -> bytes:
    """Target rare evidence types for record-level review; not a random inference sample."""
    groups: dict[str, list[dict]] = defaultdict(list)
    for record in rows:
        if not record["record_id"].isdigit():
            raise ValueError("Unexpected CIViC evidence ID")
        url_prefix = "https://civicdb.org/links/evidence_items/"
        if (not record["record_url"].startswith(url_prefix)
                or record["record_url"][len(url_prefix):] != record["record_id"]):
            raise ValueError("Unexpected CIViC record URL or ID mismatch")
        groups[record["evidence_type"]].append(record)
    quotas = {"Predictive": 8, "Prognostic": 2, "Oncogenic": 2}
    if set(groups) != set(quotas):
        raise ValueError("Evidence types changed: review sample quotas must be revisited")
    selected = []
    for category in sorted(quotas):
        candidates = groups[category]
        if len(candidates) < quotas[category]:
            raise ValueError(f"Insufficient {category} items for sample")
        candidates.sort(key=lambda row: sha(
            f"{manifest_sha}:{row['record_id']}".encode("ascii")))
        selected.extend(candidates[:quotas[category]])
    selected.sort(key=lambda row: (row["evidence_type"], int(row["record_id"])))
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(SAMPLE_FIELDS), lineterminator="\n")
    writer.writeheader()
    for item in selected:
        writer.writerow({key: item[key] for key in SAMPLE_FIELDS})
    return stream.getvalue().encode("utf-8")


def git_revision(root: Path) -> str:
    if not (root / ".git").exists():
        raise ValueError("Build from a Git checkout to retain exact source revision")
    head = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    inputs = ["catalog", "tasks", "results", "reviews", "schemas",
              "review_packets", "scripts/build_review_packet.py", "scripts/validate.py"]
    dirty = subprocess.check_output(
        ["git", "status", "--porcelain", "--untracked-files=all", "--", *inputs],
        cwd=root, text=True).strip()
    if dirty:
        raise ValueError("Review inputs differ from the checked-out commit; commit or clean before building")
    return head


def make_packet(root: Path, revision: str) -> dict[str, bytes]:
    if len(revision) != 40 or any(ch not in "0123456789abcdef" for ch in revision):
        raise ValueError("Full lower-case Git SHA required")
    errors = validate_repository(root)
    if errors:
        raise ValueError("Repository validation failed:\n" + "\n".join(errors))
    path = root / RESULT
    packet: dict[str, bytes] = {}
    for name in RAW_FILES:
        packet[name] = checked_public_file(path / name)
    for name in DOC_FILES:
        packet[name] = checked_public_file(root / DOCS / name)
    manifest = json.loads(packet["manifest.json"])
    coverage = json.loads(packet["coverage.json"])
    if (manifest.get("task_id"), manifest.get("run_id")) != (TASK, RUN):
        raise ValueError("Wrong reviewed research task or run")
    if manifest.get("findings") != [] or coverage.get("status") != "unreviewed":
        raise ValueError("This kit only supports the unreviewed, nonclinical snapshot")
    artifacts = {a["path"]: a["sha256"] for a in manifest["artifacts"]}
    if set(artifacts) != {"coverage.json", "civic_evidence_ids.csv", "cbio_protein_changes.csv"}:
        raise ValueError("Unexpected manifest artifacts: review allowlist must be updated")
    for name, digest in artifacts.items():
        if sha(packet[name]) != digest:
            raise ValueError(f"Artifact checksum mismatch: {name}")
    civic = csv_records(packet["civic_evidence_ids.csv"],
                        ["source", "record_id", "molecular_profile", "disease",
                         "doid", "evidence_type", "record_url"])
    cbio = csv_records(packet["cbio_protein_changes.csv"],
                       ["protein_change_as_reported", "mutation_call_rows"])
    if len(civic) != coverage["civic"]["egfr_accepted_evidence_exact_luad"]:
        raise ValueError("CIViC evidence rows differ from coverage summary")
    if len({r["record_id"] for r in civic}) != len(civic) or any(r["source"] != "civic" for r in civic):
        raise ValueError("Duplicate or mislabeled CIViC evidence")
    if dict(Counter(r["evidence_type"] for r in civic)) != coverage["civic"]["evidence_type_counts_exact_luad"]:
        raise ValueError("CIViC category counts differ from coverage summary")
    if sum(int(r["mutation_call_rows"]) for r in cbio) != coverage["cbioportal"]["egfr_mutation_call_rows"]:
        raise ValueError("cBioPortal aggregate rows differ from coverage summary")
    gdc = coverage["gdc"]["egfr_occurrences"]
    if gdc["status"] == "retrieved":
        count = gdc["somatic_mutation_occurrence_rows"]
        distinct = gdc["distinct_case_ids_with_recorded_occurrence"]
        if not (type(count) is int and type(distinct) is int and 0 < distinct <= count <= gdc["limit"]):
            raise ValueError("GDC summary has inconsistent occurrence and case totals")
    packet["evidence_review_sample.csv"] = review_sample(civic, sha(packet["manifest.json"]))
    file_index = [{"name": name, "bytes": len(data), "sha256": sha(data)}
                  for name, data in sorted(packet.items())]
    meta = {
        "schema_version": 1,
        "kind": "external-scientific-review-inputs-not-clinical-validation",
        "task_id": TASK,
        "run_id": RUN,
        "source_revision": revision,
        "original_manifest_sha256": sha(packet["manifest.json"]),
        "source_response_sha256": coverage["source_content_sha256"],
        "qa_sample": {
            "selection": "SHA256 ordering of original manifest hash plus evidence ID; 8 predictive, 2 prognostic, 2 oncogenic; nonrepresentative",
            "size": 12,
        },
        "files": file_index,
    }
    packet["packet-index.json"] = (
        json.dumps(meta, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    ).encode("utf-8")
    packet["SHA256SUMS"] = "".join(
        f"{sha(data)}  {name}\n" for name, data in sorted(packet.items())
    ).encode("ascii")
    return packet


def save_packet(out: Path, files: dict[str, bytes], archive: Path | None = None) -> None:
    if out.exists():
        raise ValueError("Review output directory already exists; choose a new path")
    if archive is not None and archive.exists():
        raise ValueError("Review archive already exists; choose a new path")
    out.parent.mkdir(parents=True, exist_ok=True)
    if archive:
        archive.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out.parent) as temp:
        temp_root = Path(temp)
        stage = temp_root / "packet"
        stage.mkdir()
        for name, content in files.items():
            if not name or Path(name).name != name:
                raise ValueError("Unsafe archive member")
            (stage / name).write_bytes(content)
        if archive:
            zip_path = temp_root / "review.zip"
            with zipfile.ZipFile(zip_path, "w") as z:
                for name, content in sorted(files.items()):
                    member = zipfile.ZipInfo(name, (1980, 1, 1, 0, 0, 0))
                    member.compress_type = zipfile.ZIP_DEFLATED
                    member.external_attr = 0o644 << 16
                    z.writestr(member, content)
            zip_path.replace(archive)
        stage.replace(out)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--archive", type=Path, help="Optional deterministic ZIP destination")
    args = parser.parse_args()
    try:
        revision = git_revision(ROOT)
        packet = make_packet(ROOT, revision)
        save_packet(args.output, packet, args.archive)
    except (ValueError, OSError, subprocess.CalledProcessError, UnicodeError) as exc:
        parser.exit(1, f"Review packet failed: {exc}\n")
    print(f"Generated {len(packet)} files for original manifest "
          f"{json.loads(packet['packet-index.json'])['original_manifest_sha256']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
