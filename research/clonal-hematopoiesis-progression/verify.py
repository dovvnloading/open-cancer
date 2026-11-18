"""Integrity and source transcription checks; no clinical validity claim.

Offline checks use only stdlib and committed aggregates. Optional --online
fetches a specifically identified public model-summary workbook and the
separate aggregate coefficient CSV in memory. Only selected shared model
coefficients are compared; unselected fields are not output or persisted.
It does not enumerate repositories or download the dataset's bulk archive.
All textual file digests normalize CRLF and CR to LF before UTF-8 SHA-256.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
from urllib.request import urlopen
import xml.etree.ElementTree as ET
import zipfile

import analysis

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
TEST_PATH = "tests/test_research_clonal_hematopoiesis_progression.py"
MANIFEST_NAME = "manifest.json"
AML_SOURCE = "https://ndownloader.figshare.com/files/28917537"


def normalized_bytes(content: bytes | str) -> bytes:
    text = content.decode("utf-8") if isinstance(content, bytes) else content
    return text.replace("\r\n", "\n").replace("\r", "\n").encode("utf-8")


def digest(content: bytes | str) -> str:
    return hashlib.sha256(normalized_bytes(content)).hexdigest()


def tracked_payloads(root: Path = ROOT) -> dict[str, Path]:
    repo = root.parents[1]
    paths = {p.relative_to(repo).as_posix(): p for p in root.rglob("*")
             if p.is_file() and p.name != MANIFEST_NAME and "__pycache__" not in p.parts}
    paths[TEST_PATH] = repo / TEST_PATH
    return dict(sorted(paths.items()))


def write_manifest(root: Path = ROOT) -> dict:
    payload = {"format": "utf8-lf-sha256-v1",
               "meaning": "integrity of local textual artifacts; not source truth or clinical validity",
               "base_main_commit": "e6558b46e1f7c3ee59f590e0b93ceae48cea2f1c",
               "prespec_commit": "7409aa2",
               "files": {name: digest(path.read_bytes())
                         for name, path in tracked_payloads(root).items()}}
    (root / MANIFEST_NAME).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                                      encoding="utf-8", newline="\n")
    return payload


def verify_manifest(root: Path = ROOT) -> list[str]:
    manifest = json.loads((root / MANIFEST_NAME).read_text(encoding="utf-8"))
    if manifest["format"] != "utf8-lf-sha256-v1":
        raise AssertionError("unknown digest convention")
    actual = {name: digest(path.read_bytes()) for name, path in tracked_payloads(root).items()}
    if actual != manifest["files"]:
        missing = set(actual) ^ set(manifest["files"])
        changed = [name for name in set(actual) & set(manifest["files"])
                   if actual[name] != manifest["files"][name]]
        raise AssertionError(f"manifest mismatch: paths={sorted(missing)}, changed={sorted(changed)}")
    return list(actual)


def verify_provenance(root: Path = ROOT) -> dict[str, int]:
    required = ("source_url", "locator", "access_date", "reuse_basis")
    counts = {}
    for path in sorted((root / "data").glob("*.csv")):
        rows = analysis.read_csv(path)
        if not rows:
            raise AssertionError(f"empty evidence file: {path.name}")
        for row in rows:
            for key in required:
                if not row.get(key, "").strip():
                    raise AssertionError(f"missing {key}: {path.name}")
            if row["access_date"] != "2026-11-17":
                raise AssertionError("source snapshot uses the user-facing access date")
            if not row["source_url"].startswith("https://"):
                raise AssertionError("source URLs must be explicit")
            if path.name != "literature.csv" and "CC BY 4.0" not in row["reuse_basis"]:
                raise AssertionError("numeric evidence requires the documented redistribution basis")
        counts[path.name] = len(rows)
    return counts


def verify_derived_outcomes(root: Path = ROOT) -> int:
    count = 0
    for row in analysis.read_csv(root / "data" / "outcomes.csv"):
        if row["raw_log_estimate"]:
            beta, se = float(row["raw_log_estimate"]), float(row["log_standard_error"])
            expected = {"estimate": math.exp(beta), "lower": math.exp(beta - 1.96 * se),
                        "upper": math.exp(beta + 1.96 * se)}
            for key, value in expected.items():
                if not math.isclose(float(row[key]), value, rel_tol=2e-15, abs_tol=0):
                    raise AssertionError(f"wrong exponentiation {row['outcome_id']} {key}")
            count += 1
    return count


def verify_readme(root: Path = ROOT) -> dict:
    result = analysis.run(root)
    text = (root / "README.md").read_text(encoding="utf-8")
    block = analysis.BEGIN + text.split(analysis.BEGIN, 1)[1].split(analysis.END, 1)[0] + analysis.END
    if block != analysis.render_report(result):
        raise AssertionError("README quantitative section differs from CSV recomputation")
    return result


def _download(url: str) -> bytes:
    with urlopen(url, timeout=45) as response:
        return response.read()


def online_source_audit(root: Path = ROOT) -> dict:
    """Compare stored numbers with source cells/rows, retaining no source file.

    The fitness workbook contains model summaries; only specified shared
    gene/site coefficients are selected. No other columns are output or kept.
    The other download is the standalone 22-gene coefficient CSV, not the
    Figshare dataset archive or participant-level files.
    """
    fitness = analysis.read_csv(root / "data" / "fitness.csv")
    urls = {row["source_url"] for row in fitness}
    if len(urls) != 1:
        raise AssertionError("source audit expects the explicitly selected Table 6 source")
    namespace = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(io.BytesIO(_download(next(iter(urls))))) as archive:
        selected = {row[k] for row in fitness for k in ("estimate_cell", "lower_cell", "upper_cell")}
        xml = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        cells = {c.attrib["r"]: c.find("x:v", namespace).text
                 for c in xml.findall(".//x:c", namespace)
                 if c.attrib["r"] in selected and c.find("x:v", namespace) is not None}
    for row in fitness:
        for key in ("estimate", "lower", "upper"):
            upstream = float(cells[row[key + "_cell"]])
            if not math.isclose(float(row[key]), upstream, rel_tol=1e-12, abs_tol=1e-14):
                raise AssertionError(f"source mismatch: {row['fit_id']} {key}")
    text = _download(AML_SOURCE).decode("utf-8-sig")
    source_rows = {row["gene"]: row for row in csv.DictReader(io.StringIO(text))
                   if row["group"] == "Genes"}
    outcomes = [r for r in analysis.read_csv(root / "data" / "outcomes.csv")
                if r["study_id"] == "fabre2022_abelson_refit"]
    if set(source_rows) != {r["gene"] for r in outcomes}:
        raise AssertionError("gene omission/addition relative to published coefficient CSV")
    for row in outcomes:
        source = source_rows[row["gene"]]
        if float(row["raw_log_estimate"]) != float(source["coef"]):
            raise AssertionError(f"wrong raw log-HR: {row['gene']}")
        if float(row["log_standard_error"]) != float(source["sd"]):
            raise AssertionError(f"wrong raw SE: {row['gene']}")
    return {"fitness_numeric_cells_matched": 3 * len(fitness),
            "outcome_raw_numeric_cells_matched": 2 * len(outcomes),
            "source_gene_set_complete": True,
            "aml_source_lf_sha256": digest(text),
            "qualification": "source transcription only; rights and estimands require the documented audit"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--online", action="store_true", help="also check the two identified aggregate source files")
    parser.add_argument("--write-manifest", action="store_true", help="intentionally record current LF-normalized text digests")
    args = parser.parse_args()
    if args.write_manifest:
        write_manifest()
    result = {"status": "nonclinical_integrity_checks", "provenance_rows": verify_provenance(),
              "derived_outcomes_checked": verify_derived_outcomes(),
              "manifest_paths_checked": len(verify_manifest())}
    report = verify_readme()
    result["mechanistic_verdict"] = report["mechanistic_verdict"]
    result["proxy_verdict"] = report["proxy_verdict"]
    if args.online:
        result["source_audit"] = online_source_audit()
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
