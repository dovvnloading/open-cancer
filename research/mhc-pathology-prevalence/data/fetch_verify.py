#!/usr/bin/env python3
"""Derive narrow species aggregates from three pinned public UTF-8 sources.

Nonclinical, unreviewed research; no medical or veterinary advice.
Network access is opt-in with --fetch and uses GitHub raw text, never a browser
download. Raw source assets are kept in memory and are not written to this repo.
Both repositories lack an explicit licence; both final articles are CC BY-NC-ND
4.0. This implements the user's narrow factual-aggregate fallback, not a claim
that the article licences permit adapted source assets or relicense the data.
"""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal, InvalidOperation, localcontext
import hashlib
import io
import json
import math
from pathlib import Path
import urllib.request


BASE = Path(__file__).resolve().parent
STATUS = "nonclinical_unreviewed_research_no_medical_or_veterinary_advice"
EXPORT_TOLERANCE = Decimal("2e-15")
SOURCE_FILES = {
    Path("vincze.csv").as_posix(): {
        "repository": "OrsolyaVincze/VinczeEtalNatComm_MHC-cancer",
        "commit": "2c8a5c5b5fda6fb34c86b0e5a2b6a58dcc85cc3e",
        "path": Path("SourceData_FINAL.csv").as_posix(),
        "git_blob_sha1": "8a937a2c8d4f1749ae33c1b12d2ed52e895852f9",
        "sha256_lf": "b9d6dedc9447e9f99b275cbfbdb54e3b0e318d08c12f0384572a65101d057744",
    },
    Path("compton.csv").as_posix(): {
        "repository": "zacharycompton/cancerAcrossVertebrates",
        "commit": "4819e0933358c5fc7eab483a0c05161ea8c074b7",
        "path": Path("species-cancer-prevalence-data.csv").as_posix(),
        "git_blob_sha1": "d5e18150557a80b2dff1f53779182abd99e7d5ec",
        "sha256_lf": "70f8e2f6675e4c42dcf87acb0b826bb7c4ea32fed27a46900c1cfdf0bce70390",
    },
    Path("dictionary.txt").as_posix(): {
        "repository": "zacharycompton/cancerAcrossVertebrates",
        "commit": "4819e0933358c5fc7eab483a0c05161ea8c074b7",
        "path": Path("data-description-species-cancer-prevalence-data.txt").as_posix(),
        "git_blob_sha1": "293705259ec5d0e4fabad000e183e162057e14c5",
        "sha256_lf": "39eec5bd11bd14b09ddbf6f3248322869711dd4968e853dafe0a238724a36140",
    },
}
AGGREGATE_FIELDS = (
    "species", "order", "log_body_mass_kg", "log_adult_life_expectancy_days",
    "mammal_prey", "log1p_mhc_i", "records", "neoplasia", "malignant",
)
RECONCILIATION_FIELDS = (
    "normalized_species", "vincze_species", "compton_species", "compton_class",
    "vincze_order", "compton_order", "status", "reason",
)
AUDIT_FIELDS = (
    "species", "normalized_species", "neoplasia_residual", "neoplasia_tolerance",
    "malignant_residual", "malignancyknown_residual", "malignancy_tolerance",
)
PREDICTOR_FIELDS = (
    "Body mass (kg)", "Adult life expectancy", "Consuming mammal prey",
    "MHC I gene family size",
)


def normalize_newlines(text):
    """Change line endings only, without trimming or adding a final newline."""
    return text.replace("\r\n", "\n").replace("\r", "\n")


def normalized_sha256(text):
    return hashlib.sha256(normalize_newlines(text).encode("utf-8")).hexdigest()


def normalize_species(name):
    """Exact normalized string identity; no synonyms or taxonomic inference."""
    return " ".join(name.replace("_", " ").split()).casefold()


def load_sources(source_dir=None, fetch=False):
    """Verify LF-normalized hashes; preserve no raw files in the repository.

    The Git blob hashes are provenance metadata from the GitHub API. Content
    validation uses only the LF-normalized SHA-256 so a local CRLF/CR copy is
    equivalent on Windows and Linux.
    """
    if bool(source_dir is not None) == bool(fetch):
        raise ValueError("Choose exactly one of source_dir or fetch")
    texts = {}
    for name, source in SOURCE_FILES.items():
        if fetch:
            url = ("https://raw.githubusercontent.com/" + source["repository"]
                   + "/" + source["commit"] + "/" + source["path"])
            request = urllib.request.Request(
                url, headers={"User-Agent": "open-cancer-pathology-audit/1"})
            with urllib.request.urlopen(request, timeout=30) as response:
                raw = response.read(2_000_001)
            if len(raw) > 2_000_000:
                raise ValueError("Unexpectedly large UTF-8 source: " + name)
        else:
            raw = (Path(source_dir) / name).read_bytes()
        text = normalize_newlines(raw.decode("utf-8"))
        if normalized_sha256(text) != source["sha256_lf"]:
            raise ValueError("Pinned normalized source hash mismatch: " + name)
        texts[Path(name).as_posix()] = text
    return texts


def parse_csv(text, label):
    records = []
    for raw in csv.DictReader(io.StringIO(normalize_newlines(text))):
        if None in raw or any(value is None for value in raw.values()):
            raise ValueError("Malformed source CSV: " + label)
        records.append({key.strip(): value.strip() for key, value in raw.items()})
    indexed = {}
    for row in records:
        key = normalize_species(row["Species"])
        if not key or key in indexed:
            raise ValueError("Empty or duplicate normalized species in " + label)
        indexed[key] = row
    return indexed


def decimal_number(token, field):
    try:
        value = Decimal(token)
    except InvalidOperation as error:
        raise ValueError("Invalid numeric source value: " + field) from error
    if not value.is_finite():
        raise ValueError("Nonfinite numeric source value: " + field)
    return value


def integer(token, field, minimum=0):
    value = decimal_number(token, field)
    if value != value.to_integral_value() or value < minimum:
        raise ValueError("Invalid integer source value: " + field)
    return int(value)


def rounding_tolerance(token):
    value = decimal_number(token, "printed prevalence")
    if not 0 <= value <= 1:
        raise ValueError("Printed prevalence outside [0,1]")
    # A literal endpoint represents an observed endpoint, not rounding to an
    # integer probability. This prevents a spurious tolerance of 0.5 at zero.
    half_unit = (Decimal(0) if value in (0, 1)
                 else Decimal(1).scaleb(value.as_tuple().exponent) / 2)
    return half_unit + EXPORT_TOLERANCE


def endpoint_audit_row(row):
    """Derived arithmetic discrepancies, not an upstream dictionary rewrite."""
    n = integer(row["RecordsWithDenominators"], "RecordsWithDenominators", 1)
    k = integer(row["NeoplasiaWithDenominators"], "NeoplasiaWithDenominators")
    m = integer(row["Malignant"], "Malignant")
    known = integer(row["MalignancyKnown"], "MalignancyKnown")
    if not 0 <= m <= k <= n:
        raise ValueError("Invalid neoplasia/malignancy counts: " + row["Species"])
    with localcontext() as context:
        context.prec = 60
        neo = decimal_number(row["NeoplasiaPrevalence"], "NeoplasiaPrevalence")
        malignancy = decimal_number(row["MalignancyPrevalence"], "MalignancyPrevalence")
        return {
            "species": row["Species"],
            "normalized_species": normalize_species(row["Species"]),
            "neoplasia_residual": str(neo - Decimal(k) / n),
            "neoplasia_tolerance": str(rounding_tolerance(row["NeoplasiaPrevalence"])),
            "malignant_residual": str(malignancy - Decimal(m) / n),
            "malignancyknown_residual": str(malignancy - Decimal(known) / n),
            "malignancy_tolerance": str(rounding_tolerance(row["MalignancyPrevalence"])),
        }


def csv_content(rows, fieldnames):
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def summarize(aggregates, reconciliation, endpoint_audit):
    """All numerical prefit summaries can be recovered from committed CSVs."""
    selected_keys = {normalize_species(row["species"]) for row in aggregates}
    panel_audit = [row for row in endpoint_audit if row["normalized_species"] in selected_keys]
    def mismatches(rows, residual, tolerance):
        return sum(abs(Decimal(row[residual])) > Decimal(row[tolerance]) for row in rows)
    orders = sorted({row["order"] for row in aggregates})
    return {
        "vincze_species": sum(bool(row["vincze_species"]) for row in reconciliation),
        "compton_species": sum(bool(row["compton_species"]) for row in reconciliation),
        "normalized_matches": sum(bool(row["vincze_species"] and row["compton_species"])
                                  for row in reconciliation),
        "union_species": len(reconciliation),
        "primary_species": len(aggregates),
        "orders": len(orders),
        "species_per_order": {order: sum(row["order"] == order for row in aggregates)
                              for order in orders},
        "records": sum(int(row["records"]) for row in aggregates),
        "neoplasia": sum(int(row["neoplasia"]) for row in aggregates),
        "malignant": sum(int(row["malignant"]) for row in aggregates),
        "zero_neoplasia_species": sum(int(row["neoplasia"]) == 0 for row in aggregates),
        "zero_malignant_species": sum(int(row["malignant"]) == 0 for row in aggregates),
        "matched_exclusions": [row["vincze_species"] for row in reconciliation
                               if row["vincze_species"] and row["compton_species"]
                               and row["status"] != "included"],
        "primary_classes": sorted({row["compton_class"] for row in reconciliation
                                   if row["status"] == "included"}),
        "order_label_disagreements": [
            {"species": row["vincze_species"], "vincze_order": row["vincze_order"],
             "compton_order": row["compton_order"]}
            for row in reconciliation if row["status"] == "included"
            and row["vincze_order"] != row["compton_order"]],
        "endpoint_audit": {
            "audited_compton_species": len(endpoint_audit),
            "neoplasia_prevalence_mismatches": mismatches(
                endpoint_audit, "neoplasia_residual", "neoplasia_tolerance"),
            "malignant_prevalence_mismatches": mismatches(
                endpoint_audit, "malignant_residual", "malignancy_tolerance"),
            "malignancyknown_prevalence_mismatches": mismatches(
                endpoint_audit, "malignancyknown_residual", "malignancy_tolerance"),
            "panel_malignancyknown_prevalence_mismatches": mismatches(
                panel_audit, "malignancyknown_residual", "malignancy_tolerance"),
        },
    }


def derive_outputs(source_texts):
    """Reconcile sources and return only deterministic derived output text.

    This is data inspection, not fitting. CMR, ICM and RDI do not enter
    eligibility or calculations. Malignancy availability cannot remove a row
    from the primary panel: an invalid secondary input stops this frozen
    derivation rather than silently changing primary eligibility.
    """
    for name, metadata in SOURCE_FILES.items():
        if normalized_sha256(source_texts[name]) != metadata["sha256_lf"]:
            raise ValueError("Pinned normalized source hash mismatch: " + name)
    vincze = parse_csv(source_texts[Path("vincze.csv").as_posix()], "Vincze")
    compton = parse_csv(source_texts[Path("compton.csv").as_posix()], "Compton")
    aggregates = []
    ledger = []
    for key in sorted(set(vincze) | set(compton)):
        v = vincze.get(key)
        c = compton.get(key)
        reasons = []
        if v is None:
            reasons.append("absent_from_vincze")
        if c is None:
            reasons.append("absent_from_compton")
        if v is not None:
            for field in PREDICTOR_FIELDS:
                if not v[field]:
                    reasons.append("missing:" + field)
        status = "excluded"
        if not reasons:
            mass = decimal_number(v["Body mass (kg)"], "Body mass (kg)")
            adult_life = decimal_number(v["Adult life expectancy"], "Adult life expectancy")
            mhc = integer(v["MHC I gene family size"], "MHC I gene family size")
            diet = integer(v["Consuming mammal prey"], "Consuming mammal prey")
            if mass <= 0 or adult_life <= 0 or diet not in (0, 1):
                raise ValueError("Invalid complete predictor values: " + v["Species"])
            if not v["Taxonomic order"] or v["Taxonomic order"] in ("NA", "N/A"):
                raise ValueError("Missing Vincze order: " + v["Species"])
            n = integer(c["RecordsWithDenominators"], "RecordsWithDenominators", 1)
            k = integer(c["NeoplasiaWithDenominators"], "NeoplasiaWithDenominators")
            m = integer(c["Malignant"], "Malignant")
            if not 0 <= m <= k <= n:
                raise ValueError("Invalid common-panel event counts: " + v["Species"])
            aggregates.append({
                "species": v["Species"], "order": v["Taxonomic order"],
                "log_body_mass_kg": format(math.log(float(mass)), ".17g"),
                "log_adult_life_expectancy_days": format(math.log(float(adult_life)), ".17g"),
                "mammal_prey": diet, "log1p_mhc_i": format(math.log1p(mhc), ".17g"),
                "records": n, "neoplasia": k, "malignant": m,
            })
            status = "included"
        ledger.append({
            "normalized_species": key,
            "vincze_species": v["Species"] if v else "",
            "compton_species": c["Species"] if c else "",
            "compton_class": c["Class"] if c else "",
            "vincze_order": v["Taxonomic order"] if v else "",
            "compton_order": c["Orders"] if c else "",
            "status": status, "reason": ";".join(reasons),
        })
    aggregates.sort(key=lambda row: row["species"])
    endpoint_audit = [endpoint_audit_row(compton[key]) for key in sorted(compton)]
    summary = summarize(aggregates, ledger, endpoint_audit)
    if summary["endpoint_audit"]["neoplasia_prevalence_mismatches"]:
        raise ValueError("Pinned primary prevalence arithmetic is inconsistent")
    if summary["endpoint_audit"]["malignant_prevalence_mismatches"]:
        raise ValueError("Secondary malignancy numerator is not resolved by arithmetic")
    outputs = {
        Path("species_aggregates.csv").as_posix(): csv_content(aggregates, AGGREGATE_FIELDS),
        Path("reconciliation.csv").as_posix(): csv_content(ledger, RECONCILIATION_FIELDS),
        Path("endpoint_audit.csv").as_posix(): csv_content(endpoint_audit, AUDIT_FIELDS),
    }
    sources = {}
    for name, metadata in SOURCE_FILES.items():
        item = dict(metadata)
        item["url"] = ("https://github.com/" + item["repository"] + "/blob/"
                       + item["commit"] + "/" + item["path"])
        sources[Path(name).as_posix()] = item
    audit = {
        "status": STATUS,
        "access_date": "2026-11-18",
        "access_timezone": "America/New_York",
        "sources": sources,
        "licence_audit": {
            "checked_before_source_data_read": True,
            "vincze_repository_license": None,
            "compton_repository_license": None,
            "pinned_recursive_trees_have_license_or_copying_files": False,
            "compton_readme_has_explicit_reuse_grant": False,
            "repository_metadata_license_field": "null for both on access date",
            "vincze_article": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12644480/",
            "vincze_doi": "https://doi.org/10.1038/s41467-025-65286-x",
            "compton_article": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11726020/",
            "compton_doi": "https://doi.org/10.1158/2159-8290.CD-24-0573",
            "final_article_licences": "CC BY-NC-ND 4.0 for both final articles",
            "terms": "https://creativecommons.org/licenses/by-nc-nd/4.0/",
            "repository_tree_checks": [
                "https://api.github.com/repos/OrsolyaVincze/VinczeEtalNatComm_MHC-cancer/git/trees/2c8a5c5b5fda6fb34c86b0e5a2b6a58dcc85cc3e?recursive=1",
                "https://api.github.com/repos/zacharycompton/cancerAcrossVertebrates/git/trees/4819e0933358c5fc7eab483a0c05161ea8c074b7?recursive=1",
            ],
            "decision": (
                "User-authorized narrow factual per-species aggregates and original code only; "
                "no upstream CSV, dictionary, R code, article text or individual records mirrored. "
                "No claim that article licensing grants adaptation/redistribution rights for "
                "repository assets; project MIT does not relicense third-party material."),
            "retrieval": "GitHub API/plain UTF-8 text; no browser download or download prompt",
        },
        "normalization": "Replace underscores with spaces; trim/collapse whitespace; casefold; exact equality only; reject duplicate keys",
        "order_source": "Vincze Taxonomic order for every included species",
        "adult_lifespan": "Vincze mean remaining life after sexual maturity, in days; not maximum lifespan",
        "primary_outcome": "NeoplasiaWithDenominators / RecordsWithDenominators",
        "secondary_outcome": "Malignant / RecordsWithDenominators",
        "secondary_resolution": (
            "Malignant reproduces printed MalignancyPrevalence at token-specific decimal "
            "precision. The dictionary's MalignancyKnown numerator fails. This is an "
            "arithmetic choice of the recorded-malignancy endpoint, not confirmation of "
            "animal-level provenance or a correction of the upstream dictionary. "
            "PropMalignant is unused and is never interpreted as progression."),
        "rounding_audit": "Half the last printed decimal unit plus 2e-15; literal zero/one endpoints have only the 2e-15 export allowance",
        "prefit_summary": summary,
        "derived_file_hashes_lf": {
            Path(name).as_posix(): normalized_sha256(content)
            for name, content in outputs.items()
        },
    }
    outputs[Path("source_audit.json").as_posix()] = json.dumps(
        audit, indent=2, sort_keys=True, allow_nan=False) + "\n"
    return outputs


def compare_outputs(outputs, destination=BASE):
    """Check exact content except platform-level last-bit differences in logs.

    Archived derived hashes are checked against the actual committed CSV text.
    Regenerated log predictors may differ by a few ulps across C libraries;
    all other derived cells, audit metadata and counts must be exact.
    """
    destination = Path(destination)
    observed = {}
    for name in outputs:
        path = destination / name
        if not path.is_file():
            raise ValueError("Missing derived output: " + Path(name).as_posix())
        observed[name] = normalize_newlines(path.read_bytes().decode("utf-8"))
    aggregate_name = Path("species_aggregates.csv").as_posix()
    logs = ("log_body_mass_kg", "log_adult_life_expectancy_days", "log1p_mhc_i")
    expected_rows = list(csv.DictReader(io.StringIO(outputs[aggregate_name])))
    actual_rows = list(csv.DictReader(io.StringIO(observed[aggregate_name])))
    if len(expected_rows) != len(actual_rows):
        raise ValueError("Aggregate row count differs")
    for expected, actual in zip(expected_rows, actual_rows):
        if expected.keys() != actual.keys():
            raise ValueError("Aggregate schema differs")
        for key in expected:
            if key in logs:
                if not math.isclose(float(expected[key]), float(actual[key]),
                                    rel_tol=2e-15, abs_tol=2e-15):
                    raise ValueError("Regenerated predictor differs: " + key)
            elif expected[key] != actual[key]:
                raise ValueError("Regenerated aggregate differs: " + key)
    audit_name = Path("source_audit.json").as_posix()
    expected_audit = json.loads(outputs[audit_name])
    actual_audit = json.loads(observed[audit_name])
    for name, digest in actual_audit["derived_file_hashes_lf"].items():
        if normalized_sha256(observed[name]) != digest:
            raise ValueError("Archived derived hash differs: " + Path(name).as_posix())
    expected_audit.pop("derived_file_hashes_lf")
    actual_audit.pop("derived_file_hashes_lf")
    if actual_audit != expected_audit:
        raise ValueError("Regenerated source audit differs")
    for name in outputs:
        if name not in (aggregate_name, audit_name) and outputs[name] != observed[name]:
            raise ValueError("Regenerated content differs: " + Path(name).as_posix())


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sources = parser.add_mutually_exclusive_group(required=True)
    sources.add_argument("--fetch", action="store_true", help="Opt in to pinned GitHub raw UTF-8 reads")
    sources.add_argument("--source-dir", type=Path, help="External directory containing vincze.csv, compton.csv and dictionary.txt")
    parser.add_argument("--check", action="store_true", help="Compare derived outputs without writing")
    args = parser.parse_args(argv)
    outputs = derive_outputs(load_sources(args.source_dir, args.fetch))
    if args.check:
        compare_outputs(outputs)
    else:
        for name, content in outputs.items():
            with (BASE / name).open("w", encoding="utf-8", newline="\n") as handle:
                handle.write(content)
    print(("Checked " if args.check else "Wrote ") + str(len(outputs))
          + " derived files; no upstream source assets were written.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
