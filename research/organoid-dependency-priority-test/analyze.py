"""Run the frozen basal EGFR test on local public Cell Model Passports files.

Nonclinical, unreviewed computational research. No treatment advice.
"""
import argparse
import csv
import hashlib
import io
import json
import math
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path
from statistics import NormalDist
from zipfile import BadZipFile, ZipFile

MATRIX = Path("Organoids_fitness_scores_20260723.zip") / "Organoids_fitness_scores_binary_matrix_20260723.tsv"
MUTATIONS = Path("mutations_summary_20260724.zip") / "mutations_summary_20260724.csv"
MODELS = Path("model_list_20261013.csv")
COVERAGE = Path("model_dataset_availability_20260727.csv")
WGS = "Mutation Sanger Organoid (WGS)"
GROUPS = ("BRAF_V600E", "BRAF_reference")
MISSING = {"", "na", "nan", "n/a", "null"}


class InputError(ValueError):
    """The frozen test cannot interpret an input unambiguously."""


@contextmanager
def open_text(data_dir, path):
    """Expose either a standalone file or one named ZIP member as normalized text."""
    if len(path.parts) == 1:
        with (data_dir / path).open(encoding="utf-8", newline=None) as stream:
            yield stream
    else:
        archive_name, member = path.parts
        with ZipFile(data_dir / archive_name) as archive:
            if archive.namelist().count(member) != 1:
                raise InputError(f"Missing or ambiguous input member in {archive_name}")
            with archive.open(member) as raw, io.TextIOWrapper(raw, encoding="utf-8", newline=None) as stream:
                yield stream


def normalized_hash(data_dir, path):
    digest = hashlib.sha256()
    with open_text(data_dir, path) as stream:
        for chunk in iter(lambda: stream.read(65536), ""):
            digest.update(chunk.encode("utf-8"))
    return digest.hexdigest()


def checked_rows(stream, columns):
    rows = csv.DictReader(stream)
    header = rows.fieldnames or []
    if not set(columns).issubset(header) or len(set(header)) != len(header):
        raise InputError("Missing or duplicate CSV fields")
    for row in rows:
        if None in row or any(row[key] is None for key in columns):
            raise InputError("Malformed CSV record")
        yield row


def model_index(stream, columns):
    index = {}
    for row in checked_rows(stream, ("model_id", *columns)):
        model = row["model_id"]
        if not model or model in index:
            raise InputError("Missing or duplicate annotation key")
        index[model] = {key: row[key] for key in columns}
    return index


def classify_genotype(variants):
    """Classify one WGS-covered model's reported Sanger variants."""
    braf = set()
    ras_driver = False
    for gene, protein, driver in variants:
        if driver not in {"t", "f"}:
            raise InputError("Unrecognized cancer driver flag")
        if gene in {"KRAS", "NRAS", "HRAS"} and driver == "t":
            ras_driver = True
        if gene == "BRAF":
            braf.add((protein, driver))
    if ras_driver:
        return "ras_driver"
    if not braf:
        return "BRAF_reference"
    if braf == {("p.V600E", "t")}:
        return "BRAF_V600E"
    return "other_or_ambiguous_braf"


def load_egfr(stream):
    reader = csv.reader(stream, delimiter="\t")
    expected = ("model_name", "model_id", "source", "qc_pass", "gene_id")
    headers = [next(reader, []) for _ in expected]
    if any(not row or row[0] != label for row, label in zip(headers, expected)):
        raise InputError("Unrecognized matrix header")
    width = len(headers[0])
    if width <= 3 or any(len(row) != width for row in headers) or headers[-1][:3] != ["gene_id", "symbol", "ensembl_gene_id"]:
        raise InputError("Inconsistent matrix header")
    egfr = None
    for row in reader:
        if len(row) < 2:
            raise InputError("Malformed gene record")
        if row[1] == "EGFR":
            if egfr is not None or len(row) != width:
                raise InputError("Ambiguous EGFR row or column count")
            egfr = row[3:]
    if egfr is None:
        raise InputError("Missing EGFR row")
    return headers[1][3:], headers[2][3:], headers[3][3:], egfr


def fisher_two_sided(a, n_a, b, n_b):
    """Probability ordering for two independent binomials conditioned on margins."""
    if not (0 <= a <= n_a and 0 <= b <= n_b and n_a and n_b):
        raise ValueError("Invalid binomial counts")
    events = a + b
    observed = math.comb(n_a, a) * math.comb(n_b, b)
    numerator = 0
    for x in range(max(0, events - n_b), min(n_a, events) + 1):
        weight = math.comb(n_a, x) * math.comb(n_b, events - x)
        if weight <= observed:
            numerator += weight
    return numerator / math.comb(n_a + n_b, events)


def wilson(events, total):
    if not (0 <= events <= total and total):
        raise ValueError("Invalid binomial counts")
    z = NormalDist().inv_cdf(0.975)
    p = events / total
    denom = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denom
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denom
    return max(0.0, center - half), min(1.0, center + half)


def newcombe_difference(a, n_a, b, n_b):
    """Newcombe method 10 for two independent proportions, without continuity correction."""
    lo_a, hi_a = wilson(a, n_a)
    lo_b, hi_b = wilson(b, n_b)
    p, q = a / n_a, b / n_b
    return (p - q - math.hypot(p - lo_a, hi_b - q),
            p - q + math.hypot(hi_a - p, q - lo_b))


def collect_groups(data_dir):
    with open_text(data_dir, MODELS) as stream:
        models = model_index(stream, ("model_type", "cancer_type", "patient_id"))
    with open_text(data_dir, COVERAGE) as stream:
        availability = model_index(stream, (WGS,))
    variants = defaultdict(list)
    with open_text(data_dir, MUTATIONS) as stream:
        for row in checked_rows(stream, ("model_id", "source", "gene_symbol", "protein_mutation", "cancer_driver")):
            if row["source"] == "Sanger" and row["gene_symbol"] in {"KRAS", "NRAS", "HRAS", "BRAF"}:
                variants[row["model_id"]].append((row["gene_symbol"], row["protein_mutation"], row["cancer_driver"]))
    with open_text(data_dir, MATRIX) as stream:
        model_ids, sources, quality, calls = load_egfr(stream)

    attrition = dict.fromkeys((
        "input_screen_columns", "non_sanger_or_failed_qc", "missing_model_metadata",
        "outside_crc_organoids", "missing_donor", "no_confirmed_organoid_wgs",
        "ras_driver", "other_or_ambiguous_braf", "extra_cultures_from_same_donor",
        "missing_selected_egfr_call", "evaluable_donors",
    ), 0)
    attrition["input_screen_columns"] = len(model_ids)
    eligible = {}
    for column, (model, source, qc) in enumerate(zip(model_ids, sources, quality)):
        if source == "Sanger" and qc not in {"TRUE", "FALSE"}:
            raise InputError("Unrecognized matrix QC coding")
        if source != "Sanger" or qc == "FALSE":
            attrition["non_sanger_or_failed_qc"] += 1
            continue
        if model not in models:
            attrition["missing_model_metadata"] += 1
            continue
        row = models[model]
        if (row["model_type"], row["cancer_type"]) != ("Organoid", "Colorectal Carcinoma"):
            attrition["outside_crc_organoids"] += 1
            continue
        donor = row["patient_id"].strip()
        if donor.lower() in MISSING:
            attrition["missing_donor"] += 1
            continue
        coverage = availability.get(model, {}).get(WGS)
        if coverage not in {None, "False", "True", ""}:
            raise InputError("Unrecognized mutation dataset availability flag")
        if coverage != "True":
            attrition["no_confirmed_organoid_wgs"] += 1
            continue
        group = classify_genotype(variants[model])
        if group not in GROUPS:
            attrition[group] += 1
            continue
        if model in eligible:
            raise InputError("Duplicate eligible assay columns for one model")
        eligible[model] = (donor, group, column)

    seen_donors = set()
    groups = {group: {"essential": 0, "tested": 0} for group in GROUPS}
    for model in sorted(eligible):
        donor, group, column = eligible[model]
        if donor in seen_donors:
            attrition["extra_cultures_from_same_donor"] += 1
            continue
        seen_donors.add(donor)
        call = calls[column].strip()
        if call.lower() in MISSING:
            attrition["missing_selected_egfr_call"] += 1
            continue
        if call not in {"0", "1"}:
            raise InputError("Unrecognized EGFR binary call")
        groups[group]["essential"] += int(call)
        groups[group]["tested"] += 1
        attrition["evaluable_donors"] += 1
    return groups, attrition


def analyze(data_dir: Path) -> dict:
    result = {
        "research_status": "nonclinical, unreviewed computational research",
        "access_date": "2026-11-22",
        "source_base_url": "https://cog.sanger.ac.uk/cmp/download/",
        "inputs_sha256": {},
        "attrition": None,
        "groups": None,
        "statistics": None,
        "verdict": "feasibility_failure",
    }
    try:
        for path in (MATRIX, MODELS, COVERAGE, MUTATIONS):
            result["inputs_sha256"][path.as_posix()] = normalized_hash(data_dir, path)
        groups, attrition = collect_groups(data_dir)
    except (InputError, OSError, BadZipFile, csv.Error, UnicodeError) as error:
        result["failure_reason"] = str(error) if isinstance(error, InputError) else "Unable to read a required input"
        return result
    result["groups"] = groups
    result["attrition"] = attrition
    if min(group["tested"] for group in groups.values()) < 5:
        result["failure_reason"] = "Fewer than five evaluable independent donors in a group"
        return result
    a, n_a = groups["BRAF_V600E"]["essential"], groups["BRAF_V600E"]["tested"]
    b, n_b = groups["BRAF_reference"]["essential"], groups["BRAF_reference"]["tested"]
    delta = a / n_a - b / n_b
    p = fisher_two_sided(a, n_a, b, n_b)
    result["statistics"] = {
        "risk_difference": delta,
        "risk_difference_ci95": list(newcombe_difference(a, n_a, b, n_b)),
        "fisher_p_two_sided": p,
    }
    result["verdict"] = "inconclusive"
    if p <= 0.05 and delta > 0:
        result["verdict"] = "supports_basal_prediction"
    elif p <= 0.05 and delta < 0:
        result["verdict"] = "contradicts_basal_prediction"
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True, help="Path for the aggregate result JSON")
    args = parser.parse_args()
    result = analyze(args.data_dir)
    args.out.parent.mkdir(exist_ok=True, parents=True)
    with args.out.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return 1 if result["verdict"] == "feasibility_failure" else 0


if __name__ == "__main__":
    raise SystemExit(main())
