"""Nonclinical, unreviewed DepMap analysis. Run without arguments for offline replay."""
import argparse
import csv
import hashlib
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

HERE = Path(__file__).resolve().parent
SCREEN = HERE.parent / "expression-loss-synthetic-lethality/derived/depmap24q4_extract_selected_pairs.csv"
EXTRACT = HERE / "derived/depmap24q4_x_cn.csv"
REFERENCE = ("MID1", "KDM5C", "HUWE1", "RPS4X", "G6PD")
GENES = dict(PUDP=8226, MID1=4281, KDM5C=8242, HUWE1=10075, RPS4X=6191,
             G6PD=2539, STS=412, PNPLA4=8228)
DRAW_COUNT = 19_999
MARGIN = 0.20
ELIGIBILITY = {"representative", "missing_assay_or_lineage", "missing_sex_or_donor",
               "missing_or_invalid_cn", "repeat_donor"}


def read_models(path):
    with path.open(encoding="utf-8", newline="") as stream:
        models = {}
        for row in csv.DictReader(stream):
            key = row["ModelID"]
            if not key or key in models:
                raise ValueError(f"Missing or duplicate ModelID in {path.name}")
            models[key] = row
    return models


def number(value):
    if value in ("", "NA", "NaN", None):
        return None
    value = float(value)
    return value if math.isfinite(value) else None


def project_cn(path, screen_ids):
    """Stream the full matrix once. Preserve the source values for eight named genes."""
    projected = {}
    seen = set()
    with path.open(encoding="utf-8", newline="") as stream:
        rows = csv.reader(stream)
        header = next(rows)
        names = [f"{gene} ({entrez})" for gene, entrez in GENES.items()]
        if header[0] not in ("", "ModelID") or any(header.count(name) != 1 for name in names):
            raise ValueError("Unexpected copy-number schema")
        indices = [header.index(name) for name in names]
        for row in rows:
            key = row[0]
            if not key or key in seen or len(row) != len(header):
                raise ValueError("Duplicate/missing copy-number key or malformed row")
            seen.add(key)
            if key in screen_ids:
                projected[key] = dict(zip(GENES, (row[i] for i in indices)))
    return projected


def x_reference(row):
    values = [number(row[gene]) for gene in REFERENCE]
    if any(x is None or x < 0 for x in values):
        return None
    center = median(values)
    return center if center > 0 else None


def prepare(cn_path, model_path, output):
    for path, md5 in ((cn_path, "4696565fa74dc805983bc278b4a41802"),
                      (model_path, "675210d17675f3517b0ce39a3c274f16")):
        with path.open("rb") as source:
            if hashlib.file_digest(source, "md5").hexdigest() != md5:
                raise ValueError(f"Figshare source MD5 mismatch for {path.name}")
    screen, metadata = read_models(SCREEN), read_models(model_path)
    projected = project_cn(cn_path, screen.keys())
    seen_donors = set()
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, ("ModelID", "sex", "eligibility", *GENES), lineterminator="\n")
        writer.writeheader()
        for model in sorted(screen):
            assay = screen[model]
            meta = metadata.get(model, {})
            sex = meta.get("Sex", "")
            donor = meta.get("PatientID", "").strip()
            genes = projected.get(model, dict.fromkeys(GENES, ""))
            pudp = number(genes["PUDP"])
            if number(assay["ITPA_gene_effect"]) is None or not assay["OncotreeLineage"].strip():
                status = "missing_assay_or_lineage"
            elif sex not in ("Female", "Male") or not donor:
                status = "missing_sex_or_donor"
            elif pudp is None or pudp < 0 or x_reference(genes) is None:
                status = "missing_or_invalid_cn"
            elif donor in seen_donors:
                status = "repeat_donor"
            else:
                status = "representative"
                seen_donors.add(donor)
            writer.writerow(dict(ModelID=model, sex=sex, eligibility=status, **genes))


def exposure(ratio):
    if not math.isfinite(ratio) or ratio < 0:
        raise ValueError("Invalid relative DNA dosage")
    if ratio <= 0.20:
        return "loss"
    if 0.75 <= ratio <= 1.25:
        return "retained"
    return "intermediate_or_gain"


def decision(delta, upper_p, lower_margin_p):
    if delta >= MARGIN and upper_p <= 0.05:
        return "SUPPORTED"
    if delta < MARGIN and lower_margin_p <= 0.05:
        return "FALSIFIED"
    return "INCONCLUSIVE"


def permutation(strata, observed):
    """Invert the same matched contrast under zero and the additive margin."""
    rng = random.Random(9401)
    prepared = []
    n_loss = 0
    for rows in strata:
        scores = [value for _, value in rows]
        labels = [group == "loss" for group, _ in rows]
        k = sum(labels)
        prepared.append((scores, labels, k, len(rows) - k, sum(scores)))
        n_loss += k

    high = low = 0
    for _ in range(DRAW_COUNT):
        delta = shifted = 0.0
        for scores, labels, k, m, total in prepared:
            chosen = rng.sample(range(k + m), k)
            case_sum = sum(scores[i] for i in chosen)
            overlap = sum(labels[i] for i in chosen)
            # Multiply each stratum's mean contrast by its number of loss donors.
            contribution = k * (total - case_sum) / m - case_sum
            delta += contribution
            shifted += contribution + MARGIN * (k * (k - overlap) / m - overlap)
        high += delta / n_loss >= observed - 1e-12
        low += shifted / n_loss <= observed - MARGIN + 1e-12
    return (high + 1) / (DRAW_COUNT + 1), (low + 1) / (DRAW_COUNT + 1)


def analyze(extract=EXTRACT, screen_file=SCREEN):
    screen, cn = read_models(screen_file), read_models(extract)
    if screen.keys() != cn.keys():
        raise ValueError("Extract does not cover the same ModelIDs as the frozen screen")
    eligible, groups = Counter(), Counter()
    cells = defaultdict(list)
    for key in sorted(screen):
        row = cn[key]
        if row["eligibility"] not in ELIGIBILITY:
            raise ValueError(f"Unknown eligibility for {key}")
        eligible[row["eligibility"]] += 1
        if row["eligibility"] != "representative":
            continue
        denominator = x_reference(row)
        pudp = number(row["PUDP"])
        assay = screen[key]
        score = number(assay["ITPA_gene_effect"])
        if (denominator is None or pudp is None or row["sex"] not in ("Female", "Male")
                or not assay["OncotreeLineage"].strip() or score is None):
            raise ValueError(f"Invalid representative inputs for {key}")
        group = exposure(pudp / denominator)
        groups[group] += 1
        if group != "intermediate_or_gain":
            cells[assay["OncotreeLineage"], row["sex"]].append((key, group, score))

    diagnostics = defaultdict(lambda: defaultdict(lambda: {"positive": 0, "denominator": 0}))
    matched = []
    details = []
    weighted = 0.0
    n_loss = n_retained = 0
    for (lineage, sex), models in sorted(cells.items()):
        losses = [score for _, group, score in models if group == "loss"]
        retained = [score for _, group, score in models if group == "retained"]
        if not losses or len(retained) < 2:
            continue
        contrast = sum(retained) / len(retained) - sum(losses) / len(losses)
        weighted += len(losses) * contrast
        n_loss += len(losses)
        n_retained += len(retained)
        details.append(dict(lineage=lineage, sex=sex, loss=len(losses), retained=len(retained),
                            retained_minus_loss=contrast))
        matched.append([(group, score) for _, group, score in models])
        for key, group, _ in models:
            values = {"PUDP_RNA_off": number(screen[key]["PUDP_log2tpm1"])}
            for gene in ("STS", "PNPLA4"):
                value = number(cn[key][gene])
                values[gene] = None if value is None or value < 0 else value / x_reference(cn[key])
            for name, value in values.items():
                counts = diagnostics[group][name]
                if value is not None:
                    counts["denominator"] += 1
                    counts["positive"] += int(value < 1 if name == "PUDP_RNA_off" else value <= 0.20)
    delta = weighted / n_loss if n_loss else None
    report = dict(scope="Nonclinical, unreviewed computational research", n_screen=len(screen),
                  eligibility=dict(eligible), exposure=dict(groups), strata=details,
                  matched_loss=n_loss, matched_retained=n_retained, retained_minus_loss=delta,
                  diagnostics=diagnostics, draws=0, p_zero_upper=None, p_margin_lower=None)
    if n_loss < 8 or len(details) < 2:
        report.update(verdict="INCONCLUSIVE", reason="insufficient exposure")
    else:
        upper, lower = permutation(matched, delta)
        report.update(verdict=decision(delta, upper, lower), reason="frozen contrast",
                      draws=DRAW_COUNT, p_zero_upper=upper, p_margin_lower=lower)
    return report


def normalized_sha256(path):
    return hashlib.sha256(path.read_text(encoding="utf-8").encode("utf-8")).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command")
    project = commands.add_parser("prepare", help="rebuild the small extract from original Figshare files")
    project.add_argument("copy_number", type=Path)
    project.add_argument("model_metadata", type=Path)
    project.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args.copy_number, args.model_metadata, args.output)
    else:
        report = analyze()
        repo = HERE.parents[1]
        report["sha256_lf_utf8"] = {p.relative_to(repo).as_posix(): normalized_sha256(p)
                                    for p in (SCREEN, EXTRACT, HERE / "PROTOCOL.md")}
        print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
