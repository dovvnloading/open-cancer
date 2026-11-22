"""Nonclinical, unreviewed research. Replay one frozen HS-signature test."""
import argparse
import csv
import hashlib
import json
import math
import random
import re
from collections import Counter, defaultdict
from itertools import groupby
from pathlib import Path

HERE = Path(__file__).resolve().parent
GENES = {"EXT1": 2131, "EXT2": 2132, "SLC35B2": 347734, "HS2ST1": 9653}
FIELDS = ["cell_line", "rrid", "model_id", "lineage", "ic50_nm", *GENES, "match_issue", "research_status"]
STATUS = "Nonclinical, unreviewed research"
CAP = 2500
PERMUTATIONS = 9999
SOURCE_HASHES = {
    "s2.xlsx": "4c6f547d146627fe11c34c15dc291ba0e3f11397cea9971da51f22121fce2650",
    "s4.xlsx": "45c9e40ef9243d9f19a8e073da04958e78fe744ddec976d05f2a013c7c33bf53",
    "Model.csv": "b7a0c1385e6cef30132b56aff61f1261d11e3f490490b355c430d32ee0dbdcfa",
    "CRISPRGeneEffect.csv": "3d8f3ec6dbf2db7ff834b79b508622ec0b226f3518003fe96ecf5a4fcf167e3b",
}


def digest(path):
    """SHA-256 of UTF-8 text with LF endings, or the original XLSX bytes."""
    if path.suffix == ".xlsx":
        with path.open("rb") as raw:
            return hashlib.file_digest(raw, "sha256").hexdigest()
    h = hashlib.sha256()
    with path.open(encoding="utf-8", newline=None) as stream:
        while chunk := stream.read(1024 * 1024):
            h.update(chunk.encode("utf-8"))
    return h.hexdigest()


def name_key(value):
    return re.sub("[^A-Z0-9]", "", str(value).upper())


def rrid(value):
    return str(value or "").strip().removeprefix("RRID:").strip()


def workbook_table(path):
    # Optional dependency only for rebuilding the extract from the original XLSX files.
    from openpyxl import load_workbook
    book = load_workbook(path, read_only=True, data_only=True)
    try:
        rows = iter(book.active.values)
        next(rows)  # table title
        header = next(rows)
        return header, [row for row in rows if row[0] is not None]
    finally:
        book.close()


def extract(directory, output):
    for name, expected in SOURCE_HASHES.items():
        if digest(directory / name) != expected:
            raise ValueError(f"Source checksum mismatch: {name}")
    header, cell_rows = workbook_table(directory / "s2.xlsx")
    rrid_col = header.index("RRID")
    by_name = defaultdict(list)
    for row in cell_rows:
        by_name[name_key(row[0])].append(rrid(row[rrid_col]))

    models = defaultdict(list)
    with (directory / "Model.csv").open(encoding="utf-8-sig", newline="") as f:
        for model in csv.DictReader(f):
            if cell_id := rrid(model["RRID"]):
                models[cell_id].append((model["ModelID"], model["OncotreeLineage"]))

    effects = {}
    with (directory / "CRISPRGeneEffect.csv").open(encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f)
        columns = next(reader)
        indices = [columns.index(f"{gene} ({entrez})") for gene, entrez in GENES.items()]
        for row in reader:
            if row[0] in effects:
                raise ValueError(f"Duplicate gene-effect model: {row[0]}")
            effects[row[0]] = dict(zip(GENES, (row[i] for i in indices)))

    header, drug_rows = workbook_table(directory / "s4.xlsx")
    drug_col = header.index("selpercatinib")
    occurrences = Counter(name_key(row[0]) for row in drug_rows)
    with output.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        for row in drug_rows:
            name = name_key(row[0])
            record = {"cell_line": row[0], "ic50_nm": row[drug_col], "research_status": STATUS}
            ids = by_name[name]
            if occurrences[name] != 1 or len(ids) != 1 or not ids[0]:
                record["match_issue"] = "missing or ambiguous S2 identifier"
            else:
                record["rrid"] = ids[0]
                matches = models[ids[0]]
                if len(matches) == 1:
                    model_id, lineage = matches[0]
                    record.update(model_id=model_id, lineage=lineage)
                    record.update(effects.get(model_id, {}))
                else:
                    record["match_issue"] = "ambiguous DepMap RRID" if matches else "unmatched DepMap RRID"
            writer.writerow(record)


def finite(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def capped_ic50(value):
    raw = str(value).strip()
    bound = raw.startswith(">")
    number = finite(raw.removeprefix(">"))
    if number is None or number <= 0 or (bound and number < CAP):
        return None
    return min(CAP, number)


def eligibility(records):
    model_counts = Counter(r["model_id"] for r in records if r["model_id"])
    kept = defaultdict(list)
    reasons = {}

    def exclude(row, reason):
        reasons[row["cell_line"]] = reason

    for row in records:
        name = name_key(row["cell_line"])
        ic50 = capped_ic50(row["ic50_nm"])
        deps = [finite(row[gene]) for gene in GENES]
        if name in {"LC2AD", "TT"}:
            exclude(row, "source-identified RET alteration")
        elif row["match_issue"]:
            exclude(row, row["match_issue"])
        elif not row["model_id"] or model_counts[row["model_id"]] != 1:
            exclude(row, "missing or reused model ID")
        elif not row["lineage"] or None in deps:
            exclude(row, "missing lineage or knockout measurement")
        elif ic50 is None:
            exclude(row, "unusable IC50")
        else:
            kept[row["lineage"]].append((row["model_id"], deps, ic50, row))
    groups = {}
    for lineage, models in sorted(kept.items()):
        if len(models) == 1:
            exclude(models[0][3], "singleton lineage")
        else:
            groups[lineage] = sorted(models, key=lambda m: m[0])
    return groups, reasons


def centered_ranks(values):
    size = len(values)
    result = [0.0] * size
    start = 0
    for _, tied in groupby(sorted(range(size), key=values.__getitem__), key=values.__getitem__):
        indices = list(tied)
        end = start + len(indices)
        rank = (start + (len(indices) + 1) / 2 - (size + 1) / 2) / size
        for index in indices:
            result[index] = rank
        start = end
    return result


def dot(x, y):
    return math.fsum(a * b for a, b in zip(x, y, strict=True))


def analyze(records):
    groups, reasons = eligibility(records)
    xs, ys = [], []
    uncapped = informative = 0
    for members in groups.values():
        effects = [entry[1] for entry in members]
        columns = [centered_ranks(col) for col in zip(*effects, strict=True)]
        score = [math.fsum(row) / len(GENES) for row in zip(*columns, strict=True)]
        response = [entry[2] for entry in members]
        uncapped += sum(v < CAP for v in response)
        ranks = centered_ranks(response)
        informative += len(set(score)) > 1 and len(set(ranks)) > 1
        xs.append(score)
        ys.append(ranks)

    n = sum(map(len, xs))
    den = math.sqrt(math.fsum(dot(x, x) for x in xs) * math.fsum(dot(y, y) for y in ys))
    result = {
        "scope": STATUS,
        "source_rows": len(records),
        "n": n,
        "below_cap": uncapped,
        "informative_lineages": informative,
        "lineage_counts": {k: len(v) for k, v in groups.items()},
        "excluded": dict(sorted(reasons.items())),
        "exclusion_counts": dict(Counter(reasons.values())),
    }
    failed_qc = []
    if n < 40:
        failed_qc.append("fewer than 40 models")
    if uncapped < 8:
        failed_qc.append("fewer than eight responses below the cap")
    if informative < 2:
        failed_qc.append("fewer than two informative lineages")
    if den == 0:
        failed_qc.append("zero rank-vector variance")
    if failed_qc:
        return dict(result, verdict="INCONCLUSIVE", reasons=failed_qc)

    paired = list(zip(xs, ys, strict=True))
    observed = math.fsum(dot(x, y) for x, y in paired)
    rng = random.Random(96)
    hits = 0
    for _ in range(PERMUTATIONS):
        products = []
        for score, ranks in paired:
            shuffled = ranks.copy()
            rng.shuffle(shuffled)
            products.append(dot(score, shuffled))
        hits += math.fsum(products) >= observed - 1e-12
    r = observed / den
    p = (hits + 1) / (PERMUTATIONS + 1)
    return dict(result, r=r, p=p, permutations=PERMUTATIONS, at_least_observed=hits,
                verdict="PASS" if r >= 0.20 and p <= 0.05 else "FAIL")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, help="rebuild extract from the four pinned files")
    parser.add_argument("--out", type=Path, default=HERE / "derived", help="output directory")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    panel = HERE / "derived" / "panel.csv"
    if args.source_dir:
        panel = args.out / "panel.csv"
        extract(args.source_dir, panel)
    with panel.open(encoding="utf-8", newline="") as f:
        summary = analyze(list(csv.DictReader(f)))
    summary["sha256_lf"] = {Path("derived/panel.csv").as_posix(): digest(panel),
                            Path("PROTOCOL.md").as_posix(): digest(HERE / "PROTOCOL.md")}
    with (args.out / "summary.json").open("w", encoding="utf-8", newline="\n") as f:
        json.dump(summary, f, indent=2, sort_keys=True, allow_nan=False)
        f.write("\n")


if __name__ == "__main__":
    main()
