#!/usr/bin/env python3
"""Rebuild permitted numerical aggregates from pinned public text sources.

Nonclinical, unreviewed research; no medical or veterinary advice.

The upstream repository has no declared licence. This original script does not
mirror upstream files or copy upstream analysis code. It reads the two pinned
UTF-8 text objects into memory and writes only the derived analysis tables and
provenance audit. Individual Species360 records are neither requested nor used.

Network access is opt-in with --fetch. Alternatively supply --source-dir outside
the repository. --check compares canonical contents without writing anything.
"""

import argparse
import csv
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN, localcontext
import hashlib
import io
import json
import math
from pathlib import Path
import re
import statistics
import urllib.request


SOURCE_REPOSITORY = "OrsolyaVincze/VinczeEtalNatComm_MHC-cancer"
SOURCE_COMMIT = "2c8a5c5b5fda6fb34c86b0e5a2b6a58dcc85cc3e"
PAPER_DOI = "https://doi.org/10.1038/s41467-025-65286-x"
SOURCE_FILES = {
    Path("SourceData_FINAL.csv").as_posix(): {
        "git_blob_sha1": "8a937a2c8d4f1749ae33c1b12d2ed52e895852f9",
        "sha256_lf": "b9d6dedc9447e9f99b275cbfbdb54e3b0e318d08c12f0384572a65101d057744",
    },
    Path("phylogeny.tre").as_posix(): {
        "git_blob_sha1": "64970f63fc167b7fdbadb5ebabbe3833dfb07c82",
        "sha256_lf": "23dab5f42ff4c6e736486fd54d115fff67b4b82d12a57a755f6577607b596c6a",
    },
}
EXPORT_TOLERANCE = Decimal("2e-15")
REQUIRED_FIELDS = (
    "RDI", "CMR", "Adult life expectancy", "Body mass (kg)",
    "Consuming mammal prey", "MHC I gene family size",
)
AGGREGATE_FIELDS = (
    "species", "order", "known_cause_deaths", "cmr_decimal",
    "nominal_cancer_deaths", "count_lower", "count_upper", "log_body_mass_kg",
    "log_adult_life_expectancy_days", "mammal_prey", "log1p_mhc_i",
    "phylo_fair_proportion", "phylo_block",
)


def normalize_newlines(text):
    """Normalize line endings only; retain whitespace and final-newline state."""
    return text.replace("\r\n", "\n").replace("\r", "\n")


def normalized_sha256(text):
    return hashlib.sha256(normalize_newlines(text).encode("utf-8")).hexdigest()


def git_blob_sha1(raw):
    header = ("blob " + str(len(raw)) + "\0").encode("ascii")
    return hashlib.sha1(header + raw).hexdigest()


def load_sources(source_dir=None, fetch=False):
    """Read and verify sources without saving raw material.

    Git SHA-1 is verified on the source bytes, or on their LF equivalent when a
    local text-copy operation changed only newlines. SHA-256 always covers LF
    normalized UTF-8 content. Remote reads additionally require exact Git bytes.
    """
    if bool(source_dir is not None) == bool(fetch):
        raise ValueError("Choose exactly one of source_dir or fetch")
    result = {}
    for name, expected in SOURCE_FILES.items():
        if fetch:
            url = ("https://raw.githubusercontent.com/" + SOURCE_REPOSITORY
                   + "/" + SOURCE_COMMIT + "/" + name)
            request = urllib.request.Request(
                url, headers={"User-Agent": "open-cancer-mhc-aggregate-audit/1"})
            with urllib.request.urlopen(request, timeout=30) as response:
                raw = response.read(2_000_001)
            if len(raw) > 2_000_000:
                raise ValueError("Unexpectedly large source: " + name)
        else:
            raw = (Path(source_dir) / name).read_bytes()
        text = raw.decode("utf-8")
        canonical = normalize_newlines(text)
        blob_matches = git_blob_sha1(raw) == expected["git_blob_sha1"]
        local_newline_match = (
            not fetch and git_blob_sha1(canonical.encode("utf-8"))
            == expected["git_blob_sha1"])
        if not (blob_matches or local_newline_match):
            raise ValueError("Pinned Git blob mismatch: " + name)
        if normalized_sha256(canonical) != expected["sha256_lf"]:
            raise ValueError("Pinned normalized SHA-256 mismatch: " + name)
        result[name] = canonical
    return result


def count_interval(cmr_token, denominator):
    """Nominal count and all integers compatible with stated rounding precision.

    A non-endpoint token represents a probability within half its last decimal
    unit, enlarged by the frozen probability-scale export tolerance. Both
    interval boundaries are included conservatively. Values zero and one are
    exact endpoints. This is an uncertainty set, not proof of original counts.
    """
    with localcontext() as context:
        context.prec = 60
        probability = Decimal(cmr_token)
        if not probability.is_finite() or not 0 <= probability <= 1:
            raise ValueError("Invalid CMR probability")
        if denominator < 1 or int(denominator) != denominator:
            raise ValueError("RDI must be a positive integer")
        denominator = int(denominator)
        nominal = int((probability * denominator).to_integral_value(
            rounding=ROUND_HALF_EVEN))
        if probability in (Decimal(0), Decimal(1)):
            return nominal, nominal, nominal
        last_unit = Decimal(1).scaleb(probability.as_tuple().exponent)
        half_width = last_unit / 2 + EXPORT_TOLERANCE
        lower_probability = max(Decimal(0), probability - half_width)
        upper_probability = min(Decimal(1), probability + half_width)
        lower = int((lower_probability * denominator).to_integral_value(
            rounding=ROUND_CEILING))
        upper = int((upper_probability * denominator).to_integral_value(
            rounding=ROUND_FLOOR))
        if not lower <= nominal <= upper:
            raise ValueError("CMR and RDI have no compatible nominal integer")
        return nominal, lower, upper


class Node:
    """Minimal original Newick representation; length is the incoming edge."""

    def __init__(self, label, length, children=()):
        self.label = label
        self.length = length
        self.children = tuple(children)


def parse_nexus(text):
    """Parse this pinned Nexus tree, including TRANSLATE and internal labels."""
    translation_match = re.search(r"\bTRANSLATE\s+(.*?);", text, re.I | re.S)
    tree_match = re.search(r"\bTREE\s+[^=]*=\s*(.*?;)", text, re.I | re.S)
    ntax_match = re.search(r"\bNTAX\s*=\s*(\d+)", text, re.I)
    if not translation_match or not tree_match or not ntax_match:
        raise ValueError("Expected a Nexus TRANSLATE table and one tree")
    translation = {}
    for entry in translation_match.group(1).split(","):
        key, name = entry.strip().split(None, 1)
        name = name.strip().strip("'")
        if key in translation or name in translation.values():
            raise ValueError("Duplicate translated tip")
        translation[key] = name
    if len(translation) != int(ntax_match.group(1)):
        raise ValueError("Nexus NTAX and translation count disagree")
    newick = re.sub(r"\[[^\]]*\]", "", tree_match.group(1))
    tokens = re.findall(r"[(),:;]|[^\s(),:;]+", newick)
    position = 0

    def parse_node():
        nonlocal position
        if position >= len(tokens):
            raise ValueError("Truncated Newick tree")
        children = []
        if tokens[position] == "(":
            position += 1
            children.append(parse_node())
            while tokens[position] == ",":
                position += 1
                children.append(parse_node())
            if tokens[position] != ")":
                raise ValueError("Expected Newick closing parenthesis")
            position += 1
            label = ""
            if tokens[position] not in (":", ",", ")", ";"):
                label = tokens[position]
                position += 1
        else:
            key = tokens[position]
            if key not in translation:
                raise ValueError("Untranslated tip: " + key)
            label = translation[key]
            position += 1
        length = 0.0
        if tokens[position] == ":":
            position += 1
            length = float(tokens[position])
            position += 1
        if not math.isfinite(length) or length < 0:
            raise ValueError("Invalid branch length")
        return Node(label, length, children)

    root = parse_node()
    if position != len(tokens) - 1 or tokens[position] != ";":
        raise ValueError("Unexpected content after Newick tree")

    def leaf_names(node):
        if not node.children:
            return [node.label]
        return [name for child in node.children for name in leaf_names(child)]

    names = leaf_names(root)
    if sorted(names) != sorted(translation.values()):
        raise ValueError("Newick tips and TRANSLATE table disagree")
    if root.length != 0.0:
        raise ValueError("Expected the pinned rooted tree to have zero root stem")
    return root, translation


def fair_proportion(root, selected_species):
    """Allocate retained root-tip edge lengths equally among retained descendants.

    The resulting fixed predictor-only weights are phylogenetic reweighting,
    not a phylogenetic covariance model or an adjustment for residual dependence.
    Removing other tips does not erase shared ancestral branches.
    """
    selected = set(selected_species)
    contributions = {name: [] for name in selected}
    retained_edges = []

    def visit(node, is_root=False):
        if node.children:
            descendants = [name for child in node.children for name in visit(child)]
        else:
            descendants = [node.label] if node.label in selected else []
        if descendants and not is_root:
            retained_edges.append(node.length)
            contribution = node.length / len(descendants)
            for name in descendants:
                contributions[name].append(contribution)
        return descendants

    found = visit(root, is_root=True)
    if set(found) != selected or len(found) != len(selected):
        raise ValueError("Selected species missing or duplicated in phylogeny")
    values = {name: math.fsum(parts) for name, parts in contributions.items()}
    if any(not math.isfinite(value) or value <= 0 for value in values.values()):
        raise ValueError("Invalid fair-proportion weight")
    return values, math.fsum(retained_edges)


def phylogenetic_blocks(root, species_orders):
    """Construct outcome-blind half-height clade blocks without splitting orders.

    Root depth is zero and its incoming stem is ignored. Each selected tip's
    initial group is the first root-to-tip edge reaching half the median depth.
    Union-find then merges all groups touched by one order, transitively.
    """
    depths = {}

    def measure(node, depth):
        if node.children:
            for child in node.children:
                measure(child, depth + child.length)
        elif node.label in species_orders:
            depths[node.label] = depth

    measure(root, 0.0)
    if set(depths) != set(species_orders):
        raise ValueError("Primary-panel tips missing from tree-depth audit")
    median_depth = statistics.median(depths.values())
    threshold = 0.5 * median_depth
    if not math.isfinite(threshold) or threshold <= 0:
        raise ValueError("Invalid phylogenetic block threshold")
    crossing_groups = {}

    def assign(node, depth, crossing):
        if node.children:
            for child in node.children:
                next_depth = depth + child.length
                next_crossing = crossing
                if crossing is None and depth < threshold <= next_depth:
                    next_crossing = child
                assign(child, next_depth, next_crossing)
        elif node.label in species_orders:
            if crossing is None:
                raise ValueError("Selected tip has no edge crossing block threshold")
            crossing_groups.setdefault(crossing, []).append(node.label)

    assign(root, 0.0, None)
    groups = sorted(crossing_groups.values(), key=lambda group: min(group))
    group_index = {species: index for index, group in enumerate(groups)
                   for species in group}
    parents = list(range(len(groups)))

    def find(index):
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(left, right):
        a, b = find(left), find(right)
        if a != b:
            parents[max(a, b)] = min(a, b)

    for order in sorted(set(species_orders.values())):
        indices = sorted({group_index[species] for species, taxon in species_orders.items()
                          if taxon == order})
        for index in indices[1:]:
            union(indices[0], index)
    components = {}
    for species, index in group_index.items():
        components.setdefault(find(index), []).append(species)
    sorted_components = sorted(components.values(), key=lambda component: min(component))
    blocks = {}
    sizes = {}
    orders = {}
    for index, component in enumerate(sorted_components, start=1):
        label = "clade_" + format(index, "03d")
        sizes[label] = len(component)
        orders[label] = sorted({species_orders[species] for species in component})
        blocks.update({species: label for species in component})
    if any(len({blocks[species] for species, taxon in species_orders.items() if taxon == order}) != 1
           for order in set(species_orders.values())):
        raise ValueError("A taxonomic order was split between phylogenetic blocks")
    audit = {
        "root_stem_length": root.length,
        "root_stem_ignored": True,
        "all_branch_lengths_nonnegative_and_finite": True,
        "minimum_primary_root_to_tip_depth": min(depths.values()),
        "maximum_primary_root_to_tip_depth": max(depths.values()),
        "median_primary_root_to_tip_depth": median_depth,
        "half_height_threshold": threshold,
        "initial_crossing_edge_groups": len(groups),
        "merged_block_count": len(sorted_components),
        "species_per_block": sizes,
        "orders_per_block": orders,
        "orders_are_not_split": True,
        "method": "First edge crossing half the median primary-panel tip depth; union groups sharing any order; labels ordered by smallest species identifier.",
    }
    return blocks, audit


def csv_content(rows, fieldnames):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def derive_outputs(source_texts):
    """Return deterministic output path/content pairs, with no file writes."""
    records = [
        {key.strip(): value.strip() for key, value in row.items()}
        for row in csv.DictReader(io.StringIO(source_texts["SourceData_FINAL.csv"]))
    ]
    names = [row["Species"] for row in records]
    if len(names) != len(set(names)):
        raise ValueError("Duplicate source species")
    records.sort(key=lambda row: row["Species"])
    eligible = []
    exclusions = []
    for row in records:
        missing = [field for field in REQUIRED_FIELDS if not row[field]]
        if missing:
            exclusions.append({"species": row["Species"],
                               "order": row["Taxonomic order"],
                               "reason": "missing:" + ";".join(missing)})
            continue
        eligible.append(row)
    root, translation = parse_nexus(source_texts["phylogeny.tre"])
    fp, retained_length = fair_proportion(root, [row["Species"] for row in eligible])
    blocks, block_audit = phylogenetic_blocks(
        root, {row["Species"]: row["Taxonomic order"] for row in eligible})
    aggregates = []
    residuals = []
    ambiguous = []
    for row in eligible:
        rdi_decimal = Decimal(row["RDI"])
        if rdi_decimal != int(rdi_decimal) or rdi_decimal < 1:
            raise ValueError("Noninteger RDI: " + row["Species"])
        rdi = int(rdi_decimal)
        nominal, lower, upper = count_interval(row["CMR"], rdi)
        mass = float(row["Body mass (kg)"])
        life = float(row["Adult life expectancy"])
        mhc = Decimal(row["MHC I gene family size"])
        diet = row["Consuming mammal prey"]
        if not (math.isfinite(mass) and mass > 0 and math.isfinite(life) and life > 0):
            raise ValueError("Invalid mass/lifespan: " + row["Species"])
        if not mhc.is_finite() or mhc < 0 or mhc != int(mhc):
            raise ValueError("Invalid MHC count: " + row["Species"])
        if diet not in ("0", "1"):
            raise ValueError("Invalid diet: " + row["Species"])
        aggregates.append({
            "species": row["Species"],
            "order": row["Taxonomic order"],
            "known_cause_deaths": rdi,
            "cmr_decimal": row["CMR"],
            "nominal_cancer_deaths": nominal,
            "count_lower": lower,
            "count_upper": upper,
            "log_body_mass_kg": format(math.log(mass), ".17g"),
            "log_adult_life_expectancy_days": format(math.log(life), ".17g"),
            "mammal_prey": diet,
            "log1p_mhc_i": format(math.log1p(int(mhc)), ".17g"),
            "phylo_fair_proportion": format(fp[row["Species"]], ".17g"),
            "phylo_block": blocks[row["Species"]],
        })
        residual = abs(Decimal(row["CMR"]) * rdi - Decimal(nominal))
        residuals.append((residual, row["Species"]))
        if lower != upper:
            ambiguous.append({"species": row["Species"], "cmr_decimal": row["CMR"],
                              "known_cause_deaths": rdi, "nominal": nominal,
                              "lower": lower, "upper": upper})
    order_counts = {}
    for row in aggregates:
        order_counts[row["order"]] = order_counts.get(row["order"], 0) + 1
    worst_residual, worst_species = max(residuals)
    aggregate_csv = csv_content(aggregates, AGGREGATE_FIELDS)
    exclusions_csv = csv_content(exclusions, ("species", "order", "reason"))
    audit = {
        "status": "nonclinical_unreviewed_research_no_medical_or_veterinary_advice",
        "source_repository": SOURCE_REPOSITORY,
        "source_commit": SOURCE_COMMIT,
        "paper_doi": PAPER_DOI,
        "source_files": {Path(name).as_posix(): dict(metadata)
                         for name, metadata in SOURCE_FILES.items()},
        "reuse": {
            "repository_license": None,
            "pinned_repository_license_file": False,
            "article_license": "CC-BY-NC-ND-4.0",
            "license_url": "https://creativecommons.org/licenses/by-nc-nd/4.0/",
            "scope": "Original code and derived numerical species aggregates only; no upstream assets mirrored or relicensed.",
            "individual_records": "Restricted Species360 records neither requested nor accessed.",
            "access_method": "GitHub API/plain UTF-8 text; no browser download",
        },
        "source_rows": len(records),
        "primary_species": len(aggregates),
        "excluded_species": len(exclusions),
        "orders": len(order_counts),
        "species_per_order": dict(sorted(order_counts.items())),
        "known_cause_deaths": sum(row["known_cause_deaths"] for row in aggregates),
        "nominal_cancer_deaths": sum(row["nominal_cancer_deaths"] for row in aggregates),
        "count_lower_total": sum(row["count_lower"] for row in aggregates),
        "count_upper_total": sum(row["count_upper"] for row in aggregates),
        "zero_case_species": sum(row["nominal_cancer_deaths"] == 0 for row in aggregates),
        "count_audit": {
            "method": "Half last decimal unit plus fixed probability tolerance; zero/one exact; nominal ties ROUND_HALF_EVEN; interval endpoints included.",
            "probability_export_tolerance": str(EXPORT_TOLERANCE),
            "unique_count_species": sum(row["count_lower"] == row["count_upper"]
                                        for row in aggregates),
            "ambiguous_count_species": ambiguous,
            "maximum_nominal_integer_residual": str(worst_residual),
            "maximum_residual_species": worst_species,
        },
        "phylogeny": {
            "source_format": "Nexus with TRANSLATE and rooted Newick",
            "source_tip_count": len(translation),
            "matched_primary_tips": len(fp),
            "missing_primary_tips": [],
            "retained_branch_length": retained_length,
            "sum_fair_proportion": math.fsum(fp.values()),
            "method": "Incoming edge length divided by primary-panel descendants, summed root-to-tip; phylogenetic reweighting only, not covariance adjustment.",
            "clade_blocks": block_audit,
        },
        "derived_file_hashes_lf": {
            Path("species_aggregates.csv").as_posix(): normalized_sha256(aggregate_csv),
            Path("exclusions.csv").as_posix(): normalized_sha256(exclusions_csv),
        },
    }
    return {
        Path("species_aggregates.csv").as_posix(): aggregate_csv,
        Path("exclusions.csv").as_posix(): exclusions_csv,
        Path("source_audit.json").as_posix(): json.dumps(audit, indent=2, sort_keys=True) + "\n",
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sources = parser.add_mutually_exclusive_group(required=True)
    sources.add_argument("--fetch", action="store_true", help="Opt in to plain-text HTTPS reads")
    sources.add_argument("--source-dir", type=Path, help="Local external directory holding pinned source texts")
    parser.add_argument("--check", action="store_true", help="Compare existing output without writing")
    arguments = parser.parse_args(argv)
    outputs = derive_outputs(load_sources(arguments.source_dir, arguments.fetch))
    destination = Path(__file__).resolve().parent
    for name, content in outputs.items():
        path = destination / name
        if arguments.check:
            if not path.is_file() or normalize_newlines(path.read_text(encoding="utf-8")) != content:
                raise SystemExit("Derived content differs: " + path.as_posix())
        else:
            with path.open("w", encoding="utf-8", newline="\n") as handle:
                handle.write(content)
    print(("Checked " if arguments.check else "Wrote ") + str(len(outputs))
          + " derived files; upstream files were not saved.")


if __name__ == "__main__":
    main()
