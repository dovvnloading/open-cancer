"""Offline checks on the frozen EGFR test, using invented models and calls."""
import csv
import hashlib
import importlib.util
import io
import json
import tempfile
import unittest
import zipfile
from pathlib import Path

STUDY = Path(__file__).resolve().parents[1] / "research" / "organoid-dependency-priority-test"
SPEC = importlib.util.spec_from_file_location("organoid_priority_analysis", STUDY / "analyze.py")
analysis = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(analysis)


def delimited(rows, newline, sep=","):
    output = io.StringIO(newline="")
    csv.writer(output, delimiter=sep, lineterminator=newline).writerows(rows)
    return output.getvalue().encode("utf-8")


def fixture(root, newline, *, missing_first_mutant=False, duplicate_gene=False):
    """Two groups of five plus records that should drop out at distinct gates."""
    eligible_pos = [f"m{i}" for i in range(5)]
    eligible_ref = [f"r{i}" for i in range(5)]
    discarded = ["zz_duplicate", "no_wgs", "ras", "other_braf", "missing_call", "zz_spare"]
    models = eligible_pos + eligible_ref + discarded
    donors = {"zz_duplicate": "d-m0", "zz_spare": "d-missing_call"}
    metadata = [["model_id", "model_type", "cancer_type", "patient_id"]]
    coverage = [["model_id", "Mutation Sanger Organoid (WGS)"]]
    variants = [["model_id", "gene_symbol", "protein_mutation", "cancer_driver", "source"]]
    for model in models:
        metadata.append([model, "Organoid", "Colorectal Carcinoma", donors.get(model, "d-" + model)])
        coverage.append([model, str(model != "no_wgs")])
    for model in eligible_pos + ["zz_duplicate", "ras"]:
        variants.append([model, "BRAF", "p.V600E", "t", "Sanger"])
    variants += [["ras", "NRAS", "p.Q61K", "t", "Sanger"],
                 ["other_braf", "BRAF", "p.K601E", "t", "Sanger"],
                 ["r0", "BRAF", "p.V600E", "t", "Broad"]]
    missing = {"missing_call", "m0"} if missing_first_mutant else {"missing_call"}
    calls = ["" if m in missing else str(int(m in eligible_pos)) for m in models]
    scores = [["model_name", "", ""] + models, ["model_id", "", ""] + models,
              ["source", "", ""] + ["Sanger"] * len(models),
              ["qc_pass", "", ""] + ["TRUE"] * len(models),
              ["gene_id", "symbol", "ensembl_gene_id"] + [""] * len(models),
              ["g", "EGFR", "e"] + calls]
    if duplicate_gene:
        scores.append(scores[-1])
    hashes = {}
    for filename, rows in [("model_list_20261013.csv", metadata),
                           ("model_dataset_availability_20260727.csv", coverage)]:
        (root / filename).write_bytes(delimited(rows, newline))
        hashes[filename] = hashlib.sha256(delimited(rows, "\n")).hexdigest()
    tables = [
        ("Organoids_fitness_scores_20260723.zip", "Organoids_fitness_scores_binary_matrix_20260723.tsv", scores, "\t"),
        ("mutations_summary_20260724.zip", "mutations_summary_20260724.csv", variants, ","),
    ]
    for archive, member, rows, sep in tables:
        with zipfile.ZipFile(root / archive, "w") as zf:
            zf.writestr(member, delimited(rows, newline, sep))
        key = (Path(archive) / member).as_posix()
        hashes[key] = hashlib.sha256(delimited(rows, "\n", sep)).hexdigest()
    return hashes


class OrganoidPriorityTest(unittest.TestCase):
    def test_synthetic_analysis_and_newline_normalization(self):
        for ending in ("\n", "\r\n"):
            with self.subTest(ending=repr(ending)), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                hashes = fixture(root, ending)
                result = analysis.analyze(root)
                self.assertEqual(result["inputs_sha256"], hashes)
                self.assertEqual(result["verdict"], "supports_basal_prediction")
                self.assertEqual(result["groups"]["BRAF_V600E"]["tested"], 5)
                self.assertEqual(result["groups"]["BRAF_V600E"]["essential"], 5)
                self.assertEqual(result["groups"]["BRAF_reference"]["tested"], 5)
                self.assertEqual(result["groups"]["BRAF_reference"]["essential"], 0)
                # With fixed margins, two extreme tables each have mass 1 / C(10,5).
                self.assertTrue(0.0079 < result["statistics"]["fisher_p_two_sided"] < 0.0080)

    def test_exact_probabilities_and_interval(self):
        for a, n, b, m in [(1, 10, 11, 14), (11, 14, 1, 10)]:
            self.assertTrue(0.0027 < analysis.fisher_two_sided(a, n, b, m) < 0.0029)
        self.assertTrue(0.9999 < analysis.fisher_two_sided(0, 6, 0, 8) <= 1)
        low, high = analysis.newcombe_difference(0, 5, 5, 5)
        self.assertTrue(-1 <= low < -0.999)
        self.assertTrue(-0.39 < high < -0.38)

    def test_variants_are_exact_and_exclusions_take_precedence(self):
        v600e = ("BRAF", "p.V600E", "t")
        self.assertEqual(analysis.classify_genotype([v600e, v600e]), "BRAF_V600E")
        for variant in [("BRAF", "p.V600Efs", "t"), ("BRAF", "", "f")]:
            self.assertEqual(analysis.classify_genotype([v600e, variant]), "other_or_ambiguous_braf")
        for gene in ("KRAS", "NRAS", "HRAS"):
            self.assertEqual(analysis.classify_genotype([v600e, (gene, "p.X", "t")]), "ras_driver")
        self.assertEqual(analysis.classify_genotype([("KRAS", "p.X", "f")]), "BRAF_reference")

    def test_failures_stop_inference(self):
        for scenario in ("missing_input", "duplicate_gene", "small_group"):
            with self.subTest(scenario=scenario), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                if scenario != "missing_input":
                    fixture(root, "\n", missing_first_mutant=scenario == "small_group",
                            duplicate_gene=scenario == "duplicate_gene")
                result = analysis.analyze(root)
                self.assertEqual(result["verdict"], "feasibility_failure")
                self.assertIsNone(result["statistics"])
                if scenario == "small_group":
                    self.assertEqual(result["groups"]["BRAF_V600E"]["tested"], 4)

    def test_published_aggregate(self):
        result = json.loads((STUDY / "derived" / "egfr.json").read_text(encoding="utf-8"))
        self.assertEqual(result["research_status"], "nonclinical, unreviewed computational research")
        self.assertEqual(result["groups"], {"BRAF_V600E": {"essential": 6, "tested": 8},
                                            "BRAF_reference": {"essential": 31, "tested": 35}})
        self.assertEqual(result["verdict"], "inconclusive")
        attrition = result["attrition"]
        self.assertEqual(sum(attrition.values()) - attrition["input_screen_columns"],
                         attrition["input_screen_columns"])
        self.assertEqual(attrition["evaluable_donors"], 43)
        stats = result["statistics"]
        self.assertTrue(-0.14 < stats["risk_difference"] < -0.13)
        self.assertTrue(0.30 < stats["fisher_p_two_sided"] < 0.32)
        lo, hi = stats["risk_difference_ci95"]
        self.assertTrue(-0.49 < lo < -0.48 and 0.09 < hi < 0.10)
        for key, digest in result["inputs_sha256"].items():
            self.assertEqual(Path(key).as_posix(), key)
            self.assertRegex(digest, r"^[0-9a-f]{64}$")
        self.assertNotRegex(json.dumps(result), r"(?:SIDM|SIDP|SIDS)\d+")

    def test_priority_evidence_and_search_coverage(self):
        with (STUDY / "triage.csv").open(encoding="utf-8", newline="") as stream:
            evidence = list(csv.DictReader(stream))
        with (STUDY / "derived" / "search_log.csv").open(encoding="utf-8", newline="") as stream:
            searches = list(csv.DictReader(stream))
        verdicts = {"a": "PARTIALLY KNOWN", "b": "NOT FOUND", "c": "PARTIALLY KNOWN", "d": "PUBLISHED"}
        self.assertEqual({r["claim_id"] for r in evidence}, set(verdicts))
        for row in evidence:
            self.assertEqual(row["verdict"], verdicts[row["claim_id"]])
            self.assertTrue(row["doi"].startswith("10.") and row["source_url"].startswith("https://"))
            self.assertTrue(all(row[key] for key in ("citation", "specific_finding", "limit")))
        for claim in verdicts:
            rows = [r for r in searches if claim in r["claim_id"]]
            for name in ("pubmed", "europe pmc", "biorxiv", "sanger"):
                self.assertTrue(any(name in r["source"].lower() for r in rows), (claim, name))
            self.assertTrue(all(r["date"] == "2026-11-22" and r["exact_query"] for r in rows))


if __name__ == "__main__":
    unittest.main()
