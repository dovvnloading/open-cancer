"""Offline standard-library tests of the GALAXY stop-gate audit."""
from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path, PureWindowsPath

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "research" / "ctdna-false-positive-clustering" / "analysis.py"
SPEC = importlib.util.spec_from_file_location("galaxy_denominator_audit", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


class DenominatorGateTests(unittest.TestCase):
    def setUp(self) -> None:
        self.workspace = tempfile.TemporaryDirectory()
        self.addCleanup(self.workspace.cleanup)
        self.root = Path(self.workspace.name)
        (self.root / "data").mkdir()
        source = SCRIPT.parent
        for name in ("PROTOCOL.md", "data/printed_counts.csv", "data/denominator_gate.json"):
            destination = self.root / name
            destination.write_bytes((source / name).read_bytes())

    def test_untestable_and_no_occupancy(self) -> None:
        result = module.audit(self.root)
        self.assertEqual(result["verdict"], "UNTESTABLE")
        self.assertEqual(result["blocked_gates"], ["G2", "G3", "G4"])
        self.assertIs(result["occupancy_executed"], False)
        self.assertIsNone(result["occupancy_p_value"])
        self.assertIsNone(result["within_person_concentration"])
        self.assertEqual(result["sample_remainder_unassigned"], 18)

    def test_two_different_126_groups_and_two_different_10s(self) -> None:
        table = module.audit(self.root)["patient_table"]
        self.assertEqual(table, {
            "ever_positive_recurred": 54, "ever_positive_nonrecurred": 10,
            "serial_negative_recurred": 10, "serial_negative_nonrecurred": 116,
            "nonrecurrent_total": 126,
        })
        self.assertEqual(table["serial_negative_recurred"] +
                         table["serial_negative_nonrecurred"], 126)
        self.assertEqual(table["ever_positive_nonrecurred"] +
                         table["serial_negative_nonrecurred"], 126)
        self.assertNotEqual(table["serial_negative_recurred"], 0)

    def test_exact_rationals_not_float_estimates(self) -> None:
        result = module.audit(self.root)
        self.assertEqual(Fraction(result["patient_specificity_exact"]), Fraction(58, 63))
        self.assertEqual(Fraction(result["sample_specificity_exact"]), Fraction(619, 637))
        self.assertIsInstance(result["sample_specificity_exact"], str)

    def test_only_allowed_hash_paths_and_text_normalization(self) -> None:
        result = module.audit(self.root)
        keys = list(result["source_hashes_sha256_lf"])
        self.assertEqual(keys, ["PROTOCOL.md", "data/printed_counts.csv",
                                "data/denominator_gate.json"])
        self.assertTrue(all("\\" not in key for key in keys))
        self.assertEqual(module.normalized_sha256(b"a\r\nb\rc\n"),
                         module.normalized_sha256(b"a\nb\nc\n"))
        self.assertEqual(PureWindowsPath(r"research\ctdna-false-positive-clustering\data\x.csv").as_posix(),
                         "research/ctdna-false-positive-clustering/data/x.csv")
        for file in ("PROTOCOL.md", "data/printed_counts.csv", "data/denominator_gate.json"):
            p = self.root / file
            normalized = p.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
            p.write_bytes(normalized.replace(b"\n", b"\r\n"))
        self.assertEqual(module.audit(self.root)["source_hashes_sha256_lf"],
                         result["source_hashes_sha256_lf"])

    def test_missing_or_duplicate_count_fail_closed(self) -> None:
        path = self.root / "data/printed_counts.csv"
        text = path.read_text(encoding="utf-8")
        with self.subTest("missing"):
            path.write_text(text.replace("surveillance_patients,190,", "absent,190,"), encoding="utf-8")
            with self.assertRaises(ValueError):
                module.audit(self.root)
        with self.subTest("duplicate"):
            path.write_text(text + text.splitlines(keepends=True)[1], encoding="utf-8")
            with self.assertRaises(ValueError):
                module.audit(self.root)

    def test_changed_count_or_unattributed_source_fail_closed(self) -> None:
        path = self.root / "data/printed_counts.csv"
        text = path.read_text(encoding="utf-8")
        path.write_text(text.replace("serial_negative_recurrences,10,",
                                     "serial_negative_recurrences,11,"), encoding="utf-8")
        with self.assertRaises(ValueError):
            module.audit(self.root)
        path.write_text(text.replace("CC BY 4.0", "unknown"), encoding="utf-8")
        with self.assertRaises(ValueError):
            module.audit(self.root)

    def test_gate_status_must_remain_source_documented(self) -> None:
        path = self.root / "data/denominator_gate.json"
        data = json.loads(path.read_text(encoding="utf-8"))
        data["gates"]["G2"]["status"] = "supported"
        data["gates"]["G3"]["status"] = "supported"
        data["gates"]["G4"]["status"] = "supported"
        path.write_text(json.dumps(data), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "occupancy analysis is not implemented"):
            module.audit(self.root)

    def test_no_source_records_or_person_histories(self) -> None:
        data = json.loads((self.root / "data/denominator_gate.json").read_text(encoding="utf-8"))
        self.assertEqual(set(data["gates"]), {"G1", "G2", "G3", "G4"})
        self.assertTrue(all(entry["status"] in {"supported", "not_documented"}
                            for entry in data["gates"].values()))
        self.assertNotIn("samples_by_person", data)
        self.assertNotIn("patient_id", (self.root / "data/printed_counts.csv").read_text())


if __name__ == "__main__":
    unittest.main()
