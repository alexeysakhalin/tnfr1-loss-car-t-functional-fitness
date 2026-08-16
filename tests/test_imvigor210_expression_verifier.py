from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VERIFIER_PATH = ROOT / "scripts" / "verify_imvigor210_expression.py"
CONTRACT_PATH = ROOT / "resources" / "IMvigor210_expression_semantic_contract_v1.json"
IMPACT_PATH = ROOT / "resources" / "IMvigor210_fixed6_canonicalization_impact_v1.json"
PREPARER_PATH = ROOT / "scripts" / "prepare_open_cohort_analysis_tables.py"


def load_verifier():
    specification = importlib.util.spec_from_file_location(
        "verify_imvigor210_expression_for_test", VERIFIER_PATH
    )
    assert specification is not None
    assert specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


VERIFIER = load_verifier()


def load_preparer():
    specification = importlib.util.spec_from_file_location(
        "prepare_open_cohort_analysis_tables_for_semantic_test", PREPARER_PATH
    )
    assert specification is not None
    assert specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


PREPARER = load_preparer()


def write_csv(path: Path, rows: list[list[str]], line_ending: str = "\n") -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator=line_ending)
        writer.writerows(rows)


def contract_for(path: Path) -> dict[str, object]:
    skeleton: dict[str, object] = {
        "semantic_contract_version": 1,
        "shape": {"rows": 2, "columns": 2, "values": 4},
        "numeric_constraints": {"zero_count": 1},
        "ordered_identifier_digests": {
            "row_ids_sha256": "0" * 64,
            "column_ids_sha256": "0" * 64,
        },
        "semantic_digests": {
            "required_scale": 6,
            "diagnostic_scales": [7, 8],
            "reference_sha256": {
                "fixed6": "0" * 64,
                "fixed7": "0" * 64,
                "fixed8": "0" * 64,
            },
        },
        "analysis_canonicalization": {
            "required": True,
            "scale": 6,
            "rounding": "ROUND_HALF_UP",
        },
    }
    observed = VERIFIER.scan_expression(path, skeleton)
    identifiers = skeleton["ordered_identifier_digests"]
    assert isinstance(identifiers, dict)
    identifiers["row_ids_sha256"] = observed["row_ids_sha256"]
    identifiers["column_ids_sha256"] = observed["column_ids_sha256"]
    semantic = skeleton["semantic_digests"]
    assert isinstance(semantic, dict)
    semantic["reference_sha256"] = observed["semantic_sha256"]
    return skeleton


def write_contract(path: Path, contract: dict[str, object]) -> None:
    path.write_text(json.dumps(contract) + "\n", encoding="utf-8")


class Imvigor210ExpressionVerifierTests(unittest.TestCase):
    def test_embedded_reference_framing_self_test(self) -> None:
        VERIFIER.run_self_test()

    def test_repository_contract_constants_and_privacy(self) -> None:
        contract = VERIFIER.load_contract(CONTRACT_PATH)
        self.assertEqual(
            hashlib.sha256(CONTRACT_PATH.read_bytes()).hexdigest(),
            "1b484e0c2f33565e967d395b45e30d114a4495d7fc604df5aaa84ccd54f1e503",
        )
        self.assertEqual(contract["shape"], {"rows": 31286, "columns": 348, "values": 10887528})
        self.assertEqual(contract["numeric_constraints"]["zero_count"], 4363672)
        self.assertEqual(
            contract["ordered_identifier_digests"]["row_ids_sha256"],
            "cf43c23e5333fbf88c1f957dd31e5bcb123f1a2b3e33bd928ffaafe3e799a53d",
        )
        self.assertEqual(
            contract["ordered_identifier_digests"]["column_ids_sha256"],
            "c034a832b17c78d6c7eb65526322169653c89f093ec7fc0e9ccae9f5204c8e70",
        )
        self.assertEqual(
            contract["semantic_digests"]["reference_sha256"],
            {
                "fixed6": "2c62f046f93f13ffb9034c2cf7887c322cb61144c653e31f260cbbfebdd7dcfa",
                "fixed7": "3acca75b245771e0f3385ead31a0477d8366576244726cf84055e535520ae9bd",
                "fixed8": "b61085abe13e5f42949b78618509398a11e5e4b190343b423977ce9da4ace341",
            },
        )
        self.assertFalse(
            contract["privacy"]["contract_contains_identifiers_or_expression_cells"]
        )
        self.assertEqual(
            contract["analysis_canonicalization"],
            {
                "required": True,
                "scale": 6,
                "rounding": "ROUND_HALF_UP",
                "operation": "Before identifier mapping, aggregation, ranking or selected-expression output, replace every parsed expression cell by fixed6_scaled_integer / 1000000.0.",
                "guarantee": "All CSV files accepted by the required fixed6 semantic digest enter downstream analysis as the same ordered numeric matrix.",
            },
        )

    def test_historical_canonicalization_impact_is_locked_and_nonidentifying(self) -> None:
        self.assertEqual(
            hashlib.sha256(IMPACT_PATH.read_bytes()).hexdigest(),
            "2ba5e6d5c33fcebb478944a8c3fd24c14ed38f30c545b19e68a23f38d73ddd6d",
        )
        impact = json.loads(IMPACT_PATH.read_text(encoding="utf-8"))
        comparison = impact["comparison"]
        self.assertEqual(comparison["selected_expression_rows"], 61944)
        self.assertTrue(comparison["ordered_keys_equal"])
        self.assertEqual(comparison["expr_value_rows_changed"], 59954)
        self.assertEqual(
            comparison["expr_value_max_absolute_change"],
            "4.999993601373376e-7",
        )
        self.assertEqual(comparison["rank_percentile_rows_changed"], 0)
        self.assertEqual(
            comparison["new_selected_expression_sha256"],
            "63f3af1f5cc977df95448d5475a3a7ace43e15a79a99b4836b2e0321ea818bac",
        )
        self.assertEqual(
            impact["downstream_checks"]["figure_5c_tcell_rank_scores_changed"],
            0,
        )
        self.assertEqual(
            impact["downstream_checks"][
                "supplementary_s6_nominal_or_bh_significance_changes"
            ],
            0,
        )
        self.assertFalse(impact["privacy"]["contains_sample_identifiers"])
        self.assertFalse(impact["privacy"]["contains_source_feature_identifiers"])
        self.assertTrue(impact["privacy"]["contains_prespecified_gene_symbols"])
        self.assertFalse(impact["privacy"]["contains_expression_cells"])

    def test_fixed6_accepts_harmless_text_and_sub_resolution_changes(self) -> None:
        base_rows = [
            ["", "sample-secret-A", "sample-secret-B"],
            ["feature-secret-A", "0", "1.23456744"],
            ["feature-secret-B", "2.5", "3"],
        ]
        equivalent_rows = [
            ["row_label", "sample-secret-A", "sample-secret-B"],
            ["feature-secret-A", "0.0000000", "1.23456749"],
            ["feature-secret-B", "2.50000000e0", "3.0000000"],
        ]
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            base = directory / "base.csv"
            equivalent = directory / "equivalent.csv"
            contract_path = directory / "contract.json"
            report_path = directory / "report.json"
            write_csv(base, base_rows, "\n")
            write_csv(equivalent, equivalent_rows, "\r\n")
            write_contract(contract_path, contract_for(base))

            report = VERIFIER.verify_expression(
                equivalent, contract_path, report_path
            )

            self.assertEqual(report["status"], "pass")
            self.assertTrue(report["required_compatibility"]["match"])
            self.assertFalse(report["diagnostic_reference_matches"]["fixed7"])
            self.assertEqual(
                PREPARER.canonical_imvigor_fixed6("1.23456744"),
                PREPARER.canonical_imvigor_fixed6("1.23456749"),
            )
            rendered = report_path.read_text(encoding="utf-8")
            for sensitive in (
                "sample-secret-A",
                "sample-secret-B",
                "feature-secret-A",
                "feature-secret-B",
                "1.23456749",
            ):
                self.assertNotIn(sensitive, rendered)

    def test_required_fixed6_change_fails_with_redacted_report(self) -> None:
        base_rows = [
            ["", "sample-secret-A", "sample-secret-B"],
            ["feature-secret-A", "0", "1.23456744"],
            ["feature-secret-B", "2.5", "3"],
        ]
        changed_rows = [row[:] for row in base_rows]
        changed_rows[1][2] = "1.2345681"
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            base = directory / "base.csv"
            changed = directory / "changed.csv"
            contract_path = directory / "contract.json"
            report_path = directory / "report.json"
            write_csv(base, base_rows)
            write_csv(changed, changed_rows)
            write_contract(contract_path, contract_for(base))
            with self.assertRaisesRegex(
                VERIFIER.VerificationError, "required_fixed_scale_semantics"
            ) as caught:
                VERIFIER.verify_expression(changed, contract_path, report_path)
            report = json.loads(report_path.read_text(encoding="utf-8"))
            self.assertEqual(report["status"], "fail")
            self.assertIn(
                "required_fixed_scale_semantics", report["failed_checks"]
            )
            rendered = json.dumps(report) + str(caught.exception)
            for sensitive in ("sample-secret", "feature-secret", "1.2345681"):
                self.assertNotIn(sensitive, rendered)

    def test_identifier_order_change_fails_without_identifier_disclosure(self) -> None:
        base_rows = [
            ["", "sample-secret-A", "sample-secret-B"],
            ["feature-secret-A", "0", "1"],
            ["feature-secret-B", "2", "3"],
        ]
        reordered_rows = [base_rows[0], base_rows[2], base_rows[1]]
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            base = directory / "base.csv"
            reordered = directory / "reordered.csv"
            contract_path = directory / "contract.json"
            report_path = directory / "report.json"
            write_csv(base, base_rows)
            write_csv(reordered, reordered_rows)
            write_contract(contract_path, contract_for(base))
            with self.assertRaises(VERIFIER.VerificationError) as caught:
                VERIFIER.verify_expression(reordered, contract_path, report_path)
            report = json.loads(report_path.read_text(encoding="utf-8"))
            self.assertIn("ordered_row_identifiers", report["failed_checks"])
            rendered = json.dumps(report) + str(caught.exception)
            self.assertNotIn("feature-secret-A", rendered)
            self.assertNotIn("feature-secret-B", rendered)

    def test_malformed_numeric_failure_still_writes_redacted_report(self) -> None:
        valid_rows = [
            ["", "sample-secret-A", "sample-secret-B"],
            ["feature-secret-A", "0", "1"],
            ["feature-secret-B", "2", "3"],
        ]
        malformed_rows = [row[:] for row in valid_rows]
        malformed_rows[1][2] = "patient-value-secret"
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            valid = directory / "valid.csv"
            malformed = directory / "malformed.csv"
            contract_path = directory / "contract.json"
            report_path = directory / "report.json"
            write_csv(valid, valid_rows)
            write_csv(malformed, malformed_rows)
            write_contract(contract_path, contract_for(valid))
            with self.assertRaises(VERIFIER.VerificationError) as caught:
                VERIFIER.verify_expression(malformed, contract_path, report_path)
            report_text = report_path.read_text(encoding="utf-8")
            self.assertIn("input_structure_or_numeric_validity", report_text)
            self.assertNotIn("patient-value-secret", report_text)
            self.assertNotIn("patient-value-secret", str(caught.exception))

    def test_negative_and_duplicate_inputs_are_rejected_generically(self) -> None:
        base_rows = [
            ["", "sample-secret-A", "sample-secret-B"],
            ["feature-secret-A", "0", "1"],
            ["feature-secret-B", "2", "3"],
        ]
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            base = directory / "base.csv"
            write_csv(base, base_rows)
            contract = contract_for(base)
            for filename, rows, message in (
                (
                    "negative.csv",
                    [base_rows[0], ["feature-secret-A", "-0.1", "1"], base_rows[2]],
                    "negative numeric field",
                ),
                (
                    "duplicate.csv",
                    [base_rows[0], base_rows[1], base_rows[1]],
                    "not unique",
                ),
            ):
                path = directory / filename
                write_csv(path, rows)
                with self.subTest(filename=filename):
                    with self.assertRaisesRegex(VERIFIER.VerificationError, message):
                        VERIFIER.scan_expression(path, contract)

    def test_fixed6_analysis_canonicalization_uses_half_up_and_redacts_overflow(self) -> None:
        self.assertEqual(PREPARER.canonical_imvigor_fixed6("0.0000005"), 0.000001)
        self.assertEqual(PREPARER.canonical_imvigor_fixed6("1.2345675"), 1.234568)
        self.assertEqual(PREPARER.canonical_imvigor_fixed6("0.00000049"), 0.0)
        secret = "1e999999"
        with self.assertRaises(ValueError) as caught:
            PREPARER.canonical_imvigor_fixed6(secret)
        self.assertNotIn(secret, str(caught.exception))

    def test_subresolution_swaps_have_identical_canonical_matrices_and_ranks(self) -> None:
        import numpy as np

        first = np.array(
            [[1.00000041, 2.0], [1.00000049, 3.0], [4.0, 1.0]], dtype=float
        )
        swapped = np.array(
            [[1.00000049, 2.0], [1.00000041, 3.0], [4.0, 1.0]], dtype=float
        )
        self.assertFalse(
            np.array_equal(
                PREPARER.transcriptome_rank_percentiles(first),
                PREPARER.transcriptome_rank_percentiles(swapped),
            )
        )
        canonical_first = np.array(
            [
                [PREPARER.canonical_imvigor_fixed6(value) for value in row]
                for row in first.astype(str)
            ]
        )
        canonical_swapped = np.array(
            [
                [PREPARER.canonical_imvigor_fixed6(value) for value in row]
                for row in swapped.astype(str)
            ]
        )
        self.assertTrue(np.array_equal(canonical_first, canonical_swapped))
        self.assertTrue(
            np.array_equal(
                PREPARER.transcriptome_rank_percentiles(canonical_first),
                PREPARER.transcriptome_rank_percentiles(canonical_swapped),
            )
        )

    def test_analysis_scale_and_rounding_must_match_verifier_contract(self) -> None:
        contract = VERIFIER.load_contract(CONTRACT_PATH)
        PREPARER.assert_imvigor_analysis_contract(contract)
        altered = json.loads(json.dumps(contract))
        altered["analysis_canonicalization"]["scale"] = 7
        with self.assertRaisesRegex(RuntimeError, "policies differ"):
            PREPARER.assert_imvigor_analysis_contract(altered)

    def test_column_order_and_duplicate_columns_are_rejected(self) -> None:
        base_rows = [
            ["", "sample-secret-A", "sample-secret-B"],
            ["feature-secret-A", "0", "1"],
            ["feature-secret-B", "2", "3"],
        ]
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            base = directory / "base.csv"
            contract_path = directory / "contract.json"
            write_csv(base, base_rows)
            write_contract(contract_path, contract_for(base))
            cases = {
                "reordered.csv": [
                    ["", "sample-secret-B", "sample-secret-A"],
                    base_rows[1],
                    base_rows[2],
                ],
                "duplicate.csv": [
                    ["", "sample-secret-A", "sample-secret-A"],
                    base_rows[1],
                    base_rows[2],
                ],
            }
            for filename, rows in cases.items():
                path = directory / filename
                report = directory / f"{filename}.json"
                write_csv(path, rows)
                with self.subTest(filename=filename):
                    with self.assertRaises(VERIFIER.VerificationError) as caught:
                        VERIFIER.verify_expression(path, contract_path, report)
                    rendered = report.read_text(encoding="utf-8") + str(caught.exception)
                    self.assertNotIn("sample-secret-A", rendered)
                    self.assertNotIn("sample-secret-B", rendered)

    def test_row_counts_width_and_nonfinite_values_are_rejected(self) -> None:
        base_rows = [
            ["", "sample-secret-A", "sample-secret-B"],
            ["feature-secret-A", "0", "1"],
            ["feature-secret-B", "2", "3"],
        ]
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            base = directory / "base.csv"
            contract_path = directory / "contract.json"
            write_csv(base, base_rows)
            write_contract(contract_path, contract_for(base))
            cases = {
                "too-few.csv": base_rows[:2],
                "too-many.csv": [*base_rows, ["feature-secret-C", "4", "5"]],
                "wrong-width.csv": [base_rows[0], base_rows[1][:-1], base_rows[2]],
                "nan.csv": [base_rows[0], ["feature-secret-A", "NaN", "1"], base_rows[2]],
                "inf.csv": [base_rows[0], ["feature-secret-A", "Infinity", "1"], base_rows[2]],
                "overflow.csv": [base_rows[0], ["feature-secret-A", "1e999999", "1"], base_rows[2]],
            }
            for filename, rows in cases.items():
                path = directory / filename
                report = directory / f"{filename}.json"
                write_csv(path, rows)
                with self.subTest(filename=filename):
                    with self.assertRaises(VERIFIER.VerificationError) as caught:
                        VERIFIER.verify_expression(path, contract_path, report)
                    rendered = report.read_text(encoding="utf-8") + str(caught.exception)
                    for sensitive in (
                        "sample-secret-A",
                        "feature-secret-A",
                        "1e999999",
                    ):
                        self.assertNotIn(sensitive, rendered)

    def test_zero_mask_mutation_is_rejected_even_when_fixed6_is_unchanged(self) -> None:
        base_rows = [
            ["", "sample-secret-A", "sample-secret-B"],
            ["feature-secret-A", "0", "1"],
            ["feature-secret-B", "2", "3"],
        ]
        mutated_rows = [row[:] for row in base_rows]
        mutated_rows[1][1] = "0.0000004"
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            base = directory / "base.csv"
            mutated = directory / "mutated.csv"
            contract_path = directory / "contract.json"
            report_path = directory / "report.json"
            write_csv(base, base_rows)
            write_csv(mutated, mutated_rows)
            write_contract(contract_path, contract_for(base))
            with self.assertRaises(VERIFIER.VerificationError):
                VERIFIER.verify_expression(mutated, contract_path, report_path)
            report = json.loads(report_path.read_text(encoding="utf-8"))
            self.assertEqual(
                report["semantic_sha256"]["fixed6"],
                json.loads(contract_path.read_text())["semantic_digests"][
                    "reference_sha256"
                ]["fixed6"],
            )
            self.assertIn("zero_count", report["failed_checks"])

    def test_invalid_contract_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "contract.json"
            path.write_text(
                json.dumps(
                    {
                        "semantic_contract_version": 1,
                        "shape": {"rows": 2, "columns": 2, "values": 5},
                        "numeric_constraints": {"zero_count": 0},
                        "ordered_identifier_digests": {},
                        "semantic_digests": {},
                        "analysis_canonicalization": {},
                    }
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                VERIFIER.VerificationError, "shape is inconsistent"
            ):
                VERIFIER.load_contract(path)


if __name__ == "__main__":
    unittest.main()
