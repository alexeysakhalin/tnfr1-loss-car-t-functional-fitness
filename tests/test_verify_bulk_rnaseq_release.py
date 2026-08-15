from __future__ import annotations

import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from scripts import verify_bulk_rnaseq_release as verifier


ROOT = Path(__file__).resolve().parents[1]
DERIVED_DIR = ROOT / "data" / "experimental" / "bulk_rnaseq" / "derived"


class ReleaseVerifierTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.work = Path(self.temporary_directory.name)
        self.contract = verifier.TableContract(
            generated_relative_path="example.tsv.gz",
            release_filename="example.tsv.gz",
            key_columns=("condition", "gene"),
            exact_columns=("status", "integer_count"),
            numeric_columns=("estimate", "pvalue"),
            decision_rules=(verifier.DecisionRule("pvalue", "lt", 0.05),),
        )
        self.release = pd.DataFrame(
            {
                "condition": ["control", "treated", "treated"],
                "gene": ["A", "A", "B"],
                "status": ["modelled", "modelled", "all_zero"],
                "integer_count": ["12", "14", "0"],
                "estimate": ["1", "-2", "NA"],
                "pvalue": ["0.049", "0.2", "NA"],
            }
        )

    @staticmethod
    def write_table(frame: pd.DataFrame, path: Path) -> None:
        with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
            frame.to_csv(handle, sep="\t", index=False, lineterminator="\n")

    def compare(
        self, generated: pd.DataFrame, release: pd.DataFrame | None = None
    ) -> dict[str, object]:
        generated_path = self.work / "generated.tsv.gz"
        release_path = self.work / "release.tsv.gz"
        self.write_table(generated, generated_path)
        self.write_table(self.release if release is None else release, release_path)
        return verifier.compare_table(generated_path, release_path, self.contract)

    def test_numeric_last_bit_variation_is_accepted(self) -> None:
        generated = self.release.copy()
        generated.loc[0, "estimate"] = "1.0000005"
        generated.loc[0, "pvalue"] = str(0.049 * 10**-0.00005)
        summary = self.compare(generated)
        self.assertTrue(summary["exact_schema_keys_order_text_and_na_masks"])

    def test_semantic_key_and_exact_fields_must_match(self) -> None:
        changed_key = self.release.copy()
        changed_key.loc[1, "gene"] = "C"
        with self.assertRaisesRegex(
            verifier.ReleaseVerificationError, "Exact semantic field differs"
        ):
            self.compare(changed_key)

        changed_count = self.release.copy()
        changed_count.loc[1, "integer_count"] = "15"
        with self.assertRaisesRegex(
            verifier.ReleaseVerificationError, "Exact semantic field differs"
        ):
            self.compare(changed_count)

    def test_numeric_drift_and_missingness_are_rejected(self) -> None:
        drifted = self.release.copy()
        drifted.loc[0, "estimate"] = "1.00001"
        with self.assertRaisesRegex(
            verifier.ReleaseVerificationError, "Numeric drift exceeds"
        ):
            self.compare(drifted)

        changed_missingness = self.release.copy()
        changed_missingness.loc[2, "estimate"] = "0"
        with self.assertRaisesRegex(
            verifier.ReleaseVerificationError, "Missingness differs"
        ):
            self.compare(changed_missingness)

    def test_p_values_are_compared_in_negative_log10_space(self) -> None:
        generated = self.release.copy()
        generated.loc[0, "pvalue"] = str(0.049 * 10**-0.0002)
        with self.assertRaisesRegex(
            verifier.ReleaseVerificationError, "max .delta -log10.p..=0.0001"
        ):
            self.compare(generated)

    def test_threshold_decisions_must_match_even_within_tolerance(self) -> None:
        release = self.release.copy()
        release.loc[0, "pvalue"] = "0.05"
        generated = release.copy()
        generated.loc[0, "pvalue"] = "0.0499999999"
        with self.assertRaisesRegex(
            verifier.ReleaseVerificationError, "Scientific decision differs"
        ):
            self.compare(generated, release)

    def test_manifest_ignores_only_well_formed_output_hashes(self) -> None:
        columns = ["analysis_id", "output_file", "rows", "sha256", "contrast"]
        release = pd.DataFrame(
            [["a", "a.tsv.gz", "46425", "a" * 64, '["group","A","B"]']],
            columns=columns,
        )
        generated = release.copy()
        generated.loc[0, "sha256"] = "b" * 64
        generated_path = self.work / "generated_manifest.tsv"
        release_path = self.work / "release_manifest.tsv"
        generated.to_csv(generated_path, sep="\t", index=False)
        release.to_csv(release_path, sep="\t", index=False)
        verifier.compare_manifest(generated_path, release_path)

        generated.loc[0, "contrast"] = '["group","B","A"]'
        generated.to_csv(generated_path, sep="\t", index=False)
        with self.assertRaisesRegex(
            verifier.ReleaseVerificationError, "differs outside output SHA hashes"
        ):
            verifier.compare_manifest(generated_path, release_path)

    def test_generated_manifest_hash_must_match_its_output(self) -> None:
        output = self.work / "result.tsv.gz"
        with gzip.open(output, "wt", encoding="utf-8") as handle:
            handle.write("Gene_Symbol\tvalue\nA\t1\n")
        manifest = pd.DataFrame(
            [
                {
                    "analysis_id": "example",
                    "output_file": output.name,
                    "rows": "1",
                    "sha256": verifier.sha256_file(output),
                }
            ]
        )
        with (
            mock.patch.object(verifier, "EXPECTED_ANALYSES", 1),
            mock.patch.object(verifier, "EXPECTED_RESULT_ROWS", 1),
        ):
            verifier.verify_complete_exports(self.work, manifest)
            manifest.loc[0, "sha256"] = "0" * 64
            with self.assertRaisesRegex(
                verifier.ReleaseVerificationError,
                "Manifest output SHA-256 does not match",
            ):
                verifier.verify_complete_exports(self.work, manifest)

    def test_metadata_and_environment_remain_exact(self) -> None:
        generated_metadata = self.work / "generated.json"
        release_metadata = self.work / "release.json"
        generated_metadata.write_text(json.dumps({"pydeseq2": "0.5.4"}))
        release_metadata.write_text(
            json.dumps({"pydeseq2": "0.5.4"}, indent=2) + "\n"
        )
        verifier.compare_metadata(generated_metadata, release_metadata)
        generated_metadata.write_text(json.dumps({"pydeseq2": "0.5.5"}))
        with self.assertRaisesRegex(
            verifier.ReleaseVerificationError, "run_metadata.json differs"
        ):
            verifier.compare_metadata(generated_metadata, release_metadata)

        generated_environment = self.work / "generated.txt"
        release_environment = self.work / "release.txt"
        generated_environment.write_text("numpy==2.5.2\n")
        release_environment.write_text("numpy==2.5.2\n")
        verifier.compare_environment(generated_environment, release_environment)
        generated_environment.write_text("numpy==2.5.3\n")
        with self.assertRaisesRegex(
            verifier.ReleaseVerificationError, "environment.freeze.txt differs"
        ):
            verifier.compare_environment(generated_environment, release_environment)

    def test_versioned_tables_satisfy_the_declared_contracts(self) -> None:
        for contract in verifier.TABLE_CONTRACTS:
            path = DERIVED_DIR / contract.release_filename
            with self.subTest(filename=contract.release_filename):
                summary = verifier.compare_table(path, path, contract)
                self.assertEqual(
                    summary["generated_compressed_sha256"],
                    summary["release_compressed_sha256"],
                )
                self.assertEqual(
                    summary["generated_decompressed_sha256"],
                    summary["release_decompressed_sha256"],
                )

    def test_workflow_calls_the_tested_verifier(self) -> None:
        workflow = (ROOT / ".github" / "workflows" / "bulk-rnaseq.yml").read_text()
        self.assertIn("python scripts/verify_bulk_rnaseq_release.py", workflow)
        self.assertIn("--report results/bulk_rnaseq/release_comparison_report.json", workflow)
        self.assertIn("if: always()", workflow)
        self.assertNotIn("generated.read_bytes()", workflow)

    def test_cli_writes_a_failure_report_for_artifact_upload(self) -> None:
        generated_dir = self.work / "generated"
        release_dir = self.work / "release"
        generated_dir.mkdir()
        release_dir.mkdir()
        report_path = generated_dir / "release_comparison_report.json"
        completed = subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts" / "verify_bulk_rnaseq_release.py"),
                "--generated-dir",
                str(generated_dir),
                "--release-dir",
                str(release_dir),
                "--report",
                str(report_path),
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(completed.returncode, 0)
        report = json.loads(report_path.read_text(encoding="utf-8"))
        self.assertEqual(report["status"], "failed")
        self.assertIn("error", report)


if __name__ == "__main__":
    unittest.main()
