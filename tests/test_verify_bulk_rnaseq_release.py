from __future__ import annotations

import gzip
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest
from dataclasses import replace
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
            integer_columns=("integer_count",),
            numeric_columns=(
                "base_mean",
                "estimate",
                "se",
                "wald",
                "pvalue",
                "padj",
            ),
            base_mean_column="base_mean",
            effect_column="estimate",
            standard_error_column="se",
            wald_column="wald",
            adjusted_p_column="padj",
        )
        self.release = pd.DataFrame(
            {
                "condition": ["control", "treated", "treated"],
                "gene": ["A", "B", "C"],
                "status": ["modelled", "modelled", "all_zero"],
                "integer_count": ["12", "14", "0"],
                "base_mean": ["40", "40", "0"],
                "estimate": ["2", "0.5", "NA"],
                "se": ["1", "2", "NA"],
                "wald": ["2", "0.25", "NA"],
                "pvalue": ["0.01", "0.1", "NA"],
                "padj": ["0.049", "0.2", "NA"],
            }
        )

    @staticmethod
    def write_table(frame: pd.DataFrame, path: Path) -> None:
        with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
            frame.to_csv(handle, sep="\t", index=False, lineterminator="\n")

    def compare(
        self,
        generated: pd.DataFrame,
        release: pd.DataFrame | None = None,
        contract: verifier.TableContract | None = None,
    ) -> dict[str, object]:
        generated_path = self.work / "generated.tsv.gz"
        release_path = self.work / "release.tsv.gz"
        self.write_table(generated, generated_path)
        self.write_table(self.release if release is None else release, release_path)
        return verifier.compare_table(
            generated_path,
            release_path,
            self.contract if contract is None else contract,
        )

    def test_large_numeric_drift_is_diagnostic_when_outcomes_are_unchanged(self) -> None:
        generated = self.release.copy()
        generated.loc[0, "estimate"] = "20"
        generated.loc[0, "se"] = "10"
        summary = self.compare(generated)
        delta = summary["numeric_deltas_diagnostic_only"]["estimate"]
        self.assertEqual(delta["max_absolute_delta"], 18.0)
        self.assertTrue(
            summary["exact_schema_keys_order_text_integer_fields_and_na_masks"]
        )

    def test_standalone_adjusted_p_flip_is_reported_when_deg_call_is_unchanged(
        self,
    ) -> None:
        generated = self.release.copy()
        generated.loc[1, "padj"] = "0.049"
        summary = self.compare(generated)
        flips = summary["scientific_outcomes"][
            "standalone_adjusted_p_threshold_flips"
        ]
        self.assertEqual(flips["count"], 1)
        self.assertTrue(flips["outcome_categories_unchanged"])
        self.assertEqual(flips["rows"][0]["key"]["gene"], "B")

    def test_manuscript_label_gene_lfc_guard_is_narrow_and_reported(self) -> None:
        guarded_release = (
            self.release.iloc[[0]]
            .rename(columns={"gene": "gene_symbol"})
            .reset_index(drop=True)
        )
        contract = replace(
            self.contract,
            key_columns=("condition", "gene_symbol"),
            guarded_effect_genes=("A",),
        )
        within = guarded_release.copy()
        within.loc[0, "estimate"] = "2.0005"
        within.loc[0, "wald"] = "2.0005"
        summary = self.compare(within, release=guarded_release, contract=contract)
        guard = summary["scientific_outcomes"][
            "manuscript_label_gene_effect_guard"
        ]
        self.assertEqual(guard["expected_rows"], 1)
        self.assertAlmostEqual(guard["max_absolute_delta"], 0.0005)
        self.assertTrue(guard["all_within_tolerance"])
        self.assertAlmostEqual(guard["rows"][0]["absolute_delta"], 0.0005)

        outside = guarded_release.copy()
        outside.loc[0, "estimate"] = "2.0011"
        outside.loc[0, "wald"] = "2.0011"
        with self.assertRaisesRegex(
            verifier.ReleaseVerificationError,
            "Manuscript label-gene log2FC differs by more than 0.001",
        ):
            self.compare(outside, release=guarded_release, contract=contract)

    def test_schema_keys_text_integer_fields_and_missingness_are_exact(self) -> None:
        changed_key = self.release.copy()
        changed_key.loc[1, "gene"] = "D"
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

        noncanonical_count = self.release.copy()
        noncanonical_count.loc[1, "integer_count"] = "014"
        with self.assertRaisesRegex(
            verifier.ReleaseVerificationError, "canonical non-negative integer"
        ):
            self.compare(noncanonical_count)

        changed_missingness = self.release.copy()
        changed_missingness.loc[2, "estimate"] = "0"
        with self.assertRaisesRegex(
            verifier.ReleaseVerificationError, "Missingness differs"
        ):
            self.compare(changed_missingness)

    def test_base_mean_and_lfc_threshold_masks_are_exact(self) -> None:
        changed_base = self.release.copy()
        changed_base.loc[0, "base_mean"] = "29.999"
        with self.assertRaisesRegex(
            verifier.ReleaseVerificationError, "baseMean >= 30"
        ):
            self.compare(changed_base)

        changed_lfc = self.release.copy()
        changed_lfc.loc[1, "estimate"] = "1.1"
        changed_lfc.loc[1, "se"] = "2.2"
        changed_lfc.loc[1, "wald"] = "0.5"
        with self.assertRaisesRegex(
            verifier.ReleaseVerificationError, "log2FC > 1"
        ):
            self.compare(changed_lfc)

    def test_combined_deg_category_change_is_rejected(self) -> None:
        generated = self.release.copy()
        generated.loc[0, "padj"] = "0.051"
        with self.assertRaisesRegex(
            verifier.ReleaseVerificationError, "combined DEG category"
        ):
            self.compare(generated)

    def test_wald_identity_is_checked_within_each_table(self) -> None:
        generated = self.release.copy()
        generated.loc[0, "wald"] = "1.5"
        with self.assertRaisesRegex(
            verifier.ReleaseVerificationError, "Wald statistic is not log2FC/SE"
        ):
            self.compare(generated)

    def test_numeric_domains_remain_valid(self) -> None:
        generated = self.release.copy()
        generated.loc[1, "padj"] = "1.1"
        with self.assertRaisesRegex(
            verifier.ReleaseVerificationError, r"outside \[0, 1\]"
        ):
            self.compare(generated)

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

    def test_versioned_tables_satisfy_declared_outcome_contracts(self) -> None:
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
                self.assertEqual(
                    summary["scientific_outcomes"][
                        "standalone_adjusted_p_threshold_flips"
                    ]["count"],
                    0,
                )

    def test_prespecified_interaction_genes_are_reported(self) -> None:
        contract = verifier.TABLE_CONTRACTS[2]
        path = DERIVED_DIR / contract.release_filename
        summary = verifier.compare_table(path, path, contract)
        records = summary["scientific_outcomes"][
            "prespecified_interaction_genes"
        ]
        self.assertEqual(
            [record["gene"] for record in records],
            list(verifier.PRESPECIFIED_INTERACTION_GENES),
        )

    def test_exact_check_progress_is_retained_after_table_failure(self) -> None:
        report = verifier.new_report(self.work / "generated", self.work / "release")
        with (
            mock.patch.object(verifier, "compare_manifest", return_value=pd.DataFrame()),
            mock.patch.object(verifier, "verify_complete_exports"),
            mock.patch.object(verifier, "compare_metadata"),
            mock.patch.object(verifier, "compare_environment"),
            mock.patch.object(
                verifier,
                "compare_table",
                side_effect=verifier.ReleaseVerificationError("table failed"),
            ),
        ):
            with self.assertRaisesRegex(
                verifier.ReleaseVerificationError, "table failed"
            ):
                verifier.verify_release(
                    self.work / "generated", self.work / "release", report
                )
        self.assertTrue(
            all(status == "passed" for status in report["exact_checks"].values())
        )
        first_table = verifier.TABLE_CONTRACTS[0].release_filename
        self.assertEqual(report["table_checks"][first_table], "failed")

    def test_workflow_uses_gate_then_committed_canonical_adapters(self) -> None:
        workflow = (ROOT / ".github" / "workflows" / "bulk-rnaseq.yml").read_text()
        self.assertIn("python scripts/verify_bulk_rnaseq_release.py", workflow)
        self.assertIn(
            "--report results/bulk_rnaseq/release_comparison_report.json", workflow
        )
        self.assertIn("python scripts/write_runtime_provenance.py", workflow)
        self.assertIn("if: always()", workflow)
        self.assertIn("sha256sum --check SHA256SUMS", workflow)
        self.assertIn(
            "Confirm the outcome gate and committed release inputs", workflow
        )
        self.assertNotIn("Confirm the equivalence gate", workflow)
        self.assertIn("Render figures from the committed canonical adapters", workflow)
        self.assertNotIn("generated.read_bytes()", workflow)
        self.assertNotIn("cp results/bulk_rnaseq/figure_inputs", workflow)
        for script_name in ("03_figure_1_B_C.R", "04_figure_2_B_suppl_S2D.R"):
            script = (ROOT / "R" / script_name).read_text()
            self.assertIn(
                '"data", "experimental", "bulk_rnaseq", "derived"', script
            )
            label_block = re.search(
                r"label_genes <- c\((.*?)\n\)", script, flags=re.DOTALL
            )
            self.assertIsNotNone(label_block)
            observed_label_genes = tuple(
                re.findall(r'"([A-Z0-9-]+)"', label_block.group(1))
            )
            self.assertEqual(
                observed_label_genes, verifier.MANUSCRIPT_LABEL_GENES
            )

    def test_runtime_provenance_contains_platform_and_numeric_library_details(
        self,
    ) -> None:
        output = self.work / "runtime_provenance.json"
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts" / "write_runtime_provenance.py"),
                "--output",
                str(output),
            ],
            check=True,
        )
        provenance = json.loads(output.read_text(encoding="utf-8"))
        self.assertIn("platform", provenance)
        self.assertIn("python", provenance)
        self.assertIn("zlib", provenance)
        self.assertIn("show_config", provenance["numpy"])
        self.assertIn("OPENBLAS_NUM_THREADS", provenance["thread_environment"])

    def test_cli_writes_stage_specific_failure_progress(self) -> None:
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
        self.assertEqual(
            report["exact_checks"]["manifest_excluding_output_sha256"], "failed"
        )
        self.assertEqual(
            report["exact_checks"]["generated_output_sha256_self_consistency"],
            "not_started",
        )

    def test_downloaded_ci_artifact_satisfies_v2_outcome_contract(self) -> None:
        artifact_value = os.environ.get("BULK_RNASEQ_CI_ARTIFACT_DIR")
        if not artifact_value:
            self.skipTest("BULK_RNASEQ_CI_ARTIFACT_DIR is not set")
        artifact = Path(artifact_value)
        report = verifier.new_report(artifact, DERIVED_DIR)
        summaries = verifier.verify_release(artifact, DERIVED_DIR, report)
        self.assertEqual(len(summaries), 5)
        self.assertTrue(
            all(status == "passed" for status in report["exact_checks"].values())
        )
        self.assertTrue(
            all(status == "passed" for status in report["table_checks"].values())
        )
        flip_counts = [
            summary["scientific_outcomes"][
                "standalone_adjusted_p_threshold_flips"
            ]["count"]
            for summary in summaries
        ]
        self.assertEqual(flip_counts, [8, 12, 0, 7, 2])
        first_venn = summaries[0]["scientific_outcomes"]["venn_membership"]
        second_venn = summaries[1]["scientific_outcomes"]["venn_membership"]
        self.assertEqual(
            [first_venn["sets"][name]["count"] for name in ("TNF", "IFNG", "TNF_IFNG")],
            [458, 1659, 1595],
        )
        self.assertEqual(
            [second_venn["sets"][name]["count"] for name in ("TNF", "IFNG", "TNF_IFNG")],
            [446, 1222, 844],
        )
        self.assertGreater(
            summaries[0]["numeric_deltas_diagnostic_only"][
                verifier.FIGURE_1_EFFECT
            ]["max_absolute_delta"],
            1.0,
        )
        first_label_guard = summaries[0]["scientific_outcomes"][
            "manuscript_label_gene_effect_guard"
        ]
        second_label_guard = summaries[1]["scientific_outcomes"][
            "manuscript_label_gene_effect_guard"
        ]
        self.assertEqual(first_label_guard["expected_rows"], 63)
        self.assertEqual(second_label_guard["expected_rows"], 84)
        self.assertLessEqual(
            first_label_guard["max_absolute_delta"],
            verifier.MANUSCRIPT_LABEL_LFC_ATOL,
        )
        self.assertLessEqual(
            second_label_guard["max_absolute_delta"],
            verifier.MANUSCRIPT_LABEL_LFC_ATOL,
        )
        self.assertTrue(first_label_guard["all_within_tolerance"])
        self.assertTrue(second_label_guard["all_within_tolerance"])


if __name__ == "__main__":
    unittest.main()
