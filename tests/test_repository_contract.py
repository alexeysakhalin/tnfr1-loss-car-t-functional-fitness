from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


class RepositoryContractTests(unittest.TestCase):
    def test_git_tracks_no_raw_or_sample_level_tables(self) -> None:
        completed = subprocess.run(
            ["git", "ls-files"], cwd=ROOT, check=True,
            text=True, capture_output=True,
        )
        tracked = set(completed.stdout.splitlines())
        forbidden_prefixes = (
            "data/raw/", "data/processed/", "results/",
        )
        self.assertFalse(
            [path for path in tracked if path.startswith(forbidden_prefixes)]
        )
        allowed_analysis = {
            "data/analysis/checkmate_c6_global_gene_models.tsv.gz",
            "data/analysis/checkmate_c6_group_balance.tsv",
            "data/analysis/checkmate_c6_aggregate_qc.json",
        }
        self.assertFalse(
            [
                path for path in tracked
                if path.startswith("data/analysis/") and path not in allowed_analysis
            ]
        )

    def test_source_manifest_is_complete_and_checksum_pinned(self) -> None:
        path = ROOT / "data" / "source_manifest.tsv"
        with path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        source_ids = [row["source_id"] for row in rows]
        self.assertEqual(len(source_ids), len(set(source_ids)))
        self.assertGreaterEqual(len(rows), 8)
        imvigor_rows = [
            row for row in rows if row["source_id"].startswith("imvigor210_")
        ]
        self.assertEqual(len(imvigor_rows), 3)
        self.assertEqual(
            {row["cohort_id"] for row in imvigor_rows}, {"IMvigor210_BLCA"}
        )
        for row in rows:
            self.assertEqual(len(row["sha256"]), 64)
            int(row["sha256"], 16)
            self.assertGreater(int(row["size_bytes"]), 0)
            self.assertIn(
                row["repository_policy"],
                {"local_only", "metadata_only", "mapping_resource_committed"},
            )

    def test_signature_resource_has_ten_fixed_twenty_gene_sets(self) -> None:
        signature_path = ROOT / "resources" / "CAR_T_state_signatures.csv"
        self.assertEqual(
            hashlib.sha256(signature_path.read_bytes()).hexdigest(),
            "ecb0eeb90b1a7016f600e2c825eaae2508dae379f0b7c849dfcc10e72d8112d7",
        )
        rows = read_csv(signature_path)
        by_cluster: dict[str, list[str]] = {}
        for row in rows:
            cluster = row["cluster"]
            if not cluster.startswith("C"):
                cluster = f"C{cluster}"
            by_cluster.setdefault(cluster, []).append(row["gene"])
        self.assertEqual(set(by_cluster), {f"C{i}" for i in range(10)})
        for cluster, genes in by_cluster.items():
            self.assertEqual(len(genes), 20, cluster)
            self.assertEqual(len(genes), len(set(genes)), cluster)
        nonsignificant = {
            (f"C{row['cluster']}", row["gene"])
            for row in rows
            if float(row["p_val_adj"]) >= 0.05
        }
        self.assertEqual(nonsignificant, {("C0", "KLRC1"), ("C9", "CHI3L2")})
        helper_text = (ROOT / "R" / "00_bioinfo_helpers.R").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("is.na(.data$padj) | .data$padj < 0.05", helper_text)
        self.assertIn('cox.zph(split_fit, transform = "rank")', helper_text)

    def test_checkmate_source_qc(self) -> None:
        path = ROOT / "reference_results" / "source_qc.json"
        qc = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(qc["rna_samples"], 311)
        self.assertEqual(qc["arm_counts"], {"NIVOLUMAB": 181, "EVEROLIMUS": 130})
        self.assertEqual(
            qc["trial_by_arm"]["NIVOLUMAB"],
            {"CM-009": 16, "CM-010": 45, "CM-025": 120},
        )
        self.assertEqual(qc["OS_CNSR_counts"], {"1": 231, "0": 80})
        self.assertEqual(qc["PFS_CNSR_counts"], {"1": 276, "0": 35})

    def test_primary_c6_result_and_multiplicity(self) -> None:
        rows = read_csv(
            ROOT / "reference_results" / "survival_results_nivolumab_primary.csv"
        )
        selected = [
            row
            for row in rows
            if row["endpoint"] == "OS" and row["signature"] == "C6"
        ]
        self.assertTrue(
            all(
                math.isfinite(float(result["PH_test_p"]))
                and 0 <= float(result["PH_test_p"]) <= 1
                for result in rows
            )
        )
        self.assertEqual(len(selected), 1)
        row = selected[0]
        self.assertEqual(int(row["n"]), 181)
        self.assertEqual(int(row["events"]), 123)
        self.assertAlmostEqual(float(row["HR_high_vs_low"]), 1.50408766, places=7)
        self.assertAlmostEqual(float(row["cox_p"]), 0.0259058, places=7)
        self.assertAlmostEqual(float(row["cox_BH_10_states"]), 0.166234, places=6)
        self.assertAlmostEqual(float(row["logrank_p"]), 0.0246259, places=7)
        self.assertAlmostEqual(
            float(row["logrank_BH_10_states"]), 0.208407, places=6
        )
        self.assertGreater(float(row["cox_BH_10_states"]), 0.05)

    def test_experimental_analysis_tables(self) -> None:
        manifest_path = ROOT / "data" / "experimental" / "experimental_data_manifest.tsv"
        with manifest_path.open("r", encoding="utf-8", newline="") as handle:
            manifest = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(manifest), 18)
        source_rows = [row for row in manifest if row["record_type"] == "source_workbook"]
        output_rows = [
            row for row in manifest if row["record_type"] == "canonical_analysis_table"
        ]
        self.assertEqual(len(source_rows), 8)
        self.assertEqual(len(output_rows), 10)
        self.assertEqual(
            {row["record_id"] for row in source_rows},
            {"TNF", "IFNG", "TNF_IFNG", "T6_MATCHED", "WT", "KO1", "KO2", "TCR"},
        )
        for row in output_rows:
            path = ROOT / row["filename_or_path"]
            self.assertTrue(path.is_file(), path)
            self.assertEqual(path.stat().st_size, int(row["size_bytes"]), path)
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), row["sha256"])
            if path.suffix == ".gz":
                with gzip.open(path, "rb") as handle:
                    while handle.read(1024 * 1024):
                        pass

        significant_path = (
            ROOT
            / "data"
            / "experimental"
            / "hela_cytokine_significant_differential_expression.tsv.gz"
        )
        with gzip.open(significant_path, "rt", encoding="utf-8", newline="") as handle:
            significant = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(significant), 28186)
        self.assertEqual(
            {row["condition"] for row in significant}, {"TNF", "IFNG", "TNF_IFNG"}
        )
        self.assertTrue(all(float(row["adjusted_p_value"]) < 0.05 for row in significant))
        self.assertTrue(
            all(
                math.isclose(
                    float(row["wald_statistic_treatment_vs_untreated"]),
                    float(row["log2_fold_change_treatment_vs_untreated"])
                    / float(row["lfc_se"]),
                    rel_tol=1e-12,
                    abs_tol=1e-12,
                )
                for row in significant
            )
        )
        up_sets = {
            condition: {
                row["gene_symbol"]
                for row in significant
                if row["condition"] == condition
                and float(row["base_mean"]) >= 30
                and float(row["log2_fold_change_treatment_vs_untreated"]) > 1
                and float(row["adjusted_p_value"]) < 0.05
            }
            for condition in ("TNF", "IFNG", "TNF_IFNG")
        }
        self.assertEqual(
            {condition: len(genes) for condition, genes in up_sets.items()},
            {"TNF": 458, "IFNG": 1659, "TNF_IFNG": 1595},
        )
        self.assertEqual(len(up_sets["TNF"] & up_sets["IFNG"]), 203)
        self.assertEqual(len(up_sets["TNF"] & up_sets["TNF_IFNG"]), 266)
        self.assertEqual(len(up_sets["IFNG"] & up_sets["TNF_IFNG"]), 1184)
        self.assertEqual(len(set.intersection(*up_sets.values())), 187)

        expected_cells = {"WT": 3743, "KO1": 2475, "KO2": 1831, "TCR": 5662}
        for sample, expected_n in expected_cells.items():
            path = (
                ROOT / "data" / "experimental" / "singlecell"
                / f"{sample}_targeted_counts.tsv.gz"
            )
            with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
                reader = csv.reader(handle, delimiter="\t")
                header = next(reader)
                self.assertEqual(header[0], "cell_index")
                self.assertEqual(len(header) - 1, 259)
                self.assertEqual(sum(1 for _ in reader), expected_n)

        qc = json.loads(
            (
                ROOT / "data" / "experimental"
                / "experimental_preparation_qc.json"
            ).read_text(encoding="utf-8")
        )
        cytokine_qc = qc["cytokine_differential_expression"]
        self.assertEqual(
            {condition: values["rows"] for condition, values in cytokine_qc.items()},
            {"TNF": 6632, "IFNG": 10672, "TNF_IFNG": 10882},
        )
        for values in cytokine_qc.values():
            self.assertEqual(values["source_stat_sign_reversed_rows"], values["rows"])
            self.assertEqual(values["rows_with_adjusted_p_at_least_0_05"], 0)
        self.assertEqual(
            qc["t6_matched_design"]["all_sheet"],
            {
                "matches_ordered_condition_union": True,
                "rows": 185700,
                "source_row_digest_sha256":
                    "b59ed0d497265f6802051221a9d482d1426744e7d400b03a8d80b2b84a58a85f",
            },
        )
        self.assertEqual(qc["single_cell"]["TCR"]["cells"], 5662)
        self.assertEqual(qc["single_cell"]["TCR"]["library_size_median"], 5625.0)
        self.assertEqual(qc["single_cell"]["TCR"]["detected_features_median"], 103.0)

    def test_t6_tables_preserve_valid_deseq2_missingness(self) -> None:
        expected_patterns = {
            "untreated": {
                (False, False, False, False, False): 14100,
                (True, True, True, True, True): 21898,
                (False, False, False, False, True): 10416,
                (False, False, False, True, True): 11,
            },
            "ifn": {
                (False, False, False, False, False): 19628,
                (True, True, True, True, True): 21169,
                (False, False, False, False, True): 5617,
                (False, False, False, True, True): 11,
            },
            "tnf": {
                (False, False, False, False, False): 18059,
                (True, True, True, True, True): 22060,
                (False, False, False, False, True): 6298,
                (False, False, False, True, True): 8,
            },
            "ti": {
                (False, False, False, False, False): 18484,
                (True, True, True, True, True): 21470,
                (False, False, False, False, True): 6467,
                (False, False, False, True, True): 4,
            },
        }
        columns = (
            "log2_fold_change_t6_vs_wt",
            "lfc_se",
            "wald_statistic_t6_vs_wt",
            "p_value",
            "adjusted_p_value",
        )
        for condition, expected in expected_patterns.items():
            path = (
                ROOT / "data" / "experimental"
                / f"hela_t6_vs_wt_{condition}_differential_expression.tsv.gz"
            )
            observed: dict[tuple[bool, ...], int] = {}
            with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
                rows = csv.DictReader(handle, delimiter="\t")
                for row in rows:
                    pattern = tuple(row[column] == "" for column in columns)
                    observed[pattern] = observed.get(pattern, 0) + 1
            self.assertEqual(observed, expected, condition)
            for pattern in observed:
                lfc_missing, se_missing, stat_missing, p_missing, padj_missing = pattern
                self.assertEqual(
                    (lfc_missing, se_missing, stat_missing),
                    (lfc_missing,) * 3,
                )
                self.assertFalse((not padj_missing) and p_missing)
                self.assertFalse((not p_missing) and lfc_missing)

        down_sets: dict[str, set[str]] = {}
        for condition in ("tnf", "ifn", "ti"):
            path = (
                ROOT / "data" / "experimental"
                / f"hela_t6_vs_wt_{condition}_differential_expression.tsv.gz"
            )
            with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
                rows = csv.DictReader(handle, delimiter="\t")
                down_sets[condition] = {
                    row["gene_symbol"]
                    for row in rows
                    if float(row["base_mean"]) >= 30
                    and row["log2_fold_change_t6_vs_wt"] != ""
                    and row["adjusted_p_value"] != ""
                    and float(row["log2_fold_change_t6_vs_wt"]) < -1
                    and float(row["adjusted_p_value"]) < 0.05
                }
        self.assertEqual(
            {condition: len(genes) for condition, genes in down_sets.items()},
            {"tnf": 446, "ifn": 1222, "ti": 845},
        )
        self.assertEqual(len(down_sets["tnf"] & down_sets["ifn"]), 85)
        self.assertEqual(len(down_sets["tnf"] & down_sets["ti"]), 165)
        self.assertEqual(len(down_sets["ifn"] & down_sets["ti"]), 225)
        self.assertEqual(len(set.intersection(*down_sets.values())), 61)

    def test_figure_5g_family_has_thirty_models(self) -> None:
        rows = read_csv(
            ROOT / "reference_results" / "Figure_5G_nivolumab_associations.csv"
        )
        self.assertEqual(len(rows), 30)
        self.assertEqual({row["cluster"] for row in rows}, {f"C{i}" for i in range(10)})
        self.assertEqual(
            {row["predictor"] for row in rows},
            {"ICAM1_z", "B2M_z", "combined_B2M_ICAM1_z"},
        )
        self.assertTrue(all(int(row["n"]) == 181 for row in rows))
        self.assertTrue(
            all(0 <= float(row["BH_p_30_models"]) <= 1 for row in rows)
        )
        c6_combined = [
            row for row in rows
            if row["cluster"] == "C6"
            and row["predictor"] == "combined_B2M_ICAM1_z"
        ]
        self.assertEqual(len(c6_combined), 1)
        self.assertAlmostEqual(
            float(c6_combined[0]["beta"]), 0.4357724025102622, places=12
        )
        self.assertAlmostEqual(
            float(c6_combined[0]["BH_p_30_models"]),
            1.1354309530985604e-10,
            places=20,
        )

    def test_committed_aggregate_model_contract_when_present(self) -> None:
        path = ROOT / "data" / "analysis" / "checkmate_c6_global_gene_models.tsv.gz"
        if not path.exists():
            self.skipTest("aggregate CheckMate model table has not been generated")
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(
            hashlib.sha256(path.read_bytes()).hexdigest(),
            "21bf8a3c3af290e97190a12dd48730f0fec9ca6fb65fc74a27cb7c4827e449c6",
        )
        self.assertGreaterEqual(len(rows), 1000)
        self.assertEqual(len(rows), len({row["gene"] for row in rows}))
        for row in rows:
            self.assertGreaterEqual(int(row["n"]), 163)
            self.assertLessEqual(int(row["n"]), 181)
            self.assertTrue(math.isfinite(float(row["beta"])))
            self.assertGreaterEqual(float(row["p"]), 0)
            self.assertLessEqual(float(row["p"]), 1)

    def test_reference_figure_checksums(self) -> None:
        checksum_file = ROOT / "reference_figures" / "SHA256SUMS"
        for line in checksum_file.read_text(encoding="utf-8").splitlines():
            expected, relative = line.split(maxsplit=1)
            path = ROOT / "reference_figures" / relative
            observed = hashlib.sha256(path.read_bytes()).hexdigest()
            self.assertEqual(observed, expected, relative)


if __name__ == "__main__":
    unittest.main()
