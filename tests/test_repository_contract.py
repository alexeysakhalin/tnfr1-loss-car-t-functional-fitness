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
        for row in rows:
            self.assertEqual(len(row["sha256"]), 64)
            int(row["sha256"], 16)
            self.assertGreater(int(row["size_bytes"]), 0)
            self.assertIn(
                row["repository_policy"],
                {"local_only", "metadata_only", "mapping_resource_committed"},
            )

    def test_signature_resource_has_ten_fixed_twenty_gene_sets(self) -> None:
        rows = read_csv(ROOT / "resources" / "CAR_T_state_signatures.csv")
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
        helper_text = (ROOT / "R" / "00_bioinfo_helpers.R").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("is.na(.data$padj) | .data$padj < 0.05", helper_text)

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
