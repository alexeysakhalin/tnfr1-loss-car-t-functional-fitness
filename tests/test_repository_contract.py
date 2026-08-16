from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
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
            "data/analysis/depmap_s1b_eligible_models.tsv.gz",
            "data/analysis/depmap_s1b_preparation_qc.json",
            "data/analysis/depmap_s1b_source_provenance.json",
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

    def test_author_confirmed_experimental_aliases(self) -> None:
        path = (
            ROOT / "data" / "experimental" / "experimental_sample_aliases.tsv"
        )
        with path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(
            rows,
            [
                {
                    "source_identifier": "T6",
                    "manuscript_identifier": "TNFR1-KO1",
                    "confirmation_basis": "author confirmation",
                    "confirmation_date": "2026-08-14",
                }
            ],
        )
        figure_script = (ROOT / "R" / "04_figure_2_B_suppl_S2D.R").read_text(
            encoding="utf-8"
        )
        self.assertIn("TNFR1-KO1 versus WT", figure_script)
        self.assertNotIn("Use `T6` until", (ROOT / "docs" / "FIGURE_MAP.md").read_text(
            encoding="utf-8"
        ))

    def test_depmap_expression_metadata_contract_is_checksum_pinned(self) -> None:
        render_script = (ROOT / "R" / "11_supplementary_1B.R").read_text(
            encoding="utf-8"
        )
        self.assertIn("depmap_s1b_eligible_models.tsv.gz", render_script)
        self.assertIn("depmap_s1b_preparation_qc.json", render_script)
        self.assertIn("depmap_s1b_source_provenance.json", render_script)
        self.assertIn("depmap_s1b_statistics.csv", render_script)
        self.assertNotIn("DEPMAP_EXPRESSION_RELEASE", render_script)
        self.assertNotIn("DEPMAP_MODEL_RELEASE", render_script)
        self.assertIn('release_pair_status <- "unverified"', render_script)
        self.assertIn('expression_release <- "DepMap Public 25Q2"', render_script)
        self.assertIn('model_release_identity_status <- "unverified"', render_script)
        self.assertIn('"tissue_origin_filter_applied"', render_script)
        self.assertIn("library(data.table)", render_script)
        self.assertIn("library(ggplot2)", render_script)
        for package in ("dplyr", "readr", "tidyr", "digest", "jsonlite"):
            self.assertNotIn(f"library({package})", render_script)
        self.assertIn("Supplementary_Figure_S1B.png", render_script)
        self.assertIn("Supplementary_Figure_S1B.tiff", render_script)
        self.assertNotIn("Supplementary_Figure_S1B_strip", render_script)
        self.assertIn("n = 1,591", render_script)
        self.assertIn("n = 749 (47.1%)", render_script)
        self.assertIn("n = 254 (16.0%)", render_script)
        self.assertIn("n = 423 (26.6%)", render_script)
        self.assertIn("n = 165 (10.4%)", render_script)
        self.assertIn("RIPK3_below_threshold = 1003L", render_script)
        self.assertIn("NLRP3_below_threshold = 1172L", render_script)
        self.assertIn("both_below_threshold = 749L", render_script)

        prep_script = (ROOT / "scripts" / "prepare_depmap_s1b.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("mtime=0", prep_script)
        self.assertIn('"release_pair_status": "unverified"', prep_script)
        self.assertIn('"expression_release": "DepMap Public 25Q2"', prep_script)
        self.assertIn('"model_release_identity_status": "unverified"', prep_script)
        self.assertIn('"same_release_pair": None', prep_script)
        self.assertIn("EXPECTED_EXPRESSION_SIZE = 521_526_607", prep_script)
        self.assertIn("e0326e16eb23bea1be980fce315acb36b224dedd7af6b47e0ba37e7747dbcc47", prep_script)
        self.assertIn("EXPECTED_MODEL_SIZE = 694_278", prep_script)
        self.assertIn("b096e03bfefdc2679211545ddbf1bb7878d69ffde07ae335af5b968a7883733c", prep_script)
        self.assertIn(
            "90bfdbe5c44cbb8f822e655ba7f179f3033933116285b6b2f85153b2d3d17c75",
            prep_script,
        )
        self.assertIn(
            "9dbb9de8805696c1345816ab07edd23fb4fd95e117739f3c5c3b1cf062c1233b",
            prep_script,
        )
        self.assertIn("af4472ab734ea3aec974d992b504c7e5", prep_script)

        derived_path = ROOT / "data" / "analysis" / (
            "depmap_s1b_eligible_models.tsv.gz"
        )
        self.assertEqual(derived_path.stat().st_size, 47514)
        self.assertEqual(
            hashlib.sha256(derived_path.read_bytes()).hexdigest(),
            "368ad92b085a722d3984a5355bea3109d8e5a2b29ffe563c3fd284cf8970f354",
        )
        with gzip.open(derived_path, "rt", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 1591)
        self.assertEqual(len({row["ModelID"] for row in rows}), 1591)
        self.assertEqual(len({row["ProfileID"] for row in rows}), 1591)
        self.assertFalse(
            [
                row for row in rows
                if not row["OncotreePrimaryDisease"]
                or row["OncotreePrimaryDisease"].upper() == "NON-CANCEROUS"
            ]
        )
        ripk3_low = [row["RIPK3_below_threshold"] == "TRUE" for row in rows]
        nlrp3_low = [row["NLRP3_below_threshold"] == "TRUE" for row in rows]
        for row, ripk3_flag, nlrp3_flag in zip(rows, ripk3_low, nlrp3_low):
            self.assertEqual(
                ripk3_flag, float(row["RIPK3_log2_TPM_plus_1"]) < 0.5
            )
            self.assertEqual(
                nlrp3_flag, float(row["NLRP3_log2_TPM_plus_1"]) < 0.5
            )
            expected_category = (
                "Both below threshold" if ripk3_flag and nlrp3_flag
                else "RIPK3 below threshold only" if ripk3_flag
                else "NLRP3 below threshold only" if nlrp3_flag
                else "Neither below threshold"
            )
            self.assertEqual(row["threshold_category"], expected_category)
        self.assertEqual(sum(ripk3_low), 1003)
        self.assertEqual(sum(nlrp3_low), 1172)
        self.assertEqual(sum(a and b for a, b in zip(ripk3_low, nlrp3_low)), 749)
        self.assertEqual(
            sum(a and not b for a, b in zip(ripk3_low, nlrp3_low)), 254
        )
        self.assertEqual(
            sum(not a and b for a, b in zip(ripk3_low, nlrp3_low)), 423
        )
        self.assertEqual(
            sum(not a and not b for a, b in zip(ripk3_low, nlrp3_low)), 165
        )

        qc_path = ROOT / "data" / "analysis" / "depmap_s1b_preparation_qc.json"
        provenance_path = (
            ROOT / "data" / "analysis" / "depmap_s1b_source_provenance.json"
        )
        self.assertEqual(qc_path.stat().st_size, 1583)
        self.assertEqual(
            hashlib.sha256(qc_path.read_bytes()).hexdigest(),
            "8d58a7113ca08c7ffd8297e26f3a3a29693e61e119186365c1fb04531d988b79",
        )
        self.assertEqual(provenance_path.stat().st_size, 1895)
        self.assertEqual(
            hashlib.sha256(provenance_path.read_bytes()).hexdigest(),
            "e21b63f90835e80e71c388b13e167b214773cd6a988aa43a3aaac8cd65745242",
        )
        qc = json.loads(qc_path.read_text(encoding="utf-8"))
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        self.assertEqual(qc["join_qc"]["eligible_models"], 1591)
        self.assertEqual(qc["derived_file"]["sha256"], hashlib.sha256(
            derived_path.read_bytes()
        ).hexdigest())
        self.assertEqual(provenance["release_pair_status"], "unverified")
        self.assertEqual(provenance["expression_release"], "DepMap Public 25Q2")
        self.assertEqual(
            provenance["model_release_identity_status"], "unverified"
        )
        self.assertIsNone(provenance["same_release_pair"])
        self.assertFalse(provenance["tissue_origin_filter_applied"])

        summary_path = ROOT / "reference_results" / "depmap_s1b_statistics.csv"
        self.assertEqual(summary_path.stat().st_size, 305)
        self.assertEqual(
            hashlib.sha256(summary_path.read_bytes()).hexdigest(),
            "60585de1dd22e879220d7d9da89d6cd762685b1b11af1d60b64627490e80e990",
        )
        summary_rows = read_csv(summary_path)
        self.assertEqual(
            [(row["metric"], row["n"], row["percent"]) for row in summary_rows],
            [
                ("RIPK3 below threshold", "1003", "63.0"),
                ("NLRP3 below threshold", "1172", "73.7"),
                ("Both below threshold", "749", "47.1"),
                ("RIPK3 below threshold only", "254", "16.0"),
                ("NLRP3 below threshold only", "423", "26.6"),
                ("Neither below threshold", "165", "10.4"),
            ],
        )

        documentation = (ROOT / "docs" / "DEPMAP_S1B.md").read_text(
            encoding="utf-8"
        )
        self.assertIn("one default RNA-expression profile", documentation)
        self.assertIn("RIPK3 below threshold | 1,003 | 63.0%", documentation)
        self.assertIn("NLRP3 below threshold | 1,172 | 73.7%", documentation)
        self.assertIn("Both below threshold | 749 | 47.1%", documentation)
        self.assertIn("denominator must not", documentation)
        self.assertIn('be called "human"', documentation)
        self.assertIn("Do not silently label both files", documentation)

        manifest_path = ROOT / "data" / "source_manifest.tsv"
        with manifest_path.open("r", encoding="utf-8", newline="") as handle:
            manifest = {
                row["source_id"]: row
                for row in csv.DictReader(handle, delimiter="\t")
            }
        archive = manifest["depmap_expression_archive_supplied"]
        member = manifest["depmap_expression_csv_member_supplied"]
        model = manifest["depmap_model_metadata_supplied"]
        self.assertEqual(archive["cohort_id"], "DEPMAP_PUBLIC_25Q2")
        self.assertEqual(member["cohort_id"], "DEPMAP_PUBLIC_25Q2")
        self.assertEqual(model["cohort_id"], "DEPMAP_SOURCE_PAIR")
        self.assertEqual(archive["size_bytes"], "249426032")
        self.assertEqual(
            archive["sha256"],
            "c44524c48e20f8c5c1263eb23cd55df77ceda62cfb5246babbe22cecc90c3da0",
        )
        self.assertEqual(member["size_bytes"], "538420733")
        self.assertEqual(
            member["sha256"],
            "90bfdbe5c44cbb8f822e655ba7f179f3033933116285b6b2f85153b2d3d17c75",
        )
        self.assertEqual(model["size_bytes"], "699474")
        self.assertEqual(
            model["sha256"],
            "9dbb9de8805696c1345816ab07edd23fb4fd95e117739f3c5c3b1cf062c1233b",
        )
        self.assertIn("MD5 af4472ab734ea3aec974d992b504c7e5", model["notes"])

    def test_depmap_s1b_ci_checks_publication_outputs(self) -> None:
        workflow = (
            ROOT / ".github" / "workflows" / "depmap-s1b.yml"
        ).read_text(encoding="utf-8")
        verifier = (
            ROOT / "scripts" / "verify_depmap_s1b_outputs.py"
        ).read_text(encoding="utf-8")
        self.assertIn('r-version: "4.4.3"', workflow)
        self.assertIn("setup-renv@v2", workflow)
        self.assertIn("R/11_supplementary_1B.R", workflow)
        self.assertIn("verify_depmap_s1b_outputs.py", workflow)
        self.assertIn("depmap-supplementary-figure-s1b", workflow)
        self.assertIn("(4260, 3840)", verifier)
        self.assertIn('"release_pair_status": "unverified"', verifier)
        self.assertIn('"n_models": "1591"', verifier)
        self.assertIn("verify_image_pair", verifier)

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

    def test_figure_5e_s6d_s6e_selection_is_unambiguous(self) -> None:
        script = (
            ROOT / "R" / "08_figure_5E_S6D_S6E.R"
        ).read_text(encoding="utf-8")
        normalized = " ".join(script.split())
        endpoint_selector = ".data$endpoint == .env$endpoint"
        signature_selector = ".data$signature == .env$cluster"
        guard = "if (nrow(r) != 1L)"
        subtitle = "subtitle <- sprintf("
        self.assertIn(endpoint_selector, normalized)
        self.assertIn(signature_selector, normalized)
        self.assertNotIn(".data$endpoint == endpoint", normalized)
        self.assertNotIn(".data$signature == cluster", normalized)
        self.assertIn(guard, normalized)
        self.assertLess(normalized.index(endpoint_selector), normalized.index(guard))
        self.assertLess(normalized.index(signature_selector), normalized.index(guard))
        self.assertLess(normalized.index(guard), normalized.index(subtitle))

    def test_open_cohort_numeric_parser_is_strict(self) -> None:
        path = ROOT / "scripts" / "prepare_open_cohort_analysis_tables.py"
        spec = importlib.util.spec_from_file_location(
            "prepare_open_cohort_analysis_tables_for_test", path
        )
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        for value in (None, "", " ", "NA", "n/a", "NaN", "NE", float("nan")):
            with self.subTest(missing=value):
                self.assertTrue(math.isnan(module.as_float(value)))

        for value, expected in (
            (0, 0.0),
            (12345, 12345.0),
            (12.345, 12.345),
            ("12.345", 12.345),
            (" .5 ", 0.5),
            ("-2.5e+3", -2500.0),
        ):
            with self.subTest(valid=value):
                self.assertEqual(module.as_float(value), expected)

        for value in ("12,345", "1,234", "1,234.5", "0,5"):
            with self.subTest(ambiguous=value):
                with self.assertRaisesRegex(ValueError, "comma"):
                    module.as_float(value)

        for value in (
            "not-a-number",
            "1_000",
            "inf",
            "-Infinity",
            "1e9999",
            True,
            float("inf"),
        ):
            with self.subTest(invalid=value):
                with self.assertRaises(ValueError):
                    module.as_float(value)

    def test_r05_publishes_outputs_only_after_marker_guards(self) -> None:
        source = (ROOT / "R" / "05_figure_4_AB_suppl_S5A.R").read_text(
            encoding="utf-8"
        )
        self.assertIn('FINAL_FIG_DIR <- "figures"', source)
        self.assertIn('FIG_DIR <- tempfile(pattern = "R05-staging-"', source)
        self.assertNotIn('\nFIG_DIR <- "figures"\n', source)
        invocation = "promote_staged_outputs(FIG_DIR, FINAL_FIG_DIR)"
        self.assertEqual(source.count(invocation), 1)
        last_guard = source.rfind(
            'stop("TCR marker table is empty; cluster annotations cannot be '
            'release-validated.")'
        )
        promotion = source.rfind(invocation)
        self.assertGreater(last_guard, -1)
        self.assertGreater(promotion, last_guard)
        self.assertGreater(promotion, source.rfind("saveRDS(obj"))

    def test_figure_5f_displays_the_adjusted_gene_family(self) -> None:
        source = (ROOT / "R" / "09_figure_5F_S6G.R").read_text(encoding="utf-8")
        normalized = " ".join(source.split())
        self.assertIn(
            "neglog10_BH_p = -log10(pmax(.data$BH_p, "
            ".Machine$double.xmin))",
            normalized,
        )
        self.assertIn("aes(.data$beta, .data$neglog10_BH_p)", normalized)
        self.assertIn('expression(-log[10]("BH-adjusted p"))', source)
        self.assertNotIn('expression(-log[10]("nominal p"))', source)

    def test_bulk_p_value_floor_uses_zero_absolute_tolerance(self) -> None:
        source = (
            ROOT / "scripts" / "verify_bulk_figure_outputs.py"
        ).read_text(encoding="utf-8")
        normalized = " ".join(source.split())
        self.assertIn(
            'abs_tol=0.0 if field == "adjusted_p_value_floor" else 1e-12',
            normalized,
        )

    def test_scientific_workflows_validate_relevant_pushes_to_main(self) -> None:
        for relative_path in (
            ".github/workflows/bulk-rnaseq.yml",
            ".github/workflows/depmap-s1b.yml",
        ):
            with self.subTest(workflow=relative_path):
                source = (ROOT / relative_path).read_text(encoding="utf-8")
                self.assertIn(
                    "  push:\n    branches: [main]\n    paths:\n",
                    source,
                )

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
