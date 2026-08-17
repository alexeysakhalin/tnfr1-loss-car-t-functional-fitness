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
    def test_experimental_preparer_accepts_nested_deposit_layout(self) -> None:
        source = (
            ROOT / "scripts" / "prepare_experimental_analysis_tables.py"
        ).read_text(encoding="utf-8")
        self.assertIn("input_dir.rglob(pattern)", source)
        self.assertIn("nested bulk_rnaseq/ and targeted_single_cell/", source)

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

    def test_imvigor_exports_are_reproducible_local_inputs(self) -> None:
        manifest_path = ROOT / "data" / "source_manifest.tsv"
        with manifest_path.open("r", encoding="utf-8", newline="") as handle:
            manifest = {
                row["source_id"]: row
                for row in csv.DictReader(handle, delimiter="\t")
            }

        package = manifest["imvigor210_processed_package"]
        clinical = manifest["imvigor210_clinical_export"]
        expression = manifest["imvigor210_expression_export"]
        self.assertEqual(package["expected_filename"], (
            "IMvigor210CoreBiologies_1.0.0.tar.gz"
        ))
        self.assertEqual(package["size_bytes"], "122127298")
        self.assertEqual(package["sha256"], (
            "cfdd3176d7b34de5b04fb9416bfd2b20fa4b6e238aaad5f20b048a34329ea178"
        ))
        self.assertEqual(clinical["fetch_mode"], "local_input")
        self.assertEqual(expression["fetch_mode"], "local_input")
        self.assertEqual(clinical["repository_policy"], "local_only")
        self.assertEqual(expression["repository_policy"], "local_only")
        self.assertEqual(clinical["source_url"], "")
        self.assertEqual(expression["source_url"], "")
        for row in (package, clinical, expression):
            self.assertIn("scripts/export_imvigor210_inputs.R", row["notes"])

        exporter_path = ROOT / "scripts" / "export_imvigor210_inputs.R"
        exporter = exporter_path.read_text(encoding="utf-8")
        self.assertTrue(exporter.startswith("#!/usr/bin/env Rscript\n"))
        self.assertIn("DESeq::counts(cds)", exporter)
        self.assertIn("Biobase::pData(cds)", exporter)
        self.assertIn(
            "expression_log2cpm <- log2(edgeR::cpm(count_matrix, log = FALSE) + 1)",
            exporter,
        )
        self.assertIn("utils::write.csv(", exporter)
        self.assertIn('options(scipen = 0, OutDec = ".")', exporter)
        self.assertIn('file(path, open = "wb")', exporter)
        self.assertIn("--package-tarball", exporter)
        self.assertIn("--verify-only", exporter)
        self.assertIn("--diagnostics-path", exporter)
        self.assertIn("--semantic-report-path", exporter)
        self.assertIn("--external-semantic-file-verification", exporter)
        self.assertIn(clinical["sha256"], exporter)
        self.assertIn(expression["sha256"], exporter)
        self.assertIn(package["sha256"], exporter)
        self.assertIn(
            'sweep(count_matrix, 2L, library_sizes, "/") * 1e6 + 1',
            exporter,
        )
        self.assertIn("ALL_CELL_ABSOLUTE_TOLERANCE <- 5e-13", exporter)
        self.assertIn("ALL_CELL_RELATIVE_TOLERANCE <- 5e-14", exporter)
        self.assertIn("compare_expression_matrices", exporter)
        self.assertIn("read_expression_csv", exporter)
        self.assertIn("verify_expression_file_semantically", exporter)
        self.assertIn("IMvigor210_expression_semantic_contract_v1.json", exporter)
        self.assertNotIn("EXPECTED_EXPRESSION_SEMANTICS", exporter)
        self.assertNotIn("expression_values_8dp_sha256", exporter)
        self.assertNotIn(
            'verify_file(\n  staged_paths[["expression"]], EXPECTED_EXPORTS$expression',
            exporter,
        )
        self.assertLess(
            exporter.index("observed_sha256 <- sha256_file(path)"),
            exporter.index(
                "if (is.na(observed_size) || "
                "observed_size != specification$size_bytes)"
            ),
        )
        self.assertLess(
            exporter.index("verify_file(package_tarball"),
            exporter.index("utils::untar(package_tarball"),
        )
        self.assertLess(
            exporter.index('verify_file(staged_paths[["clinical"]]'),
            exporter.index("publish_verified_exports(staged_paths"),
        )

        contract_path = (
            ROOT / "resources" / "IMvigor210_expression_semantic_contract_v1.json"
        )
        contract = json.loads(contract_path.read_text(encoding="utf-8"))
        self.assertEqual(contract["semantic_contract_version"], 1)
        self.assertEqual(contract["semantic_digests"]["required_scale"], 6)
        self.assertEqual(contract["semantic_digests"]["diagnostic_scales"], [7, 8])
        self.assertTrue(contract["analysis_canonicalization"]["required"])
        self.assertEqual(contract["analysis_canonicalization"]["scale"], 6)
        self.assertEqual(
            contract["analysis_canonicalization"]["rounding"], "ROUND_HALF_UP"
        )
        self.assertFalse(
            contract["privacy"]["contract_contains_identifiers_or_expression_cells"]
        )
        preparer = (
            ROOT / "scripts" / "prepare_open_cohort_analysis_tables.py"
        ).read_text(encoding="utf-8")
        self.assertIn("semantically_verified_imvigor_expression", preparer)
        self.assertIn("verifier.verify_expression", preparer)
        self.assertIn("canonical_imvigor_fixed6", preparer)
        self.assertIn('"expression_semantic_scale": IMVIGOR_EXPRESSION_SCALE', preparer)
        imvigor_function = preparer.split("def prepare_imvigor(", 1)[1].split(
            "\ndef workbook_header", 1
        )[0]
        self.assertNotIn(
            'verified_input(raw_dir, manifest, "imvigor210_expression_export")',
            imvigor_function,
        )
        self.assertIn(
            "[canonical_imvigor_fixed6(value) for value in row[1:]]",
            imvigor_function,
        )
        self.assertNotIn("[as_float(value) for value in row[1:]]", imvigor_function)

        data_readme = (ROOT / "data" / "README.md").read_text(encoding="utf-8")
        input_inventory = (
            ROOT / "docs" / "REPRODUCIBILITY_INPUTS.md"
        ).read_text(encoding="utf-8")
        root_readme = (ROOT / "README.md").read_text(encoding="utf-8")
        data_license = (ROOT / "DATA_LICENSE.md").read_text(encoding="utf-8")
        self.assertNotIn("should be deposited", data_readme)
        for documentation in (
            data_readme, input_inventory, root_readme, data_license
        ):
            self.assertIn("scripts/export_imvigor210_inputs.R", documentation)
            self.assertIn("Zenodo", documentation)
            self.assertIn("local-only", documentation)
        for documentation in (data_readme, input_inventory, root_readme):
            self.assertIn("--source imvigor210_processed_package", documentation)
            self.assertIn("--accept-licensed-public-downloads", documentation)
        validation_workflow = (
            ROOT / ".github" / "workflows" / "validate.yml"
        ).read_text(encoding="utf-8")
        self.assertIn(
            'list.files("scripts", pattern = "[.]R$", full.names = TRUE)',
            validation_workflow,
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
        self.assertIn('release_pair_status <- "confirmed"', render_script)
        self.assertIn('expression_release <- "DepMap Public 25Q2"', render_script)
        self.assertIn('model_release <- "DepMap Public 25Q2"', render_script)
        self.assertIn('model_release_identity_status <- "confirmed"', render_script)
        self.assertIn('same_release_pair <- "TRUE"', render_script)
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
        self.assertIn('"release_pair_status": "confirmed"', prep_script)
        self.assertIn('"expression_release": "DepMap Public 25Q2"', prep_script)
        self.assertIn('"model_release": "DepMap Public 25Q2"', prep_script)
        self.assertIn('"model_release_identity_status": "confirmed"', prep_script)
        self.assertIn('"same_release_pair": True', prep_script)
        self.assertIn("EXPECTED_EXPRESSION_ROWS = 1_739", prep_script)
        self.assertIn("EXPECTED_DEFAULT_ROWS = 1_684", prep_script)
        self.assertIn("EXPECTED_NONDEFAULT_ROWS = 55", prep_script)
        self.assertIn("EXPECTED_MODEL_ROWS = 2_132", prep_script)
        self.assertIn("EXPECTED_ELIGIBLE_MODELS = 1_591", prep_script)
        self.assertIn(
            "c44524c48e20f8c5c1263eb23cd55df77ceda62cfb5246babbe22cecc90c3da0",
            prep_script,
        )
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
        self.assertEqual(provenance_path.stat().st_size, 1984)
        self.assertEqual(
            hashlib.sha256(provenance_path.read_bytes()).hexdigest(),
            "95f8c8f11fbb43b1bf093811d110bdd1c5b5d53ef7771c797c165b1061e14816",
        )
        qc = json.loads(qc_path.read_text(encoding="utf-8"))
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        self.assertEqual(qc["join_qc"]["eligible_models"], 1591)
        self.assertEqual(qc["derived_file"]["sha256"], hashlib.sha256(
            derived_path.read_bytes()
        ).hexdigest())
        self.assertEqual(provenance["release_pair_status"], "confirmed")
        self.assertEqual(provenance["expression_release"], "DepMap Public 25Q2")
        self.assertEqual(provenance["model_release"], "DepMap Public 25Q2")
        self.assertEqual(
            provenance["model_release_identity_status"], "confirmed"
        )
        self.assertIs(provenance["same_release_pair"], True)
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
        self.assertIn("confirmed same-release", documentation)
        self.assertIn("release-specific Figshare DOI", documentation)

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
        self.assertEqual(model["cohort_id"], "DEPMAP_PUBLIC_25Q2")
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
        self.assertIn("downloaded from the DepMap Portal All Data page", model["notes"])
        self.assertIn("DepMap Public 25Q2", model["citation"])
        self.assertIn("releasename=DepMap%20Public%2025Q2", model["source_url"])

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
        self.assertIn('"release_pair_status": "confirmed"', verifier)
        self.assertIn('"model_release": "DepMap Public 25Q2"', verifier)
        self.assertIn('"same_release_pair": "TRUE"', verifier)
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

    def test_r05_reviewed_signature_concordance_contract_is_exact(self) -> None:
        contract_path = (
            ROOT / "resources" / "CAR_T_state_signature_concordance_v1.csv"
        )
        self.assertEqual(
            hashlib.sha256(contract_path.read_bytes()).hexdigest(),
            "cb41afce9dc85a74a8c027ef764001e0179c144e1debee4cf047f1f49321badc",
        )
        rows = read_csv(contract_path)
        self.assertEqual([row["cluster"] for row in rows], [f"C{i}" for i in range(10)])
        self.assertEqual(
            {
                (
                    row["contract_version"], row["r_version"],
                    row["seurat_version"], row["matrix_version"],
                )
                for row in rows
            },
            {("1", "4.4.3", "5.5.1", "1.7.2")},
        )
        expected = {
            "C0": (20, "", ""),
            "C1": (18, "FOXP3;TSPAN32", "CD7;DUSP1"),
            "C2": (17, "CD70;ICOS;IL2RA", "DUSP2;F5;STAT5A"),
            "C3": (19, "C10ORF54", "IL9R"),
            "C4": (20, "", ""),
            "C5": (19, "GZMB", "FASLG"),
            "C6": (19, "CD4", "HLA-DQA1"),
            "C7": (18, "ARL4C;IL18RAP", "LIF;TRAT1"),
            "C8": (20, "", ""),
            "C9": (20, "", ""),
        }
        observed = {
            row["cluster"]: (
                int(row["n_overlap"]), row["frozen_only"], row["current_only"]
            )
            for row in rows
        }
        self.assertEqual(observed, expected)
        for row in rows:
            self.assertEqual(int(row["n_frozen"]), 20)
            self.assertEqual(int(row["n_current"]), 20)
            self.assertEqual(
                len([gene for gene in row["frozen_only"].split(";") if gene]),
                20 - int(row["n_overlap"]),
            )
            self.assertEqual(
                len([gene for gene in row["current_only"].split(";") if gene]),
                20 - int(row["n_overlap"]),
            )

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
        self.assertIn(
            'DIAGNOSTIC_DIR <- file.path("results", "targeted-singlecell-diagnostics")',
            source,
        )
        self.assertIn("Supplementary_Table_S5_signature_membership_differences.tsv", source)
        self.assertIn("CAR_T_state_signature_concordance_v1.csv", source)
        self.assertIn("reviewed_contract_match", source)
        self.assertIn("Supplementary_Table_S5_reviewed_concordance_contract.tsv", source)
        self.assertIn("Supplementary_Table_S5_C10_markers.tsv", source)
        self.assertIn("Supplementary_Table_S5_cell_counts_all_QC_including_C10.csv", source)
        self.assertIn("Supplementary_Table_S5_cell_counts_by_cluster.csv", source)
        self.assertIn("Supplementary_Table_S5_cell_filtering_QC.tsv", source)
        diagnostic_export = source.index(
            "# Keep the aggregate C10 and QC evidence in a stable diagnostic directory"
        )
        reviewed_guard = source.index(
            "if (any(is.na(signature_concordance$reviewed_contract_match))"
        )
        self.assertLess(diagnostic_export, reviewed_guard)
        self.assertLess(
            source.index("Supplementary_Table_S5_signature_membership_differences.tsv"),
            source.index(
                "Current C0-C9 top-20 marker differences do not match the exact"
            ),
        )
        self.assertIn(
            '"C10" = "C10 small cytokine/IFN-response-high cluster"',
            source,
        )
        self.assertIn("p_clusters_clean <- DimPlot(\n  obj,", source)
        self.assertNotIn("DimPlot(\n  obj_signature_reference,", source)
        self.assertIn("# Descriptive Figure 4B summaries retain all QC-passing C0-C10 cells.", source)
        self.assertIn("md <- md_all", source)
        self.assertNotIn('dplyr::filter(cluster_short != "C10")', source)
        self.assertIn('y = "% of C0-C10 QC-passing cells"', source)
        self.assertIn("FindAllMarkers(\n  obj_signature_reference,", source)
        self.assertIn("expected_marker_clusters <- as.character(0:9)", source)
        self.assertIn('"C3" = "CD8/TRDC-associated cytotoxic state"', source)
        self.assertIn('"C4" = "Cycling T-cell state II"', source)
        self.assertIn(
            '"C0" = "Mixed CD4/KLRB1-associated activated state"',
            source,
        )
        self.assertIn(
            '"C2" = "Cytokine-expressing effector state"',
            source,
        )
        self.assertNotIn("TH9-like", source)
        self.assertNotIn("stem-like/early-memory", source)
        self.assertIn('title = "Cluster composition"', source)
        self.assertNotIn("Cycling/effector T-cell state II", source)
        self.assertIn("Exploratory_TCR_C6_signature_projection.png", source)
        self.assertIn(
            "Exploratory_TCR_C6_signature_projection_by_cluster.tsv",
            source,
        )
        self.assertIn(
            "Exploratory_TCR_C6_signature_projection_gene_coverage.tsv",
            source,
        )
        self.assertIn("c6_projection_by_cluster", source)
        self.assertIn("rank_auc_score <- function", source)
        self.assertIn('"cluster_short", "cluster_annotation", "n_cells"', source)
        self.assertIn('"cxcl13_detected_cells", "cxcl13_detection_fraction"', source)
        self.assertNotIn("CXCL13_detected_n", source)
        self.assertNotIn("CXCL13_detected_pct", source)
        self.assertNotIn("CXCL13_avg_log_normalized_expression", source)
        self.assertIn("coord_cartesian(ylim = c(0, 1))", source)
        self.assertIn('"replicate-level inference was performed.', source)
        self.assertIn("strwrap(", source)
        self.assertIn("width = 170", source)
        self.assertIn(
            'x = "Independent repeated-stimulation cluster"',
            source,
        )
        self.assertIn(
            'ggtitle("Repeated CD3/CD28-stimulation UMAP by cluster")',
            source,
        )
        self.assertNotIn("Independent TCR cluster", source)
        self.assertNotIn("across TCR clusters", source)
        self.assertIn(
            "This projection does not modify the repeated-stimulation clustering",
            source,
        )
        self.assertNotIn("openxlsx", source)
        self.assertEqual(source.count("writexl::write_xlsx("), 1)
        self.assertEqual(source.count("write_marker_workbook(\n"), 2)

        workflow = (
            ROOT / ".github" / "workflows" / "cohort-inputs-targeted-singlecell.yml"
        ).read_text(encoding="utf-8")
        self.assertIn("id: targeted_run", workflow)
        self.assertIn("steps.targeted_run.outcome", workflow)
        self.assertIn("results/targeted-singlecell-diagnostics", workflow)
        self.assertIn("any::writexl@2.0.0", workflow)
        self.assertIn('writexl = "2.0.0"', workflow)
        self.assertNotIn("any::openxlsx", workflow)
        self.assertEqual(workflow.count("python scripts/verify_xlsx_workbook.py"), 2)
        self.assertIn("openpyxl==3.1.5", workflow)
        self.assertEqual(workflow.count("--require-openpyxl"), 2)
        targeted_job = workflow.split("  targeted-singlecell:\n", 1)[1]
        self.assertIn("uses: actions/setup-python@v6", targeted_job)
        self.assertIn('python-version: "3.12.13"', targeted_job)
        self.assertEqual(
            workflow.count('"resources/CAR_T_state_signature_concordance_v1.csv"'),
            2,
        )
        self.assertGreaterEqual(workflow.count("if: always()"), 2)
        self.assertGreaterEqual(
            workflow.count("R_LIBS=/runner-temp/imvigor210-r-library"), 3
        )
        self.assertIn("Exploratory_TCR_C6_signature_projection", workflow)
        self.assertIn(
            "Exploratory_TCR_C6_signature_projection_by_cluster.tsv",
            workflow,
        )
        self.assertIn(
            "Exploratory_TCR_C6_signature_projection_gene_coverage.tsv",
            workflow,
        )
        self.assertIn("lib = library_path", workflow)
        self.assertIn("find.package(package, lib.loc = library_path)", workflow)
        install_step = workflow.split(
            "      - name: Install and assert the legacy top-level packages\n", 1
        )[1].split(
            "      - name: Recreate and scientifically verify both exports\n", 1
        )[0]
        self.assertIn("docker run --rm -i \\", install_step)
        self.assertIn("for package in DESeq Biobase edgeR; do", install_step)
        self.assertIn('package_dir="$IMVIGOR_R_LIBRARY/$package"', install_step)
        self.assertIn('test -d "$package_dir"', install_step)
        self.assertGreater(
            install_step.index("for package in DESeq Biobase edgeR; do"),
            install_step.index("\n          RSCRIPT\n"),
        )
        for relative_path in (
            '"$package_dir/DESCRIPTION"',
            '"$package_dir/NAMESPACE"',
            '"$package_dir/Meta/package.rds"',
        ):
            self.assertIn(f"test -s {relative_path}", install_step)
        self.assertIn("-name '00LOCK*'", install_step)
        self.assertIn(
            "bioconductor/bioconductor_docker:RELEASE_3_11@sha256:"
            "cbd868b0543608c917cef2003cc8f051ba15e6633bf941f35802825b1e5551ab",
            workflow,
        )
        self.assertIn("id: imvigor_export", workflow)
        self.assertIn("steps.imvigor_export.outcome", workflow)
        self.assertIn(
            "--diagnostics-path "
            "/runner-temp/imvigor210-verification/export_diagnostics.tsv",
            workflow,
        )
        self.assertIn(
            "python scripts/verify_imvigor210_expression.py \\\n",
            workflow,
        )
        self.assertIn("--external-semantic-file-verification", workflow)
        self.assertIn(
            "expression_semantic_verification.json", workflow
        )
        self.assertEqual(
            workflow.count("Validate the non-identifying semantic framing"), 1
        )
        self.assertIn(
            "path: ${{ runner.temp }}/imvigor210-verification/",
            workflow,
        )
        imvigor_upload = workflow.split(
            "      - name: Upload IMvigor210 verification report only\n", 1
        )[1].split("\n\n  targeted-singlecell:", 1)[0]
        self.assertIn("if: always()", imvigor_upload)
        self.assertNotIn("IMVIGOR_EXPORT_DIR", imvigor_upload)
        self.assertNotIn("imvigor210-exports", imvigor_upload)

    def test_imvigor_exporter_preserves_named_transaction_paths(self) -> None:
        source = (ROOT / "scripts" / "export_imvigor210_inputs.R").read_text(
            encoding="utf-8"
        )
        normalized = " ".join(source.split())
        self.assertIn(
            'temporary_paths <- stats::setNames( paste0(unname(final_paths), ".tmp"), names(final_paths) )',
            normalized,
        )
        self.assertIn("temporary_paths[[name]]", source)

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

    def test_publication_surface_is_indexed_without_empty_placeholders(self) -> None:
        for relative in (
            ".github/workflows/README.md",
            "R/README.md",
            "data/README.md",
            "docs/README.md",
            "reference_results/README.md",
            "resources/README.md",
            "scripts/README.md",
            "tests/README.md",
            "validation/README.md",
        ):
            with self.subTest(relative=relative):
                self.assertTrue((ROOT / relative).is_file())

        reference_dir = ROOT / "reference_figures"
        self.assertFalse(reference_dir.exists() and any(reference_dir.iterdir()))
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        self.assertIn(
            "The repository is organized around the experimental TNFR1/CAR-T analyses.",
            readme,
        )
        self.assertIn("## Reproduce manuscript figures", readme)
        self.assertNotIn(
            "No SQLite database or monolithic R data object is required",
            readme,
        )


if __name__ == "__main__":
    unittest.main()
