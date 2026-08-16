from __future__ import annotations

import ast
import csv
import gzip
import hashlib
import importlib.util
import json
import sys
import types
import unittest
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "experimental" / "bulk_rnaseq"
DERIVED_DIR = DATA_DIR / "derived"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class BulkRnaSeqContractTests(unittest.TestCase):
    def test_checksums(self) -> None:
        for line in (DATA_DIR / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
            expected, filename = line.split(maxsplit=1)
            self.assertEqual(sha256_file(DATA_DIR / filename), expected, filename)

    def test_derived_release_checksums_and_universes(self) -> None:
        for line in (DERIVED_DIR / "SHA256SUMS").read_text(
            encoding="utf-8"
        ).splitlines():
            expected, filename = line.split(maxsplit=1)
            self.assertEqual(sha256_file(DERIVED_DIR / filename), expected, filename)

        expected_rows = {
            "figure_1b_1c_wt_cytokine_contrasts.unfiltered.tsv.gz": 46_425 * 3,
            "figure_2b_s2d_tnfr1_ko1_vs_wt_matched_treatments.unfiltered.tsv.gz": 46_425 * 4,
            "interaction_TNFR1_KO1_vs_WT_IFNG_vs_control.unfiltered.tsv.gz": 46_425,
            "interaction_TNFR1_KO1_vs_WT_TNF_vs_control.unfiltered.tsv.gz": 46_425,
            "interaction_TNFR1_KO1_vs_WT_TNF_IFNG_vs_control.unfiltered.tsv.gz": 46_425,
        }
        for filename, row_count in expected_rows.items():
            with gzip.open(DERIVED_DIR / filename, "rt", encoding="utf-8") as handle:
                header = next(handle).rstrip("\n").split("\t")
                observed_rows = sum(1 for _ in handle)
            self.assertEqual(observed_rows, row_count, filename)
            self.assertTrue(
                {"Gene_Symbol", "baseMean", "log2FoldChange", "stat", "padj"}
                .issubset(header)
                or {"condition", "gene_symbol", "base_mean", "lfc_se"}
                .issubset(header),
                filename,
            )

        runtime = json.loads((DERIVED_DIR / "run_metadata.json").read_text())
        self.assertEqual(runtime["python"], "3.12.13")
        self.assertEqual(runtime["pydeseq2"], "0.5.4")
        self.assertEqual(runtime["n_cpus"], 2)
        self.assertTrue(runtime["include_interaction"])

    def test_gene_id_count_matrix(self) -> None:
        with gzip.open(DATA_DIR / "gene_counts.tsv.gz", "rt", encoding="utf-8") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader)
            self.assertEqual(header[0], "gene_id")
            self.assertEqual(len(header), 25)
            gene_ids = set()
            rows = 0
            for row in reader:
                rows += 1
                self.assertEqual(len(row), 25)
                self.assertNotIn(row[0], gene_ids)
                gene_ids.add(row[0])
                values = [int(value) for value in row[1:]]
                self.assertTrue(all(value >= 0 for value in values))
        self.assertEqual(rows, 46_427)
        self.assertEqual(len(gene_ids), 46_427)

    def test_annotation_and_symbol_membership(self) -> None:
        symbol_counts: Counter[str] = Counter()
        gene_ids = set()
        trnav_ids = []
        with gzip.open(
            DATA_DIR / "gene_annotations.tsv.gz", "rt", encoding="utf-8", newline=""
        ) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                gene_ids.add(row["gene_id"])
                symbol_counts[row["gene_symbol"]] += 1
                if row["gene_symbol"] == "TRNAV-CAC":
                    trnav_ids.append(row["gene_id"])
        self.assertEqual(len(gene_ids), 46_427)
        self.assertEqual(len(symbol_counts), 46_425)
        self.assertEqual(
            {symbol: count for symbol, count in symbol_counts.items() if count > 1},
            {"TRNAV-CAC": 3},
        )
        self.assertEqual(
            sorted(trnav_ids, key=int),
            ["107985614", "107985615", "107985753"],
        )

        with gzip.open(
            DATA_DIR / "gene_symbol_membership.tsv.gz",
            "rt",
            encoding="utf-8",
            newline="",
        ) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 46_425)
        duplicate = [row for row in rows if int(row["n_gene_ids"]) > 1]
        self.assertEqual(len(duplicate), 1)
        self.assertEqual(duplicate[0]["gene_symbol"], "TRNAV-CAC")
        self.assertEqual(
            duplicate[0]["gene_ids"], "107985614;107985615;107985753"
        )

    def test_sample_metadata(self) -> None:
        with (DATA_DIR / "sample_metadata.tsv").open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 24)
        self.assertEqual(len({row["sample_id"] for row in rows}), 24)
        groups = Counter((row["genotype"], row["treatment"]) for row in rows)
        self.assertEqual(len(groups), 8)
        self.assertTrue(all(count == 3 for count in groups.values()))
        self.assertTrue(all(row["treatment_duration_h"] == "48" for row in rows))
        self.assertTrue(all(row["replicate_unit"] == "independent_experiment" for row in rows))
        self.assertTrue(all(row["paired_batch_status"] == "not_confirmed" for row in rows))
        for row in rows:
            expected_tnf = "50" if row["treatment"] in {"TNF", "TNF_IFNG"} else "0"
            expected_ifng = "50" if row["treatment"] in {"IFNG", "TNF_IFNG"} else "0"
            self.assertEqual(row["tnf_ng_ml"], expected_tnf)
            self.assertEqual(row["ifng_ng_ml"], expected_ifng)
            if row["genotype"] == "TNFR1_KO1":
                self.assertEqual(row["source_clone"], "T6")

    def test_source_provenance_and_qc(self) -> None:
        with (DATA_DIR / "source_manifest.tsv").open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(
            row["filename_as_received"],
            "Expression_Profile.GRCh38.gene(20260814-181921).xlsx",
        )
        self.assertEqual(
            row["canonical_archive_filename"], "Expression_Profile.GRCh38.gene.xlsx"
        )
        self.assertEqual(int(row["size_bytes"]), 8_022_936)
        self.assertEqual(
            row["sha256"],
            "41cee68a4cd33f72268b46dc78fd0708dd62655814ae19a8a39fefd8aa5d5989",
        )

        qc = json.loads((DATA_DIR / "source_qc.json").read_text(encoding="utf-8"))
        self.assertEqual(qc["count_matrix"]["gene_id_rows"], 46_427)
        self.assertEqual(qc["count_matrix"]["unique_gene_symbols"], 46_425)
        self.assertEqual(qc["count_matrix"]["missing_counts"], 0)
        self.assertEqual(
            qc["duplicate_gene_symbols"],
            {"TRNAV-CAC": [107985614, 107985615, 107985753]},
        )
        correlations = qc["within_group_replicate_correlations"]
        self.assertEqual(len(correlations), 8)
        self.assertTrue(all(group["minimum"] > 0.95 for group in correlations.values()))

    def test_analysis_script_static_contract(self) -> None:
        path = ROOT / "scripts" / "run_bulk_rnaseq_pydeseq2.py"
        source = path.read_text(encoding="utf-8")
        ast.parse(source)
        required_fragments = (
            'REQUIRED_PYDESEQ2_VERSION = "0.5.4"',
            '["treatment", treatment, "control"]',
            '["genotype", "TNFR1_KO1", "WT"]',
            '"~genotype + treatment + genotype:treatment"',
            'dds.cond(genotype="TNFR1_KO1", treatment=treatment)',
            "contrast = np.asarray(",
            "ko_treated - wt_treated - ko_control + wt_control",
            '"all_zero_in_subset"',
            'complete.loc[all_zero, "baseMean"] = 0.0',
            '"All-zero features unexpectedly contain test statistics"',
            "len(complete) != 46_425",
            '"environment.freeze.txt"',
        )
        for fragment in required_fragments:
            self.assertIn(fragment, source)

    def test_r_figure_scripts_use_full_unfiltered_adapters(self) -> None:
        helper = (ROOT / "R" / "bulk_rnaseq_figure_helpers.R").read_text(
            encoding="utf-8"
        )
        self.assertIn("BULK_RNASEQ_EXPECTED_SYMBOLS <- 46425L", helper)
        self.assertIn("Gene-symbol universes differ between conditions", helper)
        self.assertIn("expected_statistic <-", helper)
        self.assertIn("observed_statistic <-", helper)
        self.assertIn("lfc_se must be positive", helper)
        self.assertIn("plottable_bulk_rows <- function(adapter)", helper)
        self.assertIn("plottable_modelled_symbols = nrow(plottable)", helper)
        self.assertIn("BULK_VOLCANO_P_VALUE_FLOOR <- 1e-300", helper)
        self.assertIn("BULK_VOLCANO_Y_CAP <- 300", helper)
        self.assertIn("symmetric_bulk_volcano_x_limits <- function", helper)
        self.assertIn("write_bulk_volcano_output_contract <- function", helper)
        self.assertIn("save_bulk_volcano_triptych <- function", helper)
        self.assertIn('width = 6.5, height = 3, dpi = 600', helper)
        self.assertIn('compression = "lzw"', helper)

        figure_1 = (ROOT / "R" / "03_figure_1_B_C.R").read_text(
            encoding="utf-8"
        )
        figure_2 = (ROOT / "R" / "04_figure_2_B_suppl_S2D.R").read_text(
            encoding="utf-8"
        )
        self.assertIn(
            "figure_1b_1c_wt_cytokine_contrasts.unfiltered.tsv.gz", figure_1
        )
        self.assertIn(
            "log2_fold_change_treatment_vs_untreated", figure_1
        )
        self.assertIn("wald_statistic_treatment_vs_untreated", figure_1)
        self.assertIn('figure_dir <- file.path("figures", "figure_1")', figure_1)
        self.assertIn('result_dir <- file.path("results", "figure_1")', figure_1)
        self.assertIn("Figure_1C_upregulated_overlap", figure_1)
        self.assertIn(
            'publication_condition_order <- c("TNF_IFNG", "IFNG", "TNF")',
            figure_1,
        )
        self.assertIn("Figure_1B_triptych", figure_1)
        self.assertIn("Figure_1B_volcano_output_contract.tsv", figure_1)
        for gene in (
            "ICAM1",
            "IRF1",
            "AIM2",
            "CASP1",
            "CASP4",
            "MLKL",
            "BAK1",
            "CASP7",
            "FAS",
            "CASP8",
        ):
            self.assertIn(f'"{gene}"', figure_1)
        self.assertNotIn('"CASP5"', figure_1)
        self.assertIn(
            "figure_2b_s2d_tnfr1_ko1_vs_wt_matched_treatments.unfiltered.tsv.gz",
            figure_2,
        )
        self.assertIn("log2_fold_change_ko1_vs_wt", figure_2)
        self.assertIn("wald_statistic_ko1_vs_wt", figure_2)
        self.assertIn('figure_dir <- file.path("figures", "figure_2")', figure_2)
        self.assertIn('result_dir <- file.path("results", "figure_2")', figure_2)
        self.assertIn("Supplementary_Figure_S2D_downregulated_overlap", figure_2)
        self.assertIn(
            'publication_condition_order <- c("TNF_IFNG", "IFNG", "TNF")',
            figure_2,
        )
        self.assertIn("Figure_2B_triptych", figure_2)
        self.assertIn("Figure_2B_volcano_output_contract.tsv", figure_2)
        self.assertIn("Figure_2C_ICAM1_IRF1_effects", figure_2)
        self.assertIn("wald_ci_95_lower", figure_2)
        self.assertIn("wald_ci_95_upper", figure_2)
        self.assertNotIn("geom_signif", figure_2)
        self.assertIn("figure_2b_label_genes <- c(", figure_2)
        for gene in ("ICAM1", "MLKL", "GSDME", "IRF1"):
            self.assertIn(f'"{gene}"', figure_2)
        self.assertNotIn('"CASP5"', figure_2)
        self.assertIn(
            'expected_conditions <- c("control", "IFNG", "TNF", "TNF_IFNG")',
            figure_2,
        )
        for source in (figure_1, figure_2):
            self.assertNotIn("read_excel", source)
            self.assertNotIn("readxl::", source)
            self.assertNotIn("data/rnaseq", source)
            self.assertNotIn("Gene_ID", source)
            self.assertNotIn("writexl", source)
            self.assertNotIn("write_xlsx", source)

        verifier = (
            ROOT / "scripts" / "verify_bulk_figure_outputs.py"
        ).read_text(encoding="utf-8")
        ast.parse(verifier)
        self.assertIn("PUBLICATION_CONDITION_ORDER", verifier)
        self.assertIn("points_at_displayed_y_cap", verifier)
        self.assertIn("TIFF compression is not LZW", verifier)

    def test_publication_volcano_numeric_contract(self) -> None:
        path = ROOT / "scripts" / "verify_bulk_figure_outputs.py"
        module_name = "bulk_figure_output_verifier_for_test"
        spec = importlib.util.spec_from_file_location(module_name, path)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        previous = sys.modules.get(module_name)
        sys.modules[module_name] = module
        try:
            spec.loader.exec_module(module)
            figure_1, figure_2 = module.figure_specs(ROOT)
            figure_1_rows = module.expected_contract_rows(figure_1)
            figure_2_rows = module.expected_contract_rows(figure_2)
            figure_2c_rows = module.expected_figure_2c_rows(ROOT)
        finally:
            if previous is None:
                sys.modules.pop(module_name, None)
            else:
                sys.modules[module_name] = previous

        self.assertEqual(
            [row["condition"] for row in figure_1_rows],
            ["TNF_IFNG", "IFNG", "TNF"],
        )
        self.assertEqual(
            [row["plottable_modelled_symbols"] for row in figure_1_rows],
            [20_189, 19_827, 18_555],
        )
        self.assertEqual(
            [row["points_at_displayed_y_cap"] for row in figure_1_rows],
            [285, 96, 23],
        )
        self.assertTrue(
            all(
                (row["common_x_min"], row["common_x_max"]) == (-19, 19)
                for row in figure_1_rows
            )
        )
        self.assertEqual(
            [row["plotted_label_count"] for row in figure_1_rows],
            [10, 10, 9],
        )

        self.assertEqual(
            [row["condition"] for row in figure_2_rows],
            ["TNF_IFNG", "IFNG", "TNF"],
        )
        self.assertEqual(
            [row["plottable_modelled_symbols"] for row in figure_2_rows],
            [18_661, 19_858, 17_743],
        )
        self.assertEqual(
            [row["points_at_displayed_y_cap"] for row in figure_2_rows],
            [32, 7, 37],
        )
        self.assertTrue(
            all(
                (row["common_x_min"], row["common_x_max"]) == (-12, 12)
                for row in figure_2_rows
            )
        )
        self.assertEqual(
            [row["plotted_label_count"] for row in figure_2_rows],
            [4, 4, 4],
        )
        self.assertEqual(
            [(row["condition"], row["gene_symbol"]) for row in figure_2c_rows],
            [
                ("TNF_IFNG", "ICAM1"),
                ("TNF_IFNG", "IRF1"),
                ("IFNG", "ICAM1"),
                ("IFNG", "IRF1"),
                ("TNF", "ICAM1"),
                ("TNF", "IRF1"),
            ],
        )
        self.assertEqual(
            [row["significance_code"] for row in figure_2c_rows],
            ["****", "***", "****", "****", "****", "****"],
        )
        self.assertAlmostEqual(
            figure_2c_rows[0]["log2_fold_change_ko1_vs_wt"],
            -1.9725787735602656,
        )
        self.assertAlmostEqual(
            figure_2c_rows[1]["adjusted_p_value_bh"],
            8.785896839383e-4,
        )
        self.assertAlmostEqual(
            figure_2c_rows[0]["wald_ci_95_lower"],
            -2.01701641260892,
        )

    def test_analysis_input_loading_and_factor_levels(self) -> None:
        package = types.ModuleType("pydeseq2")
        dds_module = types.ModuleType("pydeseq2.dds")
        ds_module = types.ModuleType("pydeseq2.ds")
        dds_module.DeseqDataSet = type("DeseqDataSet", (), {})
        ds_module.DeseqStats = type("DeseqStats", (), {})
        previous = {
            name: sys.modules.get(name)
            for name in ("pydeseq2", "pydeseq2.dds", "pydeseq2.ds")
        }
        sys.modules["pydeseq2"] = package
        sys.modules["pydeseq2.dds"] = dds_module
        sys.modules["pydeseq2.ds"] = ds_module
        try:
            path = ROOT / "scripts" / "run_bulk_rnaseq_pydeseq2.py"
            spec = importlib.util.spec_from_file_location("bulk_rnaseq_module", path)
            self.assertIsNotNone(spec)
            self.assertIsNotNone(spec.loader)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)

            counts, annotations, metadata = module.load_inputs(DATA_DIR)
            sample_ids = metadata["sample_id"].tolist()
            symbol_counts, membership = module.aggregate_to_gene_symbols(
                counts, annotations, sample_ids
            )
            self.assertEqual(symbol_counts.shape, (46_425, 24))
            self.assertEqual(
                membership.loc["TRNAV-CAC", "source_gene_ids"],
                "107985614;107985615;107985753",
            )

            wt_ids = metadata.loc[
                metadata["genotype"].eq("WT")
                & metadata["treatment"].isin(["control", "IFNG"]),
                "sample_id",
            ].tolist()
            wt_metadata = module.make_metadata(metadata, wt_ids)
            self.assertEqual(wt_metadata["genotype"].cat.categories.tolist(), ["WT"])
            self.assertEqual(
                wt_metadata["treatment"].cat.categories.tolist(), ["control", "IFNG"]
            )

            matched_ids = metadata.loc[
                metadata["treatment"].eq("TNF"), "sample_id"
            ].tolist()
            matched_metadata = module.make_metadata(metadata, matched_ids)
            self.assertEqual(
                matched_metadata["genotype"].cat.categories.tolist(),
                ["WT", "TNFR1_KO1"],
            )
            self.assertEqual(
                matched_metadata["treatment"].cat.categories.tolist(), ["TNF"]
            )

            example = module.pd.DataFrame(
                {
                    "Gene_Symbol": ["ICAM1"],
                    "baseMean": [100.0],
                    "log2FoldChange": [-1.0],
                    "lfcSE": [0.2],
                    "stat": [-5.0],
                    "pvalue": [1e-5],
                    "padj": [2e-4],
                }
            )
            wt_adapter = module.figure_adapter(
                example, "IFNG", "treatment_vs_untreated"
            )
            ko_adapter = module.figure_adapter(example, "IFNG", "ko1_vs_wt")
            self.assertEqual(tuple(wt_adapter.columns), module.WT_FIGURE_COLUMNS)
            self.assertEqual(tuple(ko_adapter.columns), module.KO_FIGURE_COLUMNS)
            self.assertNotIn("Gene_ID", wt_adapter.columns)
            self.assertNotIn("Gene_ID", ko_adapter.columns)
        finally:
            for name, original in previous.items():
                if original is None:
                    sys.modules.pop(name, None)
                else:
                    sys.modules[name] = original


if __name__ == "__main__":
    unittest.main()
