from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VERIFIER_PATH = ROOT / "scripts" / "verify_tcr_c6_projection_outputs.py"
SPEC = importlib.util.spec_from_file_location("tcr_c6_projection_verifier", VERIFIER_PATH)
assert SPEC is not None and SPEC.loader is not None
VERIFIER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VERIFIER)


def write_rows(
    path: Path,
    fieldnames: list[str] | tuple[str, ...],
    rows: list[dict[str, object]],
    *,
    delimiter: str,
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


class TcrC6ProjectionVerifierTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.by_cluster = self.root / "by_cluster.tsv"
        self.gene_coverage = self.root / "gene_coverage.tsv"
        self.cluster_counts = self.root / "cluster_counts.csv"
        self.signature = self.root / "signature.csv"

        counts = [10, 20, 30, 40, 50, 60]
        total = sum(counts)
        write_rows(
            self.cluster_counts,
            VERIFIER.CLUSTER_COUNT_COLUMNS,
            [
                {
                    "cluster_short": cluster,
                    "n_cells": count,
                    "fraction": count / total,
                }
                for cluster, count in zip(VERIFIER.EXPECTED_CLUSTERS, counts)
            ],
            delimiter=",",
        )
        write_rows(
            self.by_cluster,
            VERIFIER.BY_CLUSTER_COLUMNS,
            [
                {
                    "cluster_short": cluster,
                    "cluster_annotation": VERIFIER.EXPECTED_LABELS[cluster],
                    "n_cells": count,
                    "c6_full_rank_score_mean": 0.20 + index / 100,
                    "c6_full_rank_score_median": 0.20 + index / 100,
                    "c6_full_rank_score_q25": 0.15 + index / 100,
                    "c6_full_rank_score_q75": 0.25 + index / 100,
                    "c6_cycle_rank_score_mean": 0.30 + index / 100,
                    "c6_cycle_rank_score_median": 0.30 + index / 100,
                    "c6_cycle_rank_score_q25": 0.25 + index / 100,
                    "c6_cycle_rank_score_q75": 0.35 + index / 100,
                    "c6_noncycle_rank_score_mean": 0.10 + index / 100,
                    "c6_noncycle_rank_score_median": 0.10 + index / 100,
                    "c6_noncycle_rank_score_q25": 0.05 + index / 100,
                    "c6_noncycle_rank_score_q75": 0.15 + index / 100,
                    "cxcl13_detected_cells": index,
                    "cxcl13_detection_fraction": index / count,
                    "cxcl13_mean_log_normalized_expression": index / 100,
                }
                for index, (cluster, count) in enumerate(
                    zip(VERIFIER.EXPECTED_CLUSTERS, counts)
                )
            ],
            delimiter="\t",
        )
        genes = ["CXCL13"] + [f"GENE_{index}" for index in range(1, 20)]
        signature_rows = [
            {"cluster": 6, "gene": gene, "avg_log2FC": 2.0 - index / 100}
            for index, gene in enumerate(genes)
        ]
        write_rows(
            self.signature,
            ("cluster", "gene", "avg_log2FC"),
            signature_rows,
            delimiter=",",
        )
        write_rows(
            self.gene_coverage,
            VERIFIER.GENE_COVERAGE_COLUMNS,
            [
                {
                    "gene": gene,
                    "component": (
                        "cycle" if gene in VERIFIER.EXPECTED_CYCLE_GENES else "noncycle"
                    ),
                    "frozen_avg_log2FC": signature_row["avg_log2FC"],
                    "in_targeted_panel": "TRUE",
                    "used_in_score": "TRUE",
                    "exclusion_reason": "",
                }
                for gene, signature_row in zip(genes, signature_rows)
            ],
            delimiter="\t",
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def verify(self) -> None:
        VERIFIER.verify_release(
            self.by_cluster,
            self.gene_coverage,
            self.cluster_counts,
            self.signature,
        )

    def test_valid_aggregate_outputs_pass(self) -> None:
        self.verify()

    def test_extra_cell_identifier_column_is_rejected(self) -> None:
        header, rows = VERIFIER.read_delimited(self.by_cluster, "\t")
        rows[0]["cell_id"] = "forbidden-cell-barcode"
        write_rows(
            self.by_cluster,
            header + ["cell_id"],
            rows,
            delimiter="\t",
        )
        with self.assertRaisesRegex(VERIFIER.ContractError, "Unexpected columns"):
            self.verify()

    def test_cluster_count_mismatch_is_rejected(self) -> None:
        header, rows = VERIFIER.read_delimited(self.by_cluster, "\t")
        rows[0]["n_cells"] = "9"
        write_rows(self.by_cluster, header, rows, delimiter="\t")
        with self.assertRaisesRegex(VERIFIER.ContractError, "differs from count table"):
            self.verify()

    def test_rank_score_outside_unit_interval_is_rejected(self) -> None:
        header, rows = VERIFIER.read_delimited(self.by_cluster, "\t")
        rows[0]["c6_full_rank_score_mean"] = "1.01"
        write_rows(self.by_cluster, header, rows, delimiter="\t")
        with self.assertRaisesRegex(VERIFIER.ContractError, r"rank scores in \[0, 1\]"):
            self.verify()

    def test_missing_signature_gene_is_rejected(self) -> None:
        header, rows = VERIFIER.read_delimited(self.gene_coverage, "\t")
        write_rows(self.gene_coverage, header, rows[:-1], delimiter="\t")
        with self.assertRaisesRegex(VERIFIER.ContractError, "one row per frozen C6 gene"):
            self.verify()

    def test_cxcl13_exclusion_is_rejected(self) -> None:
        header, rows = VERIFIER.read_delimited(self.gene_coverage, "\t")
        rows[0]["used_in_score"] = "FALSE"
        rows[0]["exclusion_reason"] = "excluded"
        write_rows(self.gene_coverage, header, rows, delimiter="\t")
        with self.assertRaisesRegex(VERIFIER.ContractError, "must be used"):
            self.verify()

    def test_cycle_component_change_is_rejected(self) -> None:
        header, rows = VERIFIER.read_delimited(self.gene_coverage, "\t")
        rows[0]["component"] = "cycle"
        write_rows(self.gene_coverage, header, rows, delimiter="\t")
        with self.assertRaisesRegex(VERIFIER.ContractError, "component must be noncycle"):
            self.verify()


class TcrC6ProjectionRepositoryIntegrationTests(unittest.TestCase):
    def test_workflow_uses_explicit_aggregate_artifact_allowlist(self) -> None:
        workflow = (
            ROOT / ".github" / "workflows" / "cohort-inputs-targeted-singlecell.yml"
        ).read_text(encoding="utf-8")
        self.assertIn("Exploratory_TCR_C6_signature_projection", workflow)
        self.assertIn("verify_tcr_c6_projection_outputs.py", workflow)
        self.assertIn("artifact_outputs=(", workflow)
        self.assertNotIn("find figures -maxdepth 1 -type f", workflow)
        self.assertNotIn("TCR_seurat_object.rds", workflow)

    def test_projection_paths_trigger_the_targeted_workflow(self) -> None:
        workflow = (
            ROOT / ".github" / "workflows" / "cohort-inputs-targeted-singlecell.yml"
        ).read_text(encoding="utf-8")
        for path in (
            '"scripts/verify_tcr_c6_projection_outputs.py"',
            '"resources/targeted_singlecell_cluster_annotations_v1.tsv"',
            '"tests/test_tcr_c6_projection_outputs.py"',
            '"tests/test_targeted_singlecell_annotation_manifest.py"',
        ):
            self.assertEqual(workflow.count(path), 2, path)


if __name__ == "__main__":
    unittest.main()
