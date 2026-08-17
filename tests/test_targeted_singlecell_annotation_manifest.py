from __future__ import annotations

import csv
import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "resources" / "targeted_singlecell_cluster_annotations_v1.tsv"
EXPECTED_COLUMNS = (
    "manifest_version",
    "analysis_source_commit",
    "dataset",
    "cluster_id",
    "n_cells",
    "previous_release_label",
    "submission_label",
    "defining_markers",
    "annotation_method",
    "marker_source",
    "nomenclature_framework",
    "interpretive_reference",
    "reference_atlas_or_classifier",
    "assay_scope",
    "properties_not_measured",
    "transfer_status",
)
EXPECTED_LABELS = {
    "tumor_co_culture": {
        "C0": "C0 KLRB1/LGALS3-associated activated T-cell state",
        "C1": "C1 CD4/LAG3-associated activated T-cell state",
        "C2": "C2 CXCR5/IL13/CCR4-associated activated T-cell state",
        "C3": "C3 CXCR6-associated cytotoxic activated state",
        "C4": "C4 CD8/ZNF683-associated cytotoxic state",
        "C5": "C5 TRDC-high γδ-associated cytotoxic state",
        "C6": "C6 CXCL13-associated cycling T-cell state",
        "C7": "C7 cycling effector-gene-high T-cell state",
        "C8": "C8 IL9-high activated T-cell state",
        "C9": "C9 TCF7/IL7R-high early-memory-associated state",
        "C10": "C10 small cytokine/IFN-response-high cluster",
    },
    "repeated_stimulation": {
        "C0": "Mixed CD4/KLRB1-associated activated state",
        "C1": "Cycling T-cell state I",
        "C2": "Cytokine-expressing effector state",
        "C3": "CD8/TRDC-associated cytotoxic state",
        "C4": "Cycling T-cell state II",
        "C5": "CCR7/IL7R/HLA-II-associated state",
    },
}
EXPECTED_COUNTS = {
    "tumor_co_culture": {
        "C0": 1542,
        "C1": 1448,
        "C2": 903,
        "C3": 891,
        "C4": 812,
        "C5": 594,
        "C6": 516,
        "C7": 453,
        "C8": 396,
        "C9": 133,
        "C10": 33,
    },
    "repeated_stimulation": {
        "C0": 1298,
        "C1": 1034,
        "C2": 1015,
        "C3": 833,
        "C4": 806,
        "C5": 444,
    },
}


class TargetedSingleCellAnnotationManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with MANIFEST.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            cls.header = tuple(reader.fieldnames or ())
            cls.rows = list(reader)

    def test_schema_labels_and_qc_cell_counts(self) -> None:
        self.assertEqual(self.header, EXPECTED_COLUMNS)
        self.assertEqual(len(self.rows), 17)
        observed_labels: dict[str, dict[str, str]] = {}
        observed_counts: dict[str, dict[str, int]] = {}
        for row in self.rows:
            self.assertEqual(row["manifest_version"], "1")
            self.assertRegex(row["analysis_source_commit"], r"^[0-9a-f]{40}$")
            dataset = row["dataset"]
            cluster = row["cluster_id"]
            self.assertNotIn(cluster, observed_labels.setdefault(dataset, {}))
            observed_labels[dataset][cluster] = row["submission_label"]
            observed_counts.setdefault(dataset, {})[cluster] = int(row["n_cells"])
        self.assertEqual(observed_labels, EXPECTED_LABELS)
        self.assertEqual(observed_counts, EXPECTED_COUNTS)
        self.assertEqual(sum(observed_counts["tumor_co_culture"].values()), 7721)
        self.assertEqual(sum(observed_counts["repeated_stimulation"].values()), 5430)

    def test_method_and_nomenclature_roles_are_explicit(self) -> None:
        for row in self.rows:
            self.assertEqual(
                row["annotation_method"],
                "manual_post_hoc_marker_based_after_unsupervised_louvain_clustering",
            )
            self.assertEqual(row["reference_atlas_or_classifier"], "none")
            self.assertIn(
                "doi:10.1038/s41577-025-01238-2",
                row["nomenclature_framework"],
            )
            self.assertIn("259-gene targeted single-cell mRNA panel", row["assay_scope"])
            self.assertTrue(row["defining_markers"].strip())
            self.assertGreaterEqual(len(row["defining_markers"].split(";")), 5)
            self.assertTrue(row["marker_source"].strip())
            self.assertTrue(row["properties_not_measured"].strip())
            self.assertTrue(row["transfer_status"].strip())
            if row["dataset"] == "repeated_stimulation":
                self.assertIn("no mapping to tumor clusters", row["transfer_status"])

    def test_submission_labels_match_r05_and_legacy_overclaims_are_removed(self) -> None:
        r05 = (ROOT / "R" / "05_figure_4_AB_suppl_S5A.R").read_text(
            encoding="utf-8"
        )
        for row in self.rows:
            self.assertIn(f'"{row["submission_label"]}"', r05)
            if row["previous_release_label"] != row["submission_label"]:
                self.assertNotIn(f'"{row["previous_release_label"]}"', r05)
        for legacy in (
            "innate-like activated",
            "checkpoint-high",
            "type 2 helper-like",
            "cycling effector-like",
            "TH9-like",
            "stem-like/early-memory",
            "Cytokine-producing effector state",
        ):
            self.assertNotIn(legacy, r05)

    def test_c6_and_c10_scope_is_explicit(self) -> None:
        c6 = next(
            row
            for row in self.rows
            if row["dataset"] == "tumor_co_culture"
            and row["cluster_id"] == "C6"
        )
        self.assertIn("CXCL13", c6["defining_markers"].split(";"))
        self.assertIn("doi:10.1038/s43018-022-00433-7", c6["interpretive_reference"])
        self.assertIn("uniform CXCL13 positivity", c6["properties_not_measured"])

        c10 = next(
            row
            for row in self.rows
            if row["dataset"] == "tumor_co_culture"
            and row["cluster_id"] == "C10"
        )
        self.assertIn("outside historical C0-C9", c10["transfer_status"])
        self.assertIn("replicate-level enrichment", c10["properties_not_measured"])

    def test_manifest_contains_no_cell_level_records(self) -> None:
        forbidden_columns = {
            "cell",
            "cell_id",
            "cell_identifier",
            "barcode",
            "sample_id",
            "donor_id",
            "expression",
            "raw_count",
        }
        self.assertFalse(forbidden_columns.intersection(self.header))
        for row in self.rows:
            self.assertFalse(
                any(
                    re.search(r"cell[_ -]?(id|index|barcode)", value, re.IGNORECASE)
                    for value in row.values()
                )
            )


if __name__ == "__main__":
    unittest.main()

