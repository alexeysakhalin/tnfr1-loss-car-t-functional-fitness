from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VERIFIER_PATH = ROOT / "scripts" / "verify_xlsx_workbook.py"
SPEC = importlib.util.spec_from_file_location("verify_xlsx_workbook", VERIFIER_PATH)
assert SPEC is not None and SPEC.loader is not None
VERIFIER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = VERIFIER
SPEC.loader.exec_module(VERIFIER)


CONTENT_TYPES = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
  <Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
  <Override PartName="/xl/worksheets/sheet2.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
</Types>
"""

ROOT_RELS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
</Relationships>
"""

WORKBOOK = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
 xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <sheets>
    <sheet name="All_clusters_top20" sheetId="1" r:id="rId1"/>
    <sheet name="Cluster_0" sheetId="2" r:id="rId2"/>
  </sheets>
</workbook>
"""

WORKBOOK_RELS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet2.xml"/>
</Relationships>
"""


def worksheet_xml(dimension: str) -> str:
    rows = []
    for row in range(1, 4):
        cells = []
        for column in ("A", "B"):
            cells.append(
                f'<c r="{column}{row}" t="inlineStr"><is><t>{column}{row}</t></is></c>'
            )
        rows.append(f'<row r="{row}">{"".join(cells)}</row>')
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f'<dimension ref="{dimension}"/><sheetData>{"".join(rows)}</sheetData>'
        "</worksheet>"
    )


def write_fixture(
    path: Path,
    *,
    first_dimension: str = "A1:B3",
    dangling_drawing_relationship: bool = False,
) -> None:
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml", CONTENT_TYPES)
        archive.writestr("_rels/.rels", ROOT_RELS)
        archive.writestr("xl/workbook.xml", WORKBOOK)
        archive.writestr("xl/_rels/workbook.xml.rels", WORKBOOK_RELS)
        archive.writestr("xl/worksheets/sheet1.xml", worksheet_xml(first_dimension))
        archive.writestr("xl/worksheets/sheet2.xml", worksheet_xml("A1:B3"))
        if dangling_drawing_relationship:
            archive.writestr(
                "xl/worksheets/_rels/sheet1.xml.rels",
                """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/drawing" Target="../drawings/drawing1.xml"/>
</Relationships>
""",
            )


class XlsxWorkbookVerifierTests(unittest.TestCase):
    def verify(self, path: Path) -> None:
        VERIFIER.verify_marker_workbook(
            path,
            expected_clusters=["0"],
            top_n=2,
            expected_columns=2,
        )

    def test_accepts_complete_marker_workbook(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "valid.xlsx"
            write_fixture(path)
            self.verify(path)

    def test_accepts_complete_marker_workbook_in_openpyxl_normal_mode(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "valid-openpyxl.xlsx"
            write_fixture(path)
            VERIFIER.verify_marker_workbook(
                path,
                expected_clusters=["0"],
                top_n=2,
                expected_columns=2,
                require_openpyxl=True,
            )

    def test_rejects_dimension_that_hides_populated_cells(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "bad-dimension.xlsx"
            write_fixture(path, first_dimension="A1")
            with self.assertRaisesRegex(
                VERIFIER.WorkbookValidationError,
                "declares dimension 'A1'.*require 'A1:B3'",
            ):
                self.verify(path)

    def test_rejects_relationship_to_missing_drawing_part(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "dangling-relationship.xlsx"
            write_fixture(path, dangling_drawing_relationship=True)
            with self.assertRaisesRegex(
                VERIFIER.WorkbookValidationError,
                "points to missing part 'xl/drawings/drawing1.xml'",
            ):
                self.verify(path)

    def test_rejects_unexpected_sheet_contract(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "wrong-sheets.xlsx"
            write_fixture(path)
            with self.assertRaisesRegex(
                VERIFIER.WorkbookValidationError,
                "Worksheet order/name contract failed",
            ):
                VERIFIER.verify_marker_workbook(
                    path,
                    expected_clusters=["0", "1"],
                    top_n=2,
                    expected_columns=2,
                )


if __name__ == "__main__":
    unittest.main()
