"""Offline regression tests for the Step 4.1 source registry audit."""
from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

from shapely.geometry import Polygon, mapping

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src/data_sources/01_audit_step41_sources.py"
SPEC = importlib.util.spec_from_file_location("step41_audit", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class RegistryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.csv_file = self.folder / "sources.csv"
        self.rows = [
            {col: "example" for col in MODULE.EXPECTED_FIELDS},
            {col: "example" for col in MODULE.EXPECTED_FIELDS},
        ]
        for i, row in enumerate(self.rows):
            row.update(source_id=f"s{i}", priority="P0", highway99_coverage_verified="NO",
                       dataset_url=f"https://data.example.org/dataset/{i}")
        self.write_csv()

    def write_csv(self):
        with self.csv_file.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=MODULE.EXPECTED_FIELDS)
            writer.writeheader()
            writer.writerows(self.rows)

    def test_good_registry(self):
        self.assertEqual(len(MODULE.read_and_validate_catalog(self.csv_file)), 2)

    def test_duplicate_ids_rejected(self):
        self.rows[1]["source_id"] = "s0"
        self.write_csv()
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            MODULE.read_and_validate_catalog(self.csv_file)

    def test_invalid_url_rejected(self):
        self.rows[0]["dataset_url"] = "http://example.org"
        self.write_csv()
        with self.assertRaisesRegex(ValueError, "Invalid source URL"):
            MODULE.read_and_validate_catalog(self.csv_file)

    def test_invalid_priority_rejected(self):
        self.rows[0]["priority"] = "P9"
        self.write_csv()
        with self.assertRaisesRegex(ValueError, "priority"):
            MODULE.read_and_validate_catalog(self.csv_file)

    def test_missing_field_rejected(self):
        self.rows[0]["dataset_name"] = ""
        self.write_csv()
        with self.assertRaisesRegex(ValueError, "Missing"):
            MODULE.read_and_validate_catalog(self.csv_file)

    def test_invalid_coverage_flag_rejected(self):
        self.rows[0]["highway99_coverage_verified"] = "MAYBE"
        self.write_csv()
        with self.assertRaisesRegex(ValueError, "coverage"):
            MODULE.read_and_validate_catalog(self.csv_file)

    def test_wrong_columns_rejected(self):
        self.csv_file.write_text("x,y\n1,2\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "columns"):
            MODULE.read_and_validate_catalog(self.csv_file)

    def create_aoi(self, geometry=None):
        geometry = geometry or Polygon([
            (-123.2, 49.7), (-123.19, 49.7), (-123.19, 49.71), (-123.2, 49.71), (-123.2, 49.7)
        ])
        file = self.folder / "aoi.geojson"
        collection = {
            "type": "FeatureCollection",
            "features": [{"type": "Feature", "properties": {}, "geometry": mapping(geometry)}],
        }
        file.write_text(json.dumps(collection), encoding="utf-8")
        return file

    def test_aoi_is_valid(self):
        outcome = MODULE.audit_aoi(self.create_aoi())
        self.assertGreater(outcome["aoi_unique_area_km2_metric_epsg26910"], 0)
        self.assertEqual(len(outcome["aoi_bounds_wgs84_west_south_east_north"]), 4)

    def test_aoi_expected_area_success(self):
        file = self.create_aoi()
        area = MODULE.audit_aoi(file)["aoi_unique_area_km2_metric_epsg26910"]
        qa = self.folder / "qa.json"
        qa.write_text(json.dumps({"unique_terrain_data_search_area_km2": area}), encoding="utf-8")
        self.assertTrue(MODULE.audit_aoi(file, qa)["step34_area_comparison"]["matches_0_01_km2_tolerance"])

    def test_aoi_expected_area_mismatch_rejected(self):
        file = self.create_aoi()
        qa = self.folder / "qa.json"
        qa.write_text(json.dumps({"unique_terrain_data_search_area_km2": 0.01}), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            MODULE.audit_aoi(file, qa)

    def test_produce_audit_marks_coverage_unverified(self):
        file = self.create_aoi()
        outcome = MODULE.produce_audit(self.csv_file, file, None)
        self.assertEqual(outcome["sources_with_local_corridor_coverage_verified"], 0)
        self.assertIn("COVERAGE_NOT_VERIFIED", outcome["status"])


if __name__ == "__main__":
    unittest.main()
