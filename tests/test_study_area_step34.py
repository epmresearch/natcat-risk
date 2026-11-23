"""Offline tests for Step 3.4 GIS screening corridor preparation."""

from __future__ import annotations

import importlib.util
import json
import math
import tempfile
import unittest
from pathlib import Path

import geopandas as gpd
from shapely.geometry import LineString, Point

MODULE_FILE = Path(__file__).resolve().parents[1] / "src/study_area/04_prepare_spatial_corridors.py"
spec = importlib.util.spec_from_file_location("step34_corridor", MODULE_FILE)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class CorridorTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.input = self.root / "segments.geojson"
        self.qa = self.root / "step33_quality.json"
        # Two 1km adjoining projected centreline segments (UTM zone 10N).
        metric = gpd.GeoDataFrame([
            {"segment_id": "H99-0001", "segment_no": 1, "start_m": 0.0, "end_m": 1000.0,
             "length_m": 1000.0, "geometry": LineString([(480000, 5500000), (480000, 5501000)])},
            {"segment_id": "H99-0002", "segment_no": 2, "start_m": 1000.0, "end_m": 2000.0,
             "length_m": 1000.0, "geometry": LineString([(480000, 5501000), (480000, 5502000)])},
        ], crs="EPSG:26910")
        metric.to_crs("EPSG:4326").to_file(self.input, driver="GeoJSON")
        self.qa.write_text(json.dumps({"step": "3.3", "segment_count": 2,
                                       "sum_of_segment_lengths_m": 2000.0,
                                       "status": "PROVISIONAL_SEGMENTATION_REQUIRES_REVIEW"}))

    def test_valid_input(self):
        gdf, meta = module.load_segments(self.input, self.qa)
        self.assertEqual(len(gdf), 2)
        self.assertEqual(meta["step"], "3.3")

    def test_area_union_less_than_per_segment_sum(self):
        gdf, _ = module.load_segments(self.input, self.qa)
        pieces, near, terrain, values = module.prepare_footprints(gdf, 1000, 5000)
        self.assertEqual(len(pieces), 2)
        self.assertLess(near.geometry.iloc[0].area, pieces.geometry.area.sum())
        self.assertGreater(values["overlapping_area_counted_more_than_once_km2"], 0)
        self.assertLess(near.geometry.iloc[0].area, terrain.geometry.iloc[0].area)

    def test_all_line_segments_inside_near_zone(self):
        gdf, _ = module.load_segments(self.input, self.qa)
        _, near, _, _ = module.prepare_footprints(gdf, 1000, 5000)
        for line in gdf.geometry:
            self.assertTrue(near.geometry.iloc[0].covers(line))

    def test_segment_ids_and_attributes_retained(self):
        gdf, _ = module.load_segments(self.input, self.qa)
        pieces, _, _, _ = module.prepare_footprints(gdf, 1000, 5000)
        self.assertEqual(list(pieces.segment_id), ["H99-0001", "H99-0002"])
        self.assertEqual(list(pieces.start_m), [0, 1000])

    def test_invalid_buffer_widths_rejected(self):
        gdf, _ = module.load_segments(self.input, self.qa)
        for near, search in [(0, 5000), (-1, 5000), (1000, 1000), (5000, 1000), (True, 5000)]:
            with self.subTest(near=near, search=search), self.assertRaises(ValueError):
                module.prepare_footprints(gdf, near, search)

    def test_duplicate_ids_rejected(self):
        original = gpd.read_file(self.input)
        original.loc[1, "segment_id"] = "H99-0001"
        original.to_file(self.input, driver="GeoJSON")
        with self.assertRaises(ValueError):
            module.load_segments(self.input, self.qa)

    def test_disconnected_segments_rejected(self):
        g = gpd.read_file(self.input).to_crs("EPSG:26910")
        g.loc[1, "geometry"] = LineString([(481000, 5501000), (481000, 5502000)])
        g.to_crs("EPSG:4326").to_file(self.input, driver="GeoJSON")
        with self.assertRaises(ValueError):
            module.load_segments(self.input, self.qa)

    def test_chainage_discrepancy_rejected(self):
        g = gpd.read_file(self.input)
        g.loc[1, "start_m"] = 999
        g.to_file(self.input, driver="GeoJSON")
        with self.assertRaises(ValueError):
            module.load_segments(self.input, self.qa)

    def test_quality_count_mismatch_rejected(self):
        self.qa.write_text(json.dumps({"step": "3.3", "segment_count": 9,
                                       "sum_of_segment_lengths_m": 2000.0}))
        with self.assertRaises(ValueError):
            module.load_segments(self.input, self.qa)

    def test_geojson_exports_and_qa_report(self):
        report = module.build_outputs(self.input, self.qa, self.root / "out", self.root / "maps",
                                      1000, 5000, make_map=False)
        self.assertEqual(report["segment_count"], 2)
        self.assertTrue(report["qa_flags"]["exports_reopened_and_validated"])
        self.assertEqual(report["status"], "PROVISIONAL_SPATIAL_SCREENING_NOT_A_HAZARD_MODEL")
        self.assertTrue((self.root / "out/highway99_terrain_search_5000m.geojson").is_file())
        g = gpd.read_file(self.root / "out/highway99_segment_screening_zones_1000m.geojson")
        self.assertEqual(len(g), 2)
        self.assertEqual(g.crs.to_epsg(), 4326)

    def test_default_source_present_in_repo_if_available(self):
        self.assertEqual(module.DEFAULT_NEAR_M, 1000)
        self.assertEqual(module.DEFAULT_SEARCH_M, 5000)


if __name__ == "__main__":
    unittest.main()
