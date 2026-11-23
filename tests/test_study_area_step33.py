"""Offline tests for Step 3.3 road segmentation. No OSM/Overpass access."""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

import geopandas as gpd
from shapely.geometry import LineString, MultiLineString, Point

SCRIPT = Path(__file__).resolve().parents[1] / "src/study_area/03_segment_highway99.py"
spec = importlib.util.spec_from_file_location("highway99_step33", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class TestSplitLine(unittest.TestCase):
    def test_2500_m_straight_line_three_segments(self):
        line = LineString([(0, 0), (2500, 0)])
        parts = module.split_line_by_chainage(line, 1000)
        self.assertEqual([p["segment_id"] for p in parts],
                         ["H99-0001", "H99-0002", "H99-0003"])
        self.assertEqual([round(p["geometry"].length) for p in parts], [1000, 1000, 500])
        self.assertEqual(parts[-1]["end_m"], 2500)

    def test_exact_multiple_does_not_create_extra_segment(self):
        parts = module.split_line_by_chainage(LineString([(0, 0), (3000, 0)]), 1000)
        self.assertEqual(len(parts), 3)
        self.assertEqual(parts[-1]["length_m"], 1000)

    def test_line_shorter_than_one_segment(self):
        parts = module.split_line_by_chainage(LineString([(0, 0), (400, 0)]), 1000)
        self.assertEqual(len(parts), 1)
        self.assertEqual(parts[0]["length_m"], 400)

    def test_curved_polyline_keeps_original_turns(self):
        line = LineString([(0, 0), (500, 0), (500, 800), (1000, 800)])
        parts = module.split_line_by_chainage(line, 1000)
        self.assertEqual(len(parts), 2)
        self.assertTrue(Point(parts[0]["geometry"].coords[-1]).equals(
            Point(parts[1]["geometry"].coords[0])))
        self.assertAlmostEqual(sum(p["geometry"].length for p in parts), line.length, places=6)
        self.assertIn((500, 0), list(parts[0]["geometry"].coords))
        self.assertEqual(list(parts[0]["geometry"].coords)[0], (0, 0))
        self.assertEqual(list(parts[-1]["geometry"].coords)[-1], (1000, 800))

    def test_invalid_length_rejected(self):
        for value in (0, -1, 0.5, True):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    module.split_line_by_chainage(LineString([(0, 0), (10, 0)]), value)

    def test_multiline_rejected(self):
        with self.assertRaisesRegex(ValueError, "LineString"):
            module.split_line_by_chainage(
                MultiLineString([[(0, 0), (100, 0)], [(100, 0), (200, 0)]]), 1000)

    def test_empty_line_rejected(self):
        with self.assertRaises(ValueError):
            module.split_line_by_chainage(LineString(), 1000)

    def test_all_segments_max_one_kilometre(self):
        line = LineString([(0, 0), (0, 50), (50, 50), (50, 2733)])
        parts = module.split_line_by_chainage(line, 1000)
        self.assertTrue(all(x["geometry"].length <= 1000.00001 for x in parts))


class TestIntegration(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.src = self.root / "centerline.geojson"
        self.qa = self.root / "source_qa.json"
        # Synthetic line in UTM zone 10N near coastal British Columbia.
        self.metric_line = LineString([(490000, 5480000), (490000, 5481900),
                                       (490300, 5481900), (490300, 5482700)])
        gdf = gpd.GeoDataFrame({"status": ["PROVISIONAL_NOT_AUTHORITATIVE"]},
                               geometry=[self.metric_line], crs="EPSG:26910")
        gdf.to_crs("EPSG:4326").to_file(self.src, driver="GeoJSON")
        self.qa.write_text(json.dumps({
            "step": "3.2", "candidate_length_km": self.metric_line.length / 1000,
            "status": "PROVISIONAL_REQUIRES_MANUAL_VALIDATION",
        }), encoding="utf-8")

    def tearDown(self):
        self.tmp.cleanup()

    def test_input_reader_preserves_metric_length(self):
        _, geom = module.read_centerline(self.src)
        self.assertAlmostEqual(geom.length, 3000, places=3)

    def test_end_to_end_files_with_correct_qa(self):
        result = module.build_outputs(self.src, self.qa, self.root / "data",
                                      self.root / "maps", 1000, make_figure=False)
        self.assertEqual(result["segment_count"], 3)
        self.assertEqual(result["requested_segment_length_m"], 1000)
        self.assertEqual(result["max_boundary_gap_m"], 0)
        self.assertEqual(result["source_quality_status"],
                         "PROVISIONAL_REQUIRES_MANUAL_VALIDATION")
        path = self.root / "data/highway99_segments_1000m.geojson"
        self.assertTrue(path.exists())
        self.assertTrue((self.root / "data/highway99_segments_1000m.csv").exists())
        self.assertTrue((self.root / "data/highway99_segments_1000m_quality.json").exists())
        exported = gpd.read_file(path)
        self.assertEqual(len(exported), 3)
        self.assertEqual(exported["segment_id"].tolist(),
                         ["H99-0001", "H99-0002", "H99-0003"])
        self.assertAlmostEqual(exported.to_crs("EPSG:26910").length.sum(), 3000, places=3)

    def test_map_export(self):
        module.build_outputs(self.src, self.qa, self.root / "data",
                             self.root / "maps", 1000, make_figure=True)
        self.assertTrue((self.root / "maps/highway99_segments_1000m.png").exists())

    def test_mismatch_with_source_quality_rejected(self):
        self.qa.write_text(json.dumps({"step": "3.2", "candidate_length_km": 100.44}),
                           encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "differs"):
            module.read_source_quality(self.qa, 3000)

    def test_multiple_centerline_features_rejected(self):
        gdf = gpd.GeoDataFrame(
            geometry=[self.metric_line, self.metric_line], crs="EPSG:26910"
        )
        multiple = self.root / "multiple.geojson"
        gdf.to_crs("EPSG:4326").to_file(multiple, driver="GeoJSON")
        with self.assertRaisesRegex(ValueError, "exactly one"):
            module.read_centerline(multiple)

    def test_source_geometry_order_rejected_if_reversed(self):
        backward = self.root / "backward.geojson"
        gdf = gpd.GeoDataFrame(geometry=[LineString(list(self.metric_line.coords)[::-1])],
                               crs="EPSG:26910")
        gdf.to_crs("EPSG:4326").to_file(backward, driver="GeoJSON")
        with self.assertRaisesRegex(ValueError, "start is not south"):
            module.read_centerline(backward)



if __name__ == "__main__":
    unittest.main()
