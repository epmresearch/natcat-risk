"""Offline regression tests for Step 3.5 cross-stage GIS and provenance audit."""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

import geopandas as gpd
from shapely.geometry import LineString

SCRIPT = Path(__file__).resolve().parents[1] / "src/study_area/05_validate_study_area.py"
spec = importlib.util.spec_from_file_location("study_area_step35", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def make_fixture(root: Path) -> Path:
    """Generate portable tiny 3-km equivalents of the 3.2/3.3/3.4 artifacts."""
    target = root / mod.DATA_DIR
    target.mkdir(parents=True, exist_ok=True)
    lines = [LineString([(480000 + i * 1000, 5500000), (481000 + i * 1000, 5500000)])
             for i in range(3)]
    center = gpd.GeoDataFrame([{"status": "PROVISIONAL"}],
                              geometry=[LineString([(480000, 5500000), (483000, 5500000)])],
                              crs=mod.METRIC_CRS)
    center.to_crs(mod.GEO_CRS).to_file(target / "highway99_candidate_centerline.geojson", driver="GeoJSON")
    write_json(target / "highway99_candidate_quality.json", {"step": "3.2", "candidate_length_km": 3.0,
                                                               "status": "PROVISIONAL_REQUIRES_MANUAL_VALIDATION"})
    segs = gpd.GeoDataFrame([
        {"segment_id": f"H99-{i+1:04d}", "segment_no": i + 1,
         "start_m": 1000*i, "end_m": 1000*(i+1), "length_m": 1000.0}
        for i in range(3)], geometry=lines, crs=mod.METRIC_CRS)
    segs.to_crs(mod.GEO_CRS).to_file(target / "highway99_segments_1000m.geojson", driver="GeoJSON")
    segs.drop(columns="geometry").to_csv(target / "highway99_segments_1000m.csv", index=False)
    geo = target / "highway99_segments_1000m.geojson"
    s32 = target / "highway99_candidate_centerline.geojson"
    write_json(target / "highway99_segments_1000m_quality.json", {
        "step": "3.3", "segment_count": 3, "source_length_m": 3000,
        "requested_segment_length_m": 1000,
        "source_sha256": hashlib.sha256(s32.read_bytes()).hexdigest(),
    })
    zones = segs[["segment_id", "segment_no"]].copy()
    zones = gpd.GeoDataFrame(zones, geometry=segs.geometry.buffer(1000), crs=mod.METRIC_CRS)
    near_geom = zones.geometry.union_all()
    search_geom = segs.geometry.union_all().buffer(5000)
    zones.to_crs(mod.GEO_CRS).to_file(target / "highway99_segment_screening_zones_1000m.geojson", driver="GeoJSON")
    gpd.GeoDataFrame([{}], geometry=[near_geom], crs=mod.METRIC_CRS).to_crs(mod.GEO_CRS).to_file(
        target / "highway99_corridor_screening_1000m.geojson", driver="GeoJSON")
    gpd.GeoDataFrame([{}], geometry=[search_geom], crs=mod.METRIC_CRS).to_crs(mod.GEO_CRS).to_file(
        target / "highway99_terrain_search_5000m.geojson", driver="GeoJSON")
    sum_zone = float(zones.geometry.area.sum() / 1e6)
    near_area = float(near_geom.area / 1e6)
    search_area = float(search_geom.area / 1e6)
    write_json(target / "highway99_spatial_preparation_quality.json", {
        "step": "3.4", "segment_count": 3,
        "input_step33_sha256": hashlib.sha256(geo.read_bytes()).hexdigest(),
        "near_road_screening_distance_m": 1000,
        "terrain_data_search_distance_m": 5000,
        "terrain_search_polygon_part_count": 1,
        "terrain_search_smallest_part_area_m2": round(search_geom.area, 3),
        "sum_individual_segment_zone_areas_km2": round(sum_zone, 6),
        "unique_near_road_corridor_area_km2": round(near_area, 6),
        "unique_terrain_data_search_area_km2": round(search_area, 6),
        "overlapping_area_counted_more_than_once_km2": round(sum_zone - near_area, 6),
        "qa_flags": {"segment_ids_preserved": True, "near_corridor_inside_search_corridor": True,
                     "no_route_geometry_changes": True, "exports_reopened_and_validated": True},
    })
    return target


class TestStudyAreaStep35(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.data = make_fixture(self.root)

    def check_fails(self, message: str):
        with self.assertRaisesRegex(mod.ValidationError, message):
            mod.validate_inputs(self.root)

    def rewrite_geojson(self, name: str, change):
        path = self.data / name
        payload = json.loads(path.read_text(encoding="utf-8"))
        change(payload)
        write_json(path, payload)

    def test_normal_fixture_passes_with_caveats(self):
        report, _ = mod.validate_inputs(self.root)
        self.assertEqual(report["status"], "INTERNAL_QA_PASSED_WITH_CAVEATS")
        self.assertEqual(report["measurements"]["segment_count"], 3)
        self.assertEqual(report["measurements"]["terrain_search_polygon_parts"], 1)
        self.assertFalse(report["independent_official_road_validation_completed"])
        self.assertEqual(len(report["checks_passed"]), 8)

    def test_run_writes_json_without_local_absolute_paths(self):
        report = mod.run(self.root, make_map=False)
        path = self.root / mod.QA_OUT
        self.assertTrue(path.exists())
        txt = path.read_text()
        self.assertNotIn(self.tmp.name, txt)
        self.assertEqual(json.loads(txt)["status"], report["status"])

    def test_map_generation(self):
        mod.run(self.root, make_map=True)
        self.assertTrue((self.root / mod.MAP_OUT).stat().st_size > 5000)

    def test_missing_file_blocks_validation(self):
        (self.data / "highway99_terrain_search_5000m.geojson").unlink()
        self.check_fails("Missing required files")

    def test_source_sha_mismatch_blocks_validation(self):
        path = self.data / "highway99_segments_1000m_quality.json"
        qa = json.loads(path.read_text())
        qa["source_sha256"] = "0" * 64
        write_json(path, qa)
        self.check_fails("source_sha256 differs")

    def test_step34_sha_mismatch_blocks_validation(self):
        path = self.data / "highway99_spatial_preparation_quality.json"
        qa = json.loads(path.read_text())
        qa["input_step33_sha256"] = "0" * 64
        write_json(path, qa)
        self.check_fails("input_step33_sha256 differs")

    def test_nonunique_segment_ids_block_validation(self):
        self.rewrite_geojson("highway99_segments_1000m.geojson", lambda j: j["features"][1]["properties"].update(
            {"segment_id": "H99-0001"}))
        # Keep provenance consistent to test the ID gate, not only the SHA gate.
        self.sync_step34_hash()
        self.check_fails("not unique")

    def sync_step34_hash(self):
        path = self.data / "highway99_spatial_preparation_quality.json"
        qa = json.loads(path.read_text())
        qa["input_step33_sha256"] = hashlib.sha256((self.data / "highway99_segments_1000m.geojson").read_bytes()).hexdigest()
        write_json(path, qa)

    def test_wrong_chainage_blocks_validation(self):
        self.rewrite_geojson("highway99_segments_1000m.geojson", lambda j: j["features"][1]["properties"].update(
            {"start_m": 1800}))
        self.sync_step34_hash()
        self.check_fails("chainage interval")

    def test_csv_out_of_order_blocks_validation(self):
        path = self.data / "highway99_segments_1000m.csv"
        with path.open(newline="", encoding="utf-8") as f:
            records = list(csv.DictReader(f))
        with path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(list(reversed(records)))
        self.check_fails("CSV segment order")

    def test_modified_per_segment_polygon_blocks_validation(self):
        def change(j):
            ring = j["features"][0]["geometry"]["coordinates"][0]
            for xy in ring:
                xy[0] += 0.003
        self.rewrite_geojson("highway99_segment_screening_zones_1000m.geojson", change)
        self.check_fails("zone differs from reproducible buffer")

    def test_modified_search_polygon_blocks_validation(self):
        def change(j):
            ring = j["features"][0]["geometry"]["coordinates"][0]
            for xy in ring:
                xy[0] += 0.005
        self.rewrite_geojson("highway99_terrain_search_5000m.geojson", change)
        self.check_fails("buffers are not reproducible|Corridor buffers are not reproducible")

    def test_wrong_qa_area_blocks_validation(self):
        path = self.data / "highway99_spatial_preparation_quality.json"
        qa = json.loads(path.read_text())
        qa["unique_near_road_corridor_area_km2"] += 10
        write_json(path, qa)
        self.check_fails("unique_near_road_corridor_area_km2")

    def test_false_old_qa_flag_blocks_validation(self):
        path = self.data / "highway99_spatial_preparation_quality.json"
        qa = json.loads(path.read_text())
        qa["qa_flags"]["segment_ids_preserved"] = False
        write_json(path, qa)
        self.check_fails("QA flags")

    def test_wrong_segment_length_scheme_blocks_validation(self):
        path = self.data / "highway99_segments_1000m_quality.json"
        qa = json.loads(path.read_text())
        qa["requested_segment_length_m"] = 500
        write_json(path, qa)
        self.check_fails("1000 m segmentation")

    def test_wrong_screening_distance_blocks_validation(self):
        path = self.data / "highway99_spatial_preparation_quality.json"
        qa = json.loads(path.read_text())
        qa["terrain_data_search_distance_m"] = 4500
        write_json(path, qa)
        self.check_fails("1000 m and 5000 m")


if __name__ == "__main__":
    unittest.main()
