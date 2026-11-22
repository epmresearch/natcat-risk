"""Offline unit test for the Highway 99 OSM-tag matching rule."""

import importlib.util
from pathlib import Path
import unittest
from types import SimpleNamespace
from unittest.mock import Mock, patch
import geopandas as gpd
from shapely.geometry import Point

SCRIPT = Path(__file__).resolve().parents[1] / "src/study_area/02_extract_highway99_candidate.py"
spec = importlib.util.spec_from_file_location("highway99_step32", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class TestHighway99Tags(unittest.TestCase):
    def test_simple_ref(self):
        self.assertTrue(module.highway99_tagged("99", None))

    def test_route_list(self):
        self.assertTrue(module.highway99_tagged(["BC 99", "X"], None))

    def test_named_highway(self):
        self.assertTrue(module.highway99_tagged(None, "Sea to Sky Highway"))

    def test_other_roads_not_matched(self):
        self.assertFalse(module.highway99_tagged("199", "Highway 1"))
        self.assertFalse(module.highway99_tagged("99A", None))


class TestDownloadTiling(unittest.TestCase):
    def test_tiles_cover_all_waypoints(self):
        # Offline synthetic coordinates, no Overpass or geocoding involved.
        pts = gpd.GeoDataFrame(
            {"name": module.WAYPOINTS},
            geometry=gpd.points_from_xy(
                [-123.273919, -123.235330, -123.206005, -123.155861, -122.954302],
                [49.374441, 49.458750, 49.626392, 49.698074, 50.117190],
            ),
            crs="EPSG:4326",
        )
        tiles = module.build_download_tiles(pts)
        self.assertGreater(len(tiles), len(module.WAYPOINTS) - 1)
        self.assertEqual(len({key for key, _ in tiles}), len(tiles))
        for pt in pts.geometry:
            self.assertTrue(any(poly.covers(pt) for _, poly in tiles))

    def test_overlapping_adjacent_parts(self):
        pts = gpd.GeoDataFrame(
            {"name": ["A", "B"]},
            geometry=[Point(-123.2, 49.5), Point(-122.95, 50.1)],
            crs="EPSG:4326",
        )
        tiles = module.build_download_tiles(pts)
        self.assertGreater(len(tiles), 1)
        for (_, poly1), (_, poly2) in zip(tiles[:-1], tiles[1:]):
            self.assertTrue(poly1.intersects(poly2))


class TestOverpassPacing(unittest.TestCase):
    def test_30_second_interval(self):
        self.assertEqual(module.seconds_until_next_request(None, 100), 0)
        self.assertEqual(module.seconds_until_next_request(100, 105), 25)
        self.assertEqual(module.seconds_until_next_request(100, 131), 0)

    def test_retry_502_is_transient(self):
        class ResponseStatusCodeError(Exception):
            pass
        self.assertTrue(module.is_retryable_overpass_failure(
            ResponseStatusCodeError("502 Bad Gateway # capacity full")
        ))

    def test_retry_429_is_transient(self):
        class ResponseStatusCodeError(Exception):
            pass
        self.assertTrue(module.is_retryable_overpass_failure(
            ResponseStatusCodeError("429 Too Many Requests")
        ))

    def test_do_not_retry_geometry_error(self):
        self.assertFalse(module.is_retryable_overpass_failure(ValueError("invalid polygon")))

    def test_bounded_retry_settings(self):
        self.assertEqual(module.MAX_ATTEMPTS_PER_TILE, 3)
        self.assertEqual(module.RETRY_WAIT_SECONDS, (30, 60))

    def test_transient_502_retries_then_succeeds(self):
        class ResponseStatusCodeError(Exception):
            pass
        returned_graph = object()
        download = Mock(side_effect=[
            ResponseStatusCodeError("502 Bad Gateway"), returned_graph
        ])
        with patch.object(module.ox, "graph", SimpleNamespace(graph_from_polygon=download), create=True), \
             patch.object(module.time, "monotonic", side_effect=[100.0, 100.0, 131.0, 131.0]), \
             patch.object(module.time, "sleep") as sleeper:
            graph, last_finished = module.download_tile_with_pacing(
                object(), "tile01", 1, 8, None
            )
        self.assertIs(graph, returned_graph)
        self.assertEqual(last_finished, 131.0)
        self.assertEqual(download.call_count, 2)
        sleeper.assert_called_once_with(30)

    def test_failed_download_stops_after_three_attempts(self):
        class ResponseStatusCodeError(Exception):
            pass
        download = Mock(side_effect=ResponseStatusCodeError("502 Bad Gateway"))
        with patch.object(module.ox, "graph", SimpleNamespace(graph_from_polygon=download), create=True), \
             patch.object(module.time, "monotonic", side_effect=[100, 100, 131, 131, 192, 192]), \
             patch.object(module.time, "sleep") as sleeper:
            with self.assertRaisesRegex(RuntimeError, "after 3 attempts"):
                module.download_tile_with_pacing(object(), "tile01", 1, 8, None)
        self.assertEqual(download.call_count, 3)
        self.assertEqual([call.args[0] for call in sleeper.call_args_list], [30, 60])


if __name__ == "__main__":
    unittest.main()
