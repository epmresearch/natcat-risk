"""Offline tests for merging a local OSM patch without artificial links."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import networkx as nx

SCRIPT = (Path(__file__).resolve().parents[1] / "src/study_area/02d_repair_osm_gap.py")
spec = importlib.util.spec_from_file_location("osm_gap_repair", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def node(g, n, lon, lat):
    g.add_node(n, x=lon, y=lat)


class GapPatchTests(unittest.TestCase):
    def setUp(self):
        self.base = nx.MultiDiGraph()
        node(self.base, 1, -123.125, 50.04)
        node(self.base, 2, -123.123, 50.041)
        node(self.base, 3, -123.120, 50.05)
        node(self.base, 4, -123.119, 50.051)
        self.base.add_edge(1, 2, length=100)
        self.base.add_edge(3, 4, length=100)

    def test_base_disconnected(self):
        self.assertFalse(module.anchors_connected(self.base, 1, 4))

    def test_real_shared_id_patch_connects(self):
        p = nx.MultiDiGraph()
        node(p, 2, -123.123, 50.041)
        node(p, 5, -123.122, 50.045)
        node(p, 3, -123.120, 50.05)
        p.add_edge(2, 5, length=70)
        p.add_edge(5, 3, length=70)
        combined, shared = module.merge_original_osm_nodes(self.base, p)
        self.assertTrue(module.anchors_connected(combined, 1, 4))
        self.assertEqual(shared, 2)
        self.assertEqual(combined.number_of_nodes(), 5)

    def test_no_synthetic_bridge_when_patch_does_not_overlap(self):
        p = nx.MultiDiGraph()
        node(p, 6, -123.122, 50.045)
        node(p, 7, -123.121, 50.047)
        p.add_edge(6, 7)
        with self.assertRaisesRegex(RuntimeError, "no shared node"):
            module.merge_original_osm_nodes(self.base, p)

    def test_overlap_only_one_side_does_not_connect(self):
        p = nx.MultiDiGraph()
        node(p, 2, -123.123, 50.041)
        node(p, 5, -123.122, 50.045)
        p.add_edge(2, 5)
        combined, _ = module.merge_original_osm_nodes(self.base, p)
        self.assertFalse(module.anchors_connected(combined, 1, 4))

    def test_reject_mismatched_same_id_coords(self):
        p = nx.MultiDiGraph()
        node(p, 2, -122.999, 50.041)
        node(p, 3, -123.120, 50.05)
        p.add_edge(2, 3)
        with self.assertRaisesRegex(RuntimeError, "inconsistent coordinates"):
            module.merge_original_osm_nodes(self.base, p)

    def test_polygon_covers_both_ends(self):
        from shapely.geometry import Point
        poly = module.derive_gap_polygon(self.base, 2, 3, 2000)
        self.assertTrue(poly.contains(Point(-123.123, 50.041)))
        self.assertTrue(poly.contains(Point(-123.120, 50.05)))

    def test_not_retry_on_graph_issue(self):
        self.assertFalse(module.is_transient(RuntimeError("No road nodes found")))


if __name__ == "__main__":
    unittest.main()
