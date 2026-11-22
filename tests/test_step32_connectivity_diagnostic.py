"""Offline tests for the Step 3.2 diagnostic functions without OSMnx import."""
import importlib.util
from pathlib import Path
import unittest

import networkx as nx
from pyproj import Transformer

SOURCE = Path(__file__).resolve().parents[1] / "src/study_area/02c_diagnose_highway99_connectivity.py"
spec = importlib.util.spec_from_file_location("step32_connectivity", SOURCE)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class ConnectivityDiagnosticTests(unittest.TestCase):
    def setUp(self):
        self.graph = nx.MultiDiGraph()
        self.graph.add_node(1, x=-123.1, y=49.5)
        self.graph.add_node(2, x=-123.1, y=49.51)
        self.graph.add_node(3, x=-123.1, y=49.511)
        self.graph.add_node(4, x=-123.1, y=49.52)
        self.graph.add_edge(1, 2)
        self.graph.add_edge(3, 4)
        self.tf = Transformer.from_crs("EPSG:4326", "EPSG:26910", always_xy=True)

    def test_distinct_components(self):
        comps, mapping = module.component_index(self.graph)
        self.assertEqual(len(comps), 2)
        self.assertNotEqual(mapping[1], mapping[4])

    def test_nearest_pair_gap(self):
        result = module.closest_distinct_component_nodes(self.graph, [1, 2], [3, 4], self.tf)
        self.assertEqual(result["first_node"], 2)
        self.assertEqual(result["second_node"], 3)
        self.assertGreater(result["distance_m"], 100)
        self.assertLess(result["distance_m"], 130)

    def test_no_graph_mutation(self):
        count = self.graph.number_of_edges()
        module.closest_distinct_component_nodes(self.graph, [1, 2], [3, 4], self.tf)
        self.assertEqual(self.graph.number_of_edges(), count)


if __name__ == "__main__":
    unittest.main()
