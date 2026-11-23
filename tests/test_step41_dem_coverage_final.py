"""Small offline synthetic tests for independent Step 4.1 vector geometry audit."""
import importlib.util
import unittest
from pathlib import Path
from shapely.geometry import box

SCRIPT = Path(__file__).resolve().parents[1] / 'src' / 'data_sources' / '02_validate_dem_coverage.py'
spec = importlib.util.spec_from_file_location('step41_qa', SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

class TestStep41FinalQA(unittest.TestCase):
    def setUp(self):
        self.aoi = box(0, 0, 1000, 1000)

    def test_disjoint_case(self):
        a = box(0, 0, 500, 1000)
        b = box(500, 0, 1000, 1000)
        m = mod.geometry_metrics(self.aoi, a, b, a.union(b))
        self.assertAlmostEqual(m['overlap_km2'], 0)
        self.assertAlmostEqual(m['independent_union_km2'], 1)
        self.assertAlmostEqual(m['pct_combined'], 100)

    def test_overlap_case(self):
        a = box(0, 0, 700, 1000)
        b = box(300, 0, 1000, 1000)
        m = mod.geometry_metrics(self.aoi, a, b, a.union(b))
        self.assertAlmostEqual(m['overlap_km2'], 0.4)
        self.assertAlmostEqual(m['only_1m_km2'], 0.3)
        self.assertAlmostEqual(m['only_2m_km2'], 0.3)
        self.assertAlmostEqual(m['independent_union_km2'], 1.0)

    def test_detect_incorrect_combined_layer(self):
        a = box(0, 0, 700, 1000)
        b = box(300, 0, 1000, 1000)
        m = mod.geometry_metrics(self.aoi, a, b, a)
        self.assertAlmostEqual(m['independent_union_vs_saved_symmetric_difference_km2'], 0.3)

    def test_detect_coverage_outside_aoi(self):
        a = box(-100, 0, 500, 1000)
        b = box(500, 0, 1000, 1000)
        m = mod.geometry_metrics(self.aoi, a, b, a.union(b))
        self.assertAlmostEqual(m['coverage_outside_aoi_km2'], 0.1)

if __name__ == '__main__':
    unittest.main()
