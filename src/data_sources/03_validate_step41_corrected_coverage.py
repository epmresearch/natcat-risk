"""Step 4.1: Independent, read-only QA of corrected HRDEM vector coverage.

Dependencies already present in sea_to_sky_risk: geopandas, shapely.
Place at src/data_sources/03_validate_step41_corrected_coverage.py.

This examines *coverage polygons*, not elevation raster pixels or NoData.
It does not modify any GeoPackage, the QGIS project, or older QA reports.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import geopandas as gpd
from shapely.ops import unary_union

TARGET_EPSG = 26910
# Allow a few square meters for observed coordinate/overlay precision differences.
# These tolerances are far below one square kilometer; all deviations are reported.
OUTSIDE_TOLERANCE_M2 = 10.0
UNION_DIFFERENCE_TOLERANCE_M2 = 20.0
ATTRIBUTE_AREA_TOLERANCE_M2 = 2.0  # stored values rounded to 6 decimal km2
ATTRIBUTE_PERCENT_TOLERANCE = 0.0002  # percentage points

SOURCE_PATHS = {
    'study_aoi': 'data/processed/study_area/highway99_terrain_search_5000m.geojson',
    'hrdem_1m': 'data/processed/elevation/hrdem_1m_dtm_coverage_highway99.gpkg',
    'hrdem_2m_corrected': 'data/processed/elevation/hrdem_2m_dtm_coverage_highway99_corrected.gpkg',
    'combined_corrected': 'data/processed/elevation/hrdem_combined_dtm_coverage_highway99_corrected.gpkg',
}


def load_polygon(name: str, path: Path):
    if not path.is_file():
        raise FileNotFoundError(f'{name}: file not found: {path}')
    df = gpd.read_file(path)
    if df.empty:
        raise ValueError(f'{name}: file has no features')
    if df.crs is None:
        raise ValueError(f'{name}: file has no declared CRS')
    if not df.geometry.notna().all() or df.geometry.is_empty.any():
        raise ValueError(f'{name}: null or empty geometry')
    bad = df.loc[~df.geometry.is_valid]
    if not bad.empty:
        raise ValueError(f'{name}: {len(bad)} invalid geometry feature(s)')
    if not df.geometry.geom_type.isin(['Polygon', 'MultiPolygon']).all():
        raise ValueError(f'{name}: nonpolygon geometry encountered')
    metric = df.to_crs(epsg=TARGET_EPSG)
    geometry = unary_union(metric.geometry.tolist())
    if not geometry.is_valid:
        raise ValueError(f'{name}: union geometry is invalid')
    return df, metric, geometry


def execute(root: Path):
    root = root.resolve()
    report = {
        'audit': 'Step 4.1 independent corrected HRDEM vector-coverage QA',
        'timestamp_utc': datetime.now(timezone.utc).isoformat(),
        'scope': 'vector coverage polygons only; not raster valid-pixel coverage, elevation accuracy or landslide hazard',
        'working_crs_epsg': TARGET_EPSG,
        'thresholds_m2': {
            'outside_aoi': OUTSIDE_TOLERANCE_M2,
            'union_symmetric_difference': UNION_DIFFERENCE_TOLERANCE_M2,
            'saved_area_field': ATTRIBUTE_AREA_TOLERANCE_M2,
        },
        'sources': {},
        'measurements': {},
        'checks': [],
        'status': 'FAILED',
    }

    def check(label, passes, actual=None):
        rec = {'check': label, 'passed': bool(passes)}
        if actual is not None:
            rec['actual'] = actual
        report['checks'].append(rec)

    source = {}
    for key, rel in SOURCE_PATHS.items():
        path = root / rel
        original, metric, shape = load_polygon(key, path)
        source[key] = (original, metric, shape)
        report['sources'][key] = {
            'path': rel,
            'feature_count': len(original),
            'declared_crs': str(original.crs),
            'area_km2_in_epsg26910': round(shape.area / 1e6, 9),
        }
        check(f'{key}: valid polygon geometry and declared CRS', True)

    aoi = source['study_aoi'][2]
    one = source['hrdem_1m'][2]
    two = source['hrdem_2m_corrected'][2]
    saved = source['combined_corrected'][2]
    independent = unary_union([one, two])

    overlap_m2 = one.intersection(two).area
    outside_one_m2 = one.difference(aoi).area
    outside_two_m2 = two.difference(aoi).area
    outside_saved_m2 = saved.difference(aoi).area
    symdiff_m2 = saved.symmetric_difference(independent).area
    independent_inside_m2 = independent.intersection(aoi).area

    m = report['measurements']
    m.update({
        'aoi_km2': round(aoi.area / 1e6, 9),
        'one_m_km2': round(one.area / 1e6, 9),
        'two_m_corrected_km2': round(two.area / 1e6, 9),
        'one_m_two_m_overlap_km2': round(overlap_m2 / 1e6, 9),
        'independent_union_km2': round(independent.area / 1e6, 9),
        'independent_union_inside_aoi_km2': round(independent_inside_m2 / 1e6, 9),
        'saved_combined_km2': round(saved.area / 1e6, 9),
        'saved_combined_pct_of_aoi': round(100 * saved.area / aoi.area, 6),
        'saved_combined_outside_aoi_m2': round(outside_saved_m2, 6),
        'one_m_outside_aoi_m2': round(outside_one_m2, 6),
        'two_m_outside_aoi_m2': round(outside_two_m2, 6),
        'saved_vs_independent_union_symmetric_difference_m2': round(symdiff_m2, 6),
        'uncovered_aoi_km2': round(aoi.difference(saved).area / 1e6, 9),
    })

    check('corrected 2m is within AOI (<=10 m2 outside)', outside_two_m2 <= OUTSIDE_TOLERANCE_M2, round(outside_two_m2, 6))
    check('saved combined is within AOI (<=10 m2 outside)', outside_saved_m2 <= OUTSIDE_TOLERANCE_M2, round(outside_saved_m2, 6))
    check('original 1m has no material area outside AOI (<=10 m2)', outside_one_m2 <= OUTSIDE_TOLERANCE_M2, round(outside_one_m2, 6))
    check('saved combined matches independent union (<=20 m2 symmetric difference)', symdiff_m2 <= UNION_DIFFERENCE_TOLERANCE_M2, round(symdiff_m2, 6))
    check('saved combined area is <= AOI area', saved.area <= aoi.area + OUTSIDE_TOLERANCE_M2)

    saved_df = source['combined_corrected'][0]
    check('combined has exactly one feature', len(saved_df) == 1, len(saved_df))
    if len(saved_df) == 1:
        row = saved_df.iloc[0]
        for col, expected, tol in (
            ('combined_km2', saved.area / 1e6, ATTRIBUTE_AREA_TOLERANCE_M2 / 1e6),
            ('combined_pct', 100 * saved.area / aoi.area, ATTRIBUTE_PERCENT_TOLERANCE),
        ):
            if col not in saved_df.columns:
                check(f'combined attribute {col} exists', False, 'missing')
                continue
            try:
                val = float(row[col])
                delta = abs(val - expected)
                check(f'combined attribute {col} matches geometry', delta <= tol, {'saved': val, 'recomputed': round(expected, 9), 'delta': round(delta, 9)})
            except (TypeError, ValueError):
                check(f'combined attribute {col} is numeric', False, str(row[col]))

    report['status'] = 'PASSED_INTERNAL_VECTOR_COVERAGE_QA_WITH_CAVEATS' if all(c['passed'] for c in report['checks']) else 'FAILED_INTERNAL_VECTOR_COVERAGE_QA'
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[2], help='natcat-risk repository root')
    args = parser.parse_args()
    root = args.root.resolve()
    output = root / 'results/diagnostics/step4_1_corrected_coverage_qa.json'
    try:
        report = execute(root)
    except Exception as exc:
        print(f'QA ERROR: {exc}', file=sys.stderr)
        print('Check the four input GeoPackage/GeoJSON paths; no originals were changed.', file=sys.stderr)
        return 2
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print('\nSTEP 4.1 — CORRECTED HRDEM VECTOR COVERAGE QA')
    print('Status:', report['status'])
    for key, value in report['measurements'].items():
        print(f'  {key}: {value}')
    print('\nCHECKS:')
    for entry in report['checks']:
        print(f"  {'PASS' if entry['passed'] else 'FAIL'}: {entry['check']}" + (f" ({entry['actual']})" if 'actual' in entry else ''))
    print('\nReport:', output)
    print('CAVEAT: coverage polygons are not proof of valid DEM pixels or elevation accuracy.')
    return 0 if report['status'].startswith('PASSED') else 1


if __name__ == '__main__':
    sys.exit(main())
