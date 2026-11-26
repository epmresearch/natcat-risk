"""Read-only global structural and adjacent-tile seam audit for masked HRDEM 1m.

Run from the natcat-risk project root (or place this script in that root).
Uses ONLY local masked GeoTIFFs and local QA JSONs; no online access,
no files created, no terrain interpolation, no changes to original data.
This is diagnostic QA, NOT certification of vertical accuracy.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
from osgeo import gdal, osr

ROOT = Path(__file__).resolve().parent
MASK_DIR = ROOT / 'results' / 'diagnostics' / 'hrdem1m_masked_staging'
STAGE_DIR = ROOT / 'results' / 'diagnostics' / 'hrdem1m_batch_staging'
PATTERN = re.compile(r'hrdem1m_x(\d+)_y(\d+)_utm10_1km_masked\.tif$')
NODATA = -32767.0
EXPECTED = (1000, 1000)


def valid_values(data):
    return np.isfinite(data) & (data != NODATA)


def residual_at_edge(tile, neighbor, direction):
    """Signed second-difference residual crossing a 1m tile boundary."""
    if direction == 'east':
        a, av = tile['E']
        b, bv = neighbor['W']
        ok = av[:, 0] & av[:, 1] & bv[:, 0] & bv[:, 1]
        z = (b[:, 0] - a[:, 1]) - (
            (a[:, 1] - a[:, 0]) + (b[:, 1] - b[:, 0])
        ) / 2.0
    else:  # tile SOUTH, neighbor NORTH
        a, av = tile['N']
        b, bv = neighbor['S']
        ok = av[0, :] & av[1, :] & bv[0, :] & bv[1, :]
        z = (b[-1, :] - a[0, :]) - (
            (a[0, :] - a[1, :]) + (b[0, :] - b[-1, :])
        ) / 2.0
    return z[ok].astype(np.float64)


def main():
    gdal.UseExceptions()
    if not MASK_DIR.is_dir() or not STAGE_DIR.is_dir():
        raise FileNotFoundError('Missing masked/staging directory. Put the script in natcat-risk project root.')

    paths = sorted(MASK_DIR.glob('hrdem1m_x*_y*_utm10_1km_masked.tif'))
    stage_reports = list(STAGE_DIR.glob('hrdem1m_x*_y*_utm10_1km.qa.json'))
    if not paths:
        raise RuntimeError('No masked tiles found.')
    print('MASKED_TIFS_FOUND:', len(paths), flush=True)
    print('STAGING_QA_JSON_FOUND:', len(stage_reports), flush=True)

    edges = {}
    issues = []
    total_valid = 0
    total_missing_reported = 0
    for i, path in enumerate(paths, 1):
        match = PATTERN.fullmatch(path.name)
        if not match:
            continue
        x, y = map(int, match.groups())
        qa_path = path.with_suffix('.qa.json')
        if not qa_path.exists():
            issues.append(f'MISSING_QA {path.name}')
            continue
        qa = json.loads(qa_path.read_text(encoding='utf-8'))
        ds = gdal.Open(str(path), gdal.GA_ReadOnly)
        if ds is None:
            issues.append(f'UNREADABLE {path.name}')
            continue
        if (ds.RasterXSize, ds.RasterYSize) != EXPECTED or ds.RasterCount != 1:
            issues.append(f'WRONG_SIZE {path.name}')
            continue
        band = ds.GetRasterBand(1)
        if band.DataType != gdal.GDT_Float32 or band.GetNoDataValue() != NODATA:
            issues.append(f'WRONG_BAND_OR_NODATA {path.name}')
            continue
        expected_gt = (float(x), 1., 0., float(y + 1000), 0., -1.)
        if not np.allclose(ds.GetGeoTransform(), expected_gt, atol=1e-8, rtol=0):
            issues.append(f'WRONG_GEOTRANSFORM {path.name}')
            continue
        sr = osr.SpatialReference(wkt=ds.GetProjection())
        sr.AutoIdentifyEPSG()
        if sr.GetAuthorityCode(None) != '26910':
            issues.append(f'WRONG_CRS {path.name}')
            continue
        arr = band.ReadAsArray()
        mask = valid_values(arr)
        num_valid = int(mask.sum())
        total_valid += num_valid
        total_missing_reported += int(qa.get('nodata_inside_footprint', 0))
        if num_valid != qa.get('masked_valid_pixels'):
            issues.append(f'COUNT_MISMATCH {path.name}: raster={num_valid} qa={qa.get("masked_valid_pixels")}')
        if 1000000 - num_valid != qa.get('masked_nodata_pixels'):
            issues.append(f'NODATA_COUNT_MISMATCH {path.name}')
        if qa.get('status') != 'MASKED_STAGING_QA_PASSED_WITH_CAVEATS':
            issues.append(f'UNEXPECTED_QA_STATUS {path.name}: {qa.get("status")}')
        if (qa.get('tile_x'), qa.get('tile_y')) != (x, y):
            issues.append(f'QA_COORDINATE_MISMATCH {path.name}')
        edges[(x, y)] = {
            'N': (arr[:2, :].copy(), mask[:2, :].copy()),
            'S': (arr[-2:, :].copy(), mask[-2:, :].copy()),
            'W': (arr[:, :2].copy(), mask[:, :2].copy()),
            'E': (arr[:, -2:].copy(), mask[:, -2:].copy()),
        }
        ds = None
        if i % 50 == 0 or i == len(paths):
            print(f'RASTER_QA_PROGRESS: {i}/{len(paths)}', flush=True)

    results = []
    for (x, y), tile in sorted(edges.items()):
        for direction, nbr_xy in (('east', (x+1000, y)), ('north', (x, y+1000))):
            if nbr_xy not in edges:
                continue
            zz = residual_at_edge(tile, edges[nbr_xy], direction)
            if zz.size:
                abs_z = np.abs(zz)
                results.append({
                    'x': x, 'y': y, 'direction': direction,
                    'crossings': int(zz.size),
                    'p50': float(np.percentile(abs_z, 50)),
                    'p95': float(np.percentile(abs_z, 95)),
                    'max': float(abs_z.max()),
                    'signed_mean': float(zz.mean()),
                })
            else:
                results.append({'x': x, 'y': y, 'direction': direction, 'crossings': 0})

    assessable = [r for r in results if r['crossings'] >= 20]
    sparse = [r for r in results if 0 < r['crossings'] < 20]
    no_pairs = [r for r in results if r['crossings'] == 0]
    # Advisory screeners only: steep terrain, local curvature and sample count matter.
    flagged = [r for r in assessable if r['p95'] > 0.5 or abs(r['signed_mean']) > 0.2]
    worst = sorted(assessable, key=lambda r: r['p95'], reverse=True)[:12]

    print('GLOBAL_MASKED_TILES_AUDITED:', len(edges))
    print('GLOBAL_VALID_1M_PIXELS:', total_valid)
    print('GLOBAL_VALID_1M_AREA_KM2:', round(total_valid / 1e6, 6))
    print('GLOBAL_NODATA_INSIDE_FOOTPRINT_FROM_QA:', total_missing_reported)
    print('RPT_QA_COUNT_OR_METADATA_ISSUES:', len(issues))
    for item in issues[:20]:
        print('QA_ISSUE:', item)
    print('NEIGHBOR_SEAMS_FOUND:', len(results))
    print('SEAMS_ASSESSABLE_20PLUS:', len(assessable))
    print('SEAMS_SPARSE_1TO19:', len(sparse))
    print('SEAMS_NO_VALID_CROSSINGS:', len(no_pairs))
    print('ADVISORY_FLAGS_P95_GT_0_5M_OR_ABS_SIGNED_MEAN_GT_0_2M:', len(flagged))
    for item in sorted(flagged, key=lambda r: r['p95'], reverse=True)[:20]:
        print('FLAGGED_SEAM:', item)
    print('TOP_12_SEAMS_BY_P95_RESIDUAL:')
    for item in worst:
        print('  ', item)
    print('READ_ONLY_GLOBAL_QA_COMPLETE: no outputs modified')
    print('CAVEATS: second-difference seam screen is not a ground-truth accuracy test; masked footprint is AOI-clipped; sparse/no-data seams remain unassessed.')


if __name__ == '__main__':
    main()
