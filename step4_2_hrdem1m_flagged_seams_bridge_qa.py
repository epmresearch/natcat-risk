"""Diagnose 14 globally flagged HRDEM 1m tile joins against continuous VRT warps.

Usage:
    python step4_2_hrdem1m_flagged_seams_bridge_qa.py            # plan only
    python step4_2_hrdem1m_flagged_seams_bridge_qa.py --run      # run all 14

This is NOT a tile repair operation. It reads original masked staging tiles,
warps small 60m-wide diagnostic strips from the existing official two-source
VRT, and writes only NEW diagnostic strips and QA JSON to results/diagnostics.
It does not alter any 1km staging tile, mask, source, or final DEM.

A near-zero bridge residual alongside a large old-tile residual supports a
warp-extent/tile-edge effect, NOT an independent elevation-accuracy claim.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

import numpy as np
from osgeo import gdal, osr

ROOT = Path(__file__).resolve().parent
DIAG = ROOT / 'results' / 'diagnostics'
MASK = DIAG / 'hrdem1m_masked_staging'
VRT = DIAG / 'step4_2_hrdem1m_boundary_virtual.vrt'
OUTPUT = DIAG / 'hrdem1m_flagged_seam_bridge_qa'
NODATA = -32767.0
FLAGS = [
    (485000, 5490000, 'east'),
    (484000, 5488000, 'east'),
    (483000, 5503000, 'east'),
    (485000, 5489000, 'north'),
    (485000, 5474000, 'north'),
    (483000, 5467000, 'north'),
    (484000, 5505000, 'east'),
    (483000, 5470000, 'east'),
    (485000, 5472000, 'east'),
    (485000, 5475000, 'north'),
    (484000, 5487000, 'east'),
    (489000, 5503000, 'east'),
    (490000, 5536000, 'north'),
    (486000, 5506000, 'east'),
]
EXPECTED_VRT_SOURCES = {
    '/vsicurl/https://canelevation-dem.s3.ca-central-1.amazonaws.com/hrdem-mosaic-1m/2_3-mosaic-1m-dtm.tif',
    '/vsicurl/https://canelevation-dem.s3.ca-central-1.amazonaws.com/hrdem-mosaic-1m/2_4-mosaic-1m-dtm.tif',
}


def validate_vrt() -> None:
    if not VRT.is_file():
        raise FileNotFoundError(f'Missing checked two-source VRT: {VRT}')
    names = {(n.text or '').strip() for n in ET.parse(VRT).iter('SourceFilename')}
    if names != EXPECTED_VRT_SOURCES:
        raise RuntimeError(f'Unexpected VRT source URLs: {sorted(names)}')


def box_for(x: int, y: int, direction: str) -> tuple[int, int, int, int]:
    if direction == 'east':
        return (x + 970, y, x + 1030, y + 1000)
    return (x, y + 970, x + 1000, y + 1030)


def raster(path: Path, xmin: int, ymin: int, xmax: int, ymax: int) -> np.ndarray:
    ds = gdal.Open(str(path), gdal.GA_ReadOnly)
    if ds is None:
        raise RuntimeError(f'Cannot open: {path}')
    expected_size = (xmax - xmin, ymax - ymin)
    expected_gt = (float(xmin), 1., 0., float(ymax), 0., -1.)
    if (ds.RasterXSize, ds.RasterYSize) != expected_size or ds.RasterCount != 1:
        raise RuntimeError(f'Unexpected dimensions: {path}')
    if not np.allclose(ds.GetGeoTransform(), expected_gt, rtol=0, atol=1e-8):
        raise RuntimeError(f'Unexpected geotransform: {path}')
    crs = osr.SpatialReference(wkt=ds.GetProjection())
    crs.AutoIdentifyEPSG()
    if crs.GetAuthorityCode(None) != '26910':
        raise RuntimeError(f'Unexpected CRS: {path}')
    band = ds.GetRasterBand(1)
    if band.DataType != gdal.GDT_Float32 or band.GetNoDataValue() != NODATA:
        raise RuntimeError(f'Unexpected dtype / NoData: {path}')
    arr = band.ReadAsArray()
    ds = None
    return arr


def tile_path(x: int, y: int) -> Path:
    return MASK / f'hrdem1m_x{x}_y{y}_utm10_1km_masked.tif'


def combined_old(x: int, y: int, direction: str) -> np.ndarray:
    if direction == 'east':
        western = raster(tile_path(x, y), x, y, x + 1000, y + 1000)
        eastern = raster(tile_path(x + 1000, y), x + 1000, y, x + 2000, y + 1000)
        return np.concatenate((western[:, -30:], eastern[:, :30]), axis=1)
    southern = raster(tile_path(x, y), x, y, x + 1000, y + 1000)
    northern = raster(tile_path(x, y + 1000), x, y + 1000, x + 1000, y + 2000)
    return np.concatenate((northern[-30:, :], southern[:30, :]), axis=0)


def valid(a: np.ndarray) -> np.ndarray:
    return np.isfinite(a) & (a != NODATA)


def seam_residual(a: np.ndarray, direction: str, common: np.ndarray) -> np.ndarray:
    """Use the same signed second-difference as the global QA script."""
    a = a.astype(np.float64, copy=False)
    if direction == 'east':
        ok = common[:, 28:32].all(axis=1)
        z = (a[:, 30] - a[:, 29]) - 0.5 * (
            (a[:, 29] - a[:, 28]) + (a[:, 31] - a[:, 30])
        )
    else:
        ok = common[28:32, :].all(axis=0)
        z = (a[29, :] - a[30, :]) - 0.5 * (
            (a[30, :] - a[31, :]) + (a[28, :] - a[29, :])
        )
    return z[ok]


def distribution(z: np.ndarray) -> dict:
    if not z.size:
        return {'crossings': 0, 'p50_abs_m': None, 'p95_abs_m': None,
                'max_abs_m': None, 'signed_mean_m': None}
    a = np.abs(z)
    return {'crossings': int(z.size),
            'p50_abs_m': float(np.percentile(a, 50)),
            'p95_abs_m': float(np.percentile(a, 95)),
            'max_abs_m': float(a.max()),
            'signed_mean_m': float(z.mean())}


def bridge_file(x: int, y: int, direction: str) -> Path:
    return OUTPUT / f'bridge_x{x}_y{y}_{direction}_60m_utm10.tif'


def generate_bridge(path: Path, bbox: tuple[int, int, int, int]) -> None:
    if path.exists():
        print('REUSE_BRIDGE:', path.name, flush=True)
        return
    partial = path.with_name(path.stem + '.partial.tif')
    if partial.exists():
        print('DELETE_OWN_INCOMPLETE_PARTIAL:', partial.name, flush=True)
        partial.unlink()
    cmd = [
        'gdalwarp', '-q', '-of', 'GTiff', '-t_srs', 'EPSG:26910',
        '-te', *map(str, bbox), '-tr', '1', '1', '-tap', '-r', 'bilinear',
        '-et', '0', '-ovr', 'NONE', '-srcnodata', '-32767', '-dstnodata', '-32767',
        '-co', 'TILED=YES', '-co', 'COMPRESS=DEFLATE', '-co', 'PREDICTOR=3',
        '-co', 'BIGTIFF=IF_SAFER', str(VRT), str(partial),
    ]
    env = os.environ.copy()
    env.update({'GDAL_DISABLE_READDIR_ON_OPEN': 'YES',
                'GDAL_HTTP_MAX_RETRY': '3', 'GDAL_HTTP_RETRY_DELAY': '2',
                'GDAL_HTTP_TIMEOUT': '45'})
    try:
        subprocess.run(cmd, env=env, check=True)
        raster(partial, *bbox)
        os.replace(partial, path)
    finally:
        if partial.exists():
            partial.unlink()


def one(x: int, y: int, direction: str) -> dict:
    bbox = box_for(x, y, direction)
    old = combined_old(x, y, direction)
    path = bridge_file(x, y, direction)
    generate_bridge(path, bbox)
    bridge = raster(path, *bbox)
    if old.shape != bridge.shape:
        raise RuntimeError(f'Strip shape mismatch: {old.shape} vs {bridge.shape}')
    common = valid(old) & valid(bridge)
    old_residual = distribution(seam_residual(old, direction, common))
    bridge_residual = distribution(seam_residual(bridge, direction, common))
    if old_residual['crossings'] != bridge_residual['crossings']:
        raise RuntimeError('Common seam crossing count inconsistent')
    delta = np.abs(old[common].astype(np.float64) - bridge[common].astype(np.float64))
    flags = {
        'old': (old_residual['crossings'] >= 20 and (
            old_residual['p95_abs_m'] > 0.5 or abs(old_residual['signed_mean_m']) > 0.2)),
        'bridge': (bridge_residual['crossings'] >= 20 and (
            bridge_residual['p95_abs_m'] > 0.5 or abs(bridge_residual['signed_mean_m']) > 0.2)),
    }
    row = {'x': x, 'y': y, 'direction': direction,
           'bridge_file': str(path.relative_to(ROOT)),
           'common_valid_pixels': int(common.sum()),
           'old_seam': old_residual, 'bridge_seam': bridge_residual,
           'delta_abs_p50_m': float(np.percentile(delta, 50)) if delta.size else None,
           'delta_abs_p95_m': float(np.percentile(delta, 95)) if delta.size else None,
           'delta_abs_max_m': float(delta.max()) if delta.size else None,
           'old_flagged_on_common_rows': flags['old'],
           'bridge_flagged_on_common_rows': flags['bridge'],
           'method': '60m continuous VRT gdalwarp vs existing two 1km masked tiles; same resampling/target grid; common valid pixels only; diagnostic not vertical accuracy'}
    return row


def write_json_atomic(path: Path, value: dict) -> None:
    tmp = path.with_suffix(path.suffix + '.partial')
    try:
        tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', action='store_true', help='Generate and compare the 14 diagnostic bridges')
    args = parser.parse_args()
    gdal.UseExceptions()
    if not MASK.is_dir():
        raise FileNotFoundError('Missing masked staging directory; place script in natcat-risk root')
    validate_vrt()
    print('FLAGGED_SEAMS_TO_CHECK:', len(FLAGS), flush=True)
    print('OUTPUT_DIRECTORY:', OUTPUT.relative_to(ROOT), flush=True)
    print('NO_ORIGINAL_TILE_MODIFICATIONS: TRUE', flush=True)
    print('SOURCE: checked two-source official HRDEM 1m VRT', flush=True)
    if not args.run:
        print('PLAN_ONLY: use python step4_2_hrdem1m_flagged_seams_bridge_qa.py --run')
        return
    OUTPUT.mkdir(parents=True, exist_ok=True)
    results = []
    for i, (x, y, direction) in enumerate(FLAGS, 1):
        print(f'SEAM_PROGRESS: {i}/{len(FLAGS)} ({x}, {y}, {direction})', flush=True)
        result = one(x, y, direction)
        results.append(result)
        old = result['old_seam']
        new = result['bridge_seam']
        print('SEAM_COMPARISON:', f'({x},{y},{direction})',
              'crossings:', old['crossings'],
              'old_p95:', round(old['p95_abs_m'], 4) if old['p95_abs_m'] is not None else 'NONE',
              'bridge_p95:', round(new['p95_abs_m'], 4) if new['p95_abs_m'] is not None else 'NONE',
              'old_mean:', round(old['signed_mean_m'], 4) if old['signed_mean_m'] is not None else 'NONE',
              'bridge_mean:', round(new['signed_mean_m'], 4) if new['signed_mean_m'] is not None else 'NONE',
              'bridge_flagged:', result['bridge_flagged_on_common_rows'], flush=True)
        write_json_atomic(OUTPUT / f'qa_x{x}_y{y}_{direction}.json', result)
    summary = {
        'seams_tested': len(results),
        'common_assessable_20plus': sum(r['old_seam']['crossings'] >= 20 for r in results),
        'old_flagged_on_common_rows': sum(r['old_flagged_on_common_rows'] for r in results),
        'bridge_flagged_on_common_rows': sum(r['bridge_flagged_on_common_rows'] for r in results),
        'bridge_still_flagged': [
            (r['x'], r['y'], r['direction'])
            for r in results if r['bridge_flagged_on_common_rows']
        ],
        'no_completed_source_tiles_modified': True,
        'not_final_DEM_or_vertical_accuracy': True,
        'results': results,
    }
    write_json_atomic(OUTPUT / 'all_14_bridge_comparisons.json', summary)
    print('ALL_14_BRIDGE_QA_COMPLETE:', len(results), flush=True)
    print('COMMON_ASSESSABLE_20PLUS:', summary['common_assessable_20plus'], flush=True)
    print('OLD_FLAGGED_ON_COMMON_ROWS:', summary['old_flagged_on_common_rows'], flush=True)
    print('BRIDGE_STILL_FLAGGED:', summary['bridge_flagged_on_common_rows'], flush=True)
    print('BRIDGE_STILL_FLAGGED_LIST:', summary['bridge_still_flagged'], flush=True)
    print('DIAGNOSTIC_JSON:', OUTPUT / 'all_14_bridge_comparisons.json', flush=True)
    print('NO_ORIGINAL_TILE_MODIFICATIONS: TRUE', flush=True)
    print('CAVEAT: no individual or global vertical accuracy is certified', flush=True)


if __name__ == '__main__':
    try:
        main()
    except (Exception, KeyboardInterrupt) as e:
        print('ERROR:', str(e), file=sys.stderr)
        sys.exit(1)
