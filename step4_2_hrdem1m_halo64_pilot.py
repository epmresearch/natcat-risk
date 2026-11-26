"""Diagnostic test: rewarp two adjacent HRDEM 1m tiles with a 64m halo.

This does NOT repair or replace original 1km tiles. It reads the validated
2-source VRT and writes two new padded tiles into a separate diagnostic folder.
It then compares the cropped 1km cores with existing unmasked staging rasters
and with a previously generated 60m continuous VRT bridge at the shared edge.

Goal: determine whether a halo around each independent warp eliminates the
extent-dependent seam artifact observed at x=486000, y=5490000.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

import numpy as np
from osgeo import gdal, osr

ROOT = Path(__file__).resolve().parent
DIAG = ROOT / 'results' / 'diagnostics'
VRT = DIAG / 'step4_2_hrdem1m_boundary_virtual.vrt'
STAGING = DIAG / 'hrdem1m_batch_staging'
BRIDGE = DIAG / 'step4_2_hrdem1m_flagged_seam_bridge_60m.tif'
OUT = DIAG / 'hrdem1m_halo64_pilot'
HALO = 64
NODATA = -32767.0
Y = 5490000
XS = (485000, 486000)
EXPECTED_VRT = {
    '/vsicurl/https://canelevation-dem.s3.ca-central-1.amazonaws.com/hrdem-mosaic-1m/2_3-mosaic-1m-dtm.tif',
    '/vsicurl/https://canelevation-dem.s3.ca-central-1.amazonaws.com/hrdem-mosaic-1m/2_4-mosaic-1m-dtm.tif',
}


def valid(a: np.ndarray) -> np.ndarray:
    return np.isfinite(a) & (a != NODATA)


def inspect(path: Path, xmin: int, ymin: int, width: int, height: int) -> np.ndarray:
    ds = gdal.Open(str(path), gdal.GA_ReadOnly)
    if ds is None:
        raise FileNotFoundError(f'Unable to open {path}')
    if (ds.RasterCount, ds.RasterXSize, ds.RasterYSize) != (1, width, height):
        raise ValueError(f'Wrong raster dimensions: {path}')
    expected_gt = (float(xmin), 1.0, 0.0, float(ymin + height), 0.0, -1.0)
    if not np.allclose(ds.GetGeoTransform(), expected_gt, atol=1e-8, rtol=0):
        raise ValueError(f'Wrong raster grid: {path}')
    crs = osr.SpatialReference(wkt=ds.GetProjection())
    crs.AutoIdentifyEPSG()
    if crs.GetAuthorityCode(None) != '26910':
        raise ValueError(f'Wrong raster CRS: {path}')
    band = ds.GetRasterBand(1)
    if band.DataType != gdal.GDT_Float32 or band.GetNoDataValue() != NODATA:
        raise ValueError(f'Wrong raster data type or nodata: {path}')
    arr = band.ReadAsArray()
    ds = None
    return arr


def warp_with_halo(x: int) -> np.ndarray:
    outfile = OUT / f'halo{HALO}_x{x}_y{Y}_utm10.tif'
    if not outfile.is_file():
        partial = OUT / f'halo{HALO}_x{x}_y{Y}_utm10.partial.tif'
        if partial.exists():
            partial.unlink()  # Remove only this script's incomplete diagnostic.
        bbox = (x - HALO, Y - HALO, x + 1000 + HALO, Y + 1000 + HALO)
        cmd = [
            'gdalwarp', '-q', '-of', 'GTiff', '-t_srs', 'EPSG:26910',
            '-te', *[str(v) for v in bbox], '-tr', '1', '1', '-tap',
            '-r', 'bilinear', '-et', '0', '-ovr', 'NONE',
            '-srcnodata', '-32767', '-dstnodata', '-32767',
            '-co', 'TILED=YES', '-co', 'COMPRESS=DEFLATE',
            '-co', 'PREDICTOR=3', '-co', 'BIGTIFF=IF_SAFER',
            str(VRT), str(partial),
        ]
        env = os.environ.copy()
        env.update({
            'GDAL_DISABLE_READDIR_ON_OPEN': 'YES',
            'GDAL_HTTP_MAX_RETRY': '3',
            'GDAL_HTTP_RETRY_DELAY': '2',
            'GDAL_HTTP_TIMEOUT': '45',
        })
        print('WARPING_PADDED_DIAGNOSTIC:', x, Y, 'HALO_M:', HALO, flush=True)
        try:
            subprocess.run(cmd, check=True, env=env)
            inspect(partial, x - HALO, Y - HALO, 1000 + 2 * HALO, 1000 + 2 * HALO)
            os.replace(partial, outfile)
        finally:
            if partial.exists():
                partial.unlink()
    else:
        print('REUSING_PADDED_DIAGNOSTIC:', outfile.name, flush=True)
    full = inspect(outfile, x - HALO, Y - HALO, 1000 + 2 * HALO, 1000 + 2 * HALO)
    core = full[HALO:HALO + 1000, HALO:HALO + 1000].copy()
    if core.shape != (1000, 1000):
        raise RuntimeError('Unexpected core shape')
    return core


def seam_residual(strip: np.ndarray, rows: np.ndarray) -> tuple[int, float | None, float | None]:
    z = strip.astype(np.float64, copy=False)
    resid = (z[:, 30] - z[:, 29]) - 0.5 * ((z[:, 29] - z[:, 28]) + (z[:, 31] - z[:, 30]))
    v = resid[rows]
    return (int(v.size), float(np.percentile(np.abs(v), 95)) if v.size else None,
            float(v.mean()) if v.size else None)


def pct(diff: np.ndarray) -> list[float] | None:
    return [float(z) for z in np.percentile(diff, [50, 95, 100])] if diff.size else None


def main() -> None:
    gdal.UseExceptions()
    if not VRT.is_file():
        raise FileNotFoundError(f'Missing official VRT: {VRT}')
    sources = {(node.text or '').strip() for node in ET.parse(VRT).iter('SourceFilename')}
    if sources != EXPECTED_VRT:
        raise ValueError(f'Unexpected VRT sources: {sorted(sources)}')
    if not BRIDGE.is_file():
        raise FileNotFoundError(f'Missing earlier 60m diagnostic bridge: {BRIDGE}')
    OUT.mkdir(parents=True, exist_ok=True)
    old, padded, tile_reports = {}, {}, []
    for x in XS:
        old[x] = inspect(STAGING / f'hrdem1m_x{x}_y{Y}_utm10_1km.tif', x, Y, 1000, 1000)
        padded[x] = warp_with_halo(x)
        original_valid, padded_valid = valid(old[x]), valid(padded[x])
        common = original_valid & padded_valid
        diff = np.abs(old[x][common].astype(np.float64) - padded[x][common].astype(np.float64))
        result = {
            'tile': [x, Y],
            'old_valid_pixels': int(original_valid.sum()),
            'halo_valid_pixels': int(padded_valid.sum()),
            'common_valid_pixels': int(common.sum()),
            'halo_newly_valid_pixels': int((~original_valid & padded_valid).sum()),
            'halo_newly_nodata_pixels': int((original_valid & ~padded_valid).sum()),
            'old_vs_halo_abs_diff_p50_p95_max_m': pct(diff),
            'pixels_differing_over_0_5m': int((diff > 0.5).sum()),
        }
        tile_reports.append(result)
        print('HALO_TILE_COMPARISON:', result, flush=True)

    old_strip = np.concatenate((old[XS[0]][:, -30:], old[XS[1]][:, :30]), axis=1)
    halo_strip = np.concatenate((padded[XS[0]][:, -30:], padded[XS[1]][:, :30]), axis=1)
    bridge = inspect(BRIDGE, 485970, Y, 60, 1000)
    all_valid = valid(old_strip) & valid(halo_strip) & valid(bridge)
    comparable_seam_rows = all_valid[:, 28:32].all(axis=1)
    old_metrics = seam_residual(old_strip, comparable_seam_rows)
    halo_metrics = seam_residual(halo_strip, comparable_seam_rows)
    bridge_metrics = seam_residual(bridge, comparable_seam_rows)
    overlap = valid(halo_strip) & valid(bridge)
    bridge_diff = np.abs(halo_strip[overlap].astype(np.float64) - bridge[overlap].astype(np.float64))

    print('COMMON_SEAM_ROWS_ALL_THREE:', int(comparable_seam_rows.sum()))
    print('OLD_PAIR_SEAM_CROSSINGS_P95_MEAN:', old_metrics)
    print('HALO64_PAIR_SEAM_CROSSINGS_P95_MEAN:', halo_metrics)
    print('CONTINUOUS_BRIDGE_SEAM_CROSSINGS_P95_MEAN:', bridge_metrics)
    print('HALO64_MINUS_BRIDGE_COMMON_VALID_PIXELS:', int(overlap.sum()))
    print('HALO64_MINUS_BRIDGE_ABS_P50_P95_MAX_M:', pct(bridge_diff))

    summary = {
        'halo_m': HALO,
        'source': 'official_checked_two_source_1m_VRT',
        'tile_comparisons': tile_reports,
        'common_seam_rows_all_three': int(comparable_seam_rows.sum()),
        'old_pair_seam_count_p95_mean': old_metrics,
        'halo64_pair_seam_count_p95_mean': halo_metrics,
        'continuous_bridge_seam_count_p95_mean': bridge_metrics,
        'halo64_bridge_common_valid_pixels': int(overlap.sum()),
        'halo64_bridge_abs_diff_p50_p95_max_m': pct(bridge_diff),
        'notes': 'Diagnostic only; no staging/masked rasters modified, no vertical accuracy certification',
    }
    json_path = OUT / 'halo64_two_tile_comparison.json'
    tmp = json_path.with_suffix('.json.partial')
    try:
        tmp.write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
        os.replace(tmp, json_path)
    finally:
        if tmp.exists():
            tmp.unlink()
    print('DIAGNOSTIC_JSON:', json_path)
    print('ORIGINAL_TILES_MODIFIED: FALSE')


if __name__ == '__main__':
    try:
        main()
    except (Exception, KeyboardInterrupt) as e:
        print('ERROR:', str(e), file=sys.stderr)
        sys.exit(1)
