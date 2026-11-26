"""Read-only batch QA for already staged Highway 99 HRDEM 1 m tiles.

Checks original raster values against the clipped official 1m DTM footprint,
then assesses adjacent 1 km tile seams. No COG downloads, no edits, no files written.

Important: staged rasters are UNMASKED. Valid pixels outside the *AOI-clipped*
footprint are expected and must not be classified as erroneous automatically.
This is screening QA, not vertical accuracy certification.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
from scipy.ndimage import distance_transform_edt
from osgeo import gdal, ogr


STEM = re.compile(r"hrdem1m_x(\d+)_y(\d+)_utm10_1km\.tif$")
NODATA = -32767.0


def raster_mask(layer, ds, x, y):
    m = gdal.GetDriverByName("MEM").Create("", 1000, 1000, 1, gdal.GDT_Byte)
    m.SetGeoTransform(ds.GetGeoTransform())
    m.SetProjection(ds.GetProjection())
    m.GetRasterBand(1).Fill(0)
    layer.SetSpatialFilterRect(x, y, x + 1000, y + 1000)
    if gdal.RasterizeLayer(m, [1], layer, burn_values=[1]) != 0:
        raise RuntimeError("Footprint rasterization failed")
    inside = m.ReadAsArray().astype(bool)
    m = None
    return inside


def edge_residual(a, b, direction):
    """Absolute deviation of the cross-boundary slope from two nearby slopes."""
    if direction == "east":
        good = a[1][:, -2] & a[1][:, -1] & b[1][:, 0] & b[1][:, 1]
        z, w = a[0], b[0]
        residual = (w[:, 0] - z[:, -1]) - ((z[:, -1] - z[:, -2]) + (w[:, 1] - w[:, 0])) / 2
    else:  # a south; b north
        good = a[1][0, :] & a[1][1, :] & b[1][-1, :] & b[1][-2, :]
        z, w = a[0], b[0]
        residual = (w[-1, :] - z[0, :]) - ((z[0, :] - z[1, :]) + (w[-2, :] - w[-1, :])) / 2
    return np.asarray(residual[good], dtype=np.float64)


def main():
    gdal.UseExceptions()
    ogr.UseExceptions()
    root = Path(__file__).resolve().parent
    source = root / "data/processed/elevation/hrdem_1m_dtm_coverage_highway99.gpkg"
    folder = root / "results/diagnostics/hrdem1m_batch_staging"
    if not source.exists() or not folder.is_dir():
        raise FileNotFoundError("Run from project root: missing footprint or staging folder")
    vector = ogr.Open(str(source), 0)
    layer = vector.GetLayer(0)
    tiles = {}
    for path in sorted(folder.glob("hrdem1m_x*_y*_utm10_1km.tif")):
        match = STEM.fullmatch(path.name)
        if not match:
            continue
        x, y = (int(match.group(1)), int(match.group(2)))
        ds = gdal.Open(str(path), gdal.GA_ReadOnly)
        if ds.RasterCount != 1 or ds.RasterXSize != 1000 or ds.RasterYSize != 1000:
            raise ValueError(f"Invalid dimensions/bands: {path.name}")
        gt = ds.GetGeoTransform()
        if not np.allclose(gt, [x, 1, 0, y + 1000, 0, -1], atol=1e-8, rtol=0):
            raise ValueError(f"Incorrect grid: {path.name}")
        if ds.GetRasterBand(1).GetNoDataValue() != NODATA:
            raise ValueError(f"Incorrect NoData: {path.name}")
        values = ds.ReadAsArray()
        valid = np.isfinite(values) & (values != NODATA)
        fp = raster_mask(layer, ds, x, y)
        missed = (~valid) & fp
        outside_valid = valid & (~fp)
        if missed.any():
            # Distance within the footprint to its complement, at 1m pixel resolution.
            depth = distance_transform_edt(fp)
            far = int((missed & (depth > 5)).sum())
            furthest = float(depth[missed].max())
        else:
            far, furthest = 0, 0.0
        print(f"TILE {x} {y}  valid={int(valid.sum())}"
              f"  footprint={int(fp.sum())}  missing_inside={int(missed.sum())}"
              f"  missing_over_5m={far}  missing_max_depth_m={furthest:.2f}"
              f"  valid_outside_clipped_footprint={int(outside_valid.sum())}", flush=True)
        tiles[(x, y)] = (values, valid)
        ds = None

    if not tiles:
        print("NO_STAGED_TILES_FOUND")
        return

    print("STAGED_TILES_AUDITED:", len(tiles), flush=True)
    found = 0
    for (x, y), a in sorted(tiles.items()):
        for direction, neighbor in (("east", (x + 1000, y)), ("north", (x, y + 1000))):
            if neighbor not in tiles:
                continue
            z = edge_residual(a, tiles[neighbor], direction)
            found += 1
            if z.size:
                q = np.percentile(np.abs(z), [50, 95, 100])
                print(f"SEAM {x} {y} -> {neighbor[0]} {neighbor[1]}"
                      f"  direction={direction}  valid_crossings={z.size}"
                      f"  residual_abs_P50_P95_MAX_m={np.round(q, 4).tolist()}"
                      f"  signed_mean_m={float(z.mean()):.4f}", flush=True)
            else:
                print(f"SEAM {x} {y} -> {neighbor[0]} {neighbor[1]}"
                      "  valid_crossings=0  NOT_ASSESSABLE", flush=True)
    print("NEIGHBOR_SEAMS_AUDITED:", found)
    print("READ_ONLY_QA_COMPLETE: no outputs modified")
    print("CAVEATS: clipped footprint; unmasked staging; seam QA is diagnostic only")


if __name__ == "__main__":
    main()
