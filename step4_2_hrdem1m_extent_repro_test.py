"""Read-only-to-originals reproducibility test for one anomalous HRDEM 1m join.

Regenerates two 1km UNMASKED tiles (exact existing gdalwarp options) into a NEW
separate diagnostics folder, then compares them with the existing staging tiles
and the earlier continuous 60m bridge. Never overwrites original staged or
masked rasters; does not build a DEM mosaic or certify vertical accuracy.
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
DIAG = ROOT / "results" / "diagnostics"
STAGING = DIAG / "hrdem1m_batch_staging"
VRT = DIAG / "step4_2_hrdem1m_boundary_virtual.vrt"
OUTPUT = DIAG / "hrdem1m_extent_repro_qa"
BRIDGE = DIAG / "step4_2_hrdem1m_flagged_seam_bridge_60m.tif"
X1, X2, Y = 485000, 486000, 5490000
NODATA = -32767.0
EXPECTED_SOURCES = {
    "/vsicurl/https://canelevation-dem.s3.ca-central-1.amazonaws.com/hrdem-mosaic-1m/2_3-mosaic-1m-dtm.tif",
    "/vsicurl/https://canelevation-dem.s3.ca-central-1.amazonaws.com/hrdem-mosaic-1m/2_4-mosaic-1m-dtm.tif",
}


def check_source() -> None:
    if not VRT.is_file():
        raise FileNotFoundError(f"Missing original checked VRT: {VRT}")
    sources = {(e.text or "").strip() for e in ET.parse(VRT).iter("SourceFilename")}
    if sources != EXPECTED_SOURCES:
        raise ValueError(f"Unexpected source COGs: {sorted(sources)}")


def load(path: Path, xmin: int, ymin: int, width: int, height: int) -> np.ndarray:
    ds = gdal.Open(str(path), gdal.GA_ReadOnly)
    if ds is None:
        raise FileNotFoundError(f"Unable to open: {path}")
    if (ds.RasterCount, ds.RasterXSize, ds.RasterYSize) != (1, width, height):
        raise ValueError(f"Unexpected dimensions for {path.name}")
    expected_gt = (float(xmin), 1.0, 0.0, float(ymin + height), 0.0, -1.0)
    if not np.allclose(ds.GetGeoTransform(), expected_gt, rtol=0, atol=1e-8):
        raise ValueError(f"Unexpected geotransform for {path.name}")
    srs = osr.SpatialReference(wkt=ds.GetProjection())
    srs.AutoIdentifyEPSG()
    if srs.GetAuthorityCode(None) != "26910":
        raise ValueError(f"Unexpected CRS for {path.name}")
    band = ds.GetRasterBand(1)
    if band.DataType != gdal.GDT_Float32 or band.GetNoDataValue() != NODATA:
        raise ValueError(f"Unexpected band type/NoData for {path.name}")
    a = band.ReadAsArray()
    ds = None
    return a


def valid(a: np.ndarray) -> np.ndarray:
    return np.isfinite(a) & (a != NODATA)


def pcts(a: np.ndarray) -> list[float] | None:
    return np.round(np.percentile(a, [50, 95, 100]), 5).tolist() if a.size else None


def regenerate(x: int) -> Path:
    dest = OUTPUT / f"repro_hrdem1m_x{x}_y{Y}_utm10_1km.tif"
    if dest.exists():
        print("REUSING_EXISTING_DIAGNOSTIC:", dest.name, flush=True)
        load(dest, x, Y, 1000, 1000)
        return dest
    partial = OUTPUT / f"repro_hrdem1m_x{x}_y{Y}_utm10_1km.partial.tif"
    if partial.exists():
        partial.unlink()  # script-owned temporary output only
    cmd = [
        "gdalwarp", "-q", "-of", "GTiff", "-t_srs", "EPSG:26910",
        "-te", str(x), str(Y), str(x+1000), str(Y+1000),
        "-tr", "1", "1", "-tap", "-r", "bilinear", "-et", "0",
        "-ovr", "NONE", "-srcnodata", "-32767", "-dstnodata", "-32767",
        "-co", "TILED=YES", "-co", "COMPRESS=DEFLATE",
        "-co", "PREDICTOR=3", "-co", "BIGTIFF=IF_SAFER",
        str(VRT), str(partial),
    ]
    env = os.environ.copy()
    env.update({
        "GDAL_DISABLE_READDIR_ON_OPEN": "YES", "GDAL_HTTP_MAX_RETRY": "3",
        "GDAL_HTTP_RETRY_DELAY": "2", "GDAL_HTTP_TIMEOUT": "45",
    })
    try:
        print("REGENERATING_SEPARATE_DIAGNOSTIC:", x, Y, flush=True)
        subprocess.run(cmd, check=True, env=env)
        load(partial, x, Y, 1000, 1000)
        partial.replace(dest)
    finally:
        if partial.exists():
            partial.unlink()
    return dest


def seam_p95(west: np.ndarray, east: np.ndarray) -> tuple[int, float | None, float | None]:
    a0, a1 = west[:, -2], west[:, -1]
    b0, b1 = east[:, 0], east[:, 1]
    ok = valid(a0) & valid(a1) & valid(b0) & valid(b1)
    z = (b0[ok].astype(np.float64) - a1[ok]) - 0.5 * (
        (a1[ok].astype(np.float64) - a0[ok]) + (b1[ok].astype(np.float64) - b0[ok])
    )
    return int(ok.sum()), float(np.percentile(np.abs(z),95)) if z.size else None, float(z.mean()) if z.size else None


def main() -> None:
    gdal.UseExceptions()
    check_source()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    old = {}
    new = {}
    report = {"source": "checked_two_source_official_HRDEM_1m_VRT", "tile_comparisons": []}
    for x in (X1, X2):
        original = STAGING / f"hrdem1m_x{x}_y{Y}_utm10_1km.tif"
        old[x] = load(original, x, Y, 1000, 1000)
        regenerated_path = regenerate(x)
        new[x] = load(regenerated_path, x, Y, 1000, 1000)
        same_mask = np.array_equal(valid(old[x]), valid(new[x]))
        common = valid(old[x]) & valid(new[x])
        diff = np.abs(old[x][common].astype(np.float64) - new[x][common].astype(np.float64))
        exact_same = bool(same_mask and np.array_equal(old[x][common], new[x][common]))
        result = {
            "tile": [x, Y], "common_valid_pixels": int(common.sum()),
            "validity_masks_equal": bool(same_mask),
            "all_common_height_values_exactly_equal": bool(np.array_equal(old[x][common], new[x][common])),
            "fully_reproduced": exact_same,
            "abs_diff_p50_p95_max_m": pcts(diff),
            "pixels_differing_over_0_01m": int(np.count_nonzero(diff > 0.01)),
            "pixels_differing_over_0_5m": int(np.count_nonzero(diff > 0.5)),
        }
        report["tile_comparisons"].append(result)
        print("TILE_REPRO_RESULT:", result, flush=True)
    seam_old = seam_p95(old[X1], old[X2])
    seam_new = seam_p95(new[X1], new[X2])
    report["old_pair_seam_crossings_p95_mean"] = seam_old
    report["regenerated_pair_seam_crossings_p95_mean"] = seam_new
    print("OLD_PAIR_SEAM_CROSSINGS_P95_MEAN:", seam_old)
    print("REGENERATED_PAIR_SEAM_CROSSINGS_P95_MEAN:", seam_new)
    if BRIDGE.is_file():
        bridge = load(BRIDGE, 485970, Y, 60, 1000)
        seam_bridge = seam_p95(bridge[:, :30], bridge[:, 30:])
        report["existing_60m_bridge_seam_crossings_p95_mean"] = seam_bridge
        print("EXISTING_BRIDGE_SEAM_CROSSINGS_P95_MEAN:", seam_bridge)
    else:
        print("BRIDGE_NOT_PRESENT: comparison with existing bridge skipped")
    report["note"] = "Diagnostic reproduction only; no original files modified; not a vertical accuracy test"
    out_json = OUTPUT / "two_tile_warp_extent_repro_comparison.json"
    temp_json = OUTPUT / "two_tile_warp_extent_repro_comparison.json.partial"
    try:
        temp_json.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        temp_json.replace(out_json)
    finally:
        if temp_json.exists():
            temp_json.unlink()
    print("DIAGNOSTIC_JSON:", out_json)
    print("ORIGINAL_TILES_MODIFIED: FALSE")


if __name__ == "__main__":
    try:
        main()
    except (Exception, KeyboardInterrupt) as e:
        print("ERROR:", str(e), file=sys.stderr)
        sys.exit(1)
