#!/usr/bin/env python3
"""Step 4.2: pilot inter-version elevation comparison on shared 30-m AOI pixels.

Read-only inputs:
- Original archived bc_* TIFFs (never modified)
- Step 11 archived-LiDAR 30-m sample QA JSON
- Step 12 full native-pixel scan QA JSON
- Step 08 AOI-masked MRDEM raster, used ONLY for the AOI grid/mask

Outputs: two new small reports under results/diagnostics (CSV and JSON).
Does NOT create a DEM mosaic, determine which acquisition is best, verify the
DTM/DSM product type, or test against ground-truth elevation.

Run from any directory:
    python src/data_sources/13_pilot_step42_archived_lidar_overlap_elevation.py
"""
from __future__ import annotations

import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from osgeo import gdal, osr
from pyproj import CRS, Transformer


gdal.UseExceptions()
ROOT = Path(__file__).resolve().parents[2]
DIAG = ROOT / "results/diagnostics"
STEP11 = DIAG / "step4_2_archived_lidar_30m_sample_validity_qa.json"
STEP12 = DIAG / "step4_2_archived_lidar_native_pixel_all_qa.json"
REFERENCE = ROOT / "data/processed/elevation/mrdem30_dtm_highway99_aoi_utm10_masked_exact.tif"
OUT_CSV = DIAG / "step4_2_archived_lidar_overlap_elevation_pilot_by_pair.csv"
OUT_JSON = DIAG / "step4_2_archived_lidar_overlap_elevation_pilot_qa.json"
NODATA_OUT = -32767.0
PIXEL_METRES = 30.0
AOI_CENTRES_EXPECTED = 1_159_588
TARGET_CRS = CRS.from_epsg(26910)

# Deliberately paired older/newer rasters from the same named map tile.
# These are pilot comparisons, NOT an automated preference for newer years.
PAIRS = [
    ("bc_092g085_xl1m_utm10_170601_2.tif", "bc_092g085_xli1m_utm10_2018_2.tif"),
    ("bc_092g064_xl1m_utm10_2019_2.tif", "bc_092g064_xli1m_utm10_2018_2.tif"),
    ("bc_092g065_xl1m_utm10_2019_2.tif", "bc_092g065_xli1m_utm10_2018_2.tif"),
]


def reference_aoi():
    if not REFERENCE.is_file():
        raise FileNotFoundError(f"Reference AOI raster missing: {REFERENCE}")
    ds = gdal.Open(str(REFERENCE), gdal.GA_ReadOnly)
    if ds is None or ds.RasterCount != 1:
        raise RuntimeError("Could not open one-band reference AOI raster")
    srs = osr.SpatialReference()
    srs.ImportFromWkt(ds.GetProjection())
    epsg = osr.SpatialReference()
    epsg.ImportFromEPSG(26910)
    if not bool(srs.IsSame(epsg)):
        raise RuntimeError("Reference raster CRS is not EPSG:26910")
    gt = ds.GetGeoTransform()
    if not (abs(gt[1] - PIXEL_METRES) < 1e-7 and abs(gt[5] + PIXEL_METRES) < 1e-7
            and abs(gt[2]) < 1e-9 and abs(gt[4]) < 1e-9):
        raise RuntimeError(f"Reference 30-m grid unexpected: {gt}")
    band = ds.GetRasterBand(1)
    if band.GetNoDataValue() != NODATA_OUT:
        raise RuntimeError("Unexpected reference NoData value")
    pixels = band.ReadAsArray()
    inside = (np.isfinite(pixels) & (pixels != NODATA_OUT)
              & (band.GetMaskBand().ReadAsArray() != 0))
    count = int(inside.sum())
    if count != AOI_CENTRES_EXPECTED:
        raise RuntimeError(f"AOI count changed: {count}, expected {AOI_CENTRES_EXPECTED}")
    result = inside.copy(), gt
    ds = None
    return result


def get_registry():
    for p in (STEP11, STEP12):
        if not p.is_file():
            raise FileNotFoundError(f"Required earlier QA is missing: {p}")
    report11 = json.loads(STEP11.read_text(encoding="utf-8-sig"))
    report12 = json.loads(STEP12.read_text(encoding="utf-8-sig"))
    if (report11.get("candidate_tile_count") != 32 or
            report11.get("tile_samples_success") != 32 or
            report11.get("tile_samples_error") != 0):
        raise RuntimeError("Step 11 30-m candidate inventory is not a clean 32-tile scan")
    if (report12.get("status") != "FULL_PER_TILE_SCAN_COMPLETE_WITH_CAVEATS" or
            report12.get("num_succeeded") != 32 or report12.get("num_errors") != 0):
        raise RuntimeError("Step 12 complete 1-m native scan report is required")
    records = {r["filename"]: r for r in report11.get("tile_records", [])}
    if len(records) != 32:
        raise RuntimeError("Step 11 candidate records missing or duplicated")
    native_records = {r["filename"]: r for r in report12.get("records", [])}
    if len(native_records) != 32:
        raise RuntimeError("Step 12 candidate records missing or duplicated")
    expected = {name for pair in PAIRS for name in pair}
    if not expected.issubset(records) or not expected.issubset(native_records):
        raise RuntimeError("One or more pilot tiles are absent from prior scans")
    archive = Path(report11["archive_path"])
    if not archive.is_dir():
        raise FileNotFoundError(f"Archive directory does not exist: {archive}")
    return archive, records


def target_window(src, ref_gt, ncols, nrows):
    """Exact same 20-step edge-densified and 2-cell-padded window as Step 11."""
    transformer = Transformer.from_crs(CRS.from_wkt(src.GetProjection()), TARGET_CRS,
                                       always_xy=True)
    nx, ny = src.RasterXSize, src.RasterYSize
    gt = src.GetGeoTransform()
    edges = [((0, 0), (nx, 0)), ((nx, 0), (nx, ny)),
             ((nx, ny), (0, ny)), ((0, ny), (0, 0))]
    xs, ys = [], []
    for (x0, y0), (x1, y1) in edges:
        for t in np.linspace(0, 1, 21):
            col, row = x0 + t * (x1 - x0), y0 + t * (y1 - y0)
            xs.append(gt[0] + col * gt[1] + row * gt[2])
            ys.append(gt[3] + col * gt[4] + row * gt[5])
    xx, yy = transformer.transform(xs, ys)
    if not all(math.isfinite(float(v)) for v in (*xx, *yy)):
        raise RuntimeError("Non-finite source footprint coordinates")
    left, top = ref_gt[0], ref_gt[3]
    col0 = max(0, math.floor((min(xx) - left) / PIXEL_METRES) - 2)
    col1 = min(ncols, math.ceil((max(xx) - left) / PIXEL_METRES) + 2)
    row0 = max(0, math.floor((top - max(yy)) / PIXEL_METRES) - 2)
    row1 = min(nrows, math.ceil((top - min(yy)) / PIXEL_METRES) + 2)
    if col0 >= col1 or row0 >= row1:
        raise RuntimeError("Source does not intersect target grid")
    bounds = (left + col0 * PIXEL_METRES, top - row1 * PIXEL_METRES,
              left + col1 * PIXEL_METRES, top - row0 * PIXEL_METRES)
    return row0, row1, col0, col1, bounds


def sample_tile(path, aoi, ref_gt, step11_record):
    if not path.is_file():
        raise FileNotFoundError(f"Source tile not found: {path}")
    src = gdal.Open(str(path), gdal.GA_ReadOnly)
    if src is None or src.RasterCount != 1:
        raise RuntimeError(f"Could not open a single-band TIFF: {path.name}")
    if not src.GetProjection():
        raise RuntimeError(f"Missing CRS: {path.name}")
    orig_band = src.GetRasterBand(1)
    src_nodata = orig_band.GetNoDataValue()
    unit_type = (orig_band.GetUnitType() or "").strip()
    if unit_type and unit_type.lower() not in ("m", "meter", "metre", "meters", "metres"):
        raise RuntimeError(f"Non-metre vertical unit reported for {path.name}: {unit_type}")
    if src_nodata is None:
        raise RuntimeError(f"Missing source NoData definition: {path.name}")
    nrows, ncols = aoi.shape
    row0, row1, col0, col1, bounds = target_window(src, ref_gt, ncols, nrows)
    opts = gdal.WarpOptions(
        format="MEM", dstSRS="EPSG:26910", outputBounds=list(bounds),
        width=col1 - col0, height=row1 - row0,
        resampleAlg="near", errorThreshold=0.0, overviewLevel="NONE",
        outputType=gdal.GDT_Float32, dstNodata=NODATA_OUT,
        srcNodata=src_nodata, multithread=False,
    )
    warped = gdal.Warp("", src, options=opts)
    if warped is None:
        raise RuntimeError(f"Warp error: {path.name}")
    expected_gt = (bounds[0], PIXEL_METRES, 0.0, bounds[3], 0.0, -PIXEL_METRES)
    if (warped.RasterXSize != col1-col0 or warped.RasterYSize != row1-row0 or
            any(abs(v - e) > 1e-6 for v, e in zip(warped.GetGeoTransform(), expected_gt))):
        raise RuntimeError(f"Warped grid mismatch: {path.name}")
    b = warped.GetRasterBand(1)
    vals = b.ReadAsArray()
    mask = ((b.GetMaskBand().ReadAsArray() != 0) & np.isfinite(vals)
            & (vals != NODATA_OUT) & aoi[row0:row1, col0:col1])
    current_count = int(mask.sum())
    expected_count = int(step11_record["valid_samples_inside_aoi"])
    if current_count != expected_count:
        raise RuntimeError(f"30-m reproducibility FAILURE for {path.name}: "
                           f"expected {expected_count}, got {current_count}")
    # The output arrays are full-grid so all comparisons use identical coordinates.
    values_full = np.full((nrows, ncols), np.nan, dtype=np.float32)
    valid_full = np.zeros((nrows, ncols), dtype=bool)
    values_full[row0:row1, col0:col1][mask] = vals[mask]
    valid_full[row0:row1, col0:col1][mask] = True
    warped = None
    src = None
    return values_full, valid_full, current_count, unit_type


def round_number(value, digits=6):
    return round(float(value), digits)


def compare(a_name, b_name, archive, records, aoi, gt):
    a_values, a_valid, a_count, a_unit = sample_tile(archive / a_name, aoi, gt, records[a_name])
    b_values, b_valid, b_count, b_unit = sample_tile(archive / b_name, aoi, gt, records[b_name])
    shared = a_valid & b_valid
    shared_n = int(shared.sum())
    result = {
        "tile_a": a_name, "tile_b": b_name,
        "tile_a_declared_vertical_unit": a_unit or "NOT_DECLARED",
        "tile_b_declared_vertical_unit": b_unit or "NOT_DECLARED",
        "tile_a_valid_30m_aoi_samples": a_count,
        "tile_b_valid_30m_aoi_samples": b_count,
        "shared_valid_30m_aoi_samples": shared_n,
        "shared_valid_30m_sample_area_km2": round_number(shared_n * 900.0 / 1e6),
        "shared_30m_sample_pct_of_aoi": round_number(100.0 * shared_n / AOI_CENTRES_EXPECTED),
        "difference_sign": "tile_a_minus_tile_b_metres",
    }
    if shared_n:
        diffs = (a_values[shared].astype(np.float64)
                 - b_values[shared].astype(np.float64))
        abs_diffs = np.abs(diffs)
        result.update({
            "signed_mean_difference_m": round_number(np.mean(diffs)),
            "signed_median_difference_m": round_number(np.median(diffs)),
            "mean_absolute_difference_m": round_number(np.mean(abs_diffs)),
            "median_absolute_difference_m": round_number(np.median(abs_diffs)),
            "p90_absolute_difference_m": round_number(np.percentile(abs_diffs, 90)),
            "p95_absolute_difference_m": round_number(np.percentile(abs_diffs, 95)),
            "p99_absolute_difference_m": round_number(np.percentile(abs_diffs, 99)),
            "max_absolute_difference_m": round_number(np.max(abs_diffs)),
        })
        for threshold, label in [(0.1, "0p1m"), (0.5, "0p5m"),
                                 (1., "1m"), (5., "5m"), (10., "10m")]:
            count = int(np.count_nonzero(abs_diffs > threshold))
            result[f"abs_diff_gt_{label}_count"] = count
            result[f"abs_diff_gt_{label}_pct"] = round_number(100. * count / shared_n)
        result["status"] = "COMPARABLE_PILOT_OVERLAP"
    else:
        result["status"] = "NO_SHARED_VALID_SAMPLES_IN_PILOT"
    return result


def main():
    if OUT_JSON.exists() or OUT_CSV.exists():
        raise FileExistsError("Pilot outputs already exist; refusing to overwrite: "
                              f"{OUT_JSON} / {OUT_CSV}")
    archive, records = get_registry()
    aoi, gt = reference_aoi()
    rows = []
    errors = []
    print("STEP 4.2 — LIDAR INTER-VERSION OVERLAP ELEVATION PILOT")
    print("Method: same 30-m reference AOI grid; nearest-neighbour; source NoData respected")
    for idx, (a, b) in enumerate(PAIRS, 1):
        print(f"[{idx}/{len(PAIRS)}] {a} VS {b}", flush=True)
        try:
            result = compare(a, b, archive, records, aoi, gt)
            rows.append(result)
            if result["shared_valid_30m_aoi_samples"]:
                print(f"  shared={result['shared_valid_30m_aoi_samples']:,}, "
                      f"median_abs_m={result['median_absolute_difference_m']}, "
                      f"p95_abs_m={result['p95_absolute_difference_m']}, "
                      f"max_abs_m={result['max_absolute_difference_m']}, "
                      f"above_1m={result['abs_diff_gt_1m_pct']}%", flush=True)
            else:
                print("  shared=0; no comparison possible for this pair", flush=True)
        except Exception as exc:
            errors.append({"tile_a": a, "tile_b": b,
                           "error": f"{type(exc).__name__}: {exc}"})
            print(f"  ERROR: {errors[-1]['error']}", flush=True)
    status = ("PILOT_OVERLAP_COMPARISON_COMPLETE_WITH_CAVEATS" if not errors
              and all(r["shared_valid_30m_aoi_samples"] > 0 for r in rows)
              else "PILOT_OVERLAP_COMPARISON_INCOMPLETE_NEEDS_REVIEW")
    report = {
        "status": status,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_step11": str(STEP11),
        "source_step12": str(STEP12),
        "reference_aoi_grid": str(REFERENCE),
        "archive_path": str(archive),
        "pair_count_requested": len(PAIRS),
        "pair_count_completed": len(rows),
        "pair_errors": errors,
        "comparisons": rows,
        "caveats": [
            "Inter-version differences use 30-m nearest-neighbour point samples, NOT complete 1-m overlaps.",
            "Signed elevation difference is tile_a minus tile_b, in metres.",
            "Non-metre declared band units stop a comparison; absent units are labelled NOT_DECLARED, so metre assumptions need source-documentation checks.",
            "No vertical datum alignment, co-registration or independent ground truth is tested.",
            "Nominal years and tile IDs in filenames are not independently verified acquisition metadata.",
            "Differences may reflect different surveys, processing, source products, or real surface change.",
            "A completed comparison is not an accuracy PASS and does not authorize mosaic/averaging.",
            "Original archived TIFFs and project elevation files are read only.",
        ],
    }
    DIAG.mkdir(parents=True, exist_ok=True)
    if rows:
        fields = list(dict.fromkeys(k for row in rows for k in row.keys()))
        with OUT_CSV.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print("Status:", status)
    print("CSV:", OUT_CSV if rows else "not written")
    print("JSON:", OUT_JSON)
    print("CAVEAT: comparisons are screening metrics, not proof of DTM accuracy or suitability")
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
