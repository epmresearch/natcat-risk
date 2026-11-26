#!/usr/bin/env python3
"""Locate inter-version LiDAR elevation outliers on a fixed 30 m grid.

Read-only source analysis. No mosaicking and no changes to archival LiDAR.
The results are screening diagnostics, not vertical-accuracy validation.

Requires: GDAL Python bindings and NumPy in the sea_to_sky_risk environment.
Run from anywhere: python src/data_sources/14_diagnose_step42_lidar_overlap_outliers.py
"""

from __future__ import annotations

import csv
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from osgeo import gdal, osr


gdal.UseExceptions()

ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = (ROOT.parent / "old files backup" / "Climada" / "Climada_Project"
           / "Inputs" / "DEM" / "LiDAR DEMs")
REF = ROOT / "data/processed/elevation/mrdem30_dtm_highway99_aoi_utm10_masked_exact.tif"
OUT_DIR = ROOT / "results/diagnostics"
JSON_PATH = OUT_DIR / "step4_2_archived_lidar_overlap_outliers_qa.json"
CSV_PATH = OUT_DIR / "step4_2_archived_lidar_overlap_outliers_by_pair.csv"
GEOJSON_PATH = OUT_DIR / "step4_2_archived_lidar_overlap_top_outliers.geojson"
DST_NODATA = -32767.0

# Pilot pairs and shared counts observed in Step 4.2.8 (script 13).
PAIRS = [
    ("092g085_2017_vs_2018", "bc_092g085_xl1m_utm10_170601_2.tif",
     "bc_092g085_xli1m_utm10_2018_2.tif", 13687),
    ("092g064_2019_vs_2018", "bc_092g064_xl1m_utm10_2019_2.tif",
     "bc_092g064_xli1m_utm10_2018_2.tif", 7904),
    ("092g065_2019_vs_2018", "bc_092g065_xl1m_utm10_2019_2.tif",
     "bc_092g065_xli1m_utm10_2018_2.tif", 13538),
]


def ensure_safe_outputs():
    for p in (JSON_PATH, CSV_PATH, GEOJSON_PATH):
        if p.exists():
            raise FileExistsError(f"Refusing to overwrite existing diagnostic: {p}")


def describe_reference():
    if not REF.is_file():
        raise FileNotFoundError(f"Missing reference raster: {REF}")
    ds = gdal.Open(str(REF), gdal.GA_ReadOnly)
    if ds is None:
        raise RuntimeError(f"Cannot open reference: {REF}")
    srs = osr.SpatialReference(wkt=ds.GetProjection())
    if srs.GetAuthorityCode(None) != "26910":
        raise RuntimeError("Reference CRS must be EPSG:26910")
    gt = ds.GetGeoTransform()
    if not (abs(gt[1] - 30) < 1e-8 and abs(gt[5] + 30) < 1e-8
            and abs(gt[2]) < 1e-9 and abs(gt[4]) < 1e-9):
        raise RuntimeError("Expected north-up 30 m reference grid")
    b = ds.GetRasterBand(1)
    arr = b.ReadAsArray()
    nodata = b.GetNoDataValue()
    if nodata is None:
        raise RuntimeError("Reference AOI mask raster has no NoData")
    aoi = np.isfinite(arr) & (arr != nodata)
    if np.count_nonzero(aoi) != 1159588:
        raise RuntimeError("Reference AOI valid-centre count differs from Step 4.2 QA")
    xmin, ymax = gt[0], gt[3]
    xmax = xmin + ds.RasterXSize * gt[1]
    ymin = ymax + ds.RasterYSize * gt[5]
    info = {"width": ds.RasterXSize, "height": ds.RasterYSize,
            "bounds": [xmin, ymin, xmax, ymax], "gt": gt}
    ds = None
    return aoi, info


def to_wgs84_transformer():
    src = osr.SpatialReference()
    src.ImportFromEPSG(26910)
    dst = osr.SpatialReference()
    dst.ImportFromEPSG(4326)
    src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    dst.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return osr.CoordinateTransformation(src, dst)


def warp_tile_to_grid(path: Path, grid):
    if not path.is_file():
        raise FileNotFoundError(str(path))
    src = gdal.Open(str(path), gdal.GA_ReadOnly)
    if src is None:
        raise RuntimeError(f"Cannot open source: {path}")
    src_band = src.GetRasterBand(1)
    nodata = src_band.GetNoDataValue()
    if nodata is None:
        raise RuntimeError(f"Source does not declare NoData: {path.name}")
    opts = gdal.WarpOptions(
        format="MEM", dstSRS="EPSG:26910", outputBounds=grid["bounds"],
        width=grid["width"], height=grid["height"],
        resampleAlg="near", errorThreshold=0,
        overviewLevel="NONE", srcNodata=nodata, dstNodata=DST_NODATA,
        outputType=gdal.GDT_Float32,
        warpOptions=["INIT_DEST=NO_DATA"],
    )
    dst = gdal.Warp("", src, options=opts)
    if dst is None:
        raise RuntimeError(f"GDAL warp failed: {path.name}")
    actual = dst.GetGeoTransform()
    if any(abs(float(a) - float(b)) > 1e-6 for a, b in zip(actual, grid["gt"])):
        raise RuntimeError(f"Output grid differs for {path.name}")
    arr = dst.GetRasterBand(1).ReadAsArray()
    dst = None
    src = None
    if arr.shape != (grid["height"], grid["width"]):
        raise RuntimeError(f"Unexpected warped shape: {path.name}")
    return arr, float(nodata)


def summarize_pair(pair_id, file_a, file_b, expected, aoi, grid, transformer):
    arr_a, nd_a = warp_tile_to_grid(ARCHIVE / file_a, grid)
    arr_b, nd_b = warp_tile_to_grid(ARCHIVE / file_b, grid)
    va = aoi & np.isfinite(arr_a) & (arr_a != DST_NODATA)
    vb = aoi & np.isfinite(arr_b) & (arr_b != DST_NODATA)
    shared = va & vb
    count = int(np.count_nonzero(shared))
    if count == 0:
        return {"pair": pair_id, "file_a": file_a, "file_b": file_b,
                "source_nodata_a": nd_a, "source_nodata_b": nd_b,
                "shared": 0, "expected_shared_script13": expected,
                "reproduces_script13_shared": False, "status": "NO_SHARED_SAMPLES"}, []
    signed = (arr_a[shared].astype(np.float64) - arr_b[shared].astype(np.float64))
    absdiff = np.abs(signed)
    rows, cols = np.nonzero(shared)
    nlargest = min(20, count)
    largest_idx = np.argpartition(absdiff, -nlargest)[-nlargest:]
    largest_idx = largest_idx[np.argsort(-absdiff[largest_idx])]
    features = []
    for rank, idx in enumerate(largest_idx, 1):
        rr, cc = int(rows[idx]), int(cols[idx])
        x = grid["gt"][0] + (cc + 0.5) * grid["gt"][1]
        y = grid["gt"][3] + (rr + 0.5) * grid["gt"][5]
        lon, lat, _ = transformer.TransformPoint(float(x), float(y))
        features.append({"type": "Feature", "geometry": {"type": "Point", "coordinates": [lon, lat]},
                         "properties": {"pair": pair_id, "rank": rank,
                                        "absolute_difference_m": round(float(absdiff[idx]), 6),
                                        "signed_difference_a_minus_b_m": round(float(signed[idx]), 6),
                                        "elevation_a_m": round(float(arr_a[rr, cc]), 6),
                                        "elevation_b_m": round(float(arr_b[rr, cc]), 6),
                                        "utm10_easting_m": round(float(x), 3),
                                        "utm10_northing_m": round(float(y), 3),
                                        "pixel_row": rr, "pixel_col": cc}})
    rec = {
        "pair": pair_id, "file_a": file_a, "file_b": file_b,
        "source_nodata_a": nd_a, "source_nodata_b": nd_b,
        "valid_a_inside": int(np.count_nonzero(va)),
        "valid_b_inside": int(np.count_nonzero(vb)),
        "shared": count,
        "expected_shared_script13": expected,
        "reproduces_script13_shared": count == expected,
        "median_signed_a_minus_b_m": round(float(np.median(signed)), 6),
        "mean_signed_a_minus_b_m": round(float(np.mean(signed)), 6),
        "median_absolute_m": round(float(np.median(absdiff)), 6),
        "p95_absolute_m": round(float(np.percentile(absdiff, 95)), 6),
        "p99_absolute_m": round(float(np.percentile(absdiff, 99)), 6),
        "maximum_absolute_m": round(float(np.max(absdiff)), 6),
        "difference_gt_1m_count": int(np.count_nonzero(absdiff > 1)),
        "difference_gt_5m_count": int(np.count_nonzero(absdiff > 5)),
        "difference_gt_10m_count": int(np.count_nonzero(absdiff > 10)),
        "difference_gt_1m_pct": round(float(100*np.mean(absdiff > 1)), 6),
        "difference_gt_5m_pct": round(float(100*np.mean(absdiff > 5)), 6),
        "difference_gt_10m_pct": round(float(100*np.mean(absdiff > 10)), 6),
        "status": "SCREENED_NOT_ACCURACY_VALIDATED"
    }
    del arr_a, arr_b
    return rec, features


def main():
    ensure_safe_outputs()
    aoi, grid = describe_reference()
    transformer = to_wgs84_transformer()
    records, features = [], []
    print("STEP 4.2 — SPATIAL OUTLIERS OF THREE LIDAR OVERLAP PAIRS")
    for i, (pair_id, f1, f2, expected) in enumerate(PAIRS, 1):
        print(f"[{i}/3] {pair_id}", flush=True)
        rec, feats = summarize_pair(pair_id, f1, f2, expected, aoi, grid, transformer)
        records.append(rec)
        features.extend(feats)
        print(f"  shared={rec['shared']} / expected={expected}; reproduced={rec['reproduces_script13_shared']}")
        if rec["shared"]:
            print(f"  median_abs={rec['median_absolute_m']}m; p99_abs={rec['p99_absolute_m']}m; "
                  f">5m={rec['difference_gt_5m_count']}; >10m={rec['difference_gt_10m_count']}; "
                  f"max={rec['maximum_absolute_m']}m")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_status = ("THREE_PAIRS_SPATIAL_OUTLIERS_SCREENED_WITH_CAVEATS"
                      if all(r["reproduces_script13_shared"] for r in records)
                      else "COUNTS_DIFFER_FROM_SCRIPT13_NEEDS_REVIEW")
    report = {"status": summary_status, "generated_utc": datetime.now(timezone.utc).isoformat(),
              "method": "Exact 30m EPSG:26910 grid from AOI reference; GDAL nearest neighbor; "
                        "errorThreshold=0; overviewLevel=NONE; source NoData respected. "
                        "Only pixels with both tiles valid inside the reference AOI are compared.",
              "reference_mask": str(REF), "archive": str(ARCHIVE),
              "source_tiles_modified": False, "elevation_accuracy_proven": False,
              "notes": ["Different observation dates, survey methods, registration, vertical datums, "
                        "and terrain can all affect pairwise differences.",
                        "Outlier point locations are screening locations, not evidence of erroneous tiles.",
                        "This is a 30m sample-grid comparison, not a native 1m overlap analysis."],
              "pairs": records}
    with CSV_PATH.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0].keys()), extrasaction='ignore')
        writer.writeheader()
        writer.writerows(records)
    JSON_PATH.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding='utf-8')
    GEOJSON_PATH.write_text(json.dumps({"type": "FeatureCollection", "features": features},
                                      ensure_ascii=False, indent=2), encoding='utf-8')
    print(f"Status: {summary_status}")
    print(f"CSV: {CSV_PATH}")
    print(f"JSON: {JSON_PATH}")
    print(f"Outlier points (WGS84 GeoJSON): {GEOJSON_PATH}")
    print("CAVEAT: spatial outliers are not proof that any source is erroneous.")
    if summary_status != "THREE_PAIRS_SPATIAL_OUTLIERS_SCREENED_WITH_CAVEATS":
        return 2
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        print(f"ERROR: {type(e).__name__}: {e}", file=sys.stderr)
        sys.exit(1)
