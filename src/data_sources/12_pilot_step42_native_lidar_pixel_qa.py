#!/usr/bin/env python3
"""Step 4.2: pilot native-resolution valid-pixel QA for archived LiDAR tiles.

Run from anywhere, using the sea_to_sky_risk environment:
    python src/data_sources/12_pilot_step42_native_lidar_pixel_qa.py

Read-only inputs: existing AOI, Step-11 CSV/JSON, and 3 archived TIFFs.
Outputs: per-tile CSV and one JSON under results/diagnostics.
No mosaicking, resampling, reprojected raster output, or source modification.
Pixel-centre AOI inclusion is rasterized on each source raster's original grid.
The per-tile totals are NOT a de-duplicated union of 1-m source pixels.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from osgeo import gdal, ogr, osr


gdal.UseExceptions()

PROJECT = Path(__file__).resolve().parents[2]
AOI_PATH = PROJECT / "data/processed/study_area/highway99_terrain_search_5000m.geojson"
DIAGNOSTICS = PROJECT / "results/diagnostics"
STEP11_CSV = DIAGNOSTICS / "step4_2_archived_lidar_30m_sample_validity_by_tile.csv"
STEP11_JSON = DIAGNOSTICS / "step4_2_archived_lidar_30m_sample_validity_qa.json"
OUT_CSV = DIAGNOSTICS / "step4_2_archived_lidar_native_pixel_pilot_by_tile.csv"
OUT_JSON = DIAGNOSTICS / "step4_2_archived_lidar_native_pixel_pilot_qa.json"

# Different acquisition-year labels, source NoData values, and 30-m sampling behavior.
PILOT_FILES = [
    "bc_092g085_xli1m_utm10_2018_2.tif",  # High valid coverage in Step 11
    "bc_092g044_xl1m_utm10_2019_2.tif",   # Uses -9999 NoData
    "bc_092g084_xli1m_utm10_2020_2.tif", # Zero valid 30-m samples
]


def load_reference():
    with open(STEP11_JSON, encoding="utf-8-sig") as f:
        report = json.load(f)
    with open(STEP11_CSV, newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    keyed = {row["filename"]: row for row in rows}
    if len(keyed) != len(rows):
        raise RuntimeError("Duplicate filenames in Step-11 tile CSV")
    if report.get("tile_samples_error") != 0 or report.get("candidate_tile_count") != 32:
        raise RuntimeError("Expected Step-11 inventory of 32 successful candidates")
    archive = Path(report["archive_path"])
    if not archive.is_dir():
        raise FileNotFoundError(f"Archived LiDAR folder does not exist: {archive}")
    return archive, keyed


def aoi_geometry_in_srs(dst_srs):
    source_ds = ogr.Open(str(AOI_PATH), 0)
    if source_ds is None:
        raise RuntimeError(f"Cannot open AOI: {AOI_PATH}")
    layer = source_ds.GetLayer(0)
    src_srs = layer.GetSpatialRef()
    if src_srs is None:
        raise RuntimeError("AOI has no declared source CRS")
    src_srs = src_srs.Clone()
    dst_srs = dst_srs.Clone()
    src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    dst_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    transform = osr.CoordinateTransformation(src_srs, dst_srs)
    result = None
    for feature in layer:
        geometry = feature.GetGeometryRef()
        if geometry is None or geometry.IsEmpty():
            continue
        g = geometry.Clone()
        g.Transform(transform)
        result = g if result is None else result.Union(g)
    source_ds = None
    if result is None or result.IsEmpty():
        raise RuntimeError("AOI has no nonempty geometries")
    return result


def aoi_layer_in_memory(geom, srs):
    drv = ogr.GetDriverByName("Memory")
    ds = drv.CreateDataSource("aoi")
    layer = ds.CreateLayer("aoi", srs=srs, geom_type=ogr.wkbUnknown)
    feature = ogr.Feature(layer.GetLayerDefn())
    feature.SetGeometry(geom)
    if layer.CreateFeature(feature) != 0:
        raise RuntimeError("Could not create AOI feature in memory")
    feature = None
    return ds, layer


def raster_block_aoi_mask(layer, projection, src_gt, xoff, yoff, width, height):
    block_ds = gdal.GetDriverByName("MEM").Create("", width, height, 1, gdal.GDT_Byte)
    block_ds.SetProjection(projection)
    block_gt = (
        src_gt[0] + xoff * src_gt[1] + yoff * src_gt[2],
        src_gt[1],
        src_gt[2],
        src_gt[3] + xoff * src_gt[4] + yoff * src_gt[5],
        src_gt[4],
        src_gt[5],
    )
    block_ds.SetGeoTransform(block_gt)
    block_ds.GetRasterBand(1).Fill(0)
    result = gdal.RasterizeLayer(block_ds, [1], layer, burn_values=[1])
    if result != 0:
        raise RuntimeError(f"GDAL RasterizeLayer returned {result}")
    mask = block_ds.GetRasterBand(1).ReadAsArray() != 0
    block_ds = None
    return mask


def scan_native_tile(path, step11_row, chunk=1024):
    ds = gdal.Open(str(path), gdal.GA_ReadOnly)
    if ds is None:
        raise RuntimeError(f"Cannot open: {path}")
    if ds.RasterCount != 1:
        raise RuntimeError(f"Expected single-band elevation TIFF, got {ds.RasterCount}")
    source_srs = ds.GetSpatialRef()
    if source_srs is None:
        raise RuntimeError("Tile has no source CRS")
    gt = ds.GetGeoTransform()
    if not (gt[1] > 0 and gt[5] < 0 and abs(gt[2]) < 1e-10 and abs(gt[4]) < 1e-10):
        raise RuntimeError(f"Expected north-up projected raster, got geotransform {gt}")
    source_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    if not source_srs.IsProjected():
        raise RuntimeError("Expected projected metre CRS")
    native_aoi = aoi_geometry_in_srs(source_srs)
    memory_ds, aoi_layer = aoi_layer_in_memory(native_aoi, source_srs)
    left, right, bottom, top = native_aoi.GetEnvelope()
    # Crop to AOI envelope, expanded by 2 source pixels for robust boundaries.
    x0 = max(0, math.floor((left - gt[0]) / gt[1]) - 2)
    x1 = min(ds.RasterXSize, math.ceil((right - gt[0]) / gt[1]) + 2)
    y0 = max(0, math.floor((gt[3] - top) / abs(gt[5])) - 2)
    y1 = min(ds.RasterYSize, math.ceil((gt[3] - bottom) / abs(gt[5])) + 2)
    nodata = ds.GetRasterBand(1).GetNoDataValue()
    source_band = ds.GetRasterBand(1)
    validity_band = source_band.GetMaskBand()
    inside_count = 0
    valid_count = 0
    nonfinite_inside_count = 0
    blocks = 0
    intersecting_blocks = 0
    if x1 > x0 and y1 > y0:
        for yy in range(y0, y1, chunk):
            hh = min(chunk, y1 - yy)
            for xx in range(x0, x1, chunk):
                ww = min(chunk, x1 - xx)
                blocks += 1
                inside = raster_block_aoi_mask(aoi_layer, ds.GetProjection(), gt, xx, yy, ww, hh)
                n = int(np.count_nonzero(inside))
                if n == 0:
                    continue
                intersecting_blocks += 1
                inside_count += n
                sample = source_band.ReadAsArray(xx, yy, ww, hh)
                native_mask = validity_band.ReadAsArray(xx, yy, ww, hh)
                finite = np.isfinite(sample)
                valid = (native_mask != 0) & finite
                if nodata is not None and math.isfinite(nodata):
                    valid &= (sample != nodata)
                valid_count += int(np.count_nonzero(valid & inside))
                nonfinite_inside_count += int(np.count_nonzero(inside & ~finite))
    pixel_area_m2 = abs(gt[1] * gt[5])
    epsg = source_srs.GetAuthorityCode(None)
    record = {
        "filename": path.name,
        "source_epsg": epsg or "unknown",
        "source_nodata": nodata,
        "pixel_size_x_m": gt[1],
        "pixel_size_y_m": abs(gt[5]),
        "source_mask_flags": source_band.GetMaskFlags(),
        "step11_valid_30m_samples": int(step11_row["valid_samples_inside_aoi"]),
        "step11_30m_centres_in_tile_window": int(step11_row["grid_window_aoi_centres"]),
        "native_1m_aoi_pixel_centres_in_tile": inside_count,
        "native_1m_valid_pixel_centres_in_aoi": valid_count,
        "native_1m_missing_pixel_centres_in_aoi": inside_count - valid_count,
        "native_valid_pct_of_native_aoi_centres": round(100 * valid_count / inside_count, 6) if inside_count else None,
        "native_valid_km2_before_tile_overlap_dedup": round(valid_count * pixel_area_m2 / 1e6, 6),
        "native_aoi_grid_km2_before_tile_overlap_dedup": round(inside_count * pixel_area_m2 / 1e6, 6),
        "nonfinite_aoi_pixel_centres": nonfinite_inside_count,
        "examined_blocks": blocks,
        "aoi_intersecting_blocks": intersecting_blocks,
        "status": "NATIVE_PIXEL_AOI_SCAN_COMPLETE",
    }
    memory_ds = None
    ds = None
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all", action="store_true", help="After validating pilot, scan all 32 candidates")
    parser.add_argument("--chunk", type=int, default=1024, help="Pixels per source-native block edge; default 1024")
    args = parser.parse_args()
    if args.chunk < 128 or args.chunk > 2048:
        parser.error("--chunk must be between 128 and 2048")
    archive, ref_rows = load_reference()
    targets = list(ref_rows) if args.all else PILOT_FILES
    print("STEP 4.2 — ARCHIVED LIDAR NATIVE-RESOLUTION PIXEL QA")
    print(f"Mode: {'all_candidates' if args.all else 'three_tile_pilot'}")
    records = []
    errors = []
    for index, name in enumerate(targets, 1):
        print(f"[{index}/{len(targets)}] {name}", flush=True)
        try:
            if name not in ref_rows:
                raise KeyError("Not present in Step-11 candidate CSV")
            path = archive / name
            if not path.is_file():
                raise FileNotFoundError(path)
            rec = scan_native_tile(path, ref_rows[name], chunk=args.chunk)
            records.append(rec)
            print(
                "  1m AOI pixels: {native_1m_aoi_pixel_centres_in_tile:,}; "
                "valid: {native_1m_valid_pixel_centres_in_aoi:,}; "
                "valid pct: {native_valid_pct_of_native_aoi_centres}%".format(**rec),
                flush=True,
            )
        except Exception as exc:
            msg = f"{type(exc).__name__}: {exc}"
            errors.append({"filename": name, "error": msg})
            print("  ERROR: " + msg, flush=True)
    DIAGNOSTICS.mkdir(parents=True, exist_ok=True)
    if args.all:
        out_csv = DIAGNOSTICS / "step4_2_archived_lidar_native_pixel_all_by_tile.csv"
        out_json = DIAGNOSTICS / "step4_2_archived_lidar_native_pixel_all_qa.json"
    else:
        out_csv, out_json = OUT_CSV, OUT_JSON
    if records:
        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    report = {
        "status": "PILOT_COMPLETE_WITH_CAVEATS" if not errors and not args.all else (
            "FULL_PER_TILE_SCAN_COMPLETE_WITH_CAVEATS" if not errors else "SCAN_HAS_TILE_ERRORS"),
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "mode": "all_candidates" if args.all else "three_tile_pilot",
        "input_step11_csv": str(STEP11_CSV),
        "input_aoi": str(AOI_PATH),
        "archive_path": str(archive),
        "num_requested": len(targets),
        "num_succeeded": len(records),
        "num_errors": len(errors),
        "records": records,
        "errors": errors,
        "caveats": [
            "Native validity is counted at centres of original source raster pixels inside the AOI.",
            "Source tile geometries overlap: per-tile valid areas cannot be added to obtain AOI union coverage.",
            "Native-pixel coverage is not an independent DTM provenance or vertical-accuracy check.",
            "The Step-11 30-m nearest-neighbor result and native pixel counts are different sampling schemes.",
            "The source archive TIFFs are read-only; no mosaic or georeferenced output is generated.",
        ],
    }
    out_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print("Status:", report["status"])
    print("CSV:", out_csv if records else "not written (no successful tiles)")
    print("JSON:", out_json)
    print("CAVEAT: This is source-native per-tile QA, NOT a 1m union or DTM accuracy proof.")
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
