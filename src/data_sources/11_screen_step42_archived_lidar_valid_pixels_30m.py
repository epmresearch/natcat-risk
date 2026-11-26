#!/usr/bin/env python3
"""Step 4.2: sample archived bc_* LiDAR tile validity on the fixed 30-m AOI grid.

READ-ONLY with respect to archived LiDAR and existing project rasters.
Only the 32 *candidate* bc_ tiles from the verified Step 10 inventory are read.
The derived lidar_mask01.tif is ALWAYS excluded.

This is a 30-m pixel-centre SCREENING, NOT 1-m pixel coverage, DTM provenance,
vertical accuracy, or a replacement for inspecting original LiDAR validity masks.

Run from the sea_to_sky_risk environment:
    python src/data_sources/11_screen_step42_archived_lidar_valid_pixels_30m.py
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
ARCHIVE = (ROOT.parent / "old files backup" / "Climada" / "Climada_Project"
           / "Inputs" / "DEM" / "LiDAR DEMs")
SOURCE_LIST = ROOT / "results/diagnostics/step4_2_archived_lidar_footprints_without_mask_qa.json"
REFERENCE = ROOT / "data/processed/elevation/mrdem30_dtm_highway99_aoi_utm10_masked_exact.tif"
OUTDIR = ROOT / "results/diagnostics"
CSV_OUT = OUTDIR / "step4_2_archived_lidar_30m_sample_validity_by_tile.csv"
JSON_OUT = OUTDIR / "step4_2_archived_lidar_30m_sample_validity_qa.json"
TIFF_OUT = OUTDIR / "step4_2_archived_lidar_30m_sample_coverage_count.tif"
EXPECTED_AOI_PIXELS = 1_159_588
EXPECTED_CANDIDATES = 32
EXPECTED_PIXEL_M = 30.0
TARGET_CRS = CRS.from_epsg(26910)
OUT_NODATA = -32767.0


def check_reference():
    if not REFERENCE.is_file():
        raise FileNotFoundError(f"Validated Step 08 MRDEM-30 masked raster missing: {REFERENCE}")
    ds = gdal.Open(str(REFERENCE), gdal.GA_ReadOnly)
    if ds is None or ds.RasterCount != 1:
        raise ValueError("Reference must be a readable one-band GeoTIFF")
    srs = osr.SpatialReference()
    if srs.ImportFromWkt(ds.GetProjection()) != 0:
        raise ValueError("Reference has no valid CRS")
    desired = osr.SpatialReference()
    desired.ImportFromEPSG(26910)
    if not bool(srs.IsSame(desired)):
        raise ValueError("Reference must have EPSG:26910")
    gt = ds.GetGeoTransform()
    if not (abs(gt[1] - 30) < 1e-7 and abs(gt[5] + 30) < 1e-7
            and abs(gt[2]) < 1e-9 and abs(gt[4]) < 1e-9):
        raise ValueError(f"Unexpected reference 30-m grid: {gt}")
    band = ds.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    if nodata != OUT_NODATA:
        raise ValueError(f"Unexpected reference NoData: {nodata}")
    arr = band.ReadAsArray()
    if arr is None:
        raise ValueError("Cannot read reference AOI mask")
    inside = ((band.GetMaskBand().ReadAsArray() != 0)
              & np.isfinite(arr) & (arr != nodata))
    count = int(np.count_nonzero(inside))
    if count != EXPECTED_AOI_PIXELS:
        raise ValueError(f"Reference AOI pixel count has changed ({count}); expected {EXPECTED_AOI_PIXELS}")
    return ds, inside, gt


def candidate_names():
    if not SOURCE_LIST.is_file():
        raise FileNotFoundError(f"Step 10 inventory required: {SOURCE_LIST}")
    info = json.loads(SOURCE_LIST.read_text(encoding="utf-8"))
    entries = info.get("source_tile_footprints_in_aoi")
    if not isinstance(entries, list):
        raise ValueError("Step 10 inventory has no source_tile_footprints_in_aoi array")
    names = [entry.get("filename") for entry in entries]
    if (len(names) != EXPECTED_CANDIDATES or len(set(names)) != EXPECTED_CANDIDATES
            or not all(isinstance(name, str) and name.lower().startswith("bc_")
                       and name.lower().endswith((".tif", ".tiff")) for name in names)):
        raise ValueError(f"Candidate inventory changed: {len(names)}. Re-run Step 10 before screening")
    return sorted(names), info


def tile_target_window(src, ref_gt, ncols, nrows):
    """Densify source edges when projecting into EPSG:26910; pad by 2 cells."""
    src_wkt = src.GetProjection()
    if not src_wkt:
        raise ValueError("Source tile has no CRS")
    transformer = Transformer.from_crs(CRS.from_wkt(src_wkt), TARGET_CRS,
                                       always_xy=True)
    nx, ny = src.RasterXSize, src.RasterYSize
    gt = src.GetGeoTransform()
    corners = [((0, 0), (nx, 0)), ((nx, 0), (nx, ny)),
               ((nx, ny), (0, ny)), ((0, ny), (0, 0))]
    xs, ys = [], []
    for (x0, y0), (x1, y1) in corners:
        for t in np.linspace(0, 1, 21):
            col, row = x0 + t * (x1 - x0), y0 + t * (y1 - y0)
            xs.append(gt[0] + col * gt[1] + row * gt[2])
            ys.append(gt[3] + col * gt[4] + row * gt[5])
    xx, yy = transformer.transform(xs, ys)
    if not all(math.isfinite(float(v)) for v in (*xx, *yy)):
        raise ValueError("Non-finite projected tile footprint")
    x_min, x_max = min(xx), max(xx)
    y_min, y_max = min(yy), max(yy)
    left, top = ref_gt[0], ref_gt[3]
    col0 = max(0, math.floor((x_min - left) / EXPECTED_PIXEL_M) - 2)
    col1 = min(ncols, math.ceil((x_max - left) / EXPECTED_PIXEL_M) + 2)
    row0 = max(0, math.floor((top - y_max) / EXPECTED_PIXEL_M) - 2)
    row1 = min(nrows, math.ceil((top - y_min) / EXPECTED_PIXEL_M) + 2)
    if col0 >= col1 or row0 >= row1:
        return None
    bounds = (left + col0 * 30, top - row1 * 30,
              left + col1 * 30, top - row0 * 30)
    return row0, row1, col0, col1, bounds


def screen():
    for required in (ARCHIVE, SOURCE_LIST, REFERENCE):
        if not required.exists():
            raise FileNotFoundError(str(required))
    for path in (CSV_OUT, JSON_OUT, TIFF_OUT):
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite previous diagnostic output: {path}")
    names, inventory = candidate_names()
    reference, inside, ref_gt = check_reference()
    rows, cols = inside.shape
    sample_count = np.zeros((rows, cols), dtype=np.uint8)
    records, errors, cautions = [], [], []

    for index, name in enumerate(names, 1):
        path = ARCHIVE / name
        src = None
        warped = None
        try:
            if not path.is_file():
                raise FileNotFoundError(str(path))
            src = gdal.Open(str(path), gdal.GA_ReadOnly)
            if src is None or src.RasterCount != 1:
                raise ValueError("Expected a readable single-band raster")
            band = src.GetRasterBand(1)
            source_nodata = band.GetNoDataValue()
            source_mask_flags = int(band.GetMaskFlags())
            if source_nodata is None:
                cautions.append(f"{name}: missing declared source NoData; interpret samples with care")
            window = tile_target_window(src, ref_gt, cols, rows)
            if window is None:
                raise ValueError("Tile no longer overlaps target grid")
            row0, row1, col0, col1, bounds = window
            opts = gdal.WarpOptions(
                format="MEM", dstSRS="EPSG:26910",
                outputBounds=list(bounds), width=col1-col0, height=row1-row0,
                resampleAlg="near", errorThreshold=0.0, overviewLevel="NONE",
                outputType=gdal.GDT_Float32, dstNodata=OUT_NODATA,
                srcNodata=source_nodata if source_nodata is not None else None,
                multithread=False,
            )
            warped = gdal.Warp("", src, options=opts)
            if warped is None:
                raise RuntimeError("GDAL Warp to 30-m sample grid failed")
            if (warped.RasterXSize != col1-col0 or
                    warped.RasterYSize != row1-row0 or
                    any(abs(a-b) > 1e-6 for a, b in zip(
                        warped.GetGeoTransform(),
                        (bounds[0], 30.0, 0.0, bounds[3], 0.0, -30.0)))):
                raise RuntimeError("Warped sample is not aligned to reference grid")
            outband = warped.GetRasterBand(1)
            samples = outband.ReadAsArray()
            if samples is None:
                raise RuntimeError("Could not read sampled raster pixels")
            valid = ((outband.GetMaskBand().ReadAsArray() != 0)
                     & np.isfinite(samples) & (samples != OUT_NODATA))
            subset_aoi = inside[row0:row1, col0:col1]
            usable = valid & subset_aoi
            tile_valid = int(usable.sum())
            all_aoi_centres_within_window = int(subset_aoi.sum())
            current = sample_count[row0:row1, col0:col1]
            new = int(np.count_nonzero(usable & (current == 0)))
            duplicated = int(np.count_nonzero(usable & (current != 0)))
            current[usable] += 1
            records.append({
                "filename": name,
                "source_crs_epsg": CRS.from_wkt(src.GetProjection()).to_epsg(),
                "source_nodata": source_nodata,
                "source_mask_flags": source_mask_flags,
                "sampling_resolution_m": 30,
                "grid_window_aoi_centres": all_aoi_centres_within_window,
                "valid_samples_inside_aoi": tile_valid,
                "new_union_valid_samples": new,
                "overlapping_previously_valid_samples": duplicated,
                "window_cols": col1-col0,
                "window_rows": row1-row0,
                "status": "SAMPLED",
                "error": "",
            })
        except Exception as exc:
            errors.append({"filename": name, "error": f"{type(exc).__name__}: {exc}"})
            records.append({"filename": name, "status": "ERROR", "error": str(exc)})
        finally:
            warped = None
            src = None
        if index % 5 == 0 or index == len(names):
            print(f"Sampled {index}/{len(names)} archived bc_ candidate tiles")

    aoi_n = int(inside.sum())
    valid_union = int(np.count_nonzero(inside & (sample_count > 0)))
    no_sample = aoi_n - valid_union
    overlapping = int(np.count_nonzero(inside & (sample_count >= 2)))
    max_layers = int(sample_count[inside].max())
    status = ("ARCHIVED_LIDAR_30M_PIXEL_CENTRE_SCREENING_COMPLETE_WITH_CAVEATS"
              if not errors else "ARCHIVED_LIDAR_30M_SCREENING_INCOMPLETE_ERRORS")
    result = {
        "status": status,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_candidate_inventory": str(SOURCE_LIST),
        "reference_aoi_mask": str(REFERENCE),
        "archive_path": str(ARCHIVE),
        "excluded_derived_mask": "lidar_mask01.tif",
        "source_tile_filename_prefix": "bc_",
        "candidate_tile_count": len(names),
        "tile_samples_success": len(names) - len(errors),
        "tile_samples_error": len(errors),
        "source_rectangular_footprint_union_pct_from_step10": inventory.get("union_bc_tiles_only_aoi_pct"),
        "aoi_30m_pixel_centres": aoi_n,
        "union_valid_30m_nearest_samples_inside_aoi": valid_union,
        "union_valid_30m_nearest_samples_pct_inside_aoi": round(100.0 * valid_union / aoi_n, 6),
        "union_valid_30m_nearest_samples_km2": round(valid_union * 900.0 / 1e6, 6),
        "aoi_30m_centres_with_no_valid_tile_sample": no_sample,
        "aoi_no_valid_tile_sample_km2": round(no_sample * 900.0 / 1e6, 6),
        "aoi_centres_valid_in_two_or_more_candidate_tiles": overlapping,
        "maximum_valid_candidate_tiles_at_one_centre": max_layers,
        "method": ("Warp each candidate from its source pixels to a 30-m UTM 10N "
                   "window aligned to the exact MRDEM AOI grid, using nearest-neighbor "
                   "and source validity masks / NoData, with original-resolution source (overview NONE). "
                   "Union binary valid samples only at centres inside the validated AOI."),
        "notes": [
            "This is 30-m point sampling; it DOES NOT measure complete 1-m LiDAR valid-pixel coverage.",
            "Filenames bc_* are only candidate source tiles, not independently confirmed DTM products.",
            "No elevation comparison or vertical accuracy check is performed.",
            "Tile overlaps and different acquisition years are not resolved; no mosaic is constructed.",
            "Source archive and original project DEMs are never modified.",
        ],
        "cautions": cautions,
        "tile_errors": errors,
        "tile_records": records,
    }
    OUTDIR.mkdir(parents=True, exist_ok=True)
    with CSV_OUT.open("w", newline="", encoding="utf-8-sig") as fh:
        columns = ["filename", "source_crs_epsg", "source_nodata", "source_mask_flags",
                   "sampling_resolution_m", "grid_window_aoi_centres", "valid_samples_inside_aoi",
                   "new_union_valid_samples", "overlapping_previously_valid_samples",
                   "window_cols", "window_rows", "status", "error"]
        writer = csv.DictWriter(fh, columns)
        writer.writeheader()
        for record in records:
            writer.writerow(record)

    # Diagnostic raster: pixel value = number of bc_ candidate tiles with a valid
    # 30-m nearest-neighbor sample; 0 = no valid sample INSIDE AOI; 255 = outside.
    tif = gdal.GetDriverByName("GTiff").Create(
        str(TIFF_OUT), cols, rows, 1, gdal.GDT_Byte,
        options=["TILED=YES", "COMPRESS=DEFLATE", "BIGTIFF=IF_SAFER"])
    if tif is None:
        raise RuntimeError("Could not create diagnostic coverage-count raster")
    tif.SetProjection(reference.GetProjection())
    tif.SetGeoTransform(ref_gt)
    band = tif.GetRasterBand(1)
    band.SetNoDataValue(255)
    coverage = np.where(inside, sample_count, 255).astype(np.uint8)
    band.WriteArray(coverage)
    band.SetDescription("Count of valid 30-m nearest samples from candidate archived bc_ LiDAR tiles")
    tif.FlushCache()
    tif = None
    reference = None
    JSON_OUT.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")

    print("\nSTEP 4.2 — ARCHIVED LIDAR ACTUAL PIXEL SAMPLING (30m GRID)")
    for key in ("status", "candidate_tile_count", "tile_samples_success", "tile_samples_error",
                "source_rectangular_footprint_union_pct_from_step10", "aoi_30m_pixel_centres",
                "union_valid_30m_nearest_samples_inside_aoi",
                "union_valid_30m_nearest_samples_pct_inside_aoi",
                "aoi_30m_centres_with_no_valid_tile_sample",
                "aoi_centres_valid_in_two_or_more_candidate_tiles",
                "maximum_valid_candidate_tiles_at_one_centre"):
        print(f"  {key}: {result[key]}")
    if cautions:
        print("SOURCE METADATA CAUTIONS:")
        for item in cautions:
            print(" ", item)
    if errors:
        print("TILE READ/WARP ERRORS:")
        for item in errors:
            print(" ", item)
    print("CSV:", CSV_OUT)
    print("JSON:", JSON_OUT)
    print("Coverage-count raster:", TIFF_OUT)
    print("CAUTION: 30-m nearest-valid sampling is not full 1-m pixel coverage or elevation accuracy.")
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(screen())
