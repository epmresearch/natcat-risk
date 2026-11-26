"""Step 4.2: create an AOI-masked copy of the single exact-bbox MRDEM-30 warp.

The point is to apply an explicitly rasterized AOI mask WITHOUT a second
resampling/coordinate warp. The source remains unchanged; the new masked raster
must have exactly the same Float32 elevations wherever AOI pixels are included.

This verifies reproducible pixel masking and internal completeness, NOT
geodetic/vertical accuracy, landslide suitability, or original DEM quality.

Inputs:
    results/diagnostics/step4_2_mrdem30_exact_bbox.tif
    data/processed/study_area/highway99_terrain_search_5000m.geojson
Outputs (no overwriting):
    data/processed/elevation/mrdem30_dtm_highway99_aoi_utm10_masked_exact.tif
    results/diagnostics/step4_2_mrdem30_consistent_mask_qa.json
"""
from __future__ import annotations

import json
import math
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from osgeo import gdal, ogr, osr


gdal.UseExceptions()
ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "results/diagnostics/step4_2_mrdem30_exact_bbox.tif"
AOI = ROOT / "data/processed/study_area/highway99_terrain_search_5000m.geojson"
OUTPUT = ROOT / "data/processed/elevation/mrdem30_dtm_highway99_aoi_utm10_masked_exact.tif"
REPORT = ROOT / "results/diagnostics/step4_2_mrdem30_consistent_mask_qa.json"
EXPECTED_GRID = (1084, 3093)
EXPECTED_AOI_CENTRE_PIXELS = 1_159_588
EXPECTED_NODATA = -32767.0


def project_aoi_to_memory(target_srs: osr.SpatialReference):
    """Return projected OGR datasource + layer, keeping sources alive."""
    ds = ogr.Open(str(AOI), 0)
    if ds is None:
        raise RuntimeError(f"Cannot read AOI GeoJSON: {AOI}")
    source_layer = ds.GetLayer(0)
    source_srs = source_layer.GetSpatialRef()
    if source_srs is None:
        raise RuntimeError("AOI must have a declared CRS")
    source_srs = source_srs.Clone()
    source_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    transformer = osr.CoordinateTransformation(source_srs, target_srs)

    mem_vector = ogr.GetDriverByName("Memory").CreateDataSource("aoi_projected")
    mem_layer = mem_vector.CreateLayer("study_aoi", srs=target_srs, geom_type=ogr.wkbUnknown)
    feature_count = 0
    source_layer.ResetReading()
    for feature in source_layer:
        geom = feature.GetGeometryRef()
        if geom is None or geom.IsEmpty():
            continue
        projected = geom.Clone()
        if projected.Transform(transformer) != 0:
            raise RuntimeError("AOI geometry coordinate transformation failed")
        if projected.IsEmpty() or not projected.IsValid():
            raise RuntimeError("AOI geometry is empty/invalid after reprojection")
        out_feature = ogr.Feature(mem_layer.GetLayerDefn())
        out_feature.SetGeometry(projected)
        if mem_layer.CreateFeature(out_feature) != 0:
            raise RuntimeError("Cannot create reprojected AOI feature")
        out_feature = None
        feature_count += 1
    if feature_count == 0:
        raise RuntimeError("AOI contains no usable polygons")
    return ds, mem_vector, mem_layer, feature_count


def make_aoi_mask(src: gdal.Dataset, target_srs: osr.SpatialReference):
    vector_source, projected_source, layer, n_features = project_aoi_to_memory(target_srs)
    mask_ds = gdal.GetDriverByName("MEM").Create(
        "", src.RasterXSize, src.RasterYSize, 1, gdal.GDT_Byte
    )
    mask_ds.SetProjection(src.GetProjection())
    mask_ds.SetGeoTransform(src.GetGeoTransform())
    mask_ds.GetRasterBand(1).Fill(0)
    if gdal.RasterizeLayer(mask_ds, [1], layer, burn_values=[1]) != 0:
        raise RuntimeError("GDAL could not rasterize AOI on exact source grid")
    mask = mask_ds.GetRasterBand(1).ReadAsArray() != 0
    # Keep vector datasets alive until GDAL rasterization completes.
    _ = (vector_source, projected_source, layer, n_features)
    return mask, n_features


def run() -> int:
    for path in (SOURCE, AOI):
        if not path.is_file():
            raise FileNotFoundError(f"Required input is missing: {path}")
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {OUTPUT}")
    if REPORT.exists():
        raise FileExistsError(f"Refusing to overwrite existing report: {REPORT}")

    src = gdal.Open(str(SOURCE), gdal.GA_ReadOnly)
    if src is None or src.RasterCount != 1:
        raise RuntimeError("Exact-bbox source must be one readable raster band")
    band = src.GetRasterBand(1)
    if band.DataType != gdal.GDT_Float32:
        raise RuntimeError("Expected Float32 source raster")
    nodata = band.GetNoDataValue()
    if nodata is None or abs(nodata - EXPECTED_NODATA) > 1e-7:
        raise RuntimeError(f"Unexpected source NoData marker: {nodata}")
    if (src.RasterXSize, src.RasterYSize) != EXPECTED_GRID:
        raise RuntimeError("Raster size is inconsistent with prior validated grid")

    srs = osr.SpatialReference()
    if srs.ImportFromWkt(src.GetProjection()) != 0:
        raise RuntimeError("Raster CRS is missing or invalid")
    srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    expected_srs = osr.SpatialReference()
    expected_srs.ImportFromEPSG(26910)
    expected_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    if not srs.IsSame(expected_srs):
        raise RuntimeError("Source raster CRS must be EPSG:26910")
    gt = src.GetGeoTransform()
    if not (abs(gt[1] - 30) < 1e-8 and abs(gt[5] + 30) < 1e-8
            and abs(gt[2]) < 1e-10 and abs(gt[4]) < 1e-10):
        raise RuntimeError(f"Raster grid is not north-up 30m: {gt}")

    inside, feature_count = make_aoi_mask(src, srs)
    inside_n = int(np.count_nonzero(inside))
    if inside_n != EXPECTED_AOI_CENTRE_PIXELS:
        raise RuntimeError(
            f"AOI pixel-center count changed: {inside_n}; "
            f"expected {EXPECTED_AOI_CENTRE_PIXELS}. Do not silently proceed."
        )
    values = band.ReadAsArray()
    if values is None or values.shape != inside.shape:
        raise RuntimeError("Source array read failed")
    valid_source = (
        (band.GetMaskBand().ReadAsArray() != 0)
        & np.isfinite(values)
        & (values != nodata)
    )
    source_missing_inside = int(np.count_nonzero(inside & ~valid_source))
    if source_missing_inside:
        raise RuntimeError(
            f"Source has {source_missing_inside} missing pixels inside AOI. "
            "Cannot claim complete AOI coverage."
        )

    # Direct copy: do not call gdal.Warp or any resampling algorithm here.
    masked = np.full(values.shape, np.float32(nodata), dtype=np.float32)
    masked[inside] = values[inside]

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    staging = OUTPUT.with_name(OUTPUT.stem + ".partial.tif")
    if staging.exists():
        raise FileExistsError(f"Existing partial output needs review: {staging}")
    writer = None
    try:
        writer = gdal.GetDriverByName("GTiff").Create(
            str(staging), src.RasterXSize, src.RasterYSize, 1, gdal.GDT_Float32,
            options=[
                "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256",
                "COMPRESS=DEFLATE", "PREDICTOR=3", "BIGTIFF=IF_SAFER",
            ],
        )
        if writer is None:
            raise RuntimeError("Cannot create output GeoTIFF")
        writer.SetProjection(src.GetProjection())
        writer.SetGeoTransform(gt)
        writer.GetRasterBand(1).SetNoDataValue(nodata)
        writer.GetRasterBand(1).WriteArray(masked)
        writer.SetMetadata({
            "SOURCE_RASTER": "step4_2_mrdem30_exact_bbox.tif",
            "MASK_RULE": "AOI polygon GDAL RasterizeLayer default pixel-centre inclusion",
            "RESAMPLING_AT_MASK_STAGE": "NONE",
        })
        writer.FlushCache()
        writer = None

        check_ds = gdal.Open(str(staging), gdal.GA_ReadOnly)
        if check_ds is None:
            raise RuntimeError("Cannot reopen masked raster")
        check_band = check_ds.GetRasterBand(1)
        written = check_band.ReadAsArray()
        out_valid = (
            (check_band.GetMaskBand().ReadAsArray() != 0)
            & np.isfinite(written)
            & (written != nodata)
        )
        same_grid = (
            check_ds.RasterXSize == src.RasterXSize
            and check_ds.RasterYSize == src.RasterYSize
            and check_ds.GetProjection() == src.GetProjection()
            and all(abs(a - b) < 1e-8 for a, b in zip(check_ds.GetGeoTransform(), gt))
        )
        valid_inside = int(np.count_nonzero(out_valid & inside))
        missing_inside = int(np.count_nonzero(~out_valid & inside))
        valid_outside = int(np.count_nonzero(out_valid & ~inside))
        mismatch_count = int(np.count_nonzero(written[inside] != values[inside]))
        check_nodata = check_band.GetNoDataValue() == nodata
        checks = {
            "source_crs_epsg_26910_and_30m": True,
            "aoi_mask_matches_previous_qa": inside_n == EXPECTED_AOI_CENTRE_PIXELS,
            "source_has_no_missing_inside_aoi": source_missing_inside == 0,
            "output_grid_matches_source": same_grid,
            "output_nodata_is_minus_32767": check_nodata,
            "all_aoi_centres_have_valid_elevation": missing_inside == 0,
            "no_valid_pixels_outside_aoi_centres": valid_outside == 0,
            "all_inside_elevation_values_identical_to_source": mismatch_count == 0,
            "valid_pixel_count_matches_aoi_pixel_centres": valid_inside == inside_n,
        }
        if not all(checks.values()):
            raise RuntimeError("Masked raster QA failed: " + repr(checks))
        check_ds = None
        os.rename(staging, OUTPUT)  # no overwrite: OUTPUT existence checked above
    except Exception:
        writer = None
        # Keep failed partial artifact for inspection, rather than exposing it as final.
        raise

    status = "PASSED_INTERNAL_MRDEM30_CONSISTENT_MASK_QA_WITH_CAVEATS"
    report = {
        "step": "4.2 — deterministic AOI mask from single exact-bbox MRDEM-30 warp",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "source_exact_bbox": str(SOURCE),
        "aoi_polygon": str(AOI),
        "output_masked": str(OUTPUT),
        "method": "Reproject AOI to DEM CRS, rasterize onto exact grid with pixel-centre rule, and mask Float32 pixels without resampling",
        "source_aoi_feature_count": feature_count,
        "counts": {
            "raster_pixels_total": int(values.size),
            "aoi_pixel_centres": inside_n,
            "source_valid_inside_aoi": inside_n - source_missing_inside,
            "source_missing_inside_aoi": source_missing_inside,
            "output_valid_inside_aoi": valid_inside,
            "output_missing_inside_aoi": missing_inside,
            "output_valid_outside_aoi": valid_outside,
            "inside_values_differing_from_source": mismatch_count,
        },
        "checks": checks,
        "caveats": [
            "The selected reference is the uncut exact-bbox warp; its absolute elevation and vertical datum are not ground-truth validated.",
            "100% pixel-centre coverage follows from valid exact-bbox pixels and an explicit deterministic mask; it is not independent vertical accuracy evidence.",
            "Rasterized AOI pixel-centre area differs slightly from continuous polygon area at 30m resolution.",
            "For slope and other neighborhood terrain derivatives, use the uncut/buffered DEM as input and mask derivatives afterward to avoid AOI boundary artifacts.",
            "The underlying cause of different bilinear values from separately warped cutline and bbox inputs has not been conclusively established.",
        ],
    }
    REPORT.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("STEP 4.2 — CONSISTENT MRDEM-30 AOI MASK (NO RESAMPLING)")
    print("Status:", status)
    print("  aoi_pixel_centres:", inside_n)
    print("  source_missing_inside_aoi:", source_missing_inside)
    print("  output_valid_inside_aoi:", valid_inside)
    print("  output_missing_inside_aoi:", missing_inside)
    print("  output_valid_outside_aoi:", valid_outside)
    print("  elevations_changed_inside_aoi:", mismatch_count)
    print("CHECKS:")
    for name, passed in checks.items():
        print(f"  {'PASS' if passed else 'FAIL'}: {name}")
    print("Output:", OUTPUT)
    print("Report:", REPORT)
    print("CAVEAT: internal pixel-mask consistency does not establish absolute elevation accuracy.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(run())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
