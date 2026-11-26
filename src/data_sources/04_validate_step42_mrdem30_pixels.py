"""Step 4.2: independent pixel-mask QA for the local MRDEM-30 cutline raster.

This script checks validity INSIDE the true study-area polygon, rather than
interpreting GDAL's rectangular-raster STATISTICS_VALID_PERCENT as AOI coverage.

Read-only inputs; writes one JSON diagnostic report. Uses GDAL/OGR and NumPy.
"""

from __future__ import annotations

import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from osgeo import gdal, ogr, osr


gdal.UseExceptions()

REPO = Path(__file__).resolve().parents[2]
DEM_PATH = REPO / "data/processed/elevation/mrdem30_dtm_highway99_aoi_utm10.tif"
AOI_PATH = REPO / "data/processed/study_area/highway99_terrain_search_5000m.geojson"
REPORT_PATH = REPO / "results/diagnostics/step4_2_mrdem30_pixel_qa.json"
TARGET_EPSG = 26910
EXPECTED_PIXEL_METRES = 30.0
EXPECTED_NODATA = -32767.0


def projected_aoi_layer(aoi_path: Path, target_srs: osr.SpatialReference):
    """Return an in-memory OGR datasource and its AOI layer in target CRS."""
    source_ds = ogr.Open(str(aoi_path))
    if source_ds is None:
        raise RuntimeError(f"Cannot open AOI: {aoi_path}")
    source_layer = source_ds.GetLayer(0)
    source_srs = source_layer.GetSpatialRef()
    if source_srs is None:
        raise RuntimeError("AOI GeoJSON has no declared spatial reference system")

    source_srs = source_srs.Clone()
    source_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    target_srs = target_srs.Clone()
    target_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    transform = osr.CoordinateTransformation(source_srs, target_srs)

    memory_driver = ogr.GetDriverByName("Memory")
    memory_ds = memory_driver.CreateDataSource("aoi_in_raster_crs")
    aoi_layer = memory_ds.CreateLayer("study_aoi", srs=target_srs, geom_type=ogr.wkbUnknown)

    source_layer.ResetReading()
    for feature in source_layer:
        geometry = feature.GetGeometryRef()
        if geometry is None:
            continue
        clone = geometry.Clone()
        result = clone.Transform(transform)
        if result != 0:
            raise RuntimeError(f"AOI coordinate transformation failed: code {result}")
        if clone.IsEmpty() or not clone.IsValid():
            raise RuntimeError("AOI has an empty or invalid polygon after reprojection")
        new_feature = ogr.Feature(aoi_layer.GetLayerDefn())
        new_feature.SetGeometry(clone)
        if aoi_layer.CreateFeature(new_feature) != 0:
            raise RuntimeError("Could not write projected AOI geometry")
        new_feature = None

    if aoi_layer.GetFeatureCount() == 0:
        raise RuntimeError("Projected AOI layer is empty")
    return memory_ds, aoi_layer, source_ds


def run() -> int:
    for p in (DEM_PATH, AOI_PATH):
        if not p.is_file():
            raise FileNotFoundError(f"Required file not found: {p}")

    dem = gdal.Open(str(DEM_PATH), gdal.GA_ReadOnly)
    if dem is None or dem.RasterCount != 1:
        raise RuntimeError("DEM must be a readable, single-band raster")

    geo = dem.GetGeoTransform()
    srs = osr.SpatialReference()
    if srs.ImportFromWkt(dem.GetProjection()) != 0:
        raise RuntimeError("DEM projection WKT is missing or invalid")
    target = osr.SpatialReference()
    target.ImportFromEPSG(TARGET_EPSG)
    srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    target.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    projected_ds, projected_layer, source_ds = projected_aoi_layer(AOI_PATH, srs)

    # Use the EXACT grid of the DEM; pixels are selected by their centres.
    mask_ds = gdal.GetDriverByName("MEM").Create("", dem.RasterXSize, dem.RasterYSize, 1, gdal.GDT_Byte)
    mask_ds.SetGeoTransform(geo)
    mask_ds.SetProjection(dem.GetProjection())
    mask_ds.GetRasterBand(1).Fill(0)
    burn_result = gdal.RasterizeLayer(mask_ds, [1], projected_layer, burn_values=[1])
    if burn_result != 0:
        raise RuntimeError(f"Rasterizing AOI failed with GDAL error {burn_result}")
    in_aoi = mask_ds.GetRasterBand(1).ReadAsArray() == 1

    band = dem.GetRasterBand(1)
    data = band.ReadAsArray()
    nodata = band.GetNoDataValue()
    # The GDAL mask accounts for nodata; also explicitly exclude NaN and Inf.
    valid = (band.GetMaskBand().ReadAsArray() != 0) & np.isfinite(data)
    if nodata is not None and not math.isnan(nodata):
        valid &= data != nodata

    inside_total = int(np.count_nonzero(in_aoi))
    inside_valid = int(np.count_nonzero(valid & in_aoi))
    inside_missing = inside_total - inside_valid
    outside_valid = int(np.count_nonzero(valid & ~in_aoi))
    total_valid = int(np.count_nonzero(valid))
    total_pixels = int(data.size)
    valid_values = data[valid & in_aoi]

    pixel_area_m2 = abs(geo[1] * geo[5] - geo[2] * geo[4])
    pixel_size_ok = (abs(geo[1] - EXPECTED_PIXEL_METRES) < 1e-7
                     and abs(geo[5] + EXPECTED_PIXEL_METRES) < 1e-7
                     and abs(geo[2]) < 1e-10 and abs(geo[4]) < 1e-10)
    crs_ok = bool(srs.IsSame(target))
    nodata_ok = nodata is not None and abs(nodata - EXPECTED_NODATA) < 1e-7

    checks = {
        "raster_crs_is_epsg_26910": crs_ok,
        "raster_pixel_size_is_30m_and_north_up": pixel_size_ok,
        "nodata_is_minus_32767": nodata_ok,
        "aoi_contains_raster_pixel_centres": inside_total > 0,
        "aoi_has_valid_elevation_pixels": inside_valid > 0,
        "all_aoi_pixel_centres_have_valid_elevation": inside_missing == 0,
        "no_valid_pixels_outside_aoi_pixel_centres": outside_valid == 0,
    }

    if not all(checks.values()):
        status = "NEEDS_REVIEW_MRDEM30_AOI_PIXEL_QA"
    else:
        status = "PASSED_INTERNAL_MRDEM30_PIXEL_MASK_QA_WITH_CAVEATS"

    report = {
        "step": "4.2 — MRDEM-30 local GeoTIFF pixel-mask QA",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "inputs": {"dem": str(DEM_PATH), "aoi": str(AOI_PATH)},
        "raster": {
            "width_pixels": dem.RasterXSize,
            "height_pixels": dem.RasterYSize,
            "epsg_expected": TARGET_EPSG,
            "pixel_size_m": [geo[1], geo[5]],
            "nodata_value": nodata,
            "pixel_area_m2": pixel_area_m2,
        },
        "pixel_counts": {
            "total_rectangle": total_pixels,
            "valid_rectangle": total_valid,
            "aoi_pixel_centres": inside_total,
            "valid_inside_aoi": inside_valid,
            "missing_inside_aoi": inside_missing,
            "valid_outside_aoi_pixel_centres": outside_valid,
        },
        "coverage": {
            "rectangle_valid_pct": round(100 * total_valid / total_pixels, 6),
            "aoi_valid_pct": round(100 * inside_valid / inside_total, 6) if inside_total else None,
            "aoi_missing_km2_on_30m_grid": round(inside_missing * pixel_area_m2 / 1e6, 9),
            "aoi_rasterized_km2_pixel_centres": round(inside_total * pixel_area_m2 / 1e6, 9),
            "valid_outside_aoi_km2_pixel_centres": round(outside_valid * pixel_area_m2 / 1e6, 9),
        },
        "elevation_inside_aoi_m": {
            "min": float(np.min(valid_values)) if valid_values.size else None,
            "median": float(np.median(valid_values)) if valid_values.size else None,
            "max": float(np.max(valid_values)) if valid_values.size else None,
        },
        "checks": checks,
        "caveats": [
            "AOI rasterization uses pixel-centre inclusion; boundary pixels can differ from exact vector area.",
            "This checks pixel masking/coverage, not vertical elevation accuracy, terrain artefacts, or landslide suitability.",
            "Values slightly below 0 m near coastal areas are not automatically errors and need spatial assessment.",
        ],
    }

    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print("STEP 4.2 — MRDEM-30 LOCAL RASTER PIXEL QA")
    print("Status:", status)
    for key, value in report["pixel_counts"].items():
        print(f"  {key}: {value}")
    for key, value in report["coverage"].items():
        print(f"  {key}: {value}")
    print("CHECKS:")
    for name, passed in checks.items():
        print(f"  {'PASS' if passed else 'REVIEW'}: {name}")
    print("Report:", REPORT_PATH)
    print("CAVEAT: pixel coverage is not a proof of elevation accuracy.")

    # Keep in-memory data sources alive until calculations are finished.
    _ = (projected_ds, source_ds)
    return 0 if status.startswith("PASSED_") else 2


if __name__ == "__main__":
    try:
        sys.exit(run())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
