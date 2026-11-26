"""Step 4.2 — compare clipped and uncut MRDEM-30 on their identical 30 m grid.

The comparison tests whether previously observed near-AOI-boundary NoData in the
clipped product is also NoData in the independent *uncut warped raster*.

Read-only inputs: clipped GeoTIFF, uncut diagnostic GeoTIFF, study AOI GeoJSON.
Writes only results/diagnostics/step4_2_mrdem30_clipped_vs_uncut_qa.json.
Does not fill gaps, edit rasters, or assert vertical elevation accuracy.
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
ROOT = Path(__file__).resolve().parents[2]
CLIPPED = ROOT / "data/processed/elevation/mrdem30_dtm_highway99_aoi_utm10.tif"
UNCUT = ROOT / "results/diagnostics/step4_2_mrdem30_uncut_bbox_diagnostic.tif"
AOI = ROOT / "data/processed/study_area/highway99_terrain_search_5000m.geojson"
REPORT = ROOT / "results/diagnostics/step4_2_mrdem30_clipped_vs_uncut_qa.json"


def raster_validity(ds: gdal.Dataset) -> tuple[np.ndarray, np.ndarray]:
    band = ds.GetRasterBand(1)
    data = band.ReadAsArray()
    nodata = band.GetNoDataValue()
    valid = np.isfinite(data) & (band.GetMaskBand().ReadAsArray() != 0)
    if nodata is not None and not math.isnan(nodata):
        valid &= data != nodata
    return data, valid


def rasterized_aoi_mask(ds: gdal.Dataset) -> np.ndarray:
    source = ogr.Open(str(AOI), 0)
    if source is None:
        raise RuntimeError(f"Cannot read AOI: {AOI}")
    source_layer = source.GetLayer(0)
    src_srs = source_layer.GetSpatialRef()
    if src_srs is None:
        raise RuntimeError("AOI CRS is missing")
    src_srs = src_srs.Clone()
    src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    target_srs = osr.SpatialReference()
    if target_srs.ImportFromWkt(ds.GetProjection()) != 0:
        raise RuntimeError("Raster CRS invalid")
    target_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    transform = osr.CoordinateTransformation(src_srs, target_srs)

    memory_source = ogr.GetDriverByName("Memory").CreateDataSource("aoimask")
    layer = memory_source.CreateLayer("aoi", srs=target_srs, geom_type=ogr.wkbUnknown)
    source_layer.ResetReading()
    for f in source_layer:
        geom = f.GetGeometryRef()
        if geom is None or geom.IsEmpty():
            continue
        transformed = geom.Clone()
        if transformed.Transform(transform) != 0 or transformed.IsEmpty() or not transformed.IsValid():
            raise RuntimeError("AOI projection failed or is invalid")
        out_feature = ogr.Feature(layer.GetLayerDefn())
        out_feature.SetGeometry(transformed)
        if layer.CreateFeature(out_feature) != 0:
            raise RuntimeError("Failed to create AOI mask feature")
        out_feature = None
    if layer.GetFeatureCount() == 0:
        raise RuntimeError("No AOI polygon to rasterize")

    mem_raster = gdal.GetDriverByName("MEM").Create("", ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Byte)
    mem_raster.SetGeoTransform(ds.GetGeoTransform())
    mem_raster.SetProjection(ds.GetProjection())
    mem_raster.GetRasterBand(1).Fill(0)
    if gdal.RasterizeLayer(mem_raster, [1], layer, burn_values=[1]) != 0:
        raise RuntimeError("AOI rasterization failed")
    mask = mem_raster.GetRasterBand(1).ReadAsArray() != 0
    _ = (source, memory_source)
    return mask


def run() -> int:
    for path in (CLIPPED, UNCUT, AOI):
        if not path.is_file():
            raise FileNotFoundError(f"Required input not found: {path}")
    a = gdal.Open(str(CLIPPED), gdal.GA_ReadOnly)
    b = gdal.Open(str(UNCUT), gdal.GA_ReadOnly)
    if a is None or b is None or a.RasterCount != 1 or b.RasterCount != 1:
        raise RuntimeError("Both raster files must be readable single-band datasets")
    sa, sb = osr.SpatialReference(), osr.SpatialReference()
    if sa.ImportFromWkt(a.GetProjection()) != 0 or sb.ImportFromWkt(b.GetProjection()) != 0:
        raise RuntimeError("Missing CRS")
    sa.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    sb.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    grid_checks = {
        "same_dimensions": a.RasterXSize == b.RasterXSize and a.RasterYSize == b.RasterYSize,
        "same_crs": bool(sa.IsSame(sb)),
        "same_geotransform": all(abs(x-y) < 1e-7 for x, y in zip(a.GetGeoTransform(), b.GetGeoTransform())),
        "same_nodata": a.GetRasterBand(1).GetNoDataValue() == b.GetRasterBand(1).GetNoDataValue(),
    }
    if not all(grid_checks.values()):
        raise RuntimeError(f"Raster grids/metadata differ, comparison unsafe: {grid_checks}")

    clipped, clip_valid = raster_validity(a)
    uncut, uncut_valid = raster_validity(b)
    inside = rasterized_aoi_mask(a)

    clip_missing_inside = inside & ~clip_valid
    clip_valid_outside = ~inside & clip_valid
    both_valid = clip_valid & uncut_valid
    only_clip_valid = clip_valid & ~uncut_valid
    missing_inside_uncut_valid = clip_missing_inside & uncut_valid
    missing_inside_uncut_missing = clip_missing_inside & ~uncut_valid
    clipped_valid_outside_uncut_valid = clip_valid_outside & uncut_valid
    clipped_valid_outside_uncut_missing = clip_valid_outside & ~uncut_valid

    differences = np.abs(clipped[both_valid].astype(np.float64) - uncut[both_valid].astype(np.float64))
    max_abs_difference_m = float(np.max(differences)) if differences.size else None
    num_different_over_1mm = int(np.count_nonzero(differences > 0.001))
    recovered = int(np.count_nonzero(missing_inside_uncut_valid))
    original_missing = int(np.count_nonzero(clip_missing_inside))
    missing_in_both = int(np.count_nonzero(missing_inside_uncut_missing))

    metrics = {
        "aoi_pixel_centres": int(np.count_nonzero(inside)),
        "clipped_valid_inside_aoi": int(np.count_nonzero(inside & clip_valid)),
        "clipped_missing_inside_aoi": original_missing,
        "of_clipped_missing_inside_valid_in_uncut": recovered,
        "of_clipped_missing_inside_also_missing_in_uncut": missing_in_both,
        "uncut_valid_inside_aoi": int(np.count_nonzero(inside & uncut_valid)),
        "uncut_missing_inside_aoi": int(np.count_nonzero(inside & ~uncut_valid)),
        "clipped_valid_outside_aoi": int(np.count_nonzero(clip_valid_outside)),
        "of_clipped_valid_outside_valid_in_uncut": int(np.count_nonzero(clipped_valid_outside_uncut_valid)),
        "of_clipped_valid_outside_missing_in_uncut": int(np.count_nonzero(clipped_valid_outside_uncut_missing)),
        "clipped_valid_but_uncut_missing_anywhere": int(np.count_nonzero(only_clip_valid)),
        "shared_valid_pixels": int(np.count_nonzero(both_valid)),
        "shared_valid_max_absolute_elevation_difference_m": max_abs_difference_m,
        "shared_valid_pixels_differing_over_1mm": num_different_over_1mm,
    }
    baseline_reproduced = (
        metrics["aoi_pixel_centres"] == 1159588
        and metrics["clipped_missing_inside_aoi"] == 920
        and metrics["clipped_valid_outside_aoi"] == 933
    )
    if original_missing == recovered and num_different_over_1mm == 0 and metrics["clipped_valid_but_uncut_missing_anywhere"] == 0:
        status = "ALL_CLIPPED_AOI_GAPS_VALID_IN_UNCUT_REFERENCE"
    elif 0 < recovered < original_missing:
        status = "MIXED_CUTLINE_AND_UNCUT_REFERENCE_NODATA"
    else:
        status = "UNCUT_REFERENCE_COMPARISON_NEEDS_REVIEW"

    report = {
        "step": "4.2 — MRDEM-30 clipped vs uncut reference comparison",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "inputs": {"clipped": str(CLIPPED), "uncut": str(UNCUT), "aoi": str(AOI)},
        "grid_checks": grid_checks,
        "previous_qa_counts_reproduced": baseline_reproduced,
        "pixel_metrics": metrics,
        "interpretation_limits": [
            "Uncut reference is another warp of the same remote MRDEM-30 COG, not ground truth or an independent elevation product.",
            "If clipped AOI gaps are valid in the uncut raster, this supports a cutline/masking-boundary explanation, not vertical accuracy.",
            "AOI inclusion is determined by GDAL's default rasterization, not all-touched coverage.",
            "Elevation ranges and values must be reviewed separately for scientific use.",
        ],
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("STEP 4.2 — MRDEM-30 CLIPPED / UNCUT RASTER COMPARISON")
    print("Status:", status)
    print("Previous 920 / 933 mismatch counts reproduced:", "YES" if baseline_reproduced else "NO")
    for key, passed in grid_checks.items():
        print(f"  {'PASS' if passed else 'FAIL'}: {key}")
    print("PIXEL METRICS:")
    for key, value in metrics.items():
        print(f"  {key}: {value}")
    print("Report:", REPORT)
    print("CAUTION: This checks consistency/masking, not vertical DEM accuracy.")
    return 0 if baseline_reproduced and status == "ALL_CLIPPED_AOI_GAPS_VALID_IN_UNCUT_REFERENCE" else 2


if __name__ == "__main__":
    try:
        sys.exit(run())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
