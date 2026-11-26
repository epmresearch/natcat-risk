"""Step 4.2: controlled MRDEM-30 exact-warp comparison (read-only inputs).

Compare two exact reprojections of the same Canadian MRDEM-30 COG:
  - no polygon clip: -et 0 -ovr NONE -te ... -tr 30 30
  - polygon clip:    same options plus -cutline AOI
Also compare earlier default warp outputs to help isolate changed values.

Only produces results/diagnostics/step4_2_mrdem30_exact_warp_comparison.json.
No raster, source AOI, or other project file is modified. This is *not*
verification of absolute elevation accuracy, vertical datum or landslide risk.
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
FILES = {
    "original_cutline": ROOT / "data/processed/elevation/mrdem30_dtm_highway99_aoi_utm10.tif",
    "original_bbox": ROOT / "results/diagnostics/step4_2_mrdem30_uncut_bbox_diagnostic.tif",
    "exact_bbox": ROOT / "results/diagnostics/step4_2_mrdem30_exact_bbox.tif",
    "exact_cutline": ROOT / "results/diagnostics/step4_2_mrdem30_exact_cutline.tif",
}
AOI = ROOT / "data/processed/study_area/highway99_terrain_search_5000m.geojson"
REPORT = ROOT / "results/diagnostics/step4_2_mrdem30_exact_warp_comparison.json"


def read_raster(path: Path):
    ds = gdal.Open(str(path), gdal.GA_ReadOnly)
    if ds is None or ds.RasterCount != 1:
        raise RuntimeError(f"Not a readable single-band raster: {path}")
    band = ds.GetRasterBand(1)
    pixels = band.ReadAsArray()
    nodata = band.GetNoDataValue()
    valid = np.isfinite(pixels) & (band.GetMaskBand().ReadAsArray() != 0)
    if nodata is not None and not math.isnan(nodata):
        valid &= (pixels != nodata)
    srs = osr.SpatialReference()
    if srs.ImportFromWkt(ds.GetProjection()) != 0:
        raise RuntimeError(f"Missing CRS: {path}")
    srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return ds, pixels, valid, srs, nodata


def mask_from_aoi(grid_ds: gdal.Dataset, target_srs: osr.SpatialReference):
    source = ogr.Open(str(AOI), 0)
    if source is None:
        raise RuntimeError(f"Cannot open AOI: {AOI}")
    source_lyr = source.GetLayer(0)
    source_srs = source_lyr.GetSpatialRef()
    if source_srs is None:
        raise RuntimeError("AOI has no declared SRS")
    source_srs = source_srs.Clone()
    source_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    transformer = osr.CoordinateTransformation(source_srs, target_srs)

    memory_vector = ogr.GetDriverByName("Memory").CreateDataSource("aoi")
    lyr = memory_vector.CreateLayer("aoi", srs=target_srs, geom_type=ogr.wkbUnknown)
    any_features = False
    for source_feature in source_lyr:
        geometry = source_feature.GetGeometryRef()
        if geometry is None or geometry.IsEmpty():
            continue
        copied = geometry.Clone()
        if copied.Transform(transformer) != 0 or copied.IsEmpty() or not copied.IsValid():
            raise RuntimeError("AOI reprojection yielded invalid geometry")
        feature = ogr.Feature(lyr.GetLayerDefn())
        feature.SetGeometry(copied)
        if lyr.CreateFeature(feature) != 0:
            raise RuntimeError("Cannot add AOI feature to in-memory layer")
        feature = None
        any_features = True
    if not any_features:
        raise RuntimeError("AOI contains no polygon features")

    mem = gdal.GetDriverByName("MEM").Create("", grid_ds.RasterXSize, grid_ds.RasterYSize, 1, gdal.GDT_Byte)
    mem.SetGeoTransform(grid_ds.GetGeoTransform())
    mem.SetProjection(grid_ds.GetProjection())
    mem.GetRasterBand(1).Fill(0)
    if gdal.RasterizeLayer(mem, [1], lyr, burn_values=[1]) != 0:
        raise RuntimeError("Rasterization of AOI failed")
    inside = (mem.GetRasterBand(1).ReadAsArray() != 0)
    return inside


def comparison(a, av, b, bv, region=None):
    comparable = av & bv
    if region is not None:
        comparable &= region
    n = int(np.count_nonzero(comparable))
    if n == 0:
        return {"shared_valid_pixels": 0, "error": "No jointly-valid pixels in comparison area"}
    abs_delta = np.abs(a[comparable].astype(np.float64) - b[comparable].astype(np.float64))
    return {
        "shared_valid_pixels": n,
        "identical_float32_values": int(np.count_nonzero(abs_delta == 0)),
        "different_over_0_001m": int(np.count_nonzero(abs_delta > 0.001)),
        "different_over_0_1m": int(np.count_nonzero(abs_delta > 0.1)),
        "different_over_1m": int(np.count_nonzero(abs_delta > 1)),
        "different_over_5m": int(np.count_nonzero(abs_delta > 5)),
        "mean_absolute_difference_m": round(float(abs_delta.mean()), 9),
        "median_absolute_difference_m": round(float(np.median(abs_delta)), 9),
        "p95_absolute_difference_m": round(float(np.percentile(abs_delta, 95)), 9),
        "max_absolute_difference_m": round(float(abs_delta.max()), 9),
    }


def main() -> int:
    for path in (*FILES.values(), AOI):
        if not path.is_file():
            raise FileNotFoundError(f"Required file missing: {path}")

    data = {key: read_raster(path) for key, path in FILES.items()}
    reference = data["exact_bbox"]
    ref_ds, _, _, ref_srs, ref_nodata = reference
    expected = osr.SpatialReference()
    expected.ImportFromEPSG(26910)
    expected.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    if not ref_srs.IsSame(expected):
        raise RuntimeError("Exact output CRS is not EPSG:26910")

    grid_checks = {}
    for key, (ds, _arr, _valid, srs, nodata) in data.items():
        check = {
            "same_shape": (ds.RasterXSize, ds.RasterYSize) == (ref_ds.RasterXSize, ref_ds.RasterYSize),
            "same_crs": bool(srs.IsSame(ref_srs)),
            "same_geotransform": all(abs(x-y) < 1e-8 for x, y in zip(ds.GetGeoTransform(), ref_ds.GetGeoTransform())),
            "same_nodata": nodata == ref_nodata,
        }
        grid_checks[key] = check
        if not all(check.values()):
            raise RuntimeError(f"Unsafe raster grid mismatch for {key}: {check}")

    aoi_pixels = mask_from_aoi(ref_ds, ref_srs)
    counts = {"aoi_pixel_centres": int(np.count_nonzero(aoi_pixels))}
    for key, (_, _, valid, _, _) in data.items():
        counts[key] = {
            "valid_inside_aoi": int(np.count_nonzero(valid & aoi_pixels)),
            "missing_inside_aoi": int(np.count_nonzero(~valid & aoi_pixels)),
            "valid_outside_aoi": int(np.count_nonzero(valid & ~aoi_pixels)),
        }

    def pair(k1, k2, region=None):
        return comparison(data[k1][1], data[k1][2], data[k2][1], data[k2][2], region)

    pairs = {
        "exact_cutline_vs_exact_bbox_all_shared": pair("exact_cutline", "exact_bbox"),
        "exact_cutline_vs_exact_bbox_inside_aoi": pair("exact_cutline", "exact_bbox", aoi_pixels),
        "original_cutline_vs_original_bbox_all_shared": pair("original_cutline", "original_bbox"),
        "original_cutline_vs_exact_cutline_all_shared": pair("original_cutline", "exact_cutline"),
        "original_bbox_vs_exact_bbox_all_shared": pair("original_bbox", "exact_bbox"),
    }
    exact_shared = pairs["exact_cutline_vs_exact_bbox_all_shared"]
    if exact_shared.get("shared_valid_pixels", 0) == 0:
        status = "REVIEW_NO_SHARED_EXACT_PIXELS"
    elif exact_shared["different_over_0_001m"] == 0:
        status = "EXACT_PAIR_MATCHES_WITHIN_1MM_ON_SHARED_VALID_PIXELS"
    else:
        status = "EXACT_PAIR_VALUES_DIFFER_NEEDS_REVIEW"

    report = {
        "step": "4.2 — MRDEM-30 exact warp control comparison",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "inputs": {k: str(v) for k, v in FILES.items()} | {"aoi": str(AOI)},
        "method": "GDAL mask validity and pixel-centre AOI rasterization; paired absolute elevation differences",
        "grid_checks": grid_checks,
        "pixel_counts": counts,
        "comparisons": pairs,
        "limits": [
            "Both exact outputs are resampled from the same MRDEM-30 product, not independent ground truth.",
            "An identical shared-pixel value is not proof of vertical accuracy or suitable slope estimates.",
            "The cutline pixel inclusion rule may differ from the default AOI rasterization pixel-centre rule.",
            "The uncut rectangular DEM is a diagnostic; it includes valid pixels outside AOI.",
            "Outputs only a JSON report. Input rasters and AOI are not modified.",
        ],
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print("STEP 4.2 — EXACT MRDEM-30 WARP COMPARISON")
    print("Status:", status)
    print("Grid matches for all four rasters: YES")
    print("AOI pixel centres:", counts["aoi_pixel_centres"])
    for key, item in counts.items():
        if key != "aoi_pixel_centres":
            print(f"  {key}: inside_valid={item['valid_inside_aoi']}, inside_missing={item['missing_inside_aoi']}, outside_valid={item['valid_outside_aoi']}")
    for key, item in pairs.items():
        print(f"  {key}:")
        if "error" in item:
            print("    ", item["error"])
        else:
            print(f"    shared={item['shared_valid_pixels']}, differ_gt_1mm={item['different_over_0_001m']}, differ_gt_1m={item['different_over_1m']}, max_abs_diff_m={item['max_absolute_difference_m']}")
    print("Report:", REPORT)
    print("CAUTION: Reproducibility comparison only; not ground-truth elevation accuracy.")
    return 0 if status == "EXACT_PAIR_MATCHES_WITHIN_1MM_ON_SHARED_VALID_PIXELS" else 2


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
