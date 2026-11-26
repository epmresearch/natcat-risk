"""Step 4.2 — diagnose MRDEM-30/AOI pixel-mask disagreements near the boundary.

Read-only inputs: the cropped MRDEM-30 GeoTIFF and the original study AOI GeoJSON.
Only writes a JSON diagnostic report (does not rewrite, fill or resample the DEM).

Rasterizes the AOI and 30/60 m interior/exterior buffers on EXACTLY the DEM grid
using GDAL's default pixel-centre rule, then measures where the 920 missing-inside
and 933 valid-outside pixels from the first QA fall.

The result is diagnostic, NOT an elevation-accuracy assessment or QA override.
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
AOI = ROOT / "data/processed/study_area/highway99_terrain_search_5000m.geojson"
DEM = ROOT / "data/processed/elevation/mrdem30_dtm_highway99_aoi_utm10.tif"
REPORT = ROOT / "results/diagnostics/step4_2_mrdem30_boundary_diagnostics.json"
TARGET_EPSG = 26910
BUFFERS_METRES = (-60, -30, 0, 30, 60)


def reproject_and_union_aoi(path: Path, dst_srs: osr.SpatialReference):
    source = ogr.Open(str(path), 0)
    if source is None:
        raise RuntimeError(f"Cannot read AOI: {path}")
    layer = source.GetLayer(0)
    src_srs = layer.GetSpatialRef()
    if src_srs is None:
        raise RuntimeError("AOI has no declared CRS")
    src_srs = src_srs.Clone()
    src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    dst_srs = dst_srs.Clone()
    dst_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    transformer = osr.CoordinateTransformation(src_srs, dst_srs)

    geometries = []
    layer.ResetReading()
    for f in layer:
        geom = f.GetGeometryRef()
        if geom is None or geom.IsEmpty():
            continue
        candidate = geom.Clone()
        if candidate.Transform(transformer) != 0:
            raise RuntimeError("AOI reprojection failed")
        if not candidate.IsValid():
            raise RuntimeError("AOI geometry invalid after reprojection")
        geometries.append(candidate)

    if not geometries:
        raise RuntimeError("AOI has no nonempty polygon features")
    merged = geometries[0]
    for other in geometries[1:]:
        merged = merged.Union(other)
        if merged is None:
            raise RuntimeError("Failed to union multi-feature AOI")
    if merged is None or merged.IsEmpty() or not merged.IsValid():
        raise RuntimeError("Projected AOI union is empty or invalid")
    return merged, dst_srs, source


def rasterize_geometry_on_dem_grid(dem, geometry, srs):
    """Burn a single geometry using GDAL's default pixel-centre rule."""
    driver = ogr.GetDriverByName("Memory")
    vector_ds = driver.CreateDataSource("geometry_in_memory")
    vector_lyr = vector_ds.CreateLayer("mask", srs=srs, geom_type=ogr.wkbUnknown)
    feature = ogr.Feature(vector_lyr.GetLayerDefn())
    feature.SetGeometry(geometry)
    if vector_lyr.CreateFeature(feature) != 0:
        raise RuntimeError("Could not add AOI polygon to in-memory layer")
    feature = None

    mask_ds = gdal.GetDriverByName("MEM").Create(
        "", dem.RasterXSize, dem.RasterYSize, 1, gdal.GDT_Byte
    )
    if mask_ds is None:
        raise RuntimeError("Could not create in-memory raster")
    mask_ds.SetGeoTransform(dem.GetGeoTransform())
    mask_ds.SetProjection(dem.GetProjection())
    mask_ds.GetRasterBand(1).Fill(0)
    result = gdal.RasterizeLayer(mask_ds, [1], vector_lyr, burn_values=[1])
    if result != 0:
        raise RuntimeError(f"AOI rasterization error {result}")
    array = mask_ds.GetRasterBand(1).ReadAsArray() != 0
    # Keep the vector source alive until ReadAsArray completes.
    _ = vector_ds
    return array


def run():
    for input_path in (AOI, DEM):
        if not input_path.is_file():
            raise FileNotFoundError(f"Missing input: {input_path}")

    dem = gdal.Open(str(DEM), gdal.GA_ReadOnly)
    if dem is None or dem.RasterCount != 1:
        raise RuntimeError("DEM must be a readable one-band GeoTIFF")
    raster_srs = osr.SpatialReference()
    if raster_srs.ImportFromWkt(dem.GetProjection()) != 0:
        raise RuntimeError("DEM projection missing")
    raster_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    expected_srs = osr.SpatialReference()
    expected_srs.ImportFromEPSG(TARGET_EPSG)
    expected_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    if not raster_srs.IsSame(expected_srs):
        raise RuntimeError("Wrong DEM CRS; expected EPSG:26910")
    gt = dem.GetGeoTransform()
    if not (
        abs(gt[1] - 30.0) < 1e-7
        and abs(gt[5] + 30.0) < 1e-7
        and abs(gt[2]) < 1e-10
        and abs(gt[4]) < 1e-10
    ):
        raise RuntimeError("DEM raster must have a north-up 30 m grid")

    band = dem.GetRasterBand(1)
    data = band.ReadAsArray()
    nodata = band.GetNoDataValue()
    valid = (band.GetMaskBand().ReadAsArray() != 0) & np.isfinite(data)
    if nodata is not None and not math.isnan(nodata):
        valid &= data != nodata

    projected_aoi, projected_srs, source_ds = reproject_and_union_aoi(AOI, raster_srs)

    masks = {}
    for buffer_metres in BUFFERS_METRES:
        buffered = projected_aoi.Buffer(float(buffer_metres)) if buffer_metres else projected_aoi
        if buffered is None or buffered.IsEmpty():
            raise RuntimeError(f"AOI becomes empty when buffered by {buffer_metres} m")
        if not buffered.IsValid():
            raise RuntimeError(f"Invalid AOI buffer at {buffer_metres} m")
        masks[buffer_metres] = rasterize_geometry_on_dem_grid(dem, buffered, projected_srs)

    original = masks[0]
    missing_inside = ~valid & original
    valid_outside = valid & ~original
    original_missing = int(np.count_nonzero(missing_inside))
    original_outside = int(np.count_nonzero(valid_outside))

    missing_in_30m_core = int(np.count_nonzero(~valid & masks[-30]))
    missing_in_60m_core = int(np.count_nonzero(~valid & masks[-60]))
    valid_beyond_30m = int(np.count_nonzero(valid & ~masks[30]))
    valid_beyond_60m = int(np.count_nonzero(valid & ~masks[60]))

    # Classification uses actual buffer masks, not assumed pixel-distance approximations.
    metrics = {
        "reference_aoi_pixel_centres": int(np.count_nonzero(original)),
        "reference_missing_inside": original_missing,
        "reference_valid_outside": original_outside,
        "missing_inside_30m_eroded_aoi": missing_in_30m_core,
        "missing_inside_60m_eroded_aoi": missing_in_60m_core,
        "missing_within_first_30m_boundary_strip": original_missing - missing_in_30m_core,
        "missing_in_30m_to_60m_interior_strip": missing_in_30m_core - missing_in_60m_core,
        "valid_beyond_30m_expanded_aoi": valid_beyond_30m,
        "valid_beyond_60m_expanded_aoi": valid_beyond_60m,
        "valid_within_first_30m_exterior_strip": original_outside - valid_beyond_30m,
        "valid_in_30m_to_60m_exterior_strip": valid_beyond_30m - valid_beyond_60m,
        "missing_inside_60m_core_km2_30m_grid": round(missing_in_60m_core * 900 / 1e6, 9),
        "valid_beyond_60m_km2_30m_grid": round(valid_beyond_60m * 900 / 1e6, 9),
    }

    baseline_is_consistent = original_missing == 920 and original_outside == 933
    if original_missing == 0 and original_outside == 0:
        status = "MASK_MATCHES_AOI_PIXEL_CENTRES"
    elif missing_in_60m_core == 0 and valid_beyond_60m == 0:
        status = "MISMATCHES_CONFINED_TO_60M_BOUNDARY_ZONE"
    else:
        status = "INTERIOR_OR_DISTANT_PIXELS_NEED_INVESTIGATION"

    report = {
        "step": "4.2 — MRDEM-30 AOI boundary mismatch location diagnostics",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "inputs": {"aoi": str(AOI), "dem": str(DEM)},
        "method": {
            "grid_crs": "EPSG:26910",
            "pixel_size_metres": 30,
            "aoi_buffers_metres": list(BUFFERS_METRES),
            "rasterization": "GDAL default polygon pixel-centre inclusion, no ALL_TOUCHED",
        },
        "baseline_matches_previous_report_920_and_933": baseline_is_consistent,
        "metrics": metrics,
        "interpretation_limits": [
            "Buffer-zone confinement suggests, but does not prove, a cutline rasterization boundary effect.",
            "This does not certify complete source DEM coverage or elevation accuracy.",
            "No input raster or source AOI was modified. No gaps were filled.",
        ],
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print("STEP 4.2 — MRDEM-30 BOUNDARY PIXEL DIAGNOSTICS")
    print("Status:", status)
    print("Previous QA mismatch counts reproducible:", "YES" if baseline_is_consistent else "NO — REVIEW")
    for key, value in metrics.items():
        print(f"  {key}: {value}")
    print("Report:", REPORT)
    print("CAUTION: Boundary proximity is not a proof of DEM pixel correctness.")

    _ = source_ds  # keep datasource alive while using its data
    return 0 if baseline_is_consistent else 2


if __name__ == "__main__":
    try:
        sys.exit(run())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
