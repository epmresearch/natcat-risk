"""Mask ONE existing HRDEM 1 m staging tile by the verified AOI-clipped
official HRDEM 1 m DTM footprint, without resampling or touching source tiles.

Pilot only. Run from the natcat-risk project root with --tile X_MIN Y_MIN.
Creates a NEW diagnostic masked TIFF and JSON report. No network downloads.
Does not fill missing elevations, create a final mosaic, or certify accuracy.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import os

import numpy as np
from osgeo import gdal, ogr, osr

NODATA = -32767.0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tile", nargs=2, type=int, required=True,
                        metavar=("X_MIN", "Y_MIN"))
    args = parser.parse_args()
    x, y = args.tile
    gdal.UseExceptions()
    ogr.UseExceptions()

    root = Path(__file__).resolve().parent
    stem = f"hrdem1m_x{x}_y{y}_utm10_1km"
    input_path = root / "results/diagnostics/hrdem1m_batch_staging" / (stem + ".tif")
    coverage_path = root / "data/processed/elevation/hrdem_1m_dtm_coverage_highway99.gpkg"
    output_dir = root / "results/diagnostics/hrdem1m_mask_pilots"
    output_path = output_dir / (stem + "_masked.tif")
    qa_path = output_dir / (stem + "_masked.qa.json")
    temp_path = output_dir / (stem + "_masked.partial.tif")

    if not input_path.is_file() or not coverage_path.is_file():
        raise FileNotFoundError("Missing existing staging tile or clipped official footprint")
    if output_path.exists() or qa_path.exists() or temp_path.exists():
        raise FileExistsError("Pilot output already exists; refusing to overwrite: " + str(output_path))

    src = gdal.Open(str(input_path), gdal.GA_ReadOnly)
    if src.RasterCount != 1 or (src.RasterXSize, src.RasterYSize) != (1000, 1000):
        raise ValueError("Staging input must be a one-band 1000x1000 tile")
    gt = src.GetGeoTransform()
    if not np.allclose(gt, (x, 1, 0, y + 1000, 0, -1), atol=1e-8, rtol=0):
        raise ValueError("Unexpected staging geotransform")
    crs = osr.SpatialReference(wkt=src.GetProjection())
    crs.AutoIdentifyEPSG()
    if crs.GetAuthorityCode(None) != "26910":
        raise ValueError("Staging tile CRS must be EPSG:26910")
    sb = src.GetRasterBand(1)
    if sb.DataType != gdal.GDT_Float32 or sb.GetNoDataValue() != NODATA:
        raise ValueError("Staging tile must be Float32 with NoData=-32767")

    vector = ogr.Open(str(coverage_path), 0)
    layer = vector.GetLayer(0)
    vsrs = layer.GetSpatialRef()
    if vsrs is None:
        raise ValueError("Coverage CRS missing")
    vsrs.AutoIdentifyEPSG()
    if vsrs.GetAuthorityCode(None) != "26910":
        raise ValueError("Coverage must be EPSG:26910")
    layer.SetSpatialFilterRect(x, y, x + 1000, y + 1000)
    mask_ds = gdal.GetDriverByName("MEM").Create("", 1000, 1000, 1, gdal.GDT_Byte)
    mask_ds.SetGeoTransform(gt)
    mask_ds.SetProjection(src.GetProjection())
    mask_ds.GetRasterBand(1).Fill(0)
    if gdal.RasterizeLayer(mask_ds, [1], layer, burn_values=[1]) != 0:
        raise RuntimeError("Coverage rasterization failed")

    footprint = mask_ds.ReadAsArray() == 1
    elevations = sb.ReadAsArray()
    valid = np.isfinite(elevations) & (elevations != NODATA)
    kept = valid & footprint
    masked = elevations.copy()
    masked[~kept] = NODATA

    qa = {
        "status": "DIAGNOSTIC_MASK_PILOT_NOT_FINAL",
        "tile_x": x, "tile_y": y, "crs": "EPSG:26910", "pixel_size_m": 1,
        "source_tile": str(input_path.relative_to(root)),
        "mask_source": str(coverage_path.relative_to(root)),
        "mask_rule": "pixel-center rasterization of AOI-clipped 1m DTM coverage",
        "source_valid_pixels": int(valid.sum()),
        "footprint_pixels": int(footprint.sum()),
        "valid_inside_footprint": int(kept.sum()),
        "nodata_inside_footprint": int((footprint & ~valid).sum()),
        "valid_outside_footprint_removed": int((valid & ~footprint).sum()),
        "masked_valid_pixels": int(kept.sum()),
        "masked_nodata_pixels": int(masked.size - kept.sum()),
        "retained_original_pixel_values_unchanged": True,
        "caveats": [
            "Footprint is AOI-clipped; unmasked source may contain genuine elevation beyond AOI",
            "NoData inside footprint remains NoData; no interpolation or gap filling",
            "This is not a final DEM or elevation accuracy validation",
        ],
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    try:
        out = gdal.GetDriverByName("GTiff").Create(
            str(temp_path), 1000, 1000, 1, gdal.GDT_Float32,
            options=["TILED=YES", "COMPRESS=DEFLATE", "PREDICTOR=3"]
        )
        if out is None:
            raise RuntimeError("Cannot create masked pilot TIFF")
        out.SetGeoTransform(gt)
        out.SetProjection(src.GetProjection())
        out.GetRasterBand(1).SetNoDataValue(NODATA)
        out.GetRasterBand(1).WriteArray(masked)
        out.FlushCache()
        out = None

        check = gdal.Open(str(temp_path), gdal.GA_ReadOnly)
        readback = check.ReadAsArray()
        if not np.array_equal(readback, masked):
            raise RuntimeError("Output pixel values differ from expected masked array")
        check = None

        os.replace(temp_path, output_path)
        with qa_path.open("x", encoding="utf-8") as f:
            json.dump(qa, f, indent=2)
            f.write("\n")
    finally:
        if temp_path.exists():
            temp_path.unlink()

    for key in (
        "source_valid_pixels", "footprint_pixels", "valid_inside_footprint",
        "nodata_inside_footprint", "valid_outside_footprint_removed",
        "masked_valid_pixels", "masked_nodata_pixels",
    ):
        print(key.upper() + ":", qa[key])
    print("ORIGINAL_VALUES_UNCHANGED_FOR_KEPT_PIXELS: True")
    print("MASK_PILOT_QA_PASSED_WITH_CAVEATS")
    print("OUTPUT:", output_path)
    print("QA_JSON:", qa_path)
    print("NO_NETWORK_READS_AND_NO_ORIGINAL_TILES_MODIFIED")


if __name__ == "__main__":
    main()
