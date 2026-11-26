"""Mask existing official HRDEM 1 m staging tiles to the verified AOI-clipped footprint.

Run from the natcat-risk project root. By default, PLAN ONLY: read local filenames,
print counts and do not write files. Use --batch N (1..10) or --tile X Y to mask.

The original staging rasters are never modified, no network data are accessed,
no missing heights are interpolated, and no final DEM/mosaic is produced.
All outputs are diagnostic until cross-tile and full-corridor QA is complete.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import sys

import numpy as np
from osgeo import gdal, ogr, osr

NODATA = -32767.0
TILE_SIZE = 1000
PATTERN = re.compile(r"hrdem1m_x(\d+)_y(\d+)_utm10_1km\.tif$")


def validate_grid(ds, x: int, y: int, label: str):
    if ds is None:
        raise RuntimeError(f"Cannot read {label}")
    if (ds.RasterXSize, ds.RasterYSize, ds.RasterCount) != (1000, 1000, 1):
        raise ValueError(f"Unexpected dimensions/bands: {label}")
    expected = (float(x), 1.0, 0.0, float(y + 1000), 0.0, -1.0)
    if not np.allclose(ds.GetGeoTransform(), expected, atol=1e-8, rtol=0):
        raise ValueError(f"Unexpected geotransform: {label}")
    srs = osr.SpatialReference(wkt=ds.GetProjection())
    srs.AutoIdentifyEPSG()
    if srs.GetAuthorityCode(None) != "26910":
        raise ValueError(f"Expected EPSG:26910: {label}")
    band = ds.GetRasterBand(1)
    if band.DataType != gdal.GDT_Float32 or band.GetNoDataValue() != NODATA:
        raise ValueError(f"Expected Float32 with -32767 NoData: {label}")
    return band


def rasterized_footprint(coverage: Path, ds, x: int, y: int) -> np.ndarray:
    vector = ogr.Open(str(coverage), 0)
    if vector is None:
        raise RuntimeError(f"Cannot open footprint {coverage}")
    layer = vector.GetLayer(0)
    srs = layer.GetSpatialRef()
    if srs is None:
        raise RuntimeError("Footprint CRS missing")
    srs.AutoIdentifyEPSG()
    if srs.GetAuthorityCode(None) != "26910":
        raise RuntimeError("Footprint must use EPSG:26910")
    layer.SetSpatialFilterRect(x, y, x + TILE_SIZE, y + TILE_SIZE)
    mem = gdal.GetDriverByName("MEM").Create("", TILE_SIZE, TILE_SIZE, 1, gdal.GDT_Byte)
    mem.SetGeoTransform(ds.GetGeoTransform())
    mem.SetProjection(ds.GetProjection())
    mem.GetRasterBand(1).Fill(0)
    status = gdal.RasterizeLayer(mem, [1], layer, burn_values=[1])
    if status != 0:
        raise RuntimeError("Footprint rasterization failed")
    mask = mem.ReadAsArray() == 1
    return mask


def inspect_source_and_mask(source: Path, coverage: Path, x: int, y: int):
    ds = gdal.Open(str(source), gdal.GA_ReadOnly)
    band = validate_grid(ds, x, y, str(source))
    elevations = band.ReadAsArray()
    valid = np.isfinite(elevations) & (elevations != NODATA)
    inside = rasterized_footprint(coverage, ds, x, y)
    kept = valid & inside
    masked = elevations.copy()
    masked[~kept] = NODATA
    stats = {
        "source_valid_pixels": int(valid.sum()),
        "footprint_pixels": int(inside.sum()),
        "valid_inside_footprint": int(kept.sum()),
        "nodata_inside_footprint": int((inside & ~valid).sum()),
        "valid_outside_footprint_removed": int((valid & ~inside).sum()),
        "masked_valid_pixels": int(kept.sum()),
        "masked_nodata_pixels": int(masked.size - kept.sum()),
    }
    proj, geo = ds.GetProjection(), ds.GetGeoTransform()
    ds = None
    return masked, stats, proj, geo


def verify_masked_file(path: Path, masked: np.ndarray, x: int, y: int):
    ds = gdal.Open(str(path), gdal.GA_ReadOnly)
    validate_grid(ds, x, y, str(path))
    if not np.array_equal(ds.ReadAsArray(), masked):
        raise RuntimeError(f"Masked values do not exactly match expected mask: {path}")
    ds = None


def verify_existing_pair(root: Path, source: Path, coverage: Path, output_dir: Path, x: int, y: int):
    output = output_dir / (source.stem + "_masked.tif")
    qa_file = output_dir / (source.stem + "_masked.qa.json")
    expected, stats, _, _ = inspect_source_and_mask(source, coverage, x, y)
    verify_masked_file(output, expected, x, y)
    with qa_file.open(encoding="utf-8") as handle:
        previous = json.load(handle)
    if any(previous.get(key) != value for key, value in stats.items()):
        raise RuntimeError(f"Existing QA JSON does not match source and footprint: {qa_file}")
    return stats


def write_json_atomically(path: Path, payload: dict):
    tmp = path.with_name(path.name + ".partial")
    if tmp.exists():
        raise RuntimeError(f"Incomplete JSON exists; investigate before retry: {tmp}")
    try:
        with tmp.open("x", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)
            handle.write("\n")
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def run_one(root: Path, source: Path, coverage: Path, output_dir: Path, x: int, y: int):
    stem = source.stem
    output = output_dir / (stem + "_masked.tif")
    qa_file = output_dir / (stem + "_masked.qa.json")
    temp = output_dir / (stem + "_masked.partial.tif")
    if temp.exists():
        raise RuntimeError(f"Partial output exists; investigate before retry: {temp}")
    if output.exists() != qa_file.exists():
        raise RuntimeError(f"Incomplete output/QA pair; investigate: {stem}")

    expected, stats, projection, geotransform = inspect_source_and_mask(source, coverage, x, y)
    if output.exists():
        verify_masked_file(output, expected, x, y)
        with qa_file.open(encoding="utf-8") as f:
            previous = json.load(f)
        if any(previous.get(key) != value for key, value in stats.items()):
            raise RuntimeError(f"Existing QA JSON does not match source and footprint: {qa_file}")
        print(f"SKIPPED_VALID_EXISTING: {output.name} kept={stats['masked_valid_pixels']}", flush=True)
        return "SKIPPED"

    driver = gdal.GetDriverByName("GTiff")
    try:
        out = driver.Create(str(temp), TILE_SIZE, TILE_SIZE, 1, gdal.GDT_Float32,
                            options=["TILED=YES", "COMPRESS=DEFLATE", "PREDICTOR=3"])
        if out is None:
            raise RuntimeError(f"Could not create {temp}")
        out.SetGeoTransform(geotransform)
        out.SetProjection(projection)
        out.GetRasterBand(1).SetNoDataValue(NODATA)
        out.GetRasterBand(1).WriteArray(expected)
        out.FlushCache()
        out = None
        verify_masked_file(temp, expected, x, y)
        payload = {
            "status": "MASKED_STAGING_QA_PASSED_WITH_CAVEATS",
            "tile_x": x,
            "tile_y": y,
            "source_staging_tif": str(source.relative_to(root)),
            "mask_source": str(coverage.relative_to(root)),
            "mask_rule": "pixel-center rasterization of AOI-clipped official HRDEM 1m footprint",
            "crs": "EPSG:26910",
            "pixel_size_m": 1,
            "original_values_unchanged_for_retained_pixels": True,
            **stats,
            "caveats": [
                "AOI-clipped footprint, not the complete natural LiDAR footprint",
                "Missing heights inside footprint are preserved as NoData",
                "Staging only; vertical accuracy, inter-tile seam QA and final mosaic pending",
            ],
        }
        os.replace(temp, output)
        write_json_atomically(qa_file, payload)
        print(f"MASKED_QA_PASSED_WITH_CAVEATS: {output.name} kept={stats['masked_valid_pixels']} removed={stats['valid_outside_footprint_removed']} missing_inside={stats['nodata_inside_footprint']}", flush=True)
        return "SUCCESS"
    finally:
        if temp.exists():
            temp.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--batch", type=int, metavar="N", help="Mask at most N pending staging tiles (1..10)")
    modes.add_argument("--tile", type=int, nargs=2, metavar=("X_MIN", "Y_MIN"), help="Mask one existing staging tile")
    args = parser.parse_args()
    if args.batch is not None and not 1 <= args.batch <= 10:
        parser.error("--batch must be between 1 and 10")
    gdal.UseExceptions()
    ogr.UseExceptions()
    root = Path(__file__).resolve().parent
    source_dir = root / "results/diagnostics/hrdem1m_batch_staging"
    coverage = root / "data/processed/elevation/hrdem_1m_dtm_coverage_highway99.gpkg"
    output_dir = root / "results/diagnostics/hrdem1m_masked_staging"
    if not source_dir.is_dir() or not coverage.is_file():
        raise FileNotFoundError("Missing staging directory or verified AOI-clipped footprint")
    staged = []
    for p in source_dir.glob("*.tif"):
        match = PATTERN.fullmatch(p.name)
        if match:
            staged.append((int(match.group(1)), int(match.group(2)), p))
        elif not p.name.endswith(".partial.tif"):
            raise ValueError(f"Unexpected TIFF in staging directory: {p.name}")
    staged.sort(key=lambda item: (-item[1], item[0]))
    print("STAGING_FILES_AVAILABLE:", len(staged))
    print("EXISTING_MASKED_OUTPUT_FILES:", len(list(output_dir.glob("*_masked.tif"))) if output_dir.is_dir() else 0)
    print("MASK_MODE: AOI_CLIPPED_HRDEM_1M_FOOTPRINT_PIXEL_CENTERS")
    print("INPUTS: LOCAL_STAGING_ONLY_NO_NETWORK")
    if args.tile is None and args.batch is None:
        print("FIRST_5_STAGED_TILES:", [(x, y) for x, y, _ in staged[:5]])
        print("PLAN_ONLY: no files downloaded or changed")
        print("NEXT_SAFE_TEST: python step4_2_hrdem1m_mask_batch.py --batch 3")
        return
    if args.tile is not None:
        matches = [(x, y, p) for x, y, p in staged if (x, y) == tuple(args.tile)]
        if not matches:
            raise ValueError(f"Tile does not exist in local staging: {tuple(args.tile)}")
        chosen = matches
    else:
        chosen = []
        for x, y, path in staged:
            out = output_dir / (path.stem + "_masked.tif")
            qa = output_dir / (path.stem + "_masked.qa.json")
            if out.exists() != qa.exists():
                raise RuntimeError(f"Incomplete masked output/QA pair: {path.stem}")
            if out.exists():
                qa_stats = verify_existing_pair(root, path, coverage, output_dir, x, y)
                print(f"ALREADY_MASKED_VERIFIED: {out.name} kept={qa_stats['masked_valid_pixels']}")
                continue
            chosen.append((x, y, path))
            if len(chosen) == args.batch:
                break
    print("MASK_TILES_SELECTED:", len(chosen))
    if chosen:
        output_dir.mkdir(parents=True, exist_ok=True)
    counts = {"SUCCESS": 0, "SKIPPED": 0}
    for i, (x, y, path) in enumerate(chosen, 1):
        print(f"MASK_PROGRESS: {i}/{len(chosen)} tile=({x}, {y})", flush=True)
        status = run_one(root, path, coverage, output_dir, x, y)
        counts[status] += 1
    print("MASK_BATCH_FINISHED:", counts)
    print("NOTE: diagnostic masked staging only; no final mosaic or vertical accuracy certification")


if __name__ == "__main__":
    try:
        main()
    except (Exception, KeyboardInterrupt) as exc:
        print("ERROR:", exc, file=sys.stderr)
        raise SystemExit(1)
