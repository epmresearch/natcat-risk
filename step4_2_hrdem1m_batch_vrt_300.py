"""Safely stage official HRDEM 1 m DTM tiles for Highway 99 (Step 4.2).

The default action is PLAN ONLY: no network raster reads and no file writes.
Use --tile X Y for one selected candidate, or --batch N for at most N
pending candidate tiles (N <= 300). Files are staged, NOT final DEMs.

All tiles use the SAME previously checked two-source STAC COG VRT, including
north/south-only tiles. This keeps input sampling consistent across the seam.
Completed tiles are checked before being skipped; no completed TIFF is replaced.
No Old Files Backup inputs, AOI masking, blending, or final mosaic are used.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

import numpy as np
from osgeo import gdal, ogr, osr
from pyproj import Transformer


TILE_METRES = 1000
NODATA = -32767.0
GRID_TRANSFORM = (475000.0, 1000.0, 0.0, 5557000.0, 0.0, -1000.0)
GRID_WIDTH = 33
GRID_HEIGHT = 94
BASE_URL = "https://canelevation-dem.s3.ca-central-1.amazonaws.com/hrdem-mosaic-1m/"
EXPECTED_VRT_SOURCES = {
    "/vsicurl/" + BASE_URL + "2_3-mosaic-1m-dtm.tif",
    "/vsicurl/" + BASE_URL + "2_4-mosaic-1m-dtm.tif",
}


def project_root() -> Path:
    root = Path(__file__).resolve().parent
    expected = root / "data/processed/elevation/hrdem_1m_dtm_coverage_highway99.gpkg"
    if not expected.is_file():
        raise FileNotFoundError(
            "Place this script in the natcat-risk project root; missing " + str(expected)
        )
    return root


def candidate_tiles(footprint: Path) -> list[tuple[int, int, str]]:
    ds = ogr.Open(str(footprint), 0)
    if ds is None:
        raise RuntimeError("Cannot open footprint GeoPackage: " + str(footprint))
    layer = ds.GetLayer(0)
    srs = layer.GetSpatialRef()
    if srs is None:
        raise RuntimeError("Footprint CRS is missing")
    srs.AutoIdentifyEPSG()
    if srs.GetAuthorityCode(None) != "26910":
        raise RuntimeError("Footprint must be EPSG:26910")

    grid = gdal.GetDriverByName("MEM").Create(
        "", GRID_WIDTH, GRID_HEIGHT, 1, gdal.GDT_Byte
    )
    grid.SetGeoTransform(GRID_TRANSFORM)
    grid.SetProjection(srs.ExportToWkt())
    grid.GetRasterBand(1).Fill(0)
    code = gdal.RasterizeLayer(
        grid, [1], layer, burn_values=[1], options=["ALL_TOUCHED=TRUE"]
    )
    if code != 0:
        raise RuntimeError("Rasterizing candidate-tile footprint failed")

    rows, cols = np.where(grid.ReadAsArray() == 1)
    xs = 475000 + cols * TILE_METRES
    ys = 5557000 - (rows + 1) * TILE_METRES
    corner_x = np.stack([xs, xs + 1000, xs, xs + 1000], axis=1)
    corner_y = np.stack([ys, ys, ys + 1000, ys + 1000], axis=1)
    _, lat_y = Transformer.from_crs(26910, 3979, always_xy=True).transform(
        corner_x.ravel(), corner_y.ravel()
    )
    ys3979 = np.asarray(lat_y).reshape(-1, 4)
    only_south = ys3979.max(axis=1) < 499995
    only_north = ys3979.min(axis=1) > 500005
    result = []
    for idx in range(len(xs)):
        group = "2_3" if only_south[idx] else "2_4" if only_north[idx] else "BOTH"
        result.append((int(xs[idx]), int(ys[idx]), group))
    return result


def check_vrt_file(vrt: Path) -> None:
    if not vrt.is_file():
        raise FileNotFoundError("Expected previously tested VRT is missing: " + str(vrt))
    tree = ET.parse(vrt)
    names = {
        (element.text or "").strip()
        for element in tree.iter("SourceFilename")
    }
    if names != EXPECTED_VRT_SOURCES:
        raise RuntimeError(
            "VRT source URLs differ from the two checked official COG sources: "
            + repr(sorted(names))
        )


def inspect_tiff(path: Path, x: int, y: int) -> dict:
    ds = gdal.Open(str(path), gdal.GA_ReadOnly)
    if ds is None:
        raise RuntimeError("Cannot read output TIFF: " + str(path))
    if ds.RasterCount != 1 or (ds.RasterXSize, ds.RasterYSize) != (1000, 1000):
        raise RuntimeError("Wrong number of bands or tile dimensions: " + str(path))
    expected = (float(x), 1.0, 0.0, float(y + 1000), 0.0, -1.0)
    if not np.allclose(ds.GetGeoTransform(), expected, rtol=0.0, atol=1e-8):
        raise RuntimeError("Tile geotransform is not on the 1m grid: " + str(path))
    srs = osr.SpatialReference(wkt=ds.GetProjection())
    srs.AutoIdentifyEPSG()
    if srs.GetAuthorityCode(None) != "26910":
        raise RuntimeError("Output is not EPSG:26910: " + str(path))
    band = ds.GetRasterBand(1)
    if band.DataType != gdal.GDT_Float32 or band.GetNoDataValue() != NODATA:
        raise RuntimeError("Unexpected data type or NoData value: " + str(path))

    values = band.ReadAsArray()
    valid = np.isfinite(values) & (values != NODATA)
    n_valid = int(valid.sum())
    qa = {
        "qa_level": "STRUCTURAL_ONLY_NOT_VERTICAL_ACCURACY",
        "total_pixels": int(valid.size),
        "valid_pixels": n_valid,
        "nodata_pixels": int(valid.size - n_valid),
        "valid_percent": round(100.0 * n_valid / valid.size, 4),
        "min_elevation_m": float(np.min(values[valid])) if n_valid else None,
        "max_elevation_m": float(np.max(values[valid])) if n_valid else None,
        "crs": "EPSG:26910",
        "pixel_size_m": 1,
        "not_yet_done": ["AOI mask", "cross-tile seam QA", "vertical-accuracy QA", "final mosaic"],
    }
    ds = None
    return qa


def write_json_atomically(path: Path, record: dict) -> None:
    temporary = path.with_name(path.name + ".partial")
    if temporary.exists():
        temporary.unlink()  # Script-owned incomplete JSON only.
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(record, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def process_one(root: Path, vrt: Path, x: int, y: int, group: str) -> str:
    folder = root / "results/diagnostics/hrdem1m_batch_staging"
    folder.mkdir(parents=True, exist_ok=True)
    stem = f"hrdem1m_x{x}_y{y}_utm10_1km"
    output = folder / (stem + ".tif")
    sidecar = folder / (stem + ".qa.json")
    partial = folder / (stem + ".partial.tif")

    if output.exists():
        qa = inspect_tiff(output, x, y)
        if not sidecar.exists():
            write_json_atomically(sidecar, {"tile_x": x, "tile_y": y,
                "source_group": group, "source_mode": "SINGLE_TWO_SOURCE_VRT", **qa})
        print("SKIPPED_VALID_EXISTING:", output.name, "valid_percent:", qa["valid_percent"], flush=True)
        return "SKIPPED"

    if partial.exists():
        print("REMOVING_SCRIPT_OWNED_INCOMPLETE_PARTIAL:", partial.name, flush=True)
        partial.unlink()

    command = [
        "gdalwarp", "-of", "GTiff", "-t_srs", "EPSG:26910",
        "-te", str(x), str(y), str(x + 1000), str(y + 1000),
        "-tr", "1", "1", "-tap", "-r", "bilinear", "-et", "0",
        "-ovr", "NONE", "-srcnodata", "-32767", "-dstnodata", "-32767",
        "-co", "TILED=YES", "-co", "COMPRESS=DEFLATE",
        "-co", "PREDICTOR=3", "-co", "BIGTIFF=IF_SAFER",
        str(vrt), str(partial),
    ]
    env = os.environ.copy()
    env.update({
        "GDAL_DISABLE_READDIR_ON_OPEN": "YES",
        "GDAL_HTTP_MAX_RETRY": "3",
        "GDAL_HTTP_RETRY_DELAY": "2",
        "GDAL_HTTP_TIMEOUT": "45",
    })
    print("EXTRACTING:", stem, "group:", group, flush=True)
    try:
        subprocess.run(command, check=True, env=env)
        qa = inspect_tiff(partial, x, y)
        if qa["valid_pixels"] == 0:
            # Keep an auditable empty-tile marker, not a 1-million-pixel empty TIFF.
            write_json_atomically(sidecar, {
                "tile_x": x, "tile_y": y, "source_group": group,
                "status": "EMPTY_RASTER_NO_VALID_PIXELS", **qa,
            })
            print("EMPTY_TILE_NO_VALID_PIXELS:", stem, flush=True)
            return "EMPTY"
        os.replace(partial, output)
        write_json_atomically(sidecar, {
            "tile_x": x, "tile_y": y, "source_group": group,
            "status": "STRUCTURAL_QA_PASSED_WITH_CAVEATS",
            "source_mode": "SINGLE_TWO_SOURCE_VRT", **qa,
        })
        print("STRUCTURAL_QA_PASSED_WITH_CAVEATS:", output.name,
              "valid_percent:", qa["valid_percent"], flush=True)
        return "SUCCESS"
    finally:
        if partial.exists():
            partial.unlink()  # Never delete a completed output.


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--tile", nargs=2, type=int, metavar=("X_MIN", "Y_MIN"),
                      help="Extract exactly one candidate tile")
    mode.add_argument("--batch", type=int, metavar="N",
                      help="Extract at most N pending tiles (1-300); skip completed")
    args = parser.parse_args()
    if args.batch is not None and not 1 <= args.batch <= 300:
        parser.error("--batch must be between 1 and 300; run again to resume remaining tiles")

    gdal.UseExceptions()
    ogr.UseExceptions()
    root = project_root()
    footprint = root / "data/processed/elevation/hrdem_1m_dtm_coverage_highway99.gpkg"
    plan = candidate_tiles(footprint)
    print("GRID_SIZE:", GRID_HEIGHT, "rows x", GRID_WIDTH)
    print("CANDIDATE_TILES:", len(plan))
    for group in ("2_3", "2_4", "BOTH"):
        print("TILES_" + group + ":", sum(item[2] == group for item in plan))
    print("BOUNDARY_TILE_X_Y:", [(x, y) for x, y, k in plan if k == "BOTH"])
    print("ALL_TILE_INPUTS: SINGLE_TWO_SOURCE_VRT")
    print("OUTPUT_TYPE: DIAGNOSTIC_STAGING_UNMASKED_NOT_FINAL")

    if args.tile is None and args.batch is None:
        print("PLAN_ONLY: no raster files downloaded or changed")
        print("NEXT_SAFE_TEST: python step4_2_hrdem1m_batch_vrt_300.py --batch 300")
        return

    vrt = root / "results/diagnostics/step4_2_hrdem1m_boundary_virtual.vrt"
    check_vrt_file(vrt)
    lookup = {(x, y): group for x, y, group in plan}
    if args.tile is not None:
        key = (args.tile[0], args.tile[1])
        if key not in lookup:
            raise ValueError("Requested tile not in candidate grid: " + repr(key))
        todo = [(key[0], key[1], lookup[key])]
    else:
        output_dir = root / "results/diagnostics/hrdem1m_batch_staging"
        todo = []
        for x, y, group in plan:
            stem = f"hrdem1m_x{x}_y{y}_utm10_1km"
            output = output_dir / (stem + ".tif")
            qa_json = output_dir / (stem + ".qa.json")
            if output.exists():
                # Re-check completed data rather than trusting filename alone.
                inspect_tiff(output, x, y)
                if not qa_json.exists():
                    raise RuntimeError("Existing TIFF is missing QA JSON; use --tile to audit: " + str(output))
                continue
            if qa_json.exists():
                with qa_json.open(encoding="utf-8") as handle:
                    record = json.load(handle)
                if record.get("status") == "EMPTY_RASTER_NO_VALID_PIXELS":
                    continue
                raise RuntimeError("QA JSON exists without TIFF: " + str(qa_json))
            todo.append((x, y, group))
            if len(todo) == args.batch:
                break
        print("PENDING_SELECTED_FOR_BATCH:", len(todo))
    counts = {"SUCCESS": 0, "SKIPPED": 0, "EMPTY": 0}
    for i, (x, y, group) in enumerate(todo, 1):
        print(f"TILE_PROGRESS: {i}/{len(todo)}", flush=True)
        status = process_one(root, vrt, x, y, group)
        counts[status] += 1
    print("BATCH_FINISHED:", counts)
    print("NOTE: These are staging rasters; AOI masking and cross-tile QA are pending")


if __name__ == "__main__":
    try:
        main()
    except (Exception, KeyboardInterrupt) as error:
        print("ERROR:", error, file=sys.stderr)
        sys.exit(1)
