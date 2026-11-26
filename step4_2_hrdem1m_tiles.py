"""Prepare and pilot-test official HRDEM 1 m DTM tiles for Highway 99.

The default action is read-only: it prints the 1 km tile plan.
One tile can be extracted explicitly with --tile X_MIN Y_MIN.
No bulk extraction or final DEM merging is performed by this script.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from osgeo import gdal, ogr, osr
from pyproj import Transformer

GEO_TRANSFORM = (475000, 1000, 0, 5557000, 0, -1000)
GRID_WIDTH = 33
GRID_HEIGHT = 94
NO_DATA = -32767.0
BASE_URL = "https://canelevation-dem.s3.ca-central-1.amazonaws.com/hrdem-mosaic-1m/"
SOURCES = {
    "2_3": "/vsicurl/" + BASE_URL + "2_3-mosaic-1m-dtm.tif",
    "2_4": "/vsicurl/" + BASE_URL + "2_4-mosaic-1m-dtm.tif",
}


def find_project_root() -> Path:
    for folder in Path(__file__).resolve().parents:
        if (folder / "data/processed/elevation/hrdem_1m_dtm_coverage_highway99.gpkg").is_file():
            return folder
    raise FileNotFoundError(
        "Project root not found. Put this script inside the natcat-risk project folder."
    )


def create_tile_plan(coverage_path: Path) -> list[tuple[int, int, str]]:
    source = ogr.Open(str(coverage_path), 0)
    if source is None:
        raise RuntimeError(f"Cannot open the coverage file: {coverage_path}")
    layer = source.GetLayer(0)
    spatial_ref = layer.GetSpatialRef()
    if spatial_ref is None:
        raise RuntimeError("Coverage GeoPackage has no CRS.")
    spatial_ref.AutoIdentifyEPSG()
    if spatial_ref.GetAuthorityCode(None) != "26910":
        raise RuntimeError("Expected EPSG:26910 coverage; refusing to assume coordinates.")

    mask = gdal.GetDriverByName("MEM").Create(
        "", GRID_WIDTH, GRID_HEIGHT, 1, gdal.GDT_Byte
    )
    mask.SetGeoTransform(GEO_TRANSFORM)
    mask.SetProjection(spatial_ref.ExportToWkt())
    mask.GetRasterBand(1).Fill(0)
    result = gdal.RasterizeLayer(
        mask, [1], layer, burn_values=[1], options=["ALL_TOUCHED=TRUE"]
    )
    if result != 0:
        raise RuntimeError("Failed to rasterize the official coverage footprint.")

    rows, cols = np.where(mask.ReadAsArray() == 1)
    x_mins = 475000 + cols * 1000
    y_mins = 5557000 - (rows + 1) * 1000
    transformer = Transformer.from_crs("EPSG:26910", "EPSG:3979", always_xy=True)
    # Classify using all four tile corners and a conservative 5 m seam margin.
    corner_x = np.stack([x_mins, x_mins + 1000, x_mins, x_mins + 1000], axis=1)
    corner_y = np.stack([y_mins, y_mins, y_mins + 1000, y_mins + 1000], axis=1)
    _, transformed_y = transformer.transform(corner_x.ravel(), corner_y.ravel())
    y3979 = np.asarray(transformed_y).reshape(-1, 4)
    south = y3979.max(axis=1) < 499995
    north = y3979.min(axis=1) > 500005

    plan = []
    for i in range(len(rows)):
        source_name = "2_3" if south[i] else ("2_4" if north[i] else "BOTH")
        plan.append((int(x_mins[i]), int(y_mins[i]), source_name))
    return plan


def print_plan(plan: list[tuple[int, int, str]]) -> None:
    print("GRID_SIZE:", GRID_HEIGHT, "rows x", GRID_WIDTH, "columns")
    print("CANDIDATE_TILES:", len(plan))
    for name in ("2_3", "2_4", "BOTH"):
        print(f"TILES_{name}:", sum(t[2] == name for t in plan))
    print("BOUNDARY_TILE_X_Y:", [(x, y) for x, y, name in plan if name == "BOTH"])
    print("FIRST_5_TILES:", plan[:5])
    print("PLAN_ONLY: no raster files downloaded or changed")


def extract_pilot(root: Path, x_min: int, y_min: int, source_name: str) -> None:
    output_dir = root / "results" / "diagnostics" / "hrdem1m_tile_pilots"
    output = output_dir / f"hrdem1m_x{x_min}_y{y_min}_utm10_1km.tif"
    temporary = output.with_name(output.stem + ".partial.tif")
    if output.exists():
        print("SKIPPED_EXISTING:", output)
        print("Existing files are not overwritten; validate them separately before reuse.")
        return
    if temporary.exists():
        raise FileExistsError(f"Incomplete temporary file exists; inspect first: {temporary}")
    output_dir.mkdir(parents=True, exist_ok=True)

    raster_sources = (
        [SOURCES["2_3"], SOURCES["2_4"]] if source_name == "BOTH" else [SOURCES[source_name]]
    )
    command = [
        "gdalwarp", "-of", "GTiff", "-t_srs", "EPSG:26910",
        "-te", str(x_min), str(y_min), str(x_min + 1000), str(y_min + 1000),
        "-tr", "1", "1", "-tap", "-r", "bilinear", "-et", "0",
        "-ovr", "NONE", "-srcnodata", str(int(NO_DATA)),
        "-dstnodata", str(int(NO_DATA)),
        "-co", "TILED=YES", "-co", "COMPRESS=DEFLATE",
        "-co", "PREDICTOR=3", "-co", "BIGTIFF=IF_SAFER",
        *raster_sources, str(temporary),
    ]
    print("PILOT_TILE_X_Y:", x_min, y_min, flush=True)
    print("SOURCE_GROUP:", source_name, flush=True)
    print("OFFICIAL_SOURCE_COUNT:", len(raster_sources), flush=True)
    environment = os.environ.copy()
    environment["GDAL_DISABLE_READDIR_ON_OPEN"] = "YES"
    environment["GDAL_HTTP_MAX_RETRY"] = "3"
    environment["GDAL_HTTP_RETRY_DELAY"] = "2"
    subprocess.run(command, check=True, env=environment)

    dataset = gdal.Open(str(temporary))
    if dataset is None:
        raise RuntimeError(f"GDAL output could not be reopened: {temporary}")
    gt = dataset.GetGeoTransform()
    srs = osr.SpatialReference(wkt=dataset.GetProjection())
    srs.AutoIdentifyEPSG()
    expected_gt = (x_min, 1, 0, y_min + 1000, 0, -1)
    if (dataset.RasterXSize, dataset.RasterYSize) != (1000, 1000):
        raise RuntimeError("Pilot QA failed: unexpected dimensions.")
    if not np.allclose(gt, expected_gt, rtol=0, atol=1e-7):
        raise RuntimeError("Pilot QA failed: incorrect tile grid/geotransform.")
    if srs.GetAuthorityCode(None) != "26910":
        raise RuntimeError("Pilot QA failed: incorrect target CRS.")
    band = dataset.GetRasterBand(1)
    if band.GetNoDataValue() != NO_DATA:
        raise RuntimeError("Pilot QA failed: incorrect NoData value.")
    values = band.ReadAsArray()
    valid = np.isfinite(values) & (values != NO_DATA)
    valid_count = int(valid.sum())
    print("PILOT_TOTAL_PIXELS:", valid.size)
    print("PILOT_VALID_PIXELS:", valid_count)
    print("PILOT_VALID_PERCENT:", round(100 * valid_count / valid.size, 4))
    if valid_count:
        print("PILOT_ELEVATION_MIN_MAX:", float(values[valid].min()), float(values[valid].max()))
    dataset = None
    if valid_count == 0:
        temporary.unlink()
        print("NO_VALID_PIXELS: empty candidate tile not retained")
        return
    os.replace(temporary, output)
    print("PILOT_STRUCTURAL_QA_PASSED_WITH_CAVEATS")
    print("PILOT_OUTPUT:", output)
    print("NOTE: footprint/AOI masking, seam quality, and vertical accuracy are not yet verified.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tile", nargs=2, type=int, metavar=("X_MIN", "Y_MIN"),
        help="Explicitly extract ONE 1 km diagnostic tile; default is plan-only.",
    )
    args = parser.parse_args()
    gdal.UseExceptions()
    ogr.UseExceptions()
    root = find_project_root()
    coverage = root / "data/processed/elevation/hrdem_1m_dtm_coverage_highway99.gpkg"
    plan = create_tile_plan(coverage)
    print_plan(plan)
    if args.tile is None:
        return
    key = tuple(args.tile)
    candidates = {(x, y): name for x, y, name in plan}
    if key not in candidates:
        raise ValueError(f"Tile {key} is not in the coverage candidate set; no download started.")
    extract_pilot(root, key[0], key[1], candidates[key])


if __name__ == "__main__":
    try:
        main()
    except (Exception, KeyboardInterrupt) as error:
        print("ERROR:", error, file=sys.stderr)
        sys.exit(1)
