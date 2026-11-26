"""Step 4.2: audit archived LiDAR GeoTIFF *footprints*, excluding derived masks.

Metadata-only investigation. Does NOT read elevation pixels, reproject rasters,
or overwrite the originals. Run from the natcat-risk repository environment:

    python src/data_sources/10_recheck_step42_archived_lidar_footprints.py

The source file naming convention 'bc_*.tif' identifies *candidate* tiles,
not confirmed DTM products or valid pixel footprints. This audit cannot prove
coverage or elevation accuracy at raster-pixel level.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import geopandas as gpd
from osgeo import gdal
from pyproj import CRS, Transformer
from shapely.geometry import Polygon
from shapely.ops import unary_union

REPO = Path(__file__).resolve().parents[2]
ARCHIVE = (
    REPO.parent
    / "old files backup"
    / "Climada"
    / "Climada_Project"
    / "Inputs"
    / "DEM"
    / "LiDAR DEMs"
)
AOI_PATH = REPO / "data/processed/study_area/highway99_terrain_search_5000m.geojson"
REPORT = REPO / "results/diagnostics/step4_2_archived_lidar_footprints_without_mask_qa.json"
TARGET_CRS = CRS.from_epsg(26910)
BORDER_SAMPLES = 21  # denser edges for transformed raster extents


def projected_raster_rectangle(ds: gdal.Dataset) -> tuple[Polygon, str]:
    wkt = ds.GetProjectionRef()
    if not wkt:
        raise ValueError("Raster has no declared CRS")
    source_crs = CRS.from_wkt(wkt)
    transformer = Transformer.from_crs(source_crs, TARGET_CRS, always_xy=True)
    gt = ds.GetGeoTransform()
    if gt is None:
        raise ValueError("Raster has no geotransform")
    nx, ny = ds.RasterXSize, ds.RasterYSize
    if nx <= 0 or ny <= 0:
        raise ValueError("Raster is empty")

    # Densify all four edges, since a projected rectangle's curved edges
    # need not remain straight in the destination CRS.
    coords_pix = []
    edges = [((0, 0), (nx, 0)), ((nx, 0), (nx, ny)),
             ((nx, ny), (0, ny)), ((0, ny), (0, 0))]
    for (x0, y0), (x1, y1) in edges:
        for j in range(BORDER_SAMPLES - 1):
            t = j / (BORDER_SAMPLES - 1)
            coords_pix.append((x0 + t * (x1 - x0), y0 + t * (y1 - y0)))
    xs = [gt[0] + px * gt[1] + py * gt[2] for px, py in coords_pix]
    ys = [gt[3] + px * gt[4] + py * gt[5] for px, py in coords_pix]
    xx, yy = transformer.transform(xs, ys)
    if not all(math.isfinite(float(v)) for v in (*xx, *yy)):
        raise ValueError("Raster footprint transforms to non-finite coordinates")
    footprint = Polygon(zip(xx, yy))
    if not footprint.is_valid:
        footprint = footprint.buffer(0)
    if footprint.is_empty or not footprint.is_valid:
        raise ValueError("Invalid transformed raster footprint")
    auth = source_crs.to_epsg()
    return footprint, str(auth) if auth else source_crs.to_string()


def area_km2(shape) -> float:
    return round(float(shape.area) / 1_000_000, 6)


def main() -> int:
    gdal.UseExceptions()
    if not ARCHIVE.is_dir() or not AOI_PATH.is_file():
        raise FileNotFoundError(f"Check paths:\nARCHIVE={ARCHIVE}\nAOI={AOI_PATH}")

    aoi_gdf = gpd.read_file(AOI_PATH)
    if aoi_gdf.crs is None or aoi_gdf.empty:
        raise ValueError("Study AOI missing declared CRS or polygons")
    aoi = unary_union(aoi_gdf.to_crs(epsg=26910).geometry.tolist())
    if aoi.is_empty or not aoi.is_valid:
        raise ValueError("Study AOI is empty or invalid")

    paths = sorted(ARCHIVE.glob("*.tif"))  # do NOT traverse processed subfolders
    records = []
    intersections_by_name = {}
    errors = []
    for path in paths:
        try:
            ds = gdal.Open(str(path), gdal.GA_ReadOnly)
            if ds is None or ds.RasterCount < 1:
                raise ValueError("Could not open first raster band")
            footprint, crs = projected_raster_rectangle(ds)
            overlap = footprint.intersection(aoi)
            band = ds.GetRasterBand(1)
            name = path.name
            role = "source_tile_candidate" if name.lower().startswith("bc_") else "non_source_or_derived_candidate"
            record = {
                "filename": name,
                "role_by_filename_only": role,
                "crs": crs,
                "width": ds.RasterXSize,
                "height": ds.RasterYSize,
                "nodata": band.GetNoDataValue(),
                "band_type": gdal.GetDataTypeName(band.DataType),
                "bounding_rectangle_intersects_aoi": not overlap.is_empty and overlap.area > 0,
                "bounding_rectangle_aoi_intersection_km2": area_km2(overlap),
            }
            records.append(record)
            if not overlap.is_empty and overlap.area > 0:
                intersections_by_name[name] = overlap
            band = None
            ds = None
        except Exception as exc:
            errors.append({"file": path.name, "error": str(exc)})

    original_tile_names = [r["filename"] for r in records
                           if r["role_by_filename_only"] == "source_tile_candidate"]
    original_parts = [intersections_by_name[n] for n in original_tile_names
                      if n in intersections_by_name]
    all_parts = list(intersections_by_name.values())
    originals_union = unary_union(original_parts) if original_parts else Polygon()
    all_union = unary_union(all_parts) if all_parts else Polygon()
    mask_name = "lidar_mask01.tif"
    mask_record = next((r for r in records if r["filename"].lower() == mask_name), None)
    aoi_area = aoi.area
    original_area = originals_union.area
    all_area = all_union.area
    pct = lambda a: round(100 * float(a) / aoi_area, 6)

    result = {
        "status": "FOOTPRINT_SCREENING_ONLY_NOT_PIXEL_VALIDITY_QA" if not errors else "METADATA_READ_ERRORS",
        "study_area_km2": area_km2(aoi),
        "total_top_level_tiffs": len(paths),
        "metadata_read_success": len(records),
        "metadata_read_errors": errors,
        "bc_named_tile_candidates": len(original_tile_names),
        "bc_tiles_with_aoi_footprint_overlap": len(original_parts),
        "other_named_files": [r["filename"] for r in records
                              if r["role_by_filename_only"] != "source_tile_candidate"],
        "all_raster_rectangles_overlap_count": len(all_parts),
        "union_all_rectangles_aoi_km2": area_km2(all_union),
        "union_all_rectangles_aoi_pct": pct(all_area),
        "union_bc_tiles_only_aoi_km2": area_km2(originals_union),
        "union_bc_tiles_only_aoi_pct": pct(original_area),
        "area_only_covered_by_non_bc_named_footprints_km2": area_km2(all_union.difference(originals_union)),
        "not_covered_by_bc_tile_rectangles_km2": area_km2(aoi.difference(originals_union)),
        "lidar_mask01_info": mask_record,
        "source_tile_footprints_in_aoi": [r for r in records if r["role_by_filename_only"] == "source_tile_candidate" and r["bounding_rectangle_intersects_aoi"]],
        "caveats": [
            "Every rectangle uses the TIFF outer extent, not the footprint of valid elevation pixels.",
            "The bc_ filename prefix does not verify DTM type, product origin, or quality.",
            "Never interpret these percentages as actual LiDAR elevation coverage.",
            "Source raster contents, NoData distribution, and tile quality require separate QA.",
        ],
    }

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")

    print("STEP 4.2 — ARCHIVED LIDAR FOOTPRINT CHECK WITHOUT DERIVED MASKS")
    print("Status:", result["status"])
    for field in (
        "total_top_level_tiffs", "metadata_read_success", "bc_named_tile_candidates",
        "bc_tiles_with_aoi_footprint_overlap", "other_named_files", "all_raster_rectangles_overlap_count",
        "study_area_km2", "union_all_rectangles_aoi_km2", "union_all_rectangles_aoi_pct",
        "union_bc_tiles_only_aoi_km2", "union_bc_tiles_only_aoi_pct",
        "area_only_covered_by_non_bc_named_footprints_km2",
        "not_covered_by_bc_tile_rectangles_km2",
    ):
        print(f"  {field}: {result[field]}")
    if mask_record:
        print("  lidar_mask01_crs:", mask_record["crs"])
        print("  lidar_mask01_nodata:", mask_record["nodata"])
        print("  lidar_mask01_aoi_intersection_km2:", mask_record["bounding_rectangle_aoi_intersection_km2"])
    print("Report:", REPORT)
    print("CAVEAT: Rectangular footprints are NOT valid DEM-pixel coverage.")
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
