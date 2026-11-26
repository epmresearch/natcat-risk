#!/usr/bin/env python3
"""Screen archived LiDAR TIFF footprints against the Highway 99 AOI.

Metadata and geometric footprints only; no DEM pixel quality validation.
Archived TIFFs are opened read-only and never modified.
"""
from __future__ import annotations

import csv
import json
import math
from collections import Counter
from pathlib import Path

import geopandas as gpd
from osgeo import gdal
from pyproj import CRS, Transformer
from shapely.geometry import Polygon
from shapely.ops import unary_union

ROOT = Path(__file__).resolve().parents[2]
LIDAR = (ROOT.parent / "old files backup" / "Climada" / "Climada_Project"
         / "Inputs" / "DEM" / "LiDAR DEMs")
AOI = ROOT / "data/processed/study_area/highway99_terrain_search_5000m.geojson"
OUT = ROOT / "results/diagnostics"
EPSG = 26910
FIELDS = ["filename", "size_mib", "epsg", "width", "height", "pixel_x",
          "pixel_y", "dtype", "nodata", "unit", "footprint_km2",
          "overlap_aoi_km2", "intersects_aoi", "status", "error"]


def footprint(ds, transform):
    """Transform a densified raster boundary, not just its bounding box."""
    w, h = ds.RasterXSize, ds.RasterYSize
    gt = ds.GetGeoTransform()
    corners = [(0, 0), (w, 0), (w, h), (0, h), (0, 0)]
    points = []
    for (c0, r0), (c1, r1) in zip(corners, corners[1:]):
        for j in range(16):
            t = j / 16
            c, r = c0 + t * (c1 - c0), r0 + t * (r1 - r0)
            x, y = gt[0] + c * gt[1] + r * gt[2], gt[3] + c * gt[4] + r * gt[5]
            xx, yy = transform.transform(x, y)
            if not all(math.isfinite(v) for v in (xx, yy)):
                raise ValueError("Non-finite transformed coordinate")
            points.append((xx, yy))
    poly = Polygon(points)
    if not poly.is_valid or poly.is_empty or poly.area <= 0:
        raise ValueError("Invalid projected raster footprint")
    return poly


def main():
    if not LIDAR.is_dir():
        raise SystemExit(f"Archive not found: {LIDAR}")
    if not AOI.is_file():
        raise SystemExit(f"AOI not found: {AOI}")
    files = sorted([p for p in LIDAR.iterdir() if p.is_file()
                    and p.suffix.lower() in {".tif", ".tiff"}])
    if not files:
        raise SystemExit("No TIFFs in main LiDAR archive folder")
    aoi_gdf = gpd.read_file(AOI)
    if aoi_gdf.empty or aoi_gdf.crs is None:
        raise SystemExit("AOI is empty or missing its CRS")
    aoi = unary_union(list(aoi_gdf.to_crs(epsg=EPSG).geometry))
    if not aoi.is_valid or aoi.is_empty:
        raise SystemExit("Invalid AOI geometry")
    target = CRS.from_epsg(EPSG)
    gdal.UseExceptions()
    rows, intersections = [], []
    for index, path in enumerate(files, 1):
        row = dict.fromkeys(FIELDS, "")
        row["filename"] = path.name
        row["size_mib"] = round(path.stat().st_size / (1024 ** 2), 3)
        ds = None
        try:
            ds = gdal.OpenEx(str(path), gdal.OF_RASTER | gdal.OF_READONLY)
            if ds is None or ds.RasterCount == 0:
                raise ValueError("Raster could not be opened")
            wkt = ds.GetProjectionRef()
            if not wkt:
                raise ValueError("Missing raster CRS")
            source = CRS.from_wkt(wkt)
            gt = ds.GetGeoTransform()
            band = ds.GetRasterBand(1)
            row.update(epsg=source.to_epsg() or "unknown",
                       width=ds.RasterXSize, height=ds.RasterYSize,
                       pixel_x=round(math.hypot(gt[1], gt[4]), 6),
                       pixel_y=round(math.hypot(gt[2], gt[5]), 6),
                       dtype=gdal.GetDataTypeName(band.DataType),
                       nodata=band.GetNoDataValue(),
                       unit=band.GetUnitType() or "not_declared")
            poly = footprint(ds, Transformer.from_crs(source, target, always_xy=True))
            row["footprint_km2"] = round(poly.area / 1e6, 6)
            ax0, ay0, ax1, ay1 = aoi.bounds
            x0, y0, x1, y1 = poly.bounds
            overlap = (poly.intersection(aoi) if x1 >= ax0 and x0 <= ax1
                       and y1 >= ay0 and y0 <= ay1 else None)
            area = overlap.area if overlap is not None else 0.0
            row["overlap_aoi_km2"] = round(area / 1e6, 6)
            row["intersects_aoi"] = area > 0
            row["status"] = "CANDIDATE_FOOTPRINT" if area > 0 else "OUTSIDE_AOI"
            if area > 0:
                intersections.append(overlap)
        except Exception as exc:
            row["status"] = "ERROR"
            row["error"] = f"{type(exc).__name__}: {exc}"
        finally:
            ds = None
        rows.append(row)
        if index % 25 == 0 or index == len(files):
            print(f"Scanned {index}/{len(files)} files")

    candidates = [r for r in rows if r["status"] == "CANDIDATE_FOOTPRINT"]
    errors = [r for r in rows if r["status"] == "ERROR"]
    union_area = unary_union(intersections).area if intersections else 0.0
    report = {
        "status": "INVENTORY_COMPLETE_WITH_CAVEATS" if not errors else "INVENTORY_ERRORS_NEED_REVIEW",
        "archive_folder": str(LIDAR), "aoi_path": str(AOI),
        "scope": "top-level TIFFs only; derived subfolders excluded",
        "total_tiffs": len(files), "read_success": len(files) - len(errors),
        "read_errors": len(errors), "candidate_footprints": len(candidates),
        "outside_aoi": sum(r["status"] == "OUTSIDE_AOI" for r in rows),
        "aoi_km2": round(aoi.area / 1e6, 6),
        "candidate_union_aoi_km2": round(union_area / 1e6, 6),
        "candidate_union_aoi_pct": round(100 * union_area / aoi.area, 6),
        "crs_counts": dict(Counter(str(r["epsg"]) for r in rows if r["epsg"] != "")),
        "pixel_size_counts": dict(Counter(f"{r['pixel_x']}x{r['pixel_y']}" for r in rows if r["pixel_x"] != "")),
        "candidate_filenames": [r["filename"] for r in candidates],
        "errors": [{"filename": r["filename"], "error": r["error"]} for r in errors],
        "caveat": "Footprint overlap is not proof of valid DTM pixels, provenance or height accuracy."
    }
    OUT.mkdir(parents=True, exist_ok=True)
    csv_file = OUT / "step4_2_local_lidar_tile_inventory.csv"
    json_file = OUT / "step4_2_local_lidar_metadata_inventory.json"
    with csv_file.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    json_file.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print("\nSTEP 4.2 — ARCHIVED LIDAR METADATA / FOOTPRINT SCREENING")
    for key in ("status", "total_tiffs", "read_success", "read_errors",
                "candidate_footprints", "outside_aoi", "aoi_km2",
                "candidate_union_aoi_km2", "candidate_union_aoi_pct",
                "crs_counts", "pixel_size_counts"):
        print(f"  {key}: {report[key]}")
    print("CANDIDATE FILES:")
    for file in report["candidate_filenames"]:
        print(" ", file)
    if errors:
        print("ERROR FILES:", report["errors"])
    print("CSV:", csv_file)
    print("JSON:", json_file)
    print("CAUTION: This is a footprint inventory, not raster pixel validity QA.")


if __name__ == "__main__":
    main()
