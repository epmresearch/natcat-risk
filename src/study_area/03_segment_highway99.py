"""Step 3.3: segment the provisional Highway 99 centerline by metric chainage.

No network calls. Geometry is cut along the source LineString in EPSG:26910.
Input remains provisional; segments are not landslide-event footprints or routes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from shapely.geometry import LineString, Point
from shapely.ops import substring

ROOT = Path(__file__).resolve().parents[2]
METRIC_CRS = "EPSG:26910"  # NAD83 / UTM zone 10N: metres
GEOGRAPHIC_CRS = "EPSG:4326"
DEFAULT_LENGTH_M = 1000
SOURCE = ROOT / "data/processed/study_area/highway99_candidate_centerline.geojson"
SOURCE_QUALITY = ROOT / "data/processed/study_area/highway99_candidate_quality.json"
OUT_DATA = ROOT / "data/processed/study_area"
OUT_MAPS = ROOT / "results/maps"


def split_line_by_chainage(line: LineString, segment_length_m: int) -> list[dict]:
    """Cut one *projected* ordered LineString without simplifying its shape."""
    if not isinstance(line, LineString):
        raise ValueError("Expected a single ordered LineString, not a MultiLineString.")
    if line.is_empty or not line.is_valid or not math.isfinite(line.length) or line.length <= 0:
        raise ValueError("Source route must have a valid, nonempty, positive-length geometry.")
    if isinstance(segment_length_m, bool) or not isinstance(segment_length_m, int) or segment_length_m <= 0:
        raise ValueError("segment_length_m must be a positive integer number of metres.")

    total = line.length
    n = max(1, math.ceil((total - 1e-7) / segment_length_m))
    segments = []
    for i in range(n):
        start = i * segment_length_m
        end = min((i + 1) * segment_length_m, total)
        geom = substring(line, start, end)
        if not isinstance(geom, LineString) or geom.is_empty or geom.length <= 0:
            raise ValueError(f"Invalid segment geometry at index {i}.")
        if abs(geom.length - (end - start)) > 0.05:
            raise ValueError(f"Segment {i} metric length mismatch >5 cm.")
        segments.append({
            "segment_id": f"H99-{i+1:04d}",
            "segment_no": i + 1,
            "start_m": round(start, 3),
            "end_m": round(end, 3),
            "start_km": round(start / 1000, 6),
            "end_km": round(end / 1000, 6),
            "length_m": round(geom.length, 3),
            "geometry": geom,
        })

    # Explicit QA: no gaps, duplicates, loss of length, or shifted endpoints.
    if Point(segments[0]["geometry"].coords[0]).distance(Point(line.coords[0])) > 0.01:
        raise ValueError("First segment does not start at the source origin.")
    if Point(segments[-1]["geometry"].coords[-1]).distance(Point(line.coords[-1])) > 0.01:
        raise ValueError("Last segment does not end at the source destination.")
    gaps = [Point(a["geometry"].coords[-1]).distance(Point(b["geometry"].coords[0]))
            for a, b in zip(segments[:-1], segments[1:])]
    if gaps and max(gaps) > 0.01:
        raise ValueError(f"Segment boundary gap: {max(gaps):.4f} m.")
    if abs(sum(x["geometry"].length for x in segments) - total) > 0.1:
        raise ValueError("Segment total lengths do not preserve the original route length.")
    return segments


def read_centerline(source: Path) -> tuple[gpd.GeoDataFrame, LineString]:
    if not source.is_file():
        raise FileNotFoundError(f"Required Step 3.2 input not found: {source}")
    gdf = gpd.read_file(source)
    if len(gdf) != 1 or gdf.crs is None:
        raise ValueError("Expected exactly one centerline feature with a declared CRS.")
    if gdf.geometry.iloc[0] is None or gdf.geometry.iloc[0].geom_type != "LineString":
        raise ValueError("Step 3.2 centerline must contain one ordered LineString.")
    metric = gdf.to_crs(METRIC_CRS)
    geom = metric.geometry.iloc[0]
    if geom.is_empty or not geom.is_valid or geom.length <= 0:
        raise ValueError("Invalid or empty metric centerline.")
    # Step 3.2 encoded order through Horseshoe Bay -> Whistler.
    line_wgs84 = gdf.to_crs(GEOGRAPHIC_CRS).geometry.iloc[0]
    if line_wgs84.coords[0][1] >= line_wgs84.coords[-1][1]:
        raise ValueError("Unexpected route order: start is not south of end. Review source geometry.")
    return gdf, geom


def read_source_quality(path: Path, measured_length_m: float) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f"Required Step 3.2 QA metadata not found: {path}")
    qa = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(qa, dict) or qa.get("step") != "3.2":
        raise ValueError("Incorrect or unreadable Step 3.2 quality report.")
    expected_km = qa.get("candidate_length_km")
    if expected_km is None:
        raise ValueError("Step 3.2 QA metadata is missing candidate_length_km.")
    delta_m = abs(float(expected_km) * 1000 - measured_length_m)
    if delta_m > 100:
        raise ValueError(f"Centerline differs from Step 3.2 QA by {delta_m:.1f} m (>100 m).")
    return qa


def make_map(segments_metric: gpd.GeoDataFrame, out_path: Path, segment_length_m: int) -> None:
    """Draw an overview with an index colour ramp; geometry stays in metre CRS."""
    fig, ax = plt.subplots(figsize=(7.2, 11.5))
    polylines = [list(g.coords) for g in segments_metric.geometry]
    ids = segments_metric["segment_no"].to_numpy()
    lc = LineCollection(polylines, cmap="viridis", norm=Normalize(vmin=ids.min(), vmax=ids.max()),
                        linewidths=3.0, zorder=2)
    lc.set_array(ids)
    ax.add_collection(lc)
    start = Point(segments_metric.geometry.iloc[0].coords[0])
    end = Point(segments_metric.geometry.iloc[-1].coords[-1])
    ax.scatter([start.x, end.x], [start.y, end.y], s=36, color="black", zorder=4)
    ax.annotate("Horseshoe Bay", (start.x, start.y), xytext=(7, -13),
                textcoords="offset points", fontsize=9)
    ax.annotate("Whistler", (end.x, end.y), xytext=(7, 5),
                textcoords="offset points", fontsize=9)
    bounds = segments_metric.total_bounds
    dx, dy = bounds[2] - bounds[0], bounds[3] - bounds[1]
    padx, pady = max(3000, dx * 0.13), max(3000, dy * 0.045)
    ax.set_xlim(bounds[0] - padx, bounds[2] + padx)
    ax.set_ylim(bounds[1] - pady, bounds[3] + pady)
    ax.set_aspect("equal")
    ax.grid(alpha=0.18)
    ax.set_xlabel("Easting (m) — NAD83 / UTM zone 10N")
    ax.set_ylabel("Northing (m) — NAD83 / UTM zone 10N")
    ax.set_title(f"Highway 99 — {segment_length_m:,} m chainage segments\n"
                 "Horseshoe Bay to Whistler | Step 3.3", fontsize=12)
    cb = fig.colorbar(lc, ax=ax, shrink=0.64, pad=0.06)
    cb.set_label("Segment sequence (south to north)")
    fig.text(0.50, 0.025,
             "© OpenStreetMap contributors | Provisional OSM road geometry.\n"
             "Segment boundaries are analytical units, not landslide/closure boundaries.",
             ha="center", fontsize=8)
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def build_outputs(source: Path, quality: Path, out_data: Path,
                  out_maps: Path, segment_length_m: int = DEFAULT_LENGTH_M,
                  make_figure: bool = True) -> dict:
    source_gdf, route = read_centerline(source)
    source_qa = read_source_quality(quality, route.length)
    parts = split_line_by_chainage(route, segment_length_m)
    metric = gpd.GeoDataFrame(parts, crs=METRIC_CRS)
    mids = [g.interpolate(g.length / 2) for g in metric.geometry]
    mid_geo = gpd.GeoSeries(mids, crs=METRIC_CRS).to_crs(GEOGRAPHIC_CRS)
    metric["mid_lon"] = [round(float(p.x), 8) for p in mid_geo]
    metric["mid_lat"] = [round(float(p.y), 8) for p in mid_geo]
    metric["source_status"] = str(source_gdf.iloc[0].get("status", "PROVISIONAL"))

    name = f"highway99_segments_{segment_length_m}m"
    out_data.mkdir(parents=True, exist_ok=True)
    geojson_path = out_data / f"{name}.geojson"
    csv_path = out_data / f"{name}.csv"
    qa_path = out_data / f"{name}_quality.json"
    map_path = out_maps / f"{name}.png"

    metric.to_crs(GEOGRAPHIC_CRS).to_file(geojson_path, driver="GeoJSON")
    metric.drop(columns="geometry").to_csv(csv_path, index=False, float_format="%.8f")

    # Re-read the actual GeoJSON output and independently check its CRS and lengths.
    written = gpd.read_file(geojson_path)
    if written.crs is None or len(written) != len(parts):
        raise ValueError("Exported GeoJSON has incorrect feature count or missing CRS.")
    written_metric = written.to_crs(METRIC_CRS)
    max_gap = max((Point(a.coords[-1]).distance(Point(b.coords[0]))
                   for a, b in zip(written_metric.geometry[:-1], written_metric.geometry[1:])),
                  default=0.0)
    output_total_m = float(written_metric.length.sum())
    difference_m = abs(output_total_m - route.length)
    if max_gap > 0.02 or difference_m > 0.25:
        raise ValueError(f"GeoJSON round-trip QA failed: gap={max_gap}, length error={difference_m}.")
    if written["segment_id"].duplicated().any():
        raise ValueError("Duplicate segment identifiers in GeoJSON.")

    lengths = metric.geometry.length
    metadata = {
        "step": "3.3",
        "status": "PROVISIONAL_SEGMENTATION_REQUIRES_REVIEW",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_file": "data/processed/study_area/highway99_candidate_centerline.geojson",
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "source_quality_file": "data/processed/study_area/highway99_candidate_quality.json",
        "source_quality_status": source_qa.get("status"),
        "original_data_source": "OpenStreetMap via OSMnx/Overpass; © OpenStreetMap contributors",
        "route_direction": "Horseshoe Bay -> Whistler (source line order)",
        "chainage_crs": METRIC_CRS,
        "output_crs": GEOGRAPHIC_CRS,
        "requested_segment_length_m": segment_length_m,
        "source_length_m": round(route.length, 3),
        "source_quality_length_km": float(source_qa["candidate_length_km"]),
        "segment_count": len(metric),
        "full_length_segment_count": sum(abs(x - segment_length_m) < 0.1 for x in lengths),
        "last_segment_length_m": round(float(lengths.iloc[-1]), 3),
        "min_segment_length_m": round(float(lengths.min()), 3),
        "max_segment_length_m": round(float(lengths.max()), 3),
        "sum_of_segment_lengths_m": round(float(lengths.sum()), 3),
        "geojson_roundtrip_length_difference_m": round(difference_m, 6),
        "max_boundary_gap_m": round(max_gap, 6),
        "ids_unique": True,
        "input_one_ordered_linestring": True,
        "geometry_simplified": False,
        "segment_definition": "Fixed metric chainage along provisional Highway 99 centerline",
        "not_equivalent_to": ["landslide hazard zones", "slope units", "traffic-routing edges"],
        "caveats": [
            "Source Step 3.2 centerline remains provisional, not an official surveyed road centerline.",
            "One-kilometre segmentation is a first modelling resolution, not yet justified by hazard data.",
            "Results do not encode carriageways, lane count, road width, or direction of legal travel.",
            "Segmentation does not infer landslide occurrence, runout, closure, or impact.",
        ],
        "outputs": [
            f"data/processed/study_area/{name}.geojson",
            f"data/processed/study_area/{name}.csv",
            f"data/processed/study_area/{name}_quality.json",
            f"results/maps/{name}.png" if make_figure else "(map omitted)",
        ],
    }
    qa_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if make_figure:
        make_map(metric, map_path, segment_length_m)

    print(f"Source length: {route.length / 1000:.3f} km")
    print(f"Segment length target: {segment_length_m:,} m")
    print(f"Generated segments: {len(metric)}; last segment: {lengths.iloc[-1]:.2f} m")
    print(f"Sum-length deviation: {difference_m:.6f} m; maximum gap: {max_gap:.6f} m")
    print(f"Output status: {metadata['status']}")
    print("OUTPUTS:")
    for path in [geojson_path, csv_path, qa_path] + ([map_path] if make_figure else []):
        print("  ", path)
    print("Review the segment map in QGIS before committing. No OSM downloads were performed.")
    return metadata


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--segment-length-m", type=int, default=DEFAULT_LENGTH_M,
                        help="Length of each chainage segment in metres (default: 1000).")
    parser.add_argument("--source", type=Path, default=SOURCE,
                        help="Step 3.2 candidate centerline GeoJSON.")
    parser.add_argument("--source-quality", type=Path, default=SOURCE_QUALITY,
                        help="Step 3.2 quality JSON.")
    parser.add_argument("--out-data", type=Path, default=OUT_DATA)
    parser.add_argument("--out-maps", type=Path, default=OUT_MAPS)
    parser.add_argument("--skip-map", action="store_true")
    args = parser.parse_args(argv)
    build_outputs(args.source, args.source_quality, args.out_data, args.out_maps,
                  args.segment_length_m, not args.skip_map)


if __name__ == "__main__":
    main()
