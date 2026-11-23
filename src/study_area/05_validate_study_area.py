"""Step 3.5: cross-check the provisional Highway 99 study-area GIS deliverables.

This is an OFFLINE internal consistency and reproducibility audit, not an
independent survey, road-authority validation, landslide model, or risk map.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from shapely.geometry import LineString, Point

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path("data/processed/study_area")
QA_OUT = Path("results/diagnostics/highway99_study_area_step35_audit.json")
MAP_OUT = Path("results/maps/highway99_study_area_step35_audit.png")
METRIC_CRS = "EPSG:26910"
GEO_CRS = "EPSG:4326"

FILES = {
    "centerline": DATA_DIR / "highway99_candidate_centerline.geojson",
    "centerline_qa": DATA_DIR / "highway99_candidate_quality.json",
    "segments": DATA_DIR / "highway99_segments_1000m.geojson",
    "segments_csv": DATA_DIR / "highway99_segments_1000m.csv",
    "segments_qa": DATA_DIR / "highway99_segments_1000m_quality.json",
    "segment_zones": DATA_DIR / "highway99_segment_screening_zones_1000m.geojson",
    "near_corridor": DATA_DIR / "highway99_corridor_screening_1000m.geojson",
    "terrain_search": DATA_DIR / "highway99_terrain_search_5000m.geojson",
    "spatial_qa": DATA_DIR / "highway99_spatial_preparation_quality.json",
}


class ValidationError(ValueError):
    """Input is not fit for this reproducible study-area audit."""


def demand(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_qa(path: Path, expected_step: str) -> dict:
    try:
        qa = json.loads(path.read_text(encoding="utf-8"))
    except (ValueError, OSError) as exc:
        raise ValidationError(f"Cannot read QA JSON {path.name}: {exc}") from exc
    demand(isinstance(qa, dict) and qa.get("step") == expected_step,
           f"Incorrect Step {expected_step} quality metadata: {path.name}")
    return qa


def read_layer(path: Path, expected_type: str, expected_count: int | None = None) -> gpd.GeoDataFrame:
    try:
        layer = gpd.read_file(path)
    except Exception as exc:
        raise ValidationError(f"Cannot read GIS layer {path.name}: {exc}") from exc
    demand(layer.crs is not None and layer.crs.to_epsg() == 4326,
           f"{path.name} must declare geographic EPSG:4326/CRS84 coordinates.")
    demand(not layer.empty, f"{path.name} is empty.")
    if expected_count is not None:
        demand(len(layer) == expected_count,
               f"{path.name}: expected {expected_count} features, found {len(layer)}.")
    demand(not (layer.geometry.isna() | layer.geometry.is_empty).any() and layer.is_valid.all(),
           f"{path.name} contains invalid or missing geometries.")
    if expected_type == "line":
        demand(all(layer.geom_type == "LineString"), f"{path.name} must contain LineStrings.")
    else:
        demand(all(layer.geom_type.isin(["Polygon", "MultiPolygon"])),
               f"{path.name} must contain polygons.")
    return layer.to_crs(METRIC_CRS)


def approx(actual: float, target: float, tol: float, label: str) -> None:
    demand(math.isfinite(actual) and math.isfinite(target) and abs(actual - target) <= tol,
           f"{label}: observed {actual:.6f}, expected {target:.6f}, tolerance {tol:.6f}.")


def validate_inputs(root: Path) -> tuple[dict, dict]:
    """Validate actual files and return (portable audit JSON, geometries for map)."""
    root = Path(root)
    paths = {k: root / p for k, p in FILES.items()}
    missing = [str(FILES[k]).replace("\\", "/") for k, p in paths.items() if not p.is_file()]
    demand(not missing, "Missing required files:\n  " + "\n  ".join(missing))

    qa32 = read_qa(paths["centerline_qa"], "3.2")
    qa33 = read_qa(paths["segments_qa"], "3.3")
    qa34 = read_qa(paths["spatial_qa"], "3.4")
    checks: list[str] = []

    demand(qa33.get("source_sha256") == sha256(paths["centerline"]),
           "Step 3.3 source_sha256 differs from the Step 3.2 centerline file; do not mix route versions.")
    demand(qa34.get("input_step33_sha256") == sha256(paths["segments"]),
           "Step 3.4 input_step33_sha256 differs from the Step 3.3 segments GeoJSON.")
    checks.append("Step 3.2 -> Step 3.3 -> Step 3.4 SHA-256 provenance chain matches")

    demand(int(qa33.get("requested_segment_length_m", 0)) == 1000,
           "Step 3.3 must use the committed 1000 m segmentation scheme.")
    count = int(qa33.get("segment_count", 0))
    demand(count > 0 and count == int(qa34.get("segment_count", -1)),
           "Segment counts in Step 3.3 and Step 3.4 quality reports differ.")
    center = read_layer(paths["centerline"], "line", 1)
    segments = read_layer(paths["segments"], "line", count)
    zones = read_layer(paths["segment_zones"], "polygon", count)
    near = read_layer(paths["near_corridor"], "polygon", 1)
    search = read_layer(paths["terrain_search"], "polygon", 1)
    checks.append("All GIS files reopen, have declared geographic CRS and valid geometries")

    center_geom: LineString = center.geometry.iloc[0]
    center_length = float(center_geom.length)
    demand(segments["segment_id"].is_unique and segments["segment_no"].is_unique,
           "Segment IDs or segment numbers are not unique.")
    segments = segments.sort_values("segment_no").reset_index(drop=True)
    expected = list(range(1, count + 1))
    demand(list(segments["segment_no"].astype(int)) == expected, "Segment numbers must be 1..N.")
    demand(list(segments["segment_id"]) == [f"H99-{i:04d}" for i in expected],
           "Unexpected segment IDs.")
    demand("segment_id" in zones and zones["segment_id"].is_unique,
           "Per-segment zones need unique segment_id fields.")
    zones = zones.set_index("segment_id").loc[segments["segment_id"]].reset_index()
    checks.append("Sequential segment identifiers follow declared count and zone association")

    needed = {"start_m", "end_m", "length_m"}
    demand(needed.issubset(segments.columns), "Step 3.3 chainage fields are missing.")
    distances = []
    starts, ends = [], []
    for _, seg in segments.iterrows():
        measured = float(seg.geometry.length)
        approx(measured, float(seg["length_m"]), 0.2, f"{seg.segment_id} stored length")
        approx(float(seg["end_m"]) - float(seg["start_m"]), measured, 0.2,
               f"{seg.segment_id} chainage interval")
        starts.append(float(seg["start_m"]))
        ends.append(float(seg["end_m"]))
    approx(starts[0], 0.0, 0.02, "First segment chainage")
    for i in range(count - 1):
        approx(ends[i], starts[i + 1], 0.02, f"Chainage at segment {i+1}")
        distances.append(Point(segments.geometry.iloc[i].coords[-1]).distance(
            Point(segments.geometry.iloc[i + 1].coords[0])))
    max_gap = max(distances, default=0.0)
    demand(max_gap <= 0.05, f"Geometry gap between adjacent segments: {max_gap:.6f} m.")
    sum_length = float(segments.geometry.length.sum())
    approx(sum_length, center_length, 0.5, "Segment total vs centerline")
    approx(ends[-1], center_length, 0.5, "Last chainage vs centerline")
    approx(center_length, float(qa33["source_length_m"]), 0.5, "Step 3.3 QA centerline length")
    approx(center_length, float(qa32["candidate_length_km"]) * 1000, 100,
           "Step 3.2 reported rounded centerline length")

    # A single ordered polyline must coincide with the underlying unsegmented route.
    merged_coords = []
    for i, geom in enumerate(segments.geometry):
        coords = list(geom.coords)
        merged_coords.extend(coords if i == 0 else coords[1:])
    reconstructed = LineString(merged_coords)
    hausdorff = center_geom.hausdorff_distance(reconstructed)
    demand(hausdorff <= 0.2, f"Reconstructed segment path differs from centerline: {hausdorff:.4f} m.")
    demand(Point(center_geom.coords[0]).distance(Point(reconstructed.coords[0])) < 0.2
           and Point(center_geom.coords[-1]).distance(Point(reconstructed.coords[-1])) < 0.2,
           "Route endpoints/order changed in Step 3.3.")
    checks.append("Route order, segment chainage, continuous geometry and centerline overlay verified")

    with paths["segments_csv"].open(encoding="utf-8-sig", newline="") as f:
        table = list(csv.DictReader(f))
    demand(len(table) == count, "Step 3.3 CSV record count differs from GeoJSON.")
    for row, (_, seg) in zip(table, segments.iterrows()):
        demand(row.get("segment_id") == seg.segment_id,
               "CSV segment order or identifier differs from GeoJSON.")
        approx(float(row["start_m"]), float(seg.start_m), 0.001, "CSV start_m")
        approx(float(row["end_m"]), float(seg.end_m), 0.001, "CSV end_m")
    checks.append("CSV table matches GIS segment order and chainage")

    near_m, search_m = int(qa34["near_road_screening_distance_m"]), int(qa34["terrain_data_search_distance_m"])
    demand((near_m, search_m) == (1000, 5000),
           "Step 3.4 is expected to use the committed 1000 m and 5000 m screening buffers.")
    for _, seg in segments.iterrows():
        zone = zones.loc[zones["segment_id"] == seg.segment_id].geometry.iloc[0]
        expected_zone = seg.geometry.buffer(near_m)
        diff = expected_zone.symmetric_difference(zone).area
        demand(diff <= 2.0, f"{seg.segment_id} {near_m} m zone differs from reproducible buffer by {diff:.3f} m².")
    near_geom = near.geometry.iloc[0]
    search_geom = search.geometry.iloc[0]
    near_expected = zones.geometry.union_all()
    search_expected = segments.geometry.union_all().buffer(search_m)
    near_diff = near_geom.symmetric_difference(near_expected).area
    search_diff = search_geom.symmetric_difference(search_expected).area
    demand(near_diff <= 2.0 and search_diff <= 2.0,
           f"Corridor buffers are not reproducible (near XOR={near_diff:.3f}, search XOR={search_diff:.3f} m²).")
    uncovered = near_geom.difference(search_geom).area
    demand(uncovered <= 1.0, f"Near-road corridor outside terrain search: {uncovered:.3f} m².")
    checks.append("Per-segment 1 km buffers, unique 1 km union and 5 km terrain buffer reproduced")

    sum_zone_area = float(zones.geometry.area.sum()) / 1e6
    near_area = float(near_geom.area) / 1e6
    search_area = float(search_geom.area) / 1e6
    overlap = sum_zone_area - near_area
    for actual, key in (
        (sum_zone_area, "sum_individual_segment_zone_areas_km2"),
        (near_area, "unique_near_road_corridor_area_km2"),
        (search_area, "unique_terrain_data_search_area_km2"),
        (overlap, "overlapping_area_counted_more_than_once_km2"),
    ):
        approx(actual, float(qa34[key]), 0.0001, f"Step 3.4 QA {key}")
    demand(overlap >= 0, "Invalid negative overlap area.")
    demand(all(bool(x) for x in qa34.get("qa_flags", {}).values())
           and len(qa34.get("qa_flags", {})) >= 4, "Step 3.4 QA flags are missing/false.")
    checks.append("Unique areas match recorded QA; per-segment overlaps counted separately")

    parts = list(search_geom.geoms) if search_geom.geom_type == "MultiPolygon" else [search_geom]
    count_parts = len(parts)
    smallest = min(parts, key=lambda p: p.area)
    approx(float(smallest.area), float(qa34["terrain_search_smallest_part_area_m2"]), 0.1,
           "Search polygon smallest disconnected part")
    demand(count_parts == int(qa34["terrain_search_polygon_part_count"]),
           "Search polygon part count differs from Step 3.4 QA.")
    smallest_location = gpd.GeoSeries([smallest.centroid], crs=METRIC_CRS).to_crs(GEO_CRS).iloc[0]
    checks.append("Disconnected 5 km search polygon components measured and documented")

    report = {
        "step": "3.5",
        "status": "INTERNAL_QA_PASSED_WITH_CAVEATS",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "method": "Offline, reproducible cross-stage geospatial and metadata consistency audit",
        "independent_official_road_validation_completed": False,
        "hazard_or_landslide_model_completed": False,
        "readiness": "READY_FOR_TERRAIN_AND_LANDSLIDE_DATA_DISCOVERY_WITH_CAVEATS",
        "input_sha256": {str(FILES[k]).replace("\\", "/"): sha256(p) for k, p in paths.items()},
        "checks_passed": checks,
        "measurements": {
            "segment_count": count,
            "centerline_length_km": round(center_length / 1000, 6),
            "max_segment_boundary_gap_m": round(max_gap, 9),
            "reconstructed_line_hausdorff_m": round(float(hausdorff), 9),
            "near_screening_distance_m": near_m,
            "terrain_search_distance_m": search_m,
            "sum_overlapping_segment_zones_km2": round(sum_zone_area, 6),
            "unique_near_screening_area_km2": round(near_area, 6),
            "overlap_not_to_double_count_km2": round(overlap, 6),
            "unique_terrain_search_area_km2": round(search_area, 6),
            "terrain_search_polygon_parts": count_parts,
            "smallest_polygon_part_area_m2": round(smallest.area, 3),
            "smallest_polygon_part_centroid_lonlat": [round(smallest_location.x, 7), round(smallest_location.y, 7)],
        },
        "unresolved_items": [
            "The Step 3.2 OSM highway centerline remains provisional and has not been independently checked against an official BC road centerline.",
            "The 1 km segmentation and 1 km/5 km buffers are analysis/data-screening choices; conduct resolution and buffer sensitivity analyses later.",
            *(["The disconnected terrain-search polygon part was retained; inspect in QGIS and do not interpret it as a separate hazard."] if count_parts > 1 else []),
            "Per-segment buffers overlap: risk events intersecting multiple zones need explicit assignment/deduplication.",
            "No slope, landslide inventory, susceptibility, runout, road closure, loss, or CLIMADA risk outputs exist yet.",
        ],
        "attribution": "© OpenStreetMap contributors; OSM-derived provisional road geometry",
    }
    layers = {"center": center, "segments": segments, "near": near, "search": search,
              "smallest": smallest, "near_m": near_m, "search_m": search_m}
    return report, layers


def render_map(layers: dict, out: Path) -> None:
    """Overview plus a diagnostic inset around the smallest search polygon part."""
    fig, (ax, zoom) = plt.subplots(1, 2, figsize=(11, 10), gridspec_kw={"width_ratios": [2.2, 1]})
    layers["search"].plot(ax=ax, facecolor="#B9D8E8", edgecolor="#709EB6", alpha=0.50, linewidth=0.7)
    layers["near"].plot(ax=ax, facecolor="#5CB4A6", edgecolor="#2E8178", alpha=0.58, linewidth=0.6)
    layers["center"].plot(ax=ax, color="#A92B31", linewidth=1.3)
    ax.scatter([layers["smallest"].centroid.x], [layers["smallest"].centroid.y], s=28,
               c="#D58A21", edgecolor="black", linewidth=0.6, zorder=5)
    bounds = layers["search"].total_bounds
    ax.set_xlim(bounds[0]-1800, bounds[2]+1800)
    ax.set_ylim(bounds[1]-2500, bounds[3]+2500)
    ax.set_aspect("equal")
    ax.grid(alpha=0.14)
    ax.set_title("Highway 99 screening corridors\nInternal QA overview", fontsize=12)
    ax.set_xlabel("Easting (m), NAD83 / UTM 10N")
    ax.set_ylabel("Northing (m)")
    ax.legend(handles=[
        Patch(facecolor="#B9D8E8", label="5 km terrain data search"),
        Patch(facecolor="#5CB4A6", label="1 km near-road screening"),
        Patch(facecolor="#A92B31", label="Provisional OSM H99"),
        Patch(facecolor="#D58A21", label="Small disconnected part"),
    ], loc="upper left", fontsize=8)

    poly = layers["smallest"]
    zoom.fill(*poly.exterior.xy, color="#D58A21", alpha=0.50)
    zoom.plot(*poly.exterior.xy, color="#AD5D17", linewidth=1.5)
    x, y = poly.centroid.x, poly.centroid.y
    zoom.set_xlim(x - 200, x + 200)
    zoom.set_ylim(y - 200, y + 200)
    zoom.set_aspect("equal")
    zoom.grid(alpha=0.20)
    zoom.set_title(f"Smallest 5 km polygon part\n{poly.area:.1f} m² (diagnostic only)", fontsize=10)
    zoom.set_xlabel("Easting (m)")
    fig.suptitle("Step 3.5 — Cross-stage validation (not a landslide hazard map)", y=0.99, fontsize=13)
    fig.text(0.5, 0.012, "© OpenStreetMap contributors | Provisional OSM route; no official road or hazard validation",
             ha="center", fontsize=9)
    fig.tight_layout(rect=[0, 0.028, 1, 0.965])
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160)
    plt.close(fig)


def run(root: Path = ROOT, out_json: Path | None = None,
        out_map: Path | None = None, make_map: bool = True) -> dict:
    root = Path(root)
    out_json = Path(out_json) if out_json is not None else root / QA_OUT
    out_map = Path(out_map) if out_map is not None else root / MAP_OUT
    report, layers = validate_inputs(root)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if make_map:
        render_map(layers, out_map)
    print("STEP 3.5 — STUDY AREA QA: INTERNAL_QA_PASSED_WITH_CAVEATS")
    for line in report["checks_passed"]:
        print("  PASS:", line)
    m = report["measurements"]
    print(f"Segments: {m['segment_count']}; route: {m['centerline_length_km']:.3f} km")
    print(f"Near corridor: {m['unique_near_screening_area_km2']:.3f} km²; "
          f"terrain search: {m['unique_terrain_search_area_km2']:.3f} km²")
    print(f"Disconnected search parts: {m['terrain_search_polygon_parts']}; "
          f"smallest: {m['smallest_polygon_part_area_m2']:.3f} m²")
    print("IMPORTANT: this is internal QA, NOT official road validation or landslide hazard modelling.")
    print("OUTPUTS:", out_json, *([out_map] if make_map else []), sep="\n  ")
    return report


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT, help="Repository root containing data/processed/study_area")
    parser.add_argument("--out-json", type=Path, default=None)
    parser.add_argument("--out-map", type=Path, default=None)
    parser.add_argument("--skip-map", action="store_true")
    args = parser.parse_args(argv)
    run(args.root, args.out_json, args.out_map, not args.skip_map)


if __name__ == "__main__":
    main()
