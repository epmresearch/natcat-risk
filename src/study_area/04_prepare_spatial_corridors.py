"""Step 3.4: provisional GIS screening corridors around Highway 99 segments.

Offline GIS preparation ONLY. Buffered polygons are NOT landslide hazard, runout,
road-width, or impact footprints; they must not be used as a physical risk model.
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
from matplotlib.patches import Patch
from shapely.geometry import Point

ROOT = Path(__file__).resolve().parents[2]
METRIC_CRS = "EPSG:26910"  # NAD83 / UTM 10N; suitable for local metric distances
EXPORT_CRS = "EPSG:4326"  # portable GeoJSON longitude/latitude
DEFAULT_NEAR_M = 1000
DEFAULT_SEARCH_M = 5000
SOURCE = ROOT / "data/processed/study_area/highway99_segments_1000m.geojson"
SOURCE_QA = ROOT / "data/processed/study_area/highway99_segments_1000m_quality.json"
OUT_DATA = ROOT / "data/processed/study_area"
OUT_MAPS = ROOT / "results/maps"


def load_segments(source: Path, quality: Path) -> tuple[gpd.GeoDataFrame, dict]:
    """Validate segmentation provenance, IDs, chainage, and metric continuity."""
    if not source.is_file():
        raise FileNotFoundError(f"Step 3.3 GeoJSON not found: {source}")
    if not quality.is_file():
        raise FileNotFoundError(f"Step 3.3 QA JSON not found: {quality}")
    meta = json.loads(quality.read_text(encoding="utf-8"))
    if meta.get("step") != "3.3":
        raise ValueError("QA file is not from Step 3.3.")
    segments = gpd.read_file(source)
    needed = {"segment_id", "segment_no", "start_m", "end_m", "length_m"}
    if segments.empty or segments.crs is None or not needed.issubset(segments.columns):
        raise ValueError("Input must have a CRS and Step 3.3 chainage attributes.")
    if segments.geometry.isna().any() or not all(segments.geom_type == "LineString"):
        raise ValueError("Expected nonempty LineString geometries only.")
    if (segments.geometry.is_empty | ~segments.geometry.is_valid).any():
        raise ValueError("Invalid or empty segment geometry.")
    if segments["segment_id"].duplicated().any() or segments["segment_no"].duplicated().any():
        raise ValueError("Duplicate segment identifiers.")
    segments = segments.sort_values("segment_no").reset_index(drop=True)
    if list(segments["segment_no"].astype(int)) != list(range(1, len(segments) + 1)):
        raise ValueError("Segment numbers must be sequential 1..N.")
    expected_ids = [f"H99-{i:04d}" for i in range(1, len(segments) + 1)]
    if list(segments["segment_id"]) != expected_ids:
        raise ValueError("Unexpected segment IDs or sequence.")
    if int(meta.get("segment_count", -1)) != len(segments):
        raise ValueError("Input count differs from Step 3.3 quality report.")
    metric = segments.to_crs(METRIC_CRS)
    measured = metric.geometry.length.to_numpy()
    if not all(math.isfinite(float(v)) and float(v) > 0 for v in measured):
        raise ValueError("Invalid metric segment lengths.")
    length_diffs = [abs(float(r) - float(g)) for r, g in zip(metric["length_m"], measured)]
    if max(length_diffs) > 0.1:
        raise ValueError("Stored lengths disagree with segment geometries (>0.1 m).")
    chainage = metric[["start_m", "end_m"]].astype(float).to_numpy()
    if abs(chainage[0, 0]) > 0.01 or any(float(e) <= float(s) for s, e in chainage):
        raise ValueError("Invalid start or end chainage.")
    if max(abs((chainage[:, 1] - chainage[:, 0]) - measured)) > 0.1:
        raise ValueError("Chainage span disagrees with geometry length.")
    if len(metric) > 1:
        chain_gaps = [abs(float(chainage[i, 1] - chainage[i + 1, 0])) for i in range(len(metric) - 1)]
        geometry_gaps = [Point(metric.geometry.iloc[i].coords[-1]).distance(Point(metric.geometry.iloc[i + 1].coords[0]))
                         for i in range(len(metric) - 1)]
        if max(chain_gaps) > 0.01 or max(geometry_gaps) > 0.02:
            raise ValueError("Step 3.3 source contains chainage or geometry gaps.")
    if abs(float(sum(measured)) - float(meta.get("sum_of_segment_lengths_m", -1))) > 0.25:
        raise ValueError("Measured segment total differs from Step 3.3 QA (>0.25 m).")
    return metric, meta


def prepare_footprints(metric: gpd.GeoDataFrame, near_m: int, search_m: int):
    """Build overlapping per-segment proximity zones and dissolved corridor polygons."""
    for value, name in ((near_m, "near_m"), (search_m, "search_m")):
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise ValueError(f"{name} must be a positive integer number of metres.")
    if search_m <= near_m:
        raise ValueError("search_m must exceed near_m; the two corridors serve different purposes.")
    if metric.empty or metric.crs is None or str(metric.crs).upper() != METRIC_CRS:
        raise ValueError("Input segments must be nonempty and in EPSG:26910.")

    # Flat per-segment screening zones overlap deliberately near common endpoints.
    # Do NOT sum their areas as if they were disjoint exposure units.
    per_zone = metric[["segment_id", "segment_no", "start_m", "end_m", "length_m"]].copy()
    per_zone["screening_distance_m"] = near_m
    per_zone["source_status"] = "PROVISIONAL_SPATIAL_SCREENING_ONLY"
    per_zone = gpd.GeoDataFrame(per_zone, geometry=metric.geometry.buffer(near_m), crs=METRIC_CRS)
    per_zone["zone_area_m2"] = per_zone.geometry.area.round(2)

    # Dissolve the overlapping zones to represent unique area around the entire line.
    near_union = per_zone.geometry.union_all()
    # Search buffer around the same underlying road geometry, NOT a landslide runout.
    road_union = metric.geometry.union_all()
    search_union = road_union.buffer(search_m)
    near = gpd.GeoDataFrame([{"distance_m": near_m, "purpose": "NEAR_ROAD_SPATIAL_SCREENING_ONLY"}],
                            geometry=[near_union], crs=METRIC_CRS)
    search = gpd.GeoDataFrame([{"distance_m": search_m, "purpose": "TERRAIN_DATA_DISCOVERY_ONLY"}],
                              geometry=[search_union], crs=METRIC_CRS)
    for label, geometry in (("dissolved near-road corridor", near_union), ("terrain search corridor", search_union)):
        if geometry.is_empty or not geometry.is_valid or geometry.area <= 0:
            raise ValueError(f"{label} is invalid or empty.")
    uncovered_near_area = near_union.difference(search_union).area
    if uncovered_near_area > 1.0:
        raise ValueError(f"Near-road corridor extends outside terrain search corridor by {uncovered_near_area:.3f} m².")
    sum_area = float(per_zone.geometry.area.sum())
    union_area = float(near_union.area)
    if sum_area + 1.0 < union_area:
        raise ValueError("Dissolved near-road area exceeds sum of individual zone areas.")
    search_parts = list(search_union.geoms) if search_union.geom_type == "MultiPolygon" else [search_union]
    metrics = {
        "segment_count": len(per_zone),
        "terrain_search_polygon_part_count": len(search_parts),
        "terrain_search_smallest_part_area_m2": round(min(p.area for p in search_parts), 3),
        "sum_individual_segment_zone_areas_km2": round(sum_area / 1e6, 6),
        "unique_near_road_corridor_area_km2": round(union_area / 1e6, 6),
        "overlapping_area_counted_more_than_once_km2": round((sum_area - union_area) / 1e6, 6),
        "unique_terrain_data_search_area_km2": round(float(search_union.area) / 1e6, 6),
        "uncovered_near_area_m2": round(uncovered_near_area, 6),
    }
    return per_zone, near, search, metrics


def render_map(lines: gpd.GeoDataFrame, near: gpd.GeoDataFrame, search: gpd.GeoDataFrame, out: Path,
               near_m: int, search_m: int) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 11.2))
    search.plot(ax=ax, color="#B9D8E8", edgecolor="#709EB6", linewidth=0.8, alpha=0.45, zorder=1)
    near.plot(ax=ax, color="#56A8A2", edgecolor="#308E88", linewidth=0.7, alpha=0.65, zorder=2)
    lines.plot(ax=ax, color="#B63030", linewidth=1.7, zorder=3)
    a = Point(lines.geometry.iloc[0].coords[0]); b = Point(lines.geometry.iloc[-1].coords[-1])
    ax.scatter([a.x, b.x], [a.y, b.y], color="#202124", s=30, zorder=4)
    ax.annotate("Horseshoe Bay", (a.x, a.y), xytext=(6, -13), textcoords="offset points", fontsize=8)
    ax.annotate("Whistler", (b.x, b.y), xytext=(6, 5), textcoords="offset points", fontsize=8)
    bounds = search.total_bounds
    ax.set_xlim(bounds[0] - 3000, bounds[2] + 3000)
    ax.set_ylim(bounds[1] - 2500, bounds[3] + 2500)
    ax.set_aspect("equal")
    ax.grid(alpha=0.15)
    ax.set_title("Highway 99 — provisional GIS screening corridors\nHorseshoe Bay to Whistler | Step 3.4", fontsize=12)
    ax.set_xlabel("Easting (m), NAD83 / UTM zone 10N")
    ax.set_ylabel("Northing (m), NAD83 / UTM zone 10N")
    ax.legend(handles=[
        Patch(facecolor="#B9D8E8", edgecolor="#709EB6", alpha=0.65, label=f"{search_m/1000:g} km terrain data search"),
        Patch(facecolor="#56A8A2", edgecolor="#308E88", alpha=0.8, label=f"{near_m/1000:g} km near-road screening"),
        Patch(facecolor="#B63030", label="Provisional OSM Highway 99 centerline"),
    ], loc="upper left", fontsize=8)
    fig.text(0.5, 0.018,
             "© OpenStreetMap contributors | Buffers are DATA-SCREENING regions,\n"
             "NOT landslide susceptibility, initiation, runout, or impact zones.",
             ha="center", fontsize=8)
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160)
    plt.close(fig)


def build_outputs(source: Path = SOURCE, quality: Path = SOURCE_QA,
                  out_data: Path = OUT_DATA, out_maps: Path = OUT_MAPS,
                  near_m: int = DEFAULT_NEAR_M, search_m: int = DEFAULT_SEARCH_M,
                  make_map: bool = True) -> dict:
    metric, step33_qa = load_segments(source, quality)
    per_zone, near, search, metrics = prepare_footprints(metric, near_m, search_m)
    out_data.mkdir(parents=True, exist_ok=True)
    per_file = out_data / f"highway99_segment_screening_zones_{near_m}m.geojson"
    near_file = out_data / f"highway99_corridor_screening_{near_m}m.geojson"
    search_file = out_data / f"highway99_terrain_search_{search_m}m.geojson"
    qa_file = out_data / "highway99_spatial_preparation_quality.json"
    map_file = out_maps / "highway99_spatial_screening_corridors.png"
    for gdf, dest in ((per_zone, per_file), (near, near_file), (search, search_file)):
        gdf.to_crs(EXPORT_CRS).to_file(dest, driver="GeoJSON")
        check = gpd.read_file(dest)
        if len(check) != len(gdf) or check.crs is None or not check.is_valid.all():
            raise ValueError(f"GeoJSON export QA failed: {dest}")
        reprojected = check.to_crs(METRIC_CRS)
        for orig_geom, saved_geom in zip(gdf.geometry, reprojected.geometry):
            if abs(orig_geom.area - saved_geom.area) > max(0.5, orig_geom.area * 1e-5):
                raise ValueError(f"Area changed unexpectedly during GeoJSON export: {dest}")
    qa = {
        "step": "3.4",
        "status": "PROVISIONAL_SPATIAL_SCREENING_NOT_A_HAZARD_MODEL",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input_step33_geojson": "data/processed/study_area/highway99_segments_1000m.geojson",
        "input_step33_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "input_step33_quality": "data/processed/study_area/highway99_segments_1000m_quality.json",
        "step33_source_status": step33_qa.get("status"),
        "original_data_source": "OpenStreetMap via OSMnx/Overpass; © OpenStreetMap contributors",
        "metric_crs": METRIC_CRS,
        "geojson_crs": EXPORT_CRS,
        "near_road_screening_distance_m": near_m,
        "terrain_data_search_distance_m": search_m,
        "input_total_segment_length_m": round(float(metric.geometry.length.sum()), 3),
        **metrics,
        "qa_flags": {
            "segment_ids_preserved": list(per_zone["segment_id"]) == list(metric["segment_id"]),
            "near_corridor_inside_search_corridor": metrics["uncovered_near_area_m2"] <= 1.0,
            "no_route_geometry_changes": True,
            "exports_reopened_and_validated": True,
        },
        "limitations": [
            "Buffers are distance-based GIS screening regions, not landslide hazard or runout boundaries.",
            "A landslide initiated outside either buffer may still reach or disrupt the highway.",
            "Overlapping per-segment polygons MUST NOT be summed for unique area or counted as independent exposures.",
            "Future hazard points/polygons may intersect multiple adjacent segment zones: explicit assignment/deduplication is required.",
            "Buffer widths are provisional methodological parameters and require sensitivity checks.",
            "The terrain search area can include tiny disconnected GIS polygon parts; inspect part-count QA and do not infer isolated hazards.",
            "Source road geometry is a provisional OSM-derived centerline, not official surveyed geometry.",
            "No DEM, landslide inventory, susceptibility, physical damage, or CLIMADA hazard has been modelled.",
        ],
        "outputs": [
            f"data/processed/study_area/{per_file.name}",
            f"data/processed/study_area/{near_file.name}",
            f"data/processed/study_area/{search_file.name}",
            f"data/processed/study_area/{qa_file.name}",
            f"results/maps/{map_file.name}" if make_map else "(map skipped)",
        ],
    }
    qa_file.write_text(json.dumps(qa, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if make_map:
        render_map(metric, near, search, map_file, near_m, search_m)
    print("STEP 3.4 — GIS SPATIAL CORRIDORS (NOT HAZARD ZONES)")
    print(f"Input segments: {len(metric)}; near-screening: {near_m:,} m; terrain search: {search_m:,} m")
    print(f"Unique near-road corridor area: {metrics['unique_near_road_corridor_area_km2']:.3f} km²")
    print(f"Terrain data search area: {metrics['unique_terrain_data_search_area_km2']:.3f} km²")
    print("QA: segmentation preserved; overlap tracked; near-road corridor fully inside search corridor")
    print("OUTPUTS:")
    for p in [per_file, near_file, search_file, qa_file] + ([map_file] if make_map else []):
        print("  ", p)
    print("Review in QGIS with OSM basemap. Do not commit until the outputs are validated.")
    return qa


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--near-m", type=int, default=DEFAULT_NEAR_M)
    parser.add_argument("--search-m", type=int, default=DEFAULT_SEARCH_M)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--quality", type=Path, default=SOURCE_QA)
    parser.add_argument("--out-data", type=Path, default=OUT_DATA)
    parser.add_argument("--out-maps", type=Path, default=OUT_MAPS)
    parser.add_argument("--skip-map", action="store_true")
    args = parser.parse_args(argv)
    build_outputs(args.source, args.quality, args.out_data, args.out_maps,
                  args.near_m, args.search_m, not args.skip_map)


if __name__ == "__main__":
    main()
