"""Step 3.2: build a *provisional* Highway 99 road alignment from OpenStreetMap.

Run from the repository root after Step 3.1:
    python src/study_area/02_extract_highway99_candidate.py

Requirements: osmnx >= 2.1, geopandas, networkx, shapely, matplotlib.

This is a network-derived candidate alignment, NOT the authoritative BC road
centerline and NOT a traffic/closure simulation. Verify it in QGIS before use.
"""

from __future__ import annotations

from datetime import datetime, timezone
import math
import json
from pathlib import Path
import re
import time

from requests.exceptions import ConnectionError as RequestsConnectionError, Timeout as RequestsTimeout

import geopandas as gpd
import matplotlib.pyplot as plt
import networkx as nx
import osmnx as ox
from shapely.geometry import LineString

REPO = Path(__file__).resolve().parents[2]
POINTS_PATH = REPO / "data/processed/study_area/reference_places.geojson"
OUT_DATA = REPO / "data/processed/study_area"
OUT_MAPS = REPO / "results/maps"
CACHE = REPO / "data/cache/study_area"
# Smaller overlapping download tiles reduce the demand on public Overpass servers.
# Settings are operational only: these widths are not landslide hazard buffers.
MAX_TILE_AXIS_M = 16_000
TILE_BUFFER_M = 9_000
CACHE_VERSION = "v3"
CACHE_STEM = f"highway99_tiles_{MAX_TILE_AXIS_M//1000}km_{TILE_BUFFER_M//1000}km_{CACHE_VERSION}"
GRAPHML = CACHE / f"{CACHE_STEM}_combined.graphml"
TILE_CACHE = CACHE / CACHE_STEM
CRS_GEO = "EPSG:4326"
CRS_METRIC = "EPSG:26910"  # NAD83 / UTM zone 10N, metres

SEARCH_BUFFER_M = 5_000  # For map display / subsequent data discovery only.
MAX_WAYPOINT_SNAP_M = 2_500  # Avoid anchoring a town to a distant tagged segment.
MAX_LOCAL_ROAD_SNAP_M = 2_500  # Conservative fallback if H99 tags are absent nearby.
OFF_HIGHWAY_COST_MULTIPLIER = 3.0

# Courtesy pacing of consecutive *tile* calls (not a server-mandated fixed rate).
# OSMnx also checks Overpass server slot status internally when enabled.
MIN_SECONDS_BETWEEN_TILE_REQUESTS = 30
MAX_ATTEMPTS_PER_TILE = 3
RETRY_WAIT_SECONDS = (30, 60)  # Delay after failures #1 and #2.

WAYPOINTS = ["Horseshoe Bay", "Lions Bay", "Britannia Beach", "Squamish", "Whistler"]


def highway99_tagged(ref: object, name: object) -> bool:
    """Match only plausible Highway 99 labels; tolerate lists in OSM attributes."""
    def values(value: object) -> list[str]:
        if isinstance(value, (list, tuple, set)):
            return [str(v) for v in value]
        if value is None:
            return []
        return [str(value)]

    refs = values(ref)
    names = values(name)
    return any(re.search(r"(?<!\d)99(?![a-zA-Z0-9])", r) for r in refs) or any(
        "sea to sky" in n.lower() or "highway 99" in n.lower() or "hwy 99" in n.lower()
        for n in names
    )


def load_places() -> gpd.GeoDataFrame:
    if not POINTS_PATH.exists():
        raise FileNotFoundError(
            f"Missing {POINTS_PATH}. Run Step 3.1 first to create reference_places.geojson."
        )
    pts = gpd.read_file(POINTS_PATH).to_crs(CRS_GEO)
    if "name" not in pts.columns:
        raise ValueError("reference_places.geojson must contain a 'name' field.")
    if pts["name"].duplicated().any():
        raise ValueError("Duplicate reference-place names were found.")
    missing = set(WAYPOINTS) - set(pts["name"])
    if missing:
        raise ValueError(f"Missing reference places: {sorted(missing)}")
    return pts.set_index("name").loc[WAYPOINTS].reset_index()



def build_download_tiles(points: gpd.GeoDataFrame):
    """Define reproducible overlapping small-area OSM download polygons.

    Long waypoint legs are split into <=MAX_TILE_AXIS_M straight axis pieces.
    All geometry/buffering is done in metres in UTM 10N, never lon/lat degrees.
    Each piece has a TILE_BUFFER_M rounded buffer, so adjacent tile regions
    overlap around their shared endpoint and can share OpenStreetMap node IDs.
    """
    metric = points.to_crs(CRS_METRIC)
    tiles = []
    for leg_index in range(len(metric) - 1):
        start_pt = metric.geometry.iloc[leg_index]
        end_pt = metric.geometry.iloc[leg_index + 1]
        axis = LineString([start_pt, end_pt])
        n_pieces = max(1, math.ceil(axis.length / MAX_TILE_AXIS_M))
        for piece_index in range(n_pieces):
            seg = LineString([
                axis.interpolate(piece_index / n_pieces, normalized=True),
                axis.interpolate((piece_index + 1) / n_pieces, normalized=True),
            ])
            region = gpd.GeoSeries([seg.buffer(TILE_BUFFER_M)], crs=CRS_METRIC)
            region_wgs84 = region.to_crs(CRS_GEO).iloc[0]
            tile_id = f"leg{leg_index+1:02d}_part{piece_index+1:02d}"
            tiles.append((tile_id, region_wgs84))
    return tiles


def is_retryable_overpass_failure(exc: Exception) -> bool:
    """Retry only transient network/capacity errors, not geometry/code errors."""
    if isinstance(exc, (RequestsTimeout, RequestsConnectionError)):
        return True
    name = type(exc).__name__
    if name != "ResponseStatusCodeError":
        return False
    description = str(exc).lower()
    return any(term in description for term in (
        "429", "502", "503", "504", "bad gateway", "gateway timeout",
        "too many requests", "currently at capacity", "server is busy",
    ))


def seconds_until_next_request(last_finished: float | None, now: float) -> float:
    """Minimum 30-second spacing between successive OSMnx tile call completions."""
    if last_finished is None:
        return 0.0
    return max(0.0, MIN_SECONDS_BETWEEN_TILE_REQUESTS - (now - last_finished))


def download_tile_with_pacing(polygon, tile_id: str, tile_index: int,
                              tile_count: int, last_finished: float | None):
    """Serial, respectful and bounded Overpass retry; never spin indefinitely."""
    for attempt in range(1, MAX_ATTEMPTS_PER_TILE + 1):
        remaining = seconds_until_next_request(last_finished, time.monotonic())
        if remaining > 0:
            print(f"    Waiting {remaining:.0f}s to respect minimum tile spacing...", flush=True)
            time.sleep(remaining)
        print(f"    Overpass request for tile {tile_index}/{tile_count} "
              f"({tile_id}), attempt {attempt}/{MAX_ATTEMPTS_PER_TILE}", flush=True)
        try:
            # OSMnx already polls the server's slot-status endpoint if available.
            # Unsimplified graphs preserve OSM node IDs for merging tiles.
            graph = ox.graph.graph_from_polygon(
                polygon,
                network_type="drive",
                simplify=False,
                retain_all=True,
                truncate_by_edge=True,
            )
            last_finished = time.monotonic()
            return graph, last_finished
        except Exception as exc:
            last_finished = time.monotonic()
            if not is_retryable_overpass_failure(exc):
                raise RuntimeError(
                    f"Non-transient failure downloading {tile_id}; not retrying: "
                    f"{type(exc).__name__}: {exc}"
                ) from exc
            if attempt >= MAX_ATTEMPTS_PER_TILE:
                raise RuntimeError(
                    f"Overpass still unavailable after {MAX_ATTEMPTS_PER_TILE} "
                    f"attempts for {tile_id}. Existing cached tiles are safe. "
                    "Do not restart in a tight loop; use offline OSM/BC roads "
                    f"if capacity remains limited. Last error: {exc}"
                ) from exc
            retry_wait = RETRY_WAIT_SECONDS[attempt-1]
            print(f"    Transient Overpass error ({type(exc).__name__}). "
                  f"Waiting {retry_wait}s before retry...", flush=True)
            time.sleep(retry_wait)
    raise AssertionError("Unreachable retry state")


def prepare_graph(points: gpd.GeoDataFrame):
    """Download incremental OSM tiles; merge unsimplified graphs by OSM node ID.

    A failed request does not destroy the other tile caches. Re-running the
    script resumes at the first missing tile; a complete merged graph cache
    skips all further network downloads on future runs.
    """
    if GRAPHML.exists():
        print("Reusing merged OSM drive graph:", GRAPHML)
        return ox.io.load_graphml(GRAPHML)

    TILE_CACHE.mkdir(parents=True, exist_ok=True)
    tiles = build_download_tiles(points)
    print(f"Downloading or loading {len(tiles)} smaller overlapping OSM tiles")
    print(f"Tile axes <= {MAX_TILE_AXIS_M/1000:.0f} km, rounded buffer "
          f"{TILE_BUFFER_M/1000:.0f} km; these are NOT hazard boundaries.")
    graphs = []
    last_request_finished = None
    for tile_index, (tile_id, polygon) in enumerate(tiles, start=1):
        tile_file = TILE_CACHE / f"{tile_id}.graphml"
        if tile_file.exists():
            print(f"  Tile {tile_index}/{len(tiles)} ({tile_id}): cached")
            tile_graph = ox.io.load_graphml(tile_file)
        else:
            print(f"  Tile {tile_index}/{len(tiles)} ({tile_id}): downloading")
            try:
                tile_graph, last_request_finished = download_tile_with_pacing(
                    polygon, tile_id, tile_index, len(tiles), last_request_finished
                )
                if tile_graph.number_of_nodes() == 0 or tile_graph.number_of_edges() == 0:
                    raise RuntimeError(f"Tile {tile_id} returned an empty road network")
                # Save immediately: successful tiles survive a later server failure.
                ox.io.save_graphml(tile_graph, tile_file)
                print(f"    Cached {tile_graph.number_of_nodes():,} nodes and "
                      f"{tile_graph.number_of_edges():,} edges")
            except Exception as exc:
                raise RuntimeError(
                    f"OSM tile {tile_index}/{len(tiles)} ({tile_id}) failed. "
                    "Previously cached tiles remain on disk. "
                    f"Underlying error: {type(exc).__name__}: {exc}"
                ) from exc
        graphs.append(tile_graph)

    # Un-simplified graphs share original OSM node identifiers, enabling a
    # reproducible merge without artificial straight-line connections.
    graph = nx.compose_all(graphs)
    if not graph.number_of_edges():
        raise RuntimeError("Merged OSM network contains no drivable edges")
    sizes = sorted(
        (len(component) for component in nx.connected_components(graph.to_undirected(as_view=True))),
        reverse=True,
    )
    print(f"Merged graph: {graph.number_of_nodes():,} nodes, "
          f"{graph.number_of_edges():,} edges; "
          f"{len(sizes):,} connected components, top 5 sizes: {sizes[:5]}")
    ox.io.save_graphml(graph, GRAPHML)
    print(f"Merged network saved to {GRAPHML}")
    return graph


def build_candidate_route(graph, points: gpd.GeoDataFrame):
    """Prefer Highway 99-tagged roads; retain other road edges for connectivity."""
    # This is an undirected geometric trace, not a legal/time-aware driving route.
    undirected = ox.convert.to_undirected(graph)
    for u, v, k, data in undirected.edges(keys=True, data=True):
        is_h99 = highway99_tagged(data.get("ref"), data.get("name"))
        data["step32_is_h99"] = is_h99
        data["step32_cost"] = float(data.get("length", 1.0)) * (
            1.0 if is_h99 else OFF_HIGHWAY_COST_MULTIPLIER
        )

    # A center-of-town geocode is not a highway junction. Snap waypoints to
    # candidate Highway 99 edges, rather than to arbitrary nearby local roads.
    edge_records = []
    tagged_edge_records = []
    for u, v, k, data in undirected.edges(keys=True, data=True):
        geometry = data.get("geometry")
        if geometry is None:
            geometry = LineString(
                [(undirected.nodes[u]["x"], undirected.nodes[u]["y"]),
                 (undirected.nodes[v]["x"], undirected.nodes[v]["y"])]
            )
        record = {"u": u, "v": v, "geometry": geometry}
        edge_records.append(record)
        if data["step32_is_h99"]:
            tagged_edge_records.append(record)

    if not tagged_edge_records:
        raise RuntimeError(
            "No Highway 99-tagged OSM road edges were found. "
            "Inspect local OSM tags or use the authoritative BC road atlas."
        )
    all_edges = gpd.GeoDataFrame(edge_records, crs=CRS_GEO).to_crs(CRS_METRIC)
    tagged_edges = gpd.GeoDataFrame(tagged_edge_records, crs=CRS_GEO).to_crs(CRS_METRIC)
    metric_points = points.to_crs(CRS_METRIC)
    snaps = []
    for row in metric_points.itertuples():
        tagged_dists = tagged_edges.distance(row.geometry)
        tagged_dist = float(tagged_dists.min())
        local_dists = all_edges.distance(row.geometry)
        local_dist = float(local_dists.min())
        print(
            f"  Check {row.name}: nearest H99-tagged edge {tagged_dist:,.0f} m; "
            f"nearest drivable edge {local_dist:,.0f} m"
        )
        if tagged_dist <= MAX_WAYPOINT_SNAP_M:
            selected = tagged_edges.loc[tagged_dists.idxmin()]
            method = "highway99_tagged"
            selected_dist = tagged_dist
        elif local_dist <= MAX_LOCAL_ROAD_SNAP_M:
            # Do NOT invent a Highway 99 tag. A nearby drivable road is only a
            # *routing anchor* in the draft graph, not evidence that it is H99.
            selected = all_edges.loc[local_dists.idxmin()]
            method = "nearby_road_fallback_UNVERIFIED"
            selected_dist = local_dist
            print(
                f"  WARNING: {row.name} has no nearby H99-tagged edge. "
                "Using nearest drivable road as a PROVISIONAL anchor; "
                "verify the north-corridor alignment and OSM tags in QGIS."
            )
        else:
            raise RuntimeError(
                f"{row.name}: no Highway 99-tagged edge within "
                f"{MAX_WAYPOINT_SNAP_M/1000:.1f} km, and the nearest drivable "
                f"OSM road is {local_dist:,.0f} m away. The downloaded "
                "network likely misses this waypoint. Check the Step 3.1 "
                "coordinates; if correct, enlarge DOWNLOAD_BUFFER_M, delete "
                f"the stale GraphML cache at {GRAPHML}, and rerun."
            )
        candidates = (selected["u"], selected["v"])
        # Snap the anchor to an endpoint of the selected OSM road edge.
        pt_in_wgs = points.loc[points["name"].eq(row.name), "geometry"].iloc[0]
        anchor = min(
            candidates,
            key=lambda n: (undirected.nodes[n]["x"] - pt_in_wgs.x) ** 2
            + (undirected.nodes[n]["y"] - pt_in_wgs.y) ** 2,
        )
        snaps.append({
            "place": row.name,
            "node": int(anchor),
            "snap_method": method,
            "distance_to_tagged_road_m": round(tagged_dist, 1),
            "distance_to_selected_edge_m": round(selected_dist, 1),
            "distance_to_nearest_drivable_edge_m": round(local_dist, 1),
        })
        print(f"  Snap: {row.name}: {method}; {selected_dist:,.0f} m")

    route_edges = []
    coords = []
    for leg, (start, end) in enumerate(zip(snaps[:-1], snaps[1:]), start=1):
        try:
            path = nx.shortest_path(
                undirected, source=start["node"], target=end["node"],
                weight="step32_cost",
            )
        except nx.NetworkXNoPath as exc:
            raise RuntimeError(
                f"No connected road path from {start['place']} to {end['place']}. "
                "The road network has disconnected components or an incomplete corridor. "
                "Even with retain_all=True the reference anchors cannot be joined. "
                "Inspect the actual highway geometry in QGIS, adjust the download "
                "polygon/waypoints, and use a new GraphML cache filename "
                f"(current: {GRAPHML}). Do not force a straight-line connection."
            ) from exc
        for u, v in zip(path[:-1], path[1:]):
            choices = undirected.get_edge_data(u, v)
            data = min(choices.values(), key=lambda item: item["step32_cost"])
            geom = data.get("geometry")
            if geom is None:
                geom = LineString([
                    (undirected.nodes[u]["x"], undirected.nodes[u]["y"]),
                    (undirected.nodes[v]["x"], undirected.nodes[v]["y"]),
                ])
            edge_xy = list(geom.coords)
            node_u = (undirected.nodes[u]["x"], undirected.nodes[u]["y"])
            # Match edge geometry orientation to the node traversal order.
            if ((edge_xy[-1][0] - node_u[0]) ** 2 + (edge_xy[-1][1] - node_u[1]) ** 2
                    < (edge_xy[0][0] - node_u[0]) ** 2 + (edge_xy[0][1] - node_u[1]) ** 2):
                edge_xy.reverse()
            if coords and coords[-1] != edge_xy[0]:
                # A non-trivial gap indicates malformed/incorrect edge geometry.
                gap = gpd.GeoSeries(
                    [LineString([coords[-1], edge_xy[0]])], crs=CRS_GEO
                ).to_crs(CRS_METRIC).length.iloc[0]
                if gap > 100:
                    raise RuntimeError(f"Discontinuous route geometry: {gap:.1f}m gap")
            coords.extend(edge_xy if not coords else edge_xy[1:])
            route_edges.append({"leg": leg, "from_place": start["place"],
                                "to_place": end["place"],
                                "tagged_hwy99": bool(data["step32_is_h99"]),
                                "length_m": float(data.get("length", 0.0)),
                                "geometry": geom})

    line = gpd.GeoDataFrame(
        {"feature": ["candidate_highway99"],
         "status": ["PROVISIONAL_NOT_AUTHORITATIVE"],
         "method": ["OSM drive graph; highway-99 tag preferred; waypoint route"]},
        geometry=[LineString(coords)], crs=CRS_GEO,
    )
    edges = gpd.GeoDataFrame(route_edges, crs=CRS_GEO)
    length_km = float(line.to_crs(CRS_METRIC).length.iloc[0]) / 1_000
    percent_tagged = float(
        edges.loc[edges["tagged_hwy99"], "length_m"].sum()
        / max(edges["length_m"].sum(), 1) * 100
    )
    print(f"Candidate length: {length_km:.1f} km")
    print(f"Road length labelled Highway 99 in OSM: {percent_tagged:.1f}%")
    if percent_tagged < 80:
        print("REVIEW: <80% of routed distance has Highway 99 labels; inspect QGIS.")
    if any(s["snap_method"] != "highway99_tagged" for s in snaps):
        print("REVIEW REQUIRED: untagged-road waypoint fallback(s) were used. "
              "Candidate route is not validated Highway 99 geometry.")
    if not 85 <= length_km <= 175:
        print("REVIEW: route length is outside a broad plausibility band (85-175 km).")
    return line, edges, snaps, length_km, percent_tagged


def make_maps(line: gpd.GeoDataFrame, points: gpd.GeoDataFrame,
              envelope: gpd.GeoDataFrame):
    metric_line = line.to_crs(CRS_METRIC)
    metric_points = points.to_crs(CRS_METRIC)
    metric_envelope = envelope.to_crs(CRS_METRIC)

    fig, ax = plt.subplots(figsize=(8, 10))
    metric_envelope.plot(ax=ax, color="#c5ddea", alpha=0.4, edgecolor="#6b99b8")
    metric_line.plot(ax=ax, linewidth=2.7, color="#c43d28", zorder=3)
    metric_points.plot(ax=ax, color="#172a3a", markersize=30, zorder=4)
    for _, row in metric_points.iterrows():
        ax.annotate(row["name"], (row.geometry.x, row.geometry.y),
                    xytext=(6, 6), textcoords="offset points", fontsize=9)
    ax.set_title("Highway 99 corridor: OSM road-network candidate\n"
                 "Horseshoe Bay to Whistler | Step 3.2", fontsize=12)
    ax.set_xlabel("Easting (m), NAD83 / UTM zone 10N")
    ax.set_ylabel("Northing (m), NAD83 / UTM zone 10N")
    ax.set_aspect("equal")
    ax.grid(alpha=0.2)
    fig.text(0.5, 0.02, "© OpenStreetMap contributors | Candidate, not validated official centerline.\n"
             "Shading: 5-km data-search buffer, NOT a landslide hazard boundary.",
             ha="center", fontsize=8)
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    fig.savefig(OUT_MAPS / "highway99_corridor_candidate.png", dpi=190)
    plt.close(fig)

    # Regional context includes Vancouver as a point ONLY: we do not falsely
    # portray the Vancouver-Horseshoe Bay connector as an extracted road.
    van_lat, van_lon = ox.geocode("Vancouver, British Columbia, Canada")
    vancouver = gpd.GeoDataFrame(
        {"name": ["Vancouver"], "role": ["regional_reference_only"]},
        geometry=gpd.points_from_xy([van_lon], [van_lat]), crs=CRS_GEO,
    )
    regional_points = gpd.GeoDataFrame(
        {"name": ["Vancouver"] + points["name"].tolist(),
         "role": ["regional_reference_only"] + ["corridor_waypoint"] * len(points)},
        geometry=list(vancouver.geometry) + list(points.geometry), crs=CRS_GEO,
    )
    regional_points.to_file(OUT_DATA / "regional_reference_places.geojson", driver="GeoJSON")
    fig, ax = plt.subplots(figsize=(9, 10))
    metric_line.plot(ax=ax, color="#c43d28", linewidth=2.8, zorder=2)
    for _, row in regional_points.to_crs(CRS_METRIC).iterrows():
        ax.scatter(row.geometry.x, row.geometry.y,
                   c="#1675a1" if row["name"] == "Vancouver" else "#172a3a",
                   s=30, zorder=3)
        ax.annotate(row["name"], (row.geometry.x, row.geometry.y),
                    xytext=(6, 5), textcoords="offset points", fontsize=9)
    ax.set_title("Regional context: Vancouver–Whistler\n"
                 "Red line: extracted H99 candidate from Horseshoe Bay only")
    ax.set_xlabel("Easting (m), NAD83 / UTM zone 10N")
    ax.set_ylabel("Northing (m), NAD83 / UTM zone 10N")
    ax.set_aspect("equal")
    ax.grid(alpha=0.2)
    fig.text(0.5, 0.02, "© OpenStreetMap contributors | Vancouver is a reference point, NOT a routed connection.\n"
             "Transport-disruption network will be modeled later.", ha="center", fontsize=8)
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    fig.savefig(OUT_MAPS / "regional_vancouver_whistler.png", dpi=190)
    plt.close(fig)
    return regional_points


def main():
    OUT_DATA.mkdir(parents=True, exist_ok=True)
    OUT_MAPS.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    ox.settings.use_cache = True
    ox.settings.cache_folder = str(CACHE / "requests")
    ox.settings.requests_timeout = 300
    ox.settings.overpass_rate_limit = True
    # Identify the academic data-acquisition workflow to the OSM service.
    ox.settings.http_user_agent = (
        "NatCatRiskResearch/0.1 (+https://github.com/epmresearch/natcat-risk)"
    )
    print("OSMnx Overpass status-based throttling: ON")
    print("Extra minimum interval between tile requests: "
          f"{MIN_SECONDS_BETWEEN_TILE_REQUESTS}s; bounded retries: "
          f"{MAX_ATTEMPTS_PER_TILE} attempts (30s, 60s backoff)")
    # Explicitly preserve road names and route refs in graph edge attributes.
    for tag in ("ref", "name"):
        if tag not in ox.settings.useful_tags_way:
            ox.settings.useful_tags_way.append(tag)

    points = load_places()
    print("Loaded Step 3.1 points:", ", ".join(points["name"]))
    graph = prepare_graph(points)
    print(f"OSM graph: {graph.number_of_nodes():,} nodes; {graph.number_of_edges():,} edges")
    line, edges, snaps, length_km, percent_tagged = build_candidate_route(graph, points)
    envelope = gpd.GeoDataFrame(
        {"purpose": ["DATA_SEARCH_ONLY_NOT_HAZARD"], "width_m": [SEARCH_BUFFER_M]},
        geometry=line.to_crs(CRS_METRIC).buffer(SEARCH_BUFFER_M), crs=CRS_METRIC,
    ).to_crs(CRS_GEO)
    line.to_file(OUT_DATA / "highway99_candidate_centerline.geojson", driver="GeoJSON")
    edges.to_file(OUT_DATA / "highway99_candidate_edges.geojson", driver="GeoJSON")
    envelope.to_file(OUT_DATA / "highway99_data_search_buffer.geojson", driver="GeoJSON")
    region = make_maps(line, points, envelope)

    metadata = {
        "step": "3.2", "status": "PROVISIONAL_REQUIRES_MANUAL_VALIDATION",
        "utc_generated": datetime.now(timezone.utc).isoformat(),
        "source": "OpenStreetMap via OSMnx/Overpass; © OpenStreetMap contributors",
        "osmnx_version": ox.__version__, "input": str(POINTS_PATH.relative_to(REPO)),
        "input_point_names": WAYPOINTS, "analysis_crs": CRS_METRIC,
        "network_download_method": "incremental overlapping OSM tiles",
        "network_tile_axis_max_m": MAX_TILE_AXIS_M,
        "network_tile_buffer_m": TILE_BUFFER_M,
        "network_tile_count": len(build_download_tiles(points)),
        "network_cached_graph": str(GRAPHML.relative_to(REPO)),
        "search_buffer_m": SEARCH_BUFFER_M,
        "route_method": "Undirected shortest path, preferred Highway 99 OSM tags",
        "off_highway_cost_multiplier": OFF_HIGHWAY_COST_MULTIPLIER,
        "candidate_length_km": round(length_km, 2),
        "candidate_hwy99_tagged_length_percent": round(percent_tagged, 1),
        "snapped_waypoints": snaps,
        "all_anchors_tagged_highway99": all(
            s["snap_method"] == "highway99_tagged" for s in snaps
        ),
        "regional_reference_places": region["name"].tolist(),
        "limitations": [
            "Not an authoritative government-surveyed centerline.",
            "Can include connecting roads when highway refs are missing or ambiguous.",
            "Some waypoints may use unverified local-road anchors; inspect QA JSON.",
            "Undirected geometry is not a traffic/legal-access routing model.",
            "Vancouver is a geographic reference point, not a routed segment.",
            "5-km buffer is for data discovery, not landslide hazard/runout.",
            "Small OSM download tiles are operational, not physical landslide hazard buffers.",
        ],
    }
    qa_path = OUT_DATA / "highway99_candidate_quality.json"
    qa_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    print("\nOUTPUTS:")
    for path in [OUT_DATA / "highway99_candidate_centerline.geojson",
                 OUT_DATA / "highway99_candidate_edges.geojson",
                 OUT_DATA / "highway99_data_search_buffer.geojson",
                 OUT_DATA / "regional_reference_places.geojson",
                 qa_path, OUT_MAPS / "highway99_corridor_candidate.png",
                 OUT_MAPS / "regional_vancouver_whistler.png"]:
        print("  ", path)
    print("\nIMPORTANT: Review the route in QGIS against a basemap before Step 3.3.")


if __name__ == "__main__":
    main()
