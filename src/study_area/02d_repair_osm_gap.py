"""Step 3.2: repair a diagnosed OSM tile seam using a small REAL road download.

Usage from repository root (Miniforge environment active):
  python src/study_area/02d_repair_osm_gap.py

Reads the existing combined GraphML and connectivity diagnostic JSON. Downloads
ONE local OSMnx drive-network patch around the nearest Squamish/Whistler component
nodes, merges by ORIGINAL OSM node IDs, and verifies real graph connectivity.
Never creates a synthetic edge or straight-line roadway. No input caches deleted.
A successful repair reuses the original 02_extract_highway99_candidate.py to
write its normal provisional GIS outputs, after temporarily pointing it at the
separately saved repaired cache. The original script itself is not modified.
"""
from __future__ import annotations

from datetime import datetime, timezone
import importlib.util
import json
from pathlib import Path
import time

import networkx as nx
from pyproj import Transformer
from requests.exceptions import ConnectionError as RequestsConnectionError, Timeout as RequestsTimeout
from shapely.geometry import LineString
from shapely.ops import transform

REPO = Path(__file__).resolve().parents[2]
CACHE = REPO / "data/cache/study_area"
BASE = CACHE / "highway99_tiles_16km_9km_v3_combined.graphml"
PATCH = CACHE / "highway99_osm_targeted_gap_patch_4km.graphml"
REPAIRED = CACHE / "highway99_tiles_16km_9km_v3_gap_repaired.graphml"
DIAG = REPO / "results/diagnostics/step32_connectivity_report.json"
REPORT = REPO / "results/diagnostics/step32_gap_repair_report.json"
ROUTE_SCRIPT = REPO / "src/study_area/02_extract_highway99_candidate.py"
PATCH_RADIUS_M = 4000
RETRY_DELAYS_S = (30, 60)


def derive_gap_polygon(graph: nx.MultiDiGraph, node_a: int, node_b: int,
                       buffer_m: float = PATCH_RADIUS_M):
    """WGS84 query polygon buffered around measured gap, without fabricated roads."""
    if buffer_m <= 0:
        raise ValueError("Patch buffer must be positive")
    for n in (node_a, node_b):
        if n not in graph:
            raise KeyError(f"Diagnostic OSM node {n} absent from base graph")
    coords = [
        (float(graph.nodes[n]["x"]), float(graph.nodes[n]["y"]))
        for n in (node_a, node_b)
    ]
    forward = Transformer.from_crs("EPSG:4326", "EPSG:26910", always_xy=True).transform
    backward = Transformer.from_crs("EPSG:26910", "EPSG:4326", always_xy=True).transform
    metric_line = transform(forward, LineString(coords))
    return transform(backward, metric_line.buffer(buffer_m))


def merge_original_osm_nodes(base: nx.MultiDiGraph, patch: nx.MultiDiGraph):
    """Combine networks by original IDs only; reject conflicting node positions."""
    shared = set(base.nodes) & set(patch.nodes)
    if not shared:
        raise RuntimeError("OSM patch contains no shared node IDs with original graph")
    # Detect any inconsistent node IDs before composing graphs.
    for n in shared:
        dx = abs(float(base.nodes[n]["x"]) - float(patch.nodes[n]["x"]))
        dy = abs(float(base.nodes[n]["y"]) - float(patch.nodes[n]["y"]))
        if dx > 0.00002 or dy > 0.00002:
            raise RuntimeError(f"Shared OSM node {n} has inconsistent coordinates")
    combined = nx.compose(base, patch)
    return combined, len(shared)


def anchors_connected(graph, source: int, target: int) -> bool:
    """Geometry connectivity, not one-way vehicle routing validity."""
    if source not in graph or target not in graph:
        return False
    return nx.has_path(graph.to_undirected(as_view=True), source, target)


def is_transient(exc: Exception) -> bool:
    if isinstance(exc, (RequestsConnectionError, RequestsTimeout)):
        return True
    if type(exc).__name__ != "ResponseStatusCodeError":
        return False
    message = str(exc).lower()
    return any(t in message for t in (
        "429", "502", "503", "504", "timeout", "gateway", "capacity", "too many requests"
    ))


def fetch_patch(polygon, patch_file: Path):
    """One small patch, no repeated downloads if cached, bounded retry."""
    import osmnx as ox
    if patch_file.exists():
        print(f"Reusing cached local OSM patch: {patch_file}", flush=True)
        return ox.io.load_graphml(patch_file)
    for attempt in range(1, 4):
        print(f"Downloading local gap patch (attempt {attempt}/3)...", flush=True)
        try:
            result = ox.graph.graph_from_polygon(
                polygon, network_type="drive", simplify=False,
                retain_all=True, truncate_by_edge=True,
            )
            if result.number_of_nodes() == 0 or result.number_of_edges() == 0:
                raise RuntimeError("Patch query returned an empty drivable graph")
            ox.io.save_graphml(result, patch_file)
            print(f"Cached patch: {result.number_of_nodes():,} nodes, "
                  f"{result.number_of_edges():,} edges", flush=True)
            return result
        except Exception as exc:
            if not is_transient(exc) or attempt == 3:
                raise RuntimeError(f"Local OSM patch request failed: {exc}") from exc
            delay = RETRY_DELAYS_S[attempt - 1]
            print(f"Transient Overpass error; waiting {delay}s before retry.", flush=True)
            time.sleep(delay)
    raise AssertionError("Unreachable retry state")


def run_original_route_script(repaired_graph: Path):
    """Reuse existing validated output workflow, without rewriting main script."""
    spec = importlib.util.spec_from_file_location("step32_original_route", ROUTE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load original route script: {ROUTE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.GRAPHML = repaired_graph
    print("Running original Step 3.2 routing against repaired graph...", flush=True)
    module.main()


def main():
    import osmnx as ox
    for required in (BASE, DIAG, ROUTE_SCRIPT):
        if not required.is_file():
            raise SystemExit(f"Missing input: {required}")
    CACHE.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    ox.settings.use_cache = True
    ox.settings.cache_folder = str(CACHE / "requests")
    ox.settings.requests_timeout = 300
    ox.settings.overpass_rate_limit = True
    ox.settings.http_user_agent = (
        "NatCatRiskResearch/0.1 (+https://github.com/epmresearch/natcat-risk)"
    )
    for tag in ("ref", "name"):
        if tag not in ox.settings.useful_tags_way:
            ox.settings.useful_tags_way.append(tag)

    info = json.loads(DIAG.read_text(encoding="utf-8"))
    pair = info["nearest_pair_between_anchor_components"]
    first, second = int(pair["first_node"]), int(pair["second_node"])
    anchors = {row["name"]: int(row["node"]) for row in info["anchors"]}
    source, target = anchors["Squamish"], anchors["Whistler"]
    print(f"Loading the previously downloaded base graph: {BASE}", flush=True)
    base = ox.io.load_graphml(BASE)
    if anchors_connected(base, source, target):
        raise SystemExit("Source and target are already connected in base graph: no patch needed.")
    polygon = derive_gap_polygon(base, first, second, PATCH_RADIUS_M)
    print(f"Diagnosed nearest component gap: {pair['distance_m']} m", flush=True)
    print(f"Downloading only a {PATCH_RADIUS_M/1000:.0f}-km buffer "
          "around the two diagnostic nodes, not the entire corridor.", flush=True)
    patch = fetch_patch(polygon, PATCH)
    combined, shared = merge_original_osm_nodes(base, patch)
    connected = anchors_connected(combined, source, target)
    before_count = nx.number_connected_components(base.to_undirected(as_view=True))
    after_count = nx.number_connected_components(combined.to_undirected(as_view=True))
    print(f"Merged by {shared:,} SHARED real OSM node IDs", flush=True)
    print(f"Connected components: {before_count} -> {after_count}", flush=True)
    print(f"Squamish–Whistler graph connectivity: {'CONNECTED' if connected else 'STILL DISCONNECTED'}")
    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "type": "REAL_OSM_NETWORK_GAP_REPAIR_DIAGNOSTIC",
        "input_graph": str(BASE.relative_to(REPO)),
        "gap_distance_m_reported": pair["distance_m"],
        "gap_node_ids": [first, second],
        "patch_buffer_m": PATCH_RADIUS_M,
        "patch_cache": str(PATCH.relative_to(REPO)),
        "base_nodes": base.number_of_nodes(),
        "patch_nodes": patch.number_of_nodes(),
        "shared_osm_node_ids": shared,
        "merged_nodes": combined.number_of_nodes(),
        "component_count_before": before_count,
        "component_count_after": after_count,
        "anchor_names": ["Squamish", "Whistler"],
        "anchor_connected_undirected": connected,
        "geometry_rule": "Merge only original matching OSM node IDs; no synthetic edges",
        "status": "CONNECTED_PROVISIONAL" if connected else "DISCONNECTED_NEEDS_REVIEW",
    }
    REPORT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"QA report: {REPORT}")
    if not connected:
        print("DO NOT route or invent a straight connection. The local patch did "
              "not resolve connectivity; inspect the gap in QGIS or consult BC Digital Road Atlas.")
        return
    # Only write a separate merged file if the actual graph is connected.
    ox.io.save_graphml(combined, REPAIRED)
    print(f"Saved connected provisional graph: {REPAIRED}")
    try:
        run_original_route_script(REPAIRED)
    except Exception as exc:
        print("Repaired cache saved but full candidate extraction was NOT successful.")
        print(f"Route generation error: {type(exc).__name__}: {exc}")
        raise
    print("Route output generated; manual QGIS/official-map validation still required.")


if __name__ == "__main__":
    main()
