"""Diagnose disconnected components in the cached Highway 99 OSM graph.

Run (no network requests are made):
 python src/study_area/02c_diagnose_highway99_connectivity.py \
   --from-node 481392279 --to-node 3236309869

Inputs: existing combined OSMnx GraphML and Step 3.1 reference_places GeoJSON.
Diagnostic straight segments on the plot represent distances between nodes,
NOT a road, proposed connection, or safe driving route.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import networkx as nx
from pyproj import Transformer
from scipy.spatial import cKDTree

REPO = Path(__file__).resolve().parents[2]
DEFAULT_GRAPH = (REPO / "data/cache/study_area/"
                 "highway99_tiles_16km_9km_v3_combined.graphml")
OUT_DIR = REPO / "results/diagnostics"
METRIC = "EPSG:26910"


def component_index(graph):
    """Return largest-first components and node->component ID mapping."""
    components = sorted(nx.connected_components(graph.to_undirected(as_view=True)),
                        key=len, reverse=True)
    node_to_component = {node: i + 1 for i, nodes in enumerate(components)
                         for node in nodes}
    return components, node_to_component


def projected_coordinates(graph, nodes, transformer):
    """Return projected coordinates as a Nx2 array for a node collection."""
    nodes = list(nodes)
    lon = [float(graph.nodes[n]["x"]) for n in nodes]
    lat = [float(graph.nodes[n]["y"]) for n in nodes]
    east, north = transformer.transform(lon, lat)
    return np.column_stack((east, north))


def closest_distinct_component_nodes(graph, nodes_a, nodes_b, transformer):
    """Find closest nodes across two components; never add an edge."""
    nodes_a = list(nodes_a)
    nodes_b = list(nodes_b)
    if not nodes_a or not nodes_b:
        raise ValueError("Cannot compare empty components")
    xy_a = projected_coordinates(graph, nodes_a, transformer)
    xy_b = projected_coordinates(graph, nodes_b, transformer)
    distances, indices = cKDTree(xy_b).query(xy_a, k=1)
    idx_a = int(np.argmin(distances))
    idx_b = int(indices[idx_a])
    return {
        "distance_m": round(float(distances[idx_a]), 2),
        "first_node": int(nodes_a[idx_a]),
        "second_node": int(nodes_b[idx_b]),
        "first_lonlat": [float(graph.nodes[nodes_a[idx_a]]["x"]),
                         float(graph.nodes[nodes_a[idx_a]]["y"])],
        "second_lonlat": [float(graph.nodes[nodes_b[idx_b]]["x"]),
                          float(graph.nodes[nodes_b[idx_b]]["y"])],
    }


def waypoint_component_options(graph, components, transformer, points_path,
                               names=("Squamish", "Whistler"), radius_m=2500):
    """Find graph components with nodes within radius of town reference points.

    This detects alternatives to the nearest-edge snapping used by the route
    script, but does not select an anchor automatically.
    """
    if not points_path.exists():
        return {"warning": f"Reference places missing: {points_path}"}
    import geopandas as gpd
    points = gpd.read_file(points_path).to_crs("EPSG:4326")
    if "name" not in points:
        return {"warning": "No 'name' column in reference places"}
    results = {}
    for name in names:
        rows = points.loc[points["name"] == name]
        if rows.empty:
            results[name] = {"warning": "No matching reference point"}
            continue
        pt = rows.geometry.iloc[0]
        query_xy = np.array(transformer.transform(float(pt.x), float(pt.y)))
        near = []
        for comp_id, nodes in enumerate(components, 1):
            xy = projected_coordinates(graph, nodes, transformer)
            dist, idx = cKDTree(xy).query(query_xy)
            if dist <= radius_m:
                node = list(nodes)[int(idx)]
                near.append({"component_id": comp_id, "component_size": len(nodes),
                             "nearest_node": int(node),
                             "nearest_node_distance_m": round(float(dist), 1)})
        results[name] = sorted(near, key=lambda r: r["nearest_node_distance_m"])
    if isinstance(results.get("Squamish"), list) and isinstance(results.get("Whistler"), list):
        s = {x["component_id"] for x in results["Squamish"]}
        w = {x["component_id"] for x in results["Whistler"]}
        results["common_components_within_2_5km"] = sorted(s & w)
    return results


def write_diagnostic_map(graph, components, src, dst, closest, output_path, transformer):
    """Save a plot of node coverage and disconnect. No route is drawn."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 10))
    # Draw context as graph nodes only; no unverified connecting line as roadway.
    all_nodes = list(graph.nodes)
    all_xy = projected_coordinates(graph, all_nodes, transformer)
    ax.scatter(all_xy[:, 0], all_xy[:, 1], s=0.3, c="lightgrey", label="OSM graph nodes")
    for component_id, color, label in [(src, "#c0352b", "Squamish anchor component"),
                                       (dst, "#1769a8", "Whistler anchor component")]:
        xy = projected_coordinates(graph, components[component_id - 1], transformer)
        ax.scatter(xy[:, 0], xy[:, 1], s=0.6, c=color, label=label)
    first = np.array(transformer.transform(*closest["first_lonlat"]))
    second = np.array(transformer.transform(*closest["second_lonlat"]))
    ax.scatter([first[0], second[0]], [first[1], second[1]],
               c=["#c0352b", "#1769a8"], s=45, edgecolors="black", zorder=5)
    # Dotted gap annotation is explicitly *not* a road or a proposed bridge.
    ax.plot([first[0], second[0]], [first[1], second[1]], "k--", alpha=0.8,
            label="Nearest node gap (DIAGNOSTIC ONLY)")
    ax.set_aspect("equal")
    ax.set_title("Step 3.2: OSM connectivity diagnostic (NOT a route)")
    ax.set_xlabel("Easting (m), EPSG:26910")
    ax.set_ylabel("Northing (m), EPSG:26910")
    ax.legend(loc="best", fontsize=8)
    fig.text(0.5, 0.02, "© OpenStreetMap contributors. Dotted line shows only diagnostic separation, NOT a road.",
             ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", type=Path, default=DEFAULT_GRAPH)
    parser.add_argument("--from-node", type=int, default=481392279,
                        help="Squamish snapped OSM node from your latest error")
    parser.add_argument("--to-node", type=int, default=3236309869,
                        help="Whistler snapped OSM node from your latest error")
    args = parser.parse_args()
    if not args.graph.is_file():
        raise SystemExit(f"Cached graph is missing: {args.graph}. No download attempted.")

    import osmnx as ox
    print("Loading existing GraphML (OFFLINE):", args.graph, flush=True)
    graph = ox.io.load_graphml(args.graph)
    for node in (args.from_node, args.to_node):
        if node not in graph:
            raise SystemExit(f"Node ID {node} not found in this cached graph; do not infer a different node.")
    components, by_node = component_index(graph)
    src = by_node[args.from_node]
    dst = by_node[args.to_node]
    print(f"Graph: {len(graph):,} nodes, {graph.number_of_edges():,} edges, "
          f"{len(components)} weakly connected components")
    print(f"Squamish node {args.from_node}: component #{src}, size {len(components[src-1]):,}")
    print(f"Whistler node {args.to_node}: component #{dst}, size {len(components[dst-1]):,}")
    if src == dst:
        print("The two nodes share a component: inspect route graph/edge costs rather than downloading more.")
        return

    t = Transformer.from_crs("EPSG:4326", METRIC, always_xy=True)
    closest = closest_distinct_component_nodes(graph, components[src-1], components[dst-1], t)
    print("Nearest pair between these two components:")
    print(json.dumps(closest, indent=2))

    comp_ids = sorted(set(list(range(1, min(len(components), 8)+1)) + [src, dst]))
    candidates = []
    for idx, a in enumerate(comp_ids):
        for b in comp_ids[idx+1:]:
            result = closest_distinct_component_nodes(graph, components[a-1], components[b-1], t)
            candidates.append({"components": [a, b], **result})
    candidates.sort(key=lambda x: x["distance_m"])
    print("\nClosest gaps among the eight largest/anchor components (diagnostic only):")
    for result in candidates[:12]:
        print(f"  #{result['components'][0]} ↔ #{result['components'][1]}: "
              f"{result['distance_m']:,.0f}m near {result['first_lonlat']} / {result['second_lonlat']}")

    options = waypoint_component_options(
        graph, components, t,
        REPO / "data/processed/study_area/reference_places.geojson")
    print("\nComponents containing nodes within 2.5 km of Squamish/Whistler:")
    print(json.dumps(options, indent=2, ensure_ascii=False))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "type": "OFFLINE_CONNECTIVITY_DIAGNOSTIC_NOT_AN_EXTRACTED_ROUTE",
        "input_graph": str(args.graph),
        "nodes": graph.number_of_nodes(), "edges": graph.number_of_edges(),
        "component_count": len(components),
        "largest_component_sizes": [len(x) for x in components[:10]],
        "anchors": [{"name": "Squamish", "node": args.from_node,
                     "component": src},
                    {"name": "Whistler", "node": args.to_node,
                     "component": dst}],
        "nearest_pair_between_anchor_components": closest,
        "other_nearest_component_pairs": candidates[:20],
        "waypoint_anchor_component_options": options,
        "warning": "Geometrically close nodes in different components must NOT be joined without source validation."
    }
    out_json = OUT_DIR / "step32_connectivity_report.json"
    out_png = OUT_DIR / "step32_connectivity_map.png"
    out_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    write_diagnostic_map(graph, components, src, dst, closest, out_png, t)
    print("\nReport:", out_json)
    print("Map:", out_png)
    print("No download and no modification of the cached network were performed.")


if __name__ == "__main__":
    main()
