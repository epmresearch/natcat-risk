"""Create a provisional study-area map for the Sea-to-Sky corridor.

This first-pass map connects geocoded reference places with straight segments.
It DOES NOT represent the surveyed centerline of Highway 99. The accurate
road geometry will be acquired and validated in Step 3.2.

Run from the repository root:
    python src/study_area/01_create_preliminary_study_area.py
"""

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import osmnx as ox
from shapely.geometry import LineString

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_OUT = REPO_ROOT / "data" / "processed" / "study_area"
FIG_OUT = REPO_ROOT / "results" / "maps"
SEARCH_BUFFER_METERS = 10_000  # Preliminary search extent, not a landslide hazard zone.
CRS_GEOGRAPHIC = "EPSG:4326"
CRS_METRIC = "EPSG:26910"  # NAD83 / UTM zone 10N: suitable for the BC South Coast.

# Place names are geocoded when the script runs; no coordinates are hardcoded.
PLACES = [
    ("Horseshoe Bay", "Horseshoe Bay, West Vancouver, British Columbia, Canada"),
    ("Lions Bay", "Lions Bay, British Columbia, Canada"),
    ("Britannia Beach", "Britannia Beach, British Columbia, Canada"),
    ("Squamish", "Squamish, British Columbia, Canada"),
    ("Whistler", "Whistler, British Columbia, Canada"),
]


def main() -> None:
    DATA_OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)
    ox.settings.use_cache = True
    ox.settings.cache_folder = str(REPO_ROOT / "data" / "cache" / "osmnx")

    names, latitudes, longitudes = [], [], []
    print("Geocoding reference places using OpenStreetMap/Nominatim...")
    for name, query in PLACES:
        try:
            lat, lon = ox.geocode(query)  # osmnx.geocode returns (latitude, longitude).
        except Exception as exc:
            raise RuntimeError(f"Could not geocode {name}: {exc}") from exc
        print(f"  {name}: latitude={lat:.6f}, longitude={lon:.6f}")
        names.append(name)
        latitudes.append(lat)
        longitudes.append(lon)

    points = gpd.GeoDataFrame(
        {"name": names, "place_query": [query for _, query in PLACES]},
        geometry=gpd.points_from_xy(longitudes, latitudes),
        crs=CRS_GEOGRAPHIC,
    )
    approximate_axis = gpd.GeoDataFrame(
        {"description": ["Schematic waypoint connections - NOT Highway 99 geometry"]},
        geometry=[LineString(list(zip(longitudes, latitudes)))],
        crs=CRS_GEOGRAPHIC,
    )

    axis_metric = approximate_axis.to_crs(CRS_METRIC)
    draft_search_area_metric = gpd.GeoDataFrame(
        {"purpose": ["Preliminary data-search boundary only"], "buffer_m": [SEARCH_BUFFER_METERS]},
        geometry=axis_metric.geometry.buffer(SEARCH_BUFFER_METERS),
        crs=CRS_METRIC,
    )
    draft_search_area = draft_search_area_metric.to_crs(CRS_GEOGRAPHIC)

    points.to_file(DATA_OUT / "reference_places.geojson", driver="GeoJSON")
    approximate_axis.to_file(DATA_OUT / "schematic_waypoint_axis.geojson", driver="GeoJSON")
    draft_search_area.to_file(DATA_OUT / "preliminary_search_area.geojson", driver="GeoJSON")

    fig, ax = plt.subplots(figsize=(8, 11))
    draft_search_area_metric.plot(ax=ax, facecolor="lightsteelblue", edgecolor="steelblue", alpha=0.4)
    axis_metric.plot(ax=ax, color="darkorange", linewidth=2, linestyle="--")
    points_metric = points.to_crs(CRS_METRIC)
    points_metric.plot(ax=ax, color="black", markersize=28, zorder=3)
    for _, row in points_metric.iterrows():
        ax.annotate(row["name"], (row.geometry.x, row.geometry.y),
                    xytext=(7, 5), textcoords="offset points", fontsize=9)

    ax.set_title("Sea-to-Sky Study Area — Preliminary (Step 3.1)\n"
                 "Dashed line is schematic, NOT the actual Highway 99 alignment", fontsize=11)
    ax.set_xlabel("Easting (metres, NAD83 / UTM zone 10N)")
    ax.set_ylabel("Northing (metres, NAD83 / UTM zone 10N)")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.25)
    fig.text(0.5, 0.02,
             "Source: OpenStreetMap contributors (Nominatim). "
             "10-km buffer is a draft data-search boundary, not a hazard zone.",
             ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    figure_path = FIG_OUT / "study_area_preliminary.png"
    fig.savefig(figure_path, dpi=180)
    plt.close(fig)

    print("\nCreated:")
    print(f"  {DATA_OUT / 'reference_places.geojson'}")
    print(f"  {DATA_OUT / 'schematic_waypoint_axis.geojson'}")
    print(f"  {DATA_OUT / 'preliminary_search_area.geojson'}")
    print(f"  {figure_path}")
    print("\nIMPORTANT: These are provisional reference layers, not the highway centerline.")


if __name__ == "__main__":
    main()
