# Step 3.4 — Spatial screening corridors & GIS preparation (provisional)

## Research objective and precise meaning

Prepare geospatial **data-screening** regions around the **101 chainage segments** from Step 3.3. They will help locate and clip later digital elevation models (DEM), terrain derivatives, landslide inventories, geological layers, and other datasets.

**Do not interpret any buffer as an empirically validated landslide initiation zone, susceptibility zone, travel/runout distance, road impact zone, or closure boundary.** A landslide may begin well outside a screening buffer and still reach Highway 99. No hazard or risk score is generated at this stage.

## Inputs

- `data/processed/study_area/highway99_segments_1000m.geojson` — 101 one-kilometre-chainage segments; last segment shorter.
- `data/processed/study_area/highway99_segments_1000m_quality.json` — Step 3.3 QA.
- Data lineage: provisional OSM-derived Highway 99 centerline, Horseshoe Bay → Whistler; © OpenStreetMap contributors.

The script uses **no network calls** and does **not** need OSM GraphML cache files. These buffers follow the existing road geometry without resampling or modifying it.

## Coordinate systems and default buffer design

- GIS distance/area calculations: **EPSG:26910**, NAD83 / UTM zone 10N; map and QA quantities in metres and square kilometres.
- Portable GeoJSON output: **EPSG:4326**; QGIS can display layers using its project CRS, including EPSG:3857.
- **1,000 m per-segment near-road screening polygons:** one polygon for each segment; overlapping by design at boundaries and meanders.
- **1,000 m dissolved near-road corridor:** single dissolved union; count its area **once** for spatial coverage.
- **5,000 m terrain/data-discovery corridor:** dissolved buffer for searching DEM and topographic/hazard data; consistent with provisional Step 3.2 5-km data search width. This is **not** a physically justified runout radius.

These widths are **methodological placeholders**. Change them after considering source-data resolution and where initiation sites exist. Later test parameter sensitivity, and do not clip a hazard/runout assessment solely by these widths.

## Run from the repository root in the existing Miniforge environment

```bat
python -m unittest discover -s tests -p "test_study_area_step34.py"
python src/study_area/04_prepare_spatial_corridors.py
```

Optional *exploratory* alternative widths, for example 500 m near-road and 5 km search:

```bat
python src/study_area/04_prepare_spatial_corridors.py --near-m 500 --search-m 5000
```

Note: the dissolved near-corridor filename changes with width, but the shared `highway99_spatial_preparation_quality.json` and map filenames will be overwritten; if exploring alternatives, copy outputs to a different folder using `--out-data` and `--out-maps`.

## Outputs

| File | Content |
|---|---|
| `data/processed/study_area/highway99_segment_screening_zones_1000m.geojson` | 101 overlapping screening polygons, IDs/chainage, individual `zone_area_m2` |
| `data/processed/study_area/highway99_corridor_screening_1000m.geojson` | Single **dissolved** unique near-road screening area |
| `data/processed/study_area/highway99_terrain_search_5000m.geojson` | Single **dissolved** 5-km terrain data-search footprint |
| `data/processed/study_area/highway99_spatial_preparation_quality.json` | Source file hash, count, areas, overlap accounting, limitations |
| `results/maps/highway99_spatial_screening_corridors.png` | Schematic GIS overview, not a hazard map |

## Critical QA criteria

1. Count and IDs match Step 3.3 (expect 101; `H99-0001`...`H99-0101`).
2. Metric geometry lengths and chainage agree with Step 3.3 QA; no segment geometry changed.
3. GeoJSON outputs re-open correctly, retain valid polygons, and survive round-trip reprojected area tolerance.
4. The 1-km dissolved corridor is fully contained in the 5-km search corridor (numerical tolerance 1 square metre).
5. Sum of 101 per-segment polygon areas is *greater* than their dissolved union because segments overlap; **never double-count these footprints** when reporting unique study area or exposures.
6. No hazard probability, susceptibility, damage, displacement, or traffic closure is inferred.

### QGIS manual review

Open three new GeoJSON layers plus the original `highway99_segments_1000m.geojson` and OpenStreetMap base layer. Suggested order from top to bottom: road segments, 1-km segment polygons, 1-km dissolved corridor (temporarily disable if obscuring the per-segment polygons), 5-km terrain corridor, OSM base. Use transparent polygon fills. Compare the northern Squamish–Whistler transition and southern endpoint, check polygons contain road lines and there are no unexpected holes or misplaced features. Open the 101-polygon layer's attribute table to verify IDs and row count. Save screenshot and QA JSON for review.

## Risk-model integration caveats

- Later, a landslide inventory/hazard may intersect multiple neighbouring per-segment screening zones. Assign hazards by a declared rule (e.g. intersection with actual road footprint or nearest affected road point, then deduplicate); buffering alone does not establish an impact.
- Near-road polygon area does not equal paved-road asset area. Do not value it as a road.
- Direction of travel, network connectivity, bypasses and transport losses require a separate **directed network model**.
- Hazard generation, slope units, DEM processing, susceptibility, runout modelling and CLIMADA exposure/vulnerability are **later stages**.

## Version control

After successful tests and human GIS QA, commit the script, tests, docs, map, small GeoJSONs and QA JSON explicitly. `data/processed/` is ignored, so individually `git add -f` only reviewed data outputs. Do not add caches, large rasters, downloads, or the entire ignored directory. Avoid `git add -A` indiscriminately.
