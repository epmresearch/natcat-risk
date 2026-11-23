# Step 3.3 — Highway 99 chainage segmentation (provisional)

## Objective

Split the single **provisional** Highway 99 centerline produced in Step 3.2 into reproducible, equal-chainage road sections. The study corridor is **Horseshoe Bay → Whistler**; Vancouver is regional context, **not** part of the extracted centerline.

This is a GIS analysis-unit preparation step. It is **not** an estimate of slope instability, landslide runout, road width, roadway direction, travel time, or closure probability.

## Inputs and provenance

- `data/processed/study_area/highway99_candidate_centerline.geojson` (one ordered LineString)
- `data/processed/study_area/highway99_candidate_quality.json` (Step 3.2 QA; e.g. provisional 100.44 km)
- Original source: OpenStreetMap via OSMnx / Overpass, © OpenStreetMap contributors. Its accuracy is **not** equivalent to an official BC road centerline.

No graph cache, internet download, or routing calculation is required. The script does not use the Overpass API.

## Segmentation definition

- Default target: **1,000 metres** per segment; this is a provisional modelling resolution rather than a geotechnical finding. It may be revisited in sensitivity analysis (e.g., 500 m and 2,000 m) when hazard data resolution is known.
- Reproject to **EPSG:26910** (NAD83 / UTM zone 10N) for metric chainage.
- Use `shapely.ops.substring` on the **source polyline**, without simplification and without creating straight-line links between segments.
- Consecutive segments are oriented along the source line, from Horseshoe Bay to Whistler.
- Start and end chainages are measured in metres from the start of the source geometry; the last segment may be shorter than 1,000 m.
- Feature IDs use `H99-0001`, `H99-0002`, ... and are reproducible for the same input, projection and chosen length.
- GeoJSON is exported in **EPSG:4326** (as with Step 3.2); QGIS can reproject on the fly. Do not measure segment lengths from longitude/latitude degrees.

## Files generated

With 1,000 m default:

1. `data/processed/study_area/highway99_segments_1000m.geojson` — line geometry + ID, chainage, length, midpoint, provisional status.
2. `data/processed/study_area/highway99_segments_1000m.csv` — attribute table without geometry.
3. `data/processed/study_area/highway99_segments_1000m_quality.json` — provenance, feature count, topology/length checks, cautions.
4. `results/maps/highway99_segments_1000m.png` — overview plot using a metric CRS.

## Run from the repository root

```bat
python -m unittest discover -s tests -p "test_study_area_step33.py"
python src/study_area/03_segment_highway99.py
```

For a different chainage length (do **not** choose a final resolution before review):

```bat
python src/study_area/03_segment_highway99.py --segment-length-m 500
```

Run from `(sea_to_sky_risk)`; requires installed `geopandas`, `shapely`, `matplotlib`, `pandas`, `pyproj`. Expected default count is roughly **101 features** for the Step 3.2 length of approximately 100.44 km. **Use the actual generated QA report** to confirm the exact count/last length.

## Automated quality controls

- Exactly one ordered input LineString with a declared CRS.
- Input line and QA report length agree within 100 m, otherwise fail.
- Segment geometries are valid, ordered and positive-length.
- Segment boundary gaps ≤0.01 m before export; ≤0.02 m after GeoJSON round trip.
- Cumulative segment length differs from input by ≤0.10 m before export, ≤0.25 m after export.
- Segment identifiers are unique; the first and last coordinates match the source geometry.
- No gap-filling or synthetic straight-line geometry is introduced.

## Manual QGIS validation (required before Step 3.3 sign-off)

1. Open the existing Step 3.2 centerline, OSM XYZ basemap, and the new `highway99_segments_1000m.geojson`.
2. Verify the new segments sit **exactly** on the original line from Horseshoe Bay to Whistler.
3. Check a few boundaries at high zoom near Lions Bay, Squamish and the previously repaired graph gap (approximately 50.04°N, −123.12°W).
4. Review `segment_id`, `start_km`, `end_km`, `length_m` in the attribute table; last segment should be shorter than target if the length is not an exact multiple.
5. Verify there are no gaps, overlapping slivers, duplicated IDs, unintended jumps or reversed chainage.
6. Capture a QGIS screenshot of the map and send the QA JSON and terminal output for review.

## Research limitations and future validation

The segmentation resolution is an **analysis choice**, not a validated landslide footprint size. A 1-km division can obscure highly localized slope susceptibility; hazard resolution, vulnerability definitions and possible alternative lengths must be considered in later stages. Keep the original OSM-derived input and its provisional status. The segmented centerline **does not constitute a traffic routing graph** or automatically encode physical road asset valuation.

## Version control

`data/processed/` is intentionally ignored in this project. Do not `git add -f` the new GeoJSON, CSV, or QA until their contents and sizes have been reviewed. After validation, explicitly stage only vetted lightweight outputs, script, tests, docs and figure. Do not upload large GraphML/cache files or `git add -A` indiscriminately.
