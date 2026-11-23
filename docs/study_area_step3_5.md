# Step 3.5 — Provisional Highway 99 Study Area: Reproducibility and GIS QA

## Purpose and scope

This stage closes the **internal reproducibility and consistency review** of the
study-area artifacts made during Steps 3.2–3.4. It does **not** independently
validate the alignment against an authoritative road dataset, and it does not
produce landslide hazard, susceptibility, runout, road-closure, or CLIMADA risk.
The earlier Step 3.1 schematic search area is historical context, not an input
that supersedes the Step 3.2 centerline.

**Analysis region:** Highway 99, Horseshoe Bay to Whistler, British Columbia,
Canada. Vancouver remains a regional reference, not part of the 100.44 km
centerline. The centerline is derived from OpenStreetMap (© OpenStreetMap
contributors) and is still explicitly provisional.

## Prerequisites

- Activate the existing `sea_to_sky_risk` Miniforge/mamba environment (Python
  3.11, geopandas, shapely, pyproj, matplotlib, pandas/pyogrio available).
- Run from the root of the `natcat-risk` repository.
- Keep all previously committed Step 3.2/3.3/3.4 artifacts unchanged. **Do not
  regenerate or re-export their GeoJSON files merely to run this QA script:**
  the SHA-256 provenance comparisons require the original exported bytes.

The following nine files must already exist under `data/processed/study_area/`:

1. `highway99_candidate_centerline.geojson`
2. `highway99_candidate_quality.json`
3. `highway99_segments_1000m.geojson`
4. `highway99_segments_1000m.csv`
5. `highway99_segments_1000m_quality.json`
6. `highway99_segment_screening_zones_1000m.geojson`
7. `highway99_corridor_screening_1000m.geojson`
8. `highway99_terrain_search_5000m.geojson`
9. `highway99_spatial_preparation_quality.json`

No GraphML cache, live Overpass server, DEM, or network request is needed.

## Execute

From **Miniforge Prompt** in the project root:

```bat
python -m unittest discover -s tests -p "test_study_area_step35.py"
python src/study_area/05_validate_study_area.py
```

Expected tests: **15 tests, OK**. A library `DeprecationWarning` alone is not a
failed test.

## What is checked

1. **Input lineage:** Step 3.3 `source_sha256` matches the exact Step 3.2
   centerline file; Step 3.4 `input_step33_sha256` matches the exact Step 3.3
   segments GeoJSON.
2. **Spatial format:** each GeoJSON reopens with a declared geographic CRS
   (EPSG:4326 / OGC CRS84), and all geometries are nonempty and valid.
3. **Segment indexing:** unique `H99-0001` ... `H99-0101` IDs, correct order,
   and matching per-segment zone IDs. The number of segments is checked against
   both stage metadata files rather than blindly assumed.
4. **Continuous chainage:** lengths, start/end chainages and endpoints are
   consistent, and the ordered segment centerlines overlay the original
   centerline without substantial displacement.
5. **CSV consistency:** each CSV row matches the GeoJSON segment identifier
   and chainage.
6. **Polygon reproducibility:** per-segment 1,000 m polygons, the dissolved
   unique 1,000 m buffer, and the 5,000 m terrain data search buffer are
   reproduced from the stored segment geometries, allowing only small numeric
   GIS round-trip tolerances.
7. **Areas and overlaps:** individual polygon areas are summed **only as
   overlapping screening zones**; the actual unique near-road area is taken
   from the dissolved geometry. The Step 3.4 reported measurements must match.
8. **Disconnected polygon parts:** part count and smallest component area are
   verified and recorded. The tiny part is **not automatically deleted**.

The script fails with a diagnostic error if a required input is missing or a
hard internal check fails. Investigate the cause before any Git commit. It does
**not** change the source data.

## Produced files

- `results/diagnostics/highway99_study_area_step35_audit.json` — portable,
  machine-readable QA record, SHA-256 inputs, checks and unresolved limitations.
  No machine-specific absolute file paths are written to the JSON.
- `results/maps/highway99_study_area_step35_audit.png` — corridor overview and
  magnified location of the smallest 5 km polygon component.

Expected status with the current provisional OSM study-area data:

`INTERNAL_QA_PASSED_WITH_CAVEATS`

This status is NOT equivalent to independent government road validation or
finished landslide risk modelling.

## Review in QGIS

Open the existing QGIS project and overlay these layers on the OpenStreetMap
basemap, with the road centerline on top:

- `highway99_candidate_centerline.geojson`
- `highway99_corridor_screening_1000m.geojson`
- `highway99_terrain_search_5000m.geojson`

The 1 km near-road screening must be contained in the 5 km terrain discovery
area. The reported tiny disconnected part (~142.682 m²) is near the outer edge
of the 5 km polygon. The QA report contains its centroid longitude/latitude;
zoom to that location if necessary. This is a polygon geometry QA item, **not
an isolated landslide**.

## Scientific caveats to carry into the thesis

- The source road alignment and endpoint conventions have not been checked
  independently against an official BC Digital Road Atlas or similar
  authoritative centerline dataset; keep the word *provisional*.
- Segment size (1 km) and two spatial buffer sizes (1 km and 5 km) are
  provisional modelling/data-discovery choices. Later sensitivity tests are
  necessary.
- Per-segment buffers overlap markedly. Never sum them for unique area or
  double-count hazard events/exposed assets; assign events to road segments
  using a declared, deduplicated spatial rule.
- A landslide could originate beyond the 5 km search distance and still reach
  the highway. The 5 km polygon does not define the full physical hazard extent.
- This stage provides no terrain classification, slope analysis, historical
  landslide frequency, runout, damage, closures, traffic loss, or CLIMADA model.

## Completion / release gate

After 15 tests pass and the script reports `INTERNAL_QA_PASSED_WITH_CAVEATS`,
review the JSON and PNG in addition to QGIS before committing these **three
files**:

- `src/study_area/05_validate_study_area.py`
- `tests/test_study_area_step35.py`
- `docs/study_area_step3_5.md`

plus **two generated outputs**, if validated:

- `results/diagnostics/highway99_study_area_step35_audit.json`
- `results/maps/highway99_study_area_step35_audit.png`

The older source files should not be modified or re-committed by Step 3.5.
Use a distinct commit such as `test: audit provisional Highway 99 study-area GIS`.
Next, proceed to terrain and landslide data-source discovery with unresolved
limitations explicitly documented.
