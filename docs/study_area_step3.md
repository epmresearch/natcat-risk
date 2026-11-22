# Study Area — Step 3.1 (Preliminary)

## Goal
Create a preliminary spatial overview of the Highway 99 Sea-to-Sky corridor in British Columbia, Canada, before acquiring and validating the true road centerline.

## Spatial scopes
- **Physical hazard and infrastructure scope (working assumption):** the Highway 99 corridor between Horseshoe Bay and Whistler.
- **Economic and transport-disruption scope:** travel between Vancouver and Whistler, including relevant connections outside the core hazard area.
- Both scopes are subject to review as better data becomes available.

## How to run
From the repository root, inside the `sea_to_sky_risk` environment:

```bash
python src/study_area/01_create_preliminary_study_area.py
```

The script requires an internet connection to geocode five reference places through OpenStreetMap/Nominatim.

## Outputs
- `data/processed/study_area/reference_places.geojson`: geocoded reference points.
- `data/processed/study_area/schematic_waypoint_axis.geojson`: straight-line links between places **(not Highway 99)**.
- `data/processed/study_area/preliminary_search_area.geojson`: a provisional 10-km buffer around the schematic links, used for initial data discovery **(not a hazard or runout boundary)**.
- `results/maps/study_area_preliminary.png`: quick-look map for manual inspection.

## Important limitations
1. The schematic waypoint line is **not actual road geometry** and must not be used for landslide exposure, road-closure, routing, or damage calculations.
2. The 10-km buffer is a **draft** intended only for finding data; it does not represent the area affected by any landslide event.
3. Place name geocoding yields representative locations, not necessarily highway junction points.
4. The actual Highway 99 line will be acquired and checked against authoritative BC road data and/or OpenStreetMap in Step 3.2.
5. Record the dataset version, access date, coordinate system, completeness and attribution before publishing derived maps.

## Data attribution
Place reference points are derived from OpenStreetMap through the Nominatim service. Credit: © OpenStreetMap contributors. Check OSM/ODbL attribution and redistribution obligations before publication.

## Review checklist before Step 3.2
- [ ] All five reference points geocoded successfully.
- [ ] The preliminary map displays Horseshoe Bay, Lions Bay, Britannia Beach, Squamish, and Whistler in the expected order.
- [ ] GeoJSON files open in QGIS and are stored in EPSG:4326.
- [ ] The next step will replace schematic waypoint links with real highway segments.
