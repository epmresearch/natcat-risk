# Step 3.2 — targeted OSM road-network gap repair

## Why this step exists

After downloading all eight OSMnx tiles, the provisional Highway 99 network
remained disconnected between the Squamish and Whistler anchors. The **offline
connectivity diagnostic** located the nearest nodes of the two components
approximately **1,315.82 m** apart near longitude -123.12 and latitude 50.04–50.05.
A short Euclidean gap is NOT evidence that a road can be drawn between the nodes.

## Inputs

- `data/cache/study_area/highway99_tiles_16km_9km_v3_combined.graphml`
- `results/diagnostics/step32_connectivity_report.json`
- Existing `src/study_area/02_extract_highway99_candidate.py` (rate-limited v4)

## Procedure

From the repository root, with the `sea_to_sky_risk` environment activated:

```bash
python -m unittest discover -s tests -p "test_step32_gap_patch.py"
python src/study_area/02d_repair_osm_gap.py
```

The script downloads **one small, locally buffered OSM road network** surrounding
the diagnosed endpoints and caches it independently. It preserves all road-network
OSM node IDs (`simplify=False`) and composes the graph with the existing base graph
using only identical OSM node IDs. It **does not** add a synthetic edge between
nearby vertices. It verifies undirected network connectivity between the unchanged
Squamish and Whistler anchor IDs before creating a repaired GraphML file.

If connectivity is recovered, it calls the existing Step 3.2 routing code against
this *separately saved* repaired graph and generates the original candidate GIS
files and maps. The base GraphML and the eight original cached tiles are not
modified. The OSM-request code uses OSMnx rate limiting, with bounded retries
only on transient failures. If the patch fails or does not reconnect the graph,
review in QGIS or consult authoritative BC Digital Road Atlas data instead.

## Outputs

- `results/diagnostics/step32_gap_repair_report.json` (always after graph merge)
- `data/cache/study_area/highway99_osm_targeted_gap_patch_4km.graphml` (download)
- `data/cache/study_area/highway99_tiles_16km_9km_v3_gap_repaired.graphml`
  (only if Squamish–Whistler actually reconnect)
- The original Step 3.2 candidate GeoJSON and PNG outputs (only if routing succeeds)

## Important scientific caveats

- This is a **provisional OSM-derived candidate**, not the provincial authoritative
  road centerline and not a validated legally drivable route.
- The buffer is a **data-acquisition polygon**, not a landslide hazard zone.
- Undirected graph connectivity does not prove legal direction of travel.
- Validate route alignment, tagged Highway 99 percentage, and road length
  against a trusted map / provincial highway data before Step 3.3.
- Do not upload bulky GraphML caches to the shared GitHub repository.
- © OpenStreetMap contributors; adhere to ODbL attribution and data rules.
