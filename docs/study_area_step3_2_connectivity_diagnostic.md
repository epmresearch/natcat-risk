# Step 3.2 — Offline road-network connectivity diagnosis

The eight OpenStreetMap download tiles were successfully obtained, but the merged
road network is **not yet connected** between the Squamish and Whistler routing
anchors. This document explains how to diagnose the gap without downloading any
more data or constructing artificial links.

## Run

With the `sea_to_sky_risk` Conda environment active, from the repository root:

```powershell
python src/study_area/02c_diagnose_highway99_connectivity.py
```

By default, the script examines the precise OSM node IDs from the failed run:
`481392279` (Squamish) and `3236309869` (Whistler), using the existing cached
`data/cache/study_area/highway99_tiles_16km_9km_v3_combined.graphml`.
The script makes **no network requests** and preserves all cached data.

## Outputs

- `results/diagnostics/step32_connectivity_report.json`
- `results/diagnostics/step32_connectivity_map.png`

The JSON includes the connected-component membership of each anchor,
nearest coordinates of distinct disconnected components, and a list of
components with nodes close to Squamish and Whistler. The dashed line on the
image **is not a road**: it only shows straight-line separation between
nodes in different components.

## Interpretation

1. If there is a common large component within the search radius of both towns,
   route anchoring may need correction. Verify the Highway 99 labels of any
   alternative anchors before changing the routing algorithm.
2. If the two main components have a real geographic gap along Highway 99,
   review the existing tile polygon geometry and retrieve only the missing
   corridor (or compare to BC Digital Road Atlas). Do not assume a straight
   connector is valid.
3. A short distance between two disconnected nodes is *not* proof of road
   connectivity; divided highways, bridge elevations, access controls, and
   OSM topology may explain the separation.
4. Do not use the current network for transport disruption, losses, or road
   segmentation until the extracted corridor is validated.

No Git commit should be made until the diagnostic has been reviewed.
