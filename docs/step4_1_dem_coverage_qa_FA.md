# Step 4.1 — Independent DEM Coverage QA (Persian instructions)

## Purpose
Before finalizing QGIS DEM source discovery, independently audit the saved vector coverage results for HRDEM 1m and 2m. **This is not a check of actual elevation raster pixels, accuracy, landslide hazard, or road geometry.**

## Required inputs (existing project files)
- `data/processed/study_area/highway99_terrain_search_5000m.geojson`
- `data/processed/elevation/hrdem_1m_dtm_coverage_highway99.gpkg`
- `data/processed/elevation/hrdem_2m_dtm_coverage_highway99.gpkg`
- `data/processed/elevation/hrdem_combined_dtm_coverage_highway99.gpkg`
- Optional: `qgis/highway99_dem_coverage.qgz`

## Run from Miniforge Prompt (environment already active)

```bat
cd /d "F:\Projects\University of Calgary\Thesis\natcat-risk"
python -m unittest discover -s tests -p "test_step41_dem_coverage_final.py"
python src/data_sources/02_validate_dem_coverage.py
```

Report: `results/diagnostics/step4_1_dem_coverage_final_qa.json`

## QA checks
- All expected files exist and can be read; all output polygons are nonempty and valid.
- Saved 1m, 2m, and combined GeoPackages use EPSG:26910.
- Area of the AOI independently matches the previous Step 3 report.
- 1m and 2m areas independently match manual QGIS numbers.
- Independent union and overlap are calculated; the saved combined geometry is verified against that union with symmetric difference.
- Polygon coverage outside the 5km AOI is measured.
- Stored numeric area/percent fields are compared with geometry calculations where field names exist.
- Synthetic unit tests verify overlapping, disjoint, incorrect-merge, and outside-AOI cases.

## Important notes
- A `NOTE` about a missing 1m `coverage_km2` field may indicate the unusual field name seen in earlier QGIS screenshots; fix field naming before treating the output as publication-ready.
- A zero vector overlap is plausible if published 1m and 2m coverage polygons are adjacent; this check independently verifies that, rather than relying on the map's apparent colors.
- If output shows `FAILED...` or `ERROR`, do not commit or use the 88.8197% as a final result. Share the terminal output and report for diagnosis.
- `PASSED_INTERNAL_VECTOR_COVERAGE_QA_WITH_CAVEATS` means coverage polygons passed, NOT that all underlying DEM raster cells are present/valid.
- Preserve downloaded STAC GeoPackages and their URLs. Do not commit huge raw files or large derived products until reviewed.
