# Step 4.1 — Terrain, Geology and Landslide Data Source Identification

**Study area:** Sea-to-Sky Highway 99, Horseshoe Bay–Whistler, British Columbia, Canada  
**Repository:** `https://github.com/epmresearch/natcat-risk`  
**Reference GIS:** Stage 3.5 `INTERNAL_QA_PASSED_WITH_CAVEATS`  
**Prepared:** 2026-11-23 (UTC)  
**Stage 4.1 status:** *Research source shortlist prepared; local source-footprint and file-level coverage tests pending.*  

## 1. Objective and boundary

This step identifies and documents candidate datasets for later processing. **No terrain raster, landslide event file or geology layer is downloaded or certified for local Highway 99 coverage in Step 4.1.** The purpose is to decide what to inspect and retrieve next, while preserving source metadata and reproducibility.

Stage 3 source assets:

- `data/processed/study_area/highway99_candidate_centerline.geojson`: provisional OSM line; 100.443243 km.
- `data/processed/study_area/highway99_segments_1000m.geojson`: 101 chainage segments (100 × 1 km + 443.243 m).
- `data/processed/study_area/highway99_terrain_search_5000m.geojson`: screening AOI, **not** a hazard or runout envelope; 1043.616776 km² in EPSG:26910.
- `data/processed/study_area/highway99_spatial_preparation_quality.json`: source QC record. The tiny disconnected 142.682 m² part is documented, not deleted.
- `results/diagnostics/highway99_study_area_step35_audit.json`: cross-stage geospatial and SHA-256 validation provenance.

The metric CRS for local spatial analysis is `EPSG:26910` (NAD83 / UTM Zone 10N). The AOI is stored as GeoJSON in geographic coordinates. Preserve the DEM's **vertical datum**, distinguish DTM from DSM, and record datum transformations; do not assume heights in different products are interchangeable.

## 2. Recommended selection

### Priority P0: baseline DTM — NRCan MRDEM-30

Natural Resources Canada's **Medium Resolution Digital Elevation Model (MRDEM-30)** is the preferred continuous reference DTM for the initial corridor-wide terrain analysis. Its official metadata describes 30 m national coverage, both DTM and DSM products, horizontal CRS **NAD83(CSRS) Canada Atlas Lambert (EPSG:3979)**, and vertical reference **CGVD2013**. The product is composed from airborne lidar-derived NRCan data where available and processed Copernicus DEM elsewhere. The official catalogue describes MRDEM as replacing the formerly supported federal CDEM product.

- Dataset: https://open.canada.ca/data/en/dataset/18752265-bda3-498c-a4ba-9dfe68cb98da
- Access guidance: https://nrcan.github.io/CanElevation/stac-dem-mosaics/accessing-stac-with-qgis/
- Product access guide: https://download-telecharger.services.geo.ca/pub/elevation/dem_mne/MRDEM_MNEMR/CanElevation-MRDEM-Usage-Guide.pdf
- STAC search: https://datacube.services.geo.ca/stac/api/search?collections=mrdem-30
- Catalogue advertises VRT of the *DTM*: https://canelevation-dem.s3.ca-central-1.amazonaws.com/mrdem-30/mrdem-30-dtm.vrt

**Important:** the MRDEM catalogue says these assets are intended for *streaming* access. Do not download a nationwide VRT-backed raster as a single file. In Step 4.2 use QGIS/GDAL/STAC to extract or clip only the necessary tiles around the confirmed AOI, inspect pixel size and NoData, and generate a georeferenced local DTM mosaic. The coverage assertion is a producer-level national statement; we have **not** verified usable cells in every part of the local corridor.

### Priority P0: candidate landslide points — Preliminary Canadian Landslide Database v14

The **Preliminary Canadian Landslide Database** (Brideau et al., version 14.0, released 2026-05-27, DOI:10.5281/zenodo.20371365) contains **31,005 entries for Canada as a whole**, combining recorded movement events and mapped features with morphology consistent with past or ongoing movement. CSV / KMZ access is available. The record may include movement type, material, mapped point position (headscarp/source/transport/deposit), location confidence, and sometimes occurrence date, volume, trigger or references.

- Dataset and downloadable files: https://zenodo.org/records/20371365
- Version-specific DOI: https://doi.org/10.5281/zenodo.20371365

**Do not confuse the national 31,005 record count with Highway 99 local counts.** Local intersection count and date completeness are **unknown** until we ingest and filter the data. Inventory point locations are not automatically the areas of impact on the road. A feature is not necessarily a dated event. Compiling them into a frequency/return period or using them as an authoritative hazard map without independent evidence would be inappropriate. The producers emphasize spatial sampling biases and suitability for initial research/desktop review only. Verify the precise rights and reuse conditions on this *version-specific* Zenodo record before redistributing the raw file.

### Priority P1: high-resolution terrain — NRCan HRDEM

The **High Resolution Digital Elevation Model (HRDEM)** includes lidar-derived products where surveyed. Use it to assess whether higher spatial resolution is available along selected steep terrain or previously reported landslide sites, potentially to compare with MRDEM-30. HRDEM metadata and project footprints are available, but **continuous 1 m coverage along the 100 km corridor is NOT established**.

- Dataset: https://open.canada.ca/data/en/dataset/957782bf-847c-4644-a757-e383c0057995
- Access guide: https://download-telecharger.services.geo.ca/pub/elevation/dem_mne/highresolution_hauteresolution/HRDEM_Download_Instructions.pdf

For each HRDEM candidate: verify the acquisition year, DTM vs DSM, survey method, quality, boundaries, spatial coverage ratio within the AOI and any vertical datum differences. Do not artificially upsample MRDEM-30 to claim lidar-level resolution.

### Priority P1: geology — BC Bedrock Geology

The provincial **Bedrock Geology** compilation provides vector bedrock units and geological context. The official metadata refers to integration of geological map scales from **1:50,000 to 1:250,000** and download options. It is *not* a site-specific rock-mass strength/geotechnical investigation.

- Dataset: https://catalogue.data.gov.bc.ca/dataset/bedrock-geology

### Priority P1: terrain stability mapping and its coverage

Two provincial resources should be assessed together: **Terrain Mapping (TER) Project Boundaries** (where stability mapping projects exist) and **Terrain Stability Mapping (TSM) Detailed Polygons** (the mapped polygon attributes where supplied). Coverage is project-dependent, not automatically wall-to-wall across the AOI. They are complementary to slope stability interpretations, not historical landslide event counts.

- TSM polygons: https://catalogue.data.gov.bc.ca/dataset/terrain-stability-mapping-tsm-detailed-polygons-with-short-attribute-table-spatial-view
- TER project boundaries: https://catalogue.data.gov.bc.ca/dataset/terrain-mapping-ter-project-boundaries
- Broader terrain mapping attributes: https://catalogue.data.gov.bc.ca/dataset/terrain-mapping-ter-polygon-attributes

### Priority P1: independent road geometry QA — BC Digital Road Atlas

The **Digital Road Atlas (DRA), Master Partially-Attributed Roads** is a BC government network layer that exposes a `HIGHWAY_ROUTE_NUMBER` attribute. It is an independent candidate reference for the provisional OSM-derived Highway 99 centerline used in Stages 3.2–3.5. It must be *downloaded and actually compared* before any official-road validation claim.

- Dataset: https://catalogue.data.gov.bc.ca/dataset/digital-road-atlas-dra-master-partially-attributed-roads
- Service metadata example: https://delivery.maps.gov.bc.ca/arcgis/rest/services/whse/bcgw_pub_whse_basemapping/MapServer/68

### Priority P2: surficial deposits and legacy terrain

- **Quaternary Alluvium and Cover**: https://catalogue.data.gov.bc.ca/dataset/geology-quaternary-alluvium-and-cover . Check local coverage, classification and scale before treating as a proxy for material strength.
- **BC CDED 1:250,000 legacy DEM**: https://open.canada.ca/data/en/dataset/7b4fef7e-7cae-4379-97b8-62b03e9ac83d . Not preferred for corridor slope modeling compared to newer MRDEM; retain as historical reference only.

## 3. Prioritization table

| Rank | Product | First task | Local AOI coverage |
| --- | --- | --- | --- |
| P0 | MRDEM-30 DTM | Confirm data access and extract a *small test window* | National provider assertion; cells not yet inspected |
| P0 | Canadian Landslide Database v14 | Check CSV columns and filter AOI + regional context | Unknown; must calculate |
| P1 | HRDEM / lidar | Intersect footprint index with full AOI | Unknown; potentially partial |
| P1 | BC Bedrock Geology | Compare geology polygons against AOI | Not individually intersected |
| P1 | TER/TSM layers | Determine mapped-project coverage | Unknown; project-based |
| P1 | Digital Road Atlas | Compare route geometry to OSM reference | Not compared |
| P2 | Quaternary Alluvium and Cover | Inspect layer availability | Unknown |
| Reference | BC legacy CDED | Reference only, not default | Not checked |

## 4. Scientific methods and quality controls for upcoming stages

1. **AOI intersection must be measured:** `terrain_search_5000m` is an initial data-search window, not a complete terrain catchment, slope initiation domain, or runout corridor. Broad upslope source areas may extend beyond 5 km. Extend/adjust by watershed, geology or observed runout where justified.
2. **Elevation product QA:** retain source asset identifier, tile/footprint and acquisition date, native spatial resolution, horizontal CRS, vertical reference, NoData mask, min/max elevation, mosaic seams, and coverage percentage. Reproject/clip once with an explicitly logged resampling method; prevent accidental slope calculation in degrees.
3. **Inventory QA:** distinguish `event` vs `feature`, location-confidence category, geometry semantics (headscarp vs deposit), dates, duplicates and data-source bias. A point in a search polygon does **not** demonstrate that the road was reached or blocked.
4. **Mixed-resolution caution:** 30 m regional terrain may miss narrow road cuts and small initiation zones; do not claim 1 m detail where only 30 m is available. Conduct resolution sensitivity analysis when high-resolution data are accessible.
5. **No frequency invented:** a preliminary national landslide inventory is not a complete time series suitable, by itself, for estimating annual event rates or CLIMADA event frequencies. Hazard occurrence and intensities will need justified scenario/frequency construction and validation in later research stages.
6. **No GIS area double counting:** Step 3.4 per-segment 1 km buffers overlap extensively. Keep a unique union area when reporting landscape extents and define explicit incident-to-road-segment attribution later.
7. **Licences and data provenance:** retain metadata, access timestamps, the applicable dataset-specific attribution/licence, downloaded-file SHA-256 and source citations. Verify version-specific Zenodo reuse terms; never assume they match the BC/Canada open-government licences.

## 5. Reproducible offline catalog and AOI check

The package contains:

```text
src/data_sources/01_audit_step41_sources.py
tests/test_step4_1_data_sources.py
data/catalog/step4_1_source_registry.csv
docs/step4_1_data_source_assessment.md
docs/step4_1_execution_FA.md
```

Execute from the repository root (active `sea_to_sky_risk` environment):

```bat
python -m unittest discover -s tests -p "test_step4_1_data_sources.py"
python src/data_sources/01_audit_step41_sources.py
```

Generated output: `results/diagnostics/step4_1_source_discovery_audit.json`.

The audit checks local AOI validity, its bounding box and its area against Step 3.4 QA, the source registry's schema, identifiers, priorities and valid HTTPS metadata URLs. It **does not** query providers or verify coverage. The precise metadata URL, recorded license wording and source currency should be checked again before downloads or publication.

## 6. Completion conditions for Step 4.1

- [ ] Source registry and audit test pass in the actual local `natcat-risk` checkout.
- [ ] AOI EPSG, bounds and area approximately match Step 3.4 (1043.616776 km²).
- [ ] MRDEM-30 DTM access path and format assessed in QGIS or provider documentation.
- [ ] HRDEM coverage footprint inspected before selecting high-resolution tiles.
- [ ] Landslide v14 inventory format, availability and local corridor intersection plan documented; actual local count may be deferred to inventory ingestion.
- [ ] Data licenses, source versions and known sampling biases recorded.
- [ ] Limitations and methodology reviewed; stage status remains `DATA_SOURCES_SHORTLISTED_AWAITING_SPATIAL_COVERAGE_VERIFICATION` until live provider checks.

**Stage handoff:** Step 4.2 will retrieve only *tested and geographically scoped* elevation files (MRDEM-30 DTM first), derive a clipped terrain mosaic with QA and then evaluate higher-resolution supplements. GIS search buffers must not be relabelled as landslide hazard zones.

## 7. Source references (checked 2026-11-23)

1. NRCan, [MRDEM-30 catalogue](https://open.canada.ca/data/en/dataset/18752265-bda3-498c-a4ba-9dfe68cb98da); [STAC in QGIS guide](https://nrcan.github.io/CanElevation/stac-dem-mosaics/accessing-stac-with-qgis/).
2. NRCan, [HRDEM catalogue](https://open.canada.ca/data/en/dataset/957782bf-847c-4644-a757-e383c0057995); [access guide](https://download-telecharger.services.geo.ca/pub/elevation/dem_mne/highresolution_hauteresolution/HRDEM_Download_Instructions.pdf).
3. Brideau et al. (2026), [Preliminary Canadian Landslide Database v14](https://zenodo.org/records/20371365), DOI:10.5281/zenodo.20371365.
4. BC Geological Survey, [Bedrock Geology](https://catalogue.data.gov.bc.ca/dataset/bedrock-geology).
5. Government of BC, [TSM Detailed Polygons](https://catalogue.data.gov.bc.ca/dataset/terrain-stability-mapping-tsm-detailed-polygons-with-short-attribute-table-spatial-view), [TER Project Boundaries](https://catalogue.data.gov.bc.ca/dataset/terrain-mapping-ter-project-boundaries).
6. GeoBC, [Digital Road Atlas MPAR](https://catalogue.data.gov.bc.ca/dataset/digital-road-atlas-dra-master-partially-attributed-roads).
7. BC Geological Survey, [Quaternary Alluvium and Cover](https://catalogue.data.gov.bc.ca/dataset/geology-quaternary-alluvium-and-cover).
8. Government of BC, [Legacy CDED 1:250,000](https://open.canada.ca/data/en/dataset/7b4fef7e-7cae-4379-97b8-62b03e9ac83d).
