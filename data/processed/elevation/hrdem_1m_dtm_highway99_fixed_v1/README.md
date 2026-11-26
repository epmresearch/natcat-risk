# HRDEM 1 m DTM — Highway 99 — Fixed Scale v1

Technically validated research input. Not an independently validated terrain survey or landslide hazard model.

- DEM: `hrdem_1m_dtm_highway99_fixed_v1.tif`
- Spatial diagnostic register: `hrdem_1m_dtm_highway99_anomaly_register.gpkg`
- CRS: NAD83 / UTM zone 10N (EPSG:26910); 1 m grid; Float32; NoData -32767.
- Vertical reference: CGVD2013 orthometric heights, according to official product documentation. The audited source WKT is horizontal-only and the source band unit field is empty. No vertical transformation was applied.
- Source: Natural Resources Canada, HRDEM Mosaic, CanElevation Series, DTM. Open Government Licence - Canada. https://open.canada.ca/data/en/dataset/0fe65119-e96e-4a57-8bfe-9d9245fba06b
- Fixed-scale point bilinear resampling (XSCALE=YSCALE=1, exact transformer, no source overviews), pixel-center footprint masking and exact pixel mosaicking.
- Valid pixels: 464,090,223. Rasterized footprint pixels: 464,219,963. Missing inside footprint: 129,740 (99.972052% coverage). This is not whole-AOI coverage.
- All elevation values, including negative values, retained. No smoothing, gap filling, seam replacement or automatic anomaly exclusion.
- All 2,432 pixels below -5 m agree with the strict native bilinear reference within 0.001 m. This does not establish physical validity.
- Diagnostic register: 41 flagged seams, 52 insufficient-data seams, 265 tiles with internal NoData, and extreme differences between processing methods. Do not treat these features as a terrain-error or exclusion mask.
- Absolute vertical accuracy, source artifact origin and physical landslide model sensitivity remain unresolved.

`publication_manifest.json` records hashes and lineage. `qa/` contains the ten unchanged QA reports. The VRT snapshot is provenance evidence, not the final raster or a pinned copy of the remote sources. Original diagnostics and staging files remain untouched.
