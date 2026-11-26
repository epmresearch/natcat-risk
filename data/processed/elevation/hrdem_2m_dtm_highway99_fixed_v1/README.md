# Highway 99 HRDEM 2 m — fixed v1

## Status

Technically validated research input for terrain processing in the natcat-risk
project. This is a separate 2 m dataset, not the integrated 1 m / 2 m / 30 m DEM
and not a validated hazard model. No published HRDEM 1 m files were changed.

## Processing

- Source: NRCan HRDEM 2 m DTM regional mosaics 2_3 and 2_4.
- Native source CRS: EPSG:3979. Output horizontal CRS: EPSG:26910.
- Grid: 2 m, 15,500 columns by 47,000 rows; Float32, NoData -32767.
- Fixed-scale bilinear interpolation: XSCALE=1, YSCALE=1, error threshold 0,
  source overviews disabled, aligned target grid.
- Pixel-center footprint mask; 671 candidate cells, 669 extracted TIFFs and
  two empty extraction markers.
- Local mosaic assembly used direct pixel placement, without further resampling.
- Full-grid bitwise readback and independent footprint validation passed.

## Validation

- Valid cells: 116,196,821; footprint coverage:
  99.97014837%.
- NoData inside footprint: 34,697 cells.
- Valid outside footprint: 0.
- Elevation range: -29.668545 to 1752.513428 m.
- Negative cells: 5,045,476; 606 cells below -5 m.
- All 20 flagged seams reproduced numerically in continuous-source bridge
  tests and independent four-neighbor bilinear reference tests.
- The 606 deep negative cells were numerically consistent with the native
  source. They were not removed, filled or smoothed.
- The anomaly register contains 20 flagged seams, 92 insufficient-support
  boundaries, 273 tiles containing footprint NoData, 119 negative-elevation
  tiles, and 606 deep-negative pixel-center points.
- Seam lines and tile polygons represent reporting units, not the precise
  spatial extent of every affected pixel.

## Limitations

Absolute vertical accuracy, physical validity of negative source elevations,
and hazard model accuracy have NOT been validated. Numerical source fidelity
does not establish that source elevations are true ground elevations.
Coverage gaps remain NoData. The dataset is intended for research with
documented limitations and should not be treated as a hazard prediction.

## Files and reproducibility

The GeoTIFF and GeoPackage are byte-identical copies of the validated candidate
and register. The QA reports retain their original contents; their pre-publication
status flags describe the stage at which each report was generated.
The publication manifest records original relative paths, SHA-256 checksums,
configuration and the final publication status. Large staging tiles and
intermediate seam rasters remain outside this release. The provenance directory
contains the input VRT, footprint, processing configurations, production plan
and small pilot references. Original source URLs are recorded in the manifest;
remote data can change independently of this local release.
