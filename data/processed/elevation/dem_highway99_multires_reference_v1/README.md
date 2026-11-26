# Highway 99 multiresolution reference DEM — v1

## Intended use

Technically validated reference terrain elevation, with mandatory
source-provenance raster and source-transition/anomaly records.
This is not a uniformly resolved 1 m terrain model or a validated
landslide hazard product.

## Main products

- `dem_highway99_multires_reference_v1.tif`: Float32 elevation, NoData -32767.
- `dem_highway99_multires_source_provenance_v1.tif`: Byte source codes:
  0 outside AOI/NoData; 1 HRDEM 1 m; 2 HRDEM 2 m; 30 MRDEM 30 m.
- `dem_highway99_multires_transition_register.gpkg`: 50,875 cardinal
  source edges with absolute elevation jump strictly greater than 10 m.
  This is a descriptive diagnostic threshold, not an accuracy or hazard
  criterion. Other edges are not certified safe.
- `source_anomaly_registers/`: original HRDEM source anomaly registers.
- `qa/`: unchanged diagnostic reports; historical status flags refer to
  their original processing stages.
- `provenance/`: AOI, uncut MRDEM input, source publication manifests and
  processing scripts. The published HRDEM input rasters are referenced by
  path and SHA-256 in the release manifest rather than duplicated.

## Grid and coverage

NAD83 / UTM zone 10N (EPSG:26910), 33000 x 94000 cells, 1 m accounting
grid spacing, bounds (475000, 5463000, 508000, 5557000).

AOI pixel-centres: 1,043,616,657; covered 100%; no valid pixels outside AOI.
Selected source counts:

| Source | Accounting pixels | Area (km2) | Percent of AOI |
|---|---:|---:|---:|
| HRDEM 1 m | 464,090,223 | 464.090223 | 44.469415 |
| HRDEM 2 m | 698,372 | 0.698372 | 0.066918 |
| MRDEM 30 m | 578,828,062 | 578.828062 | 55.463666 |

Native cell values are copied exactly in priority order. No interpolation,
blending, smoothing, vertical correction, or negative-elevation removal
is applied at the integration stage. Source metadata documents CGVD2013,
but the local vertical realization and absolute vertical accuracy were
not independently verified.

## Critical limitations

Do not compute slope, aspect, curvature or other terrain derivatives
directly on this mixed-resolution mosaic. Use source-native neighbourhoods
with separate derivative coverage and consistency QA. Significant elevation
jumps occur at source changes and coarse native-cell boundaries. Cross-source
agreement, absolute vertical accuracy and hazard model validation remain
unaccepted or unvalidated. Negative elevations retain their source values;
the anomaly registers and QA reports must remain associated with this product.

The two release GeoTIFFs retain their validated candidate metadata, including
the historical product-status label, to preserve exact file SHA-256 checksums.
The publication manifest is authoritative for the local release status.
There has been no modification to the original published HRDEM products.
