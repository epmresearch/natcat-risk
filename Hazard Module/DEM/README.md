# DEM Preparation (LiDAR + National DEM)

A **Digital Elevation Model (DEM)** is a raster representation of the Earth’s surface elevation.  
In this project, the DEM is a key topographic input used to derive terrain variables (e.g., slope, aspect, curvature) that support **landslide susceptibility / hazard analysis**.  
Because terrain controls runoff, slope stability, and potential failure zones, a high‑quality DEM is essential for reliable downstream modeling.

---

## Data Sources

### 1) LiDAR DEM (Primary)
We use **LiDAR-derived DEM tiles** as our primary elevation dataset, downloaded from the Government of British Columbia LiDAR portal:  
`https://lidar.gov.bc.ca`

- **Count:** 139 tiles  
- **Resolution:** ~1–2 m (high detail)  
- **Limitation:** coverage is **incomplete** and does not fully cover our study area.

### 2) National DEM (Gap-Filling / Complementary)
To fill parts of the study area that are not covered by LiDAR, we use **national DEM tiles** published for all of Canada via:  
`https://open.canada.ca`

- **Count:** 128 tiles  
- **Resolution:** ~10–30 m (coarser than LiDAR)  
- **Role:** used only to **fill LiDAR gaps**, not to replace LiDAR where available.

---

## QGIS Processing Workflow

All processing was performed in **QGIS**. The target projection for this project is:

- **CRS:** EPSG:3005 (NAD83 / BC Albers)

### A) LiDAR Tiles (139 files)

1. **Reproject LiDAR layers to EPSG:3005**
   - QGIS: *Raster → Projections → Warp (Reproject)*  
   - Output CRS: **EPSG:3005**

2. **Merge all LiDAR tiles into one raster**
   - QGIS: *Raster → Miscellaneous → Merge*  
   - Output: a single merged LiDAR DEM (e.g., `LiDAR_Merged_3005.tif`)

<p align="center">
  <img src="RAW%20DEM%20LIDAR.png" width="700">
</p>


---

### B) National DEM Tiles (128 files)

1. **Reproject National DEM layers to EPSG:3005**
   - QGIS: *Raster → Projections → Warp (Reproject)*  
   - Output CRS: **EPSG:3005**

2. **Merge all national tiles into one raster**
   - QGIS: *Raster → Miscellaneous → Merge*  
   - Output: a single national DEM (e.g., `National_DEM_Merged_3005.tif`)

3. **Clip merged national DEM to the study area**
   - Because national DEM covers a very large extent, we clip it to our study boundary.
   - QGIS: *Raster → Extraction → Clip Raster by Mask Layer* (or by extent)
   - Output: clipped national DEM (e.g., `National_DEM_3005_Clipped.tif`)

<p align="center">
  <img src="National%20DEMs.png" width="700">
</p>

---

## Create the Final “Single DEM” (LiDAR priority + National gap fill)

After preparing both rasters (LiDAR merged + National clipped), we created one final DEM where:

- LiDAR values are used wherever available (higher resolution / preferred).
- National DEM values fill only the **missing / uncovered** LiDAR areas.

### Raster Calculator Step

In QGIS Raster Calculator we used the following expression:

`max("LiDAR_Merged@1", "National_DEM@1")`

This produces a single raster with **LiDAR priority** (assuming LiDAR NoData is correctly handled).

---

## Outputs (Expected)

- `LiDAR_Merged_3005.tif`  
- `National_DEM_3005_Clipped.tif`  
- `DEM_Final_3005.tif` (final combined DEM used in the project)

<p align="center">
  <img src="Output.png" width="700">
</p>

