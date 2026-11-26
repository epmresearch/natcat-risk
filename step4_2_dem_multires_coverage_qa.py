"""Read-only coverage audit for HRDEM 1 m / 2 m and MRDEM 30 m.

Run from the project root. Native pixel support is mapped onto a 1 m
pixel-centre accounting grid without interpolating elevation. Two MRDEM
fallback cases are compared. No integrated raster is created.
"""
from pathlib import Path
from hashlib import sha256
import json
import math
import tempfile
from datetime import datetime, timezone

import numpy as np
from osgeo import gdal, ogr, osr

ROOT = Path(__file__).resolve().parent
DIAG = ROOT / "results" / "diagnostics"
ELEV = ROOT / "data" / "processed" / "elevation"
NODATA = -32767.0
STRIP_ROWS = 500
METHOD = "MULTIRES_NATIVE_CELL_COVERAGE_QA_V1"
INPUTS = {
    "hrdem1m": ELEV / "hrdem_1m_dtm_highway99_fixed_v1" / "hrdem_1m_dtm_highway99_fixed_v1.tif",
    "hrdem2m": ELEV / "hrdem_2m_dtm_highway99_fixed_v1" / "hrdem_2m_dtm_highway99_fixed_v1.tif",
    "mrdem30_masked": ELEV / "mrdem30_dtm_highway99_aoi_utm10_masked_exact.tif",
    "mrdem30_uncut": DIAG / "step4_2_mrdem30_exact_bbox.tif",
    "mrdem30_qa": DIAG / "step4_2_mrdem30_consistent_mask_qa.json",
    "aoi": ROOT / "data" / "processed" / "study_area" / "highway99_terrain_search_5000m.geojson",
}
EXPECTED = {
    "hrdem1m": ((31000, 93000), (475000, 1, 0, 5557000, 0, -1),
                 "eedd95b99c4c0a56067a1c15dcce413cd8ef64b56d68b7145ec159a94643ffb6"),
    "hrdem2m": ((15500, 47000), (475000, 2, 0, 5557000, 0, -2),
                 "5ce4ffa9e13af91a72d19a799ac8dce2dfe1d1b9dea184eeb3908670aa27d6ae"),
    "mrdem30_masked": ((1084, 3093), (475410, 30, 0, 5556570, 0, -30),
                       "b8085012ab7d85a1f2028098345c579e20ee025f24b2189bd300fb2d90c1520d"),
    "mrdem30_uncut": ((1084, 3093), (475410, 30, 0, 5556570, 0, -30), None),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    h = sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def init():
    gdal.UseExceptions()
    ogr.UseExceptions()
    osr.UseExceptions()
    gdal.SetConfigOption("GDAL_PAM_ENABLED", "NO")
    gdal.SetCacheMax(128 * 1024 * 1024)


def projected_aoi(path):
    ds = ogr.Open(str(path), 0)
    require(ds is not None and ds.GetLayerCount() == 1, "Unexpected AOI datasource")
    layer = ds.GetLayer(0)
    require(layer.GetFeatureCount() == 1, "Expected one AOI feature")
    source = layer.GetSpatialRef().Clone()
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(4326)
    source.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    wgs84.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    require(source.IsSame(wgs84), "AOI CRS differs from expected WGS84")
    target = osr.SpatialReference()
    target.ImportFromEPSG(26910)
    target.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    feature = layer.GetNextFeature()
    geometry = feature.GetGeometryRef().Clone()
    require(not geometry.IsEmpty() and geometry.IsValid(), "Invalid AOI geometry")
    require(geometry.Transform(osr.CoordinateTransformation(source, target)) == 0,
            "AOI coordinate transformation failed")
    require(not geometry.IsEmpty() and geometry.IsValid(), "Invalid projected AOI geometry")
    mem = ogr.GetDriverByName("Memory").CreateDataSource("")
    ml = mem.CreateLayer("aoi", target, geometry.GetGeometryType())
    feat = ogr.Feature(ml.GetLayerDefn())
    feat.SetGeometry(geometry)
    require(ml.CreateFeature(feat) == 0, "In-memory AOI creation failed")
    return mem, ml, geometry.GetEnvelope(), geometry.GetArea(), target


def mask_grid(layer, gt, width, height):
    raster = gdal.GetDriverByName("MEM").Create("", width, height, 1, gdal.GDT_Byte)
    raster.SetGeoTransform(gt)
    raster.SetProjection(layer.GetSpatialRef().ExportToWkt())
    raster.GetRasterBand(1).Fill(0)
    layer.ResetReading()
    require(gdal.RasterizeLayer(raster, [1], layer, burn_values=[1],
                               options=["ALL_TOUCHED=FALSE"]) == 0,
            "AOI rasterization failed")
    return raster.ReadAsArray() == 1


def open_rasters(hashes, target):
    rasters = {}
    for name, (size, gt, expected_hash) in EXPECTED.items():
        if expected_hash:
            require(hashes[name] == expected_hash, f"Input hash mismatch: {name}")
        ds = gdal.Open(str(INPUTS[name]), gdal.GA_ReadOnly)
        require(ds is not None, f"Could not open {name}")
        require((ds.RasterXSize, ds.RasterYSize) == size and ds.RasterCount == 1,
                f"Unexpected raster size or band count: {name}")
        require(ds.GetGeoTransform() == gt, f"Unexpected grid: {name}")
        srs = ds.GetSpatialRef()
        require(srs is not None, f"Missing CRS: {name}")
        srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        require(srs.IsSame(target), f"Unexpected CRS: {name}")
        b = ds.GetRasterBand(1)
        require(b.DataType == gdal.GDT_Float32 and b.GetNoDataValue() == NODATA,
                f"Unexpected dtype or NoData: {name}")
        rasters[name] = ds
    return rasters


def read_valid(ds, x, y, w, h):
    band = ds.GetRasterBand(1)
    a = band.ReadAsArray(x, y, w, h)
    require(a.shape == (h, w), "Raster window shape mismatch")
    require(np.isfinite(a).all(), "Nonfinite raster values encountered")
    v = a != NODATA
    # An explicit mask must not contradict the NoData convention.
    if band.GetMaskFlags() != gdal.GMF_NODATA:
        m = band.GetMaskBand().ReadAsArray(x, y, w, h) != 0
        require(np.array_equal(m, v), "GDAL validity mask differs from numeric NoData")
    return v


def raster_bounds(ds):
    x, r, rot1, y, rot2, neg = ds.GetGeoTransform()
    require(rot1 == rot2 == 0 and r == -neg and r in (1, 2, 30)
            and x == int(x) and y == int(y), "Unsupported accounting grid")
    x, y, r = int(x), int(y), int(r)
    return x, y - ds.RasterYSize * r, x + ds.RasterXSize * r, y


def native_support(ds, bounds):
    """Boolean support of the containing native cell for each 1 m centre.

    No elevation interpolation, extrapolation, partial cell weighting,
    or GDAL overview use. All input cell boundaries are integer metres.
    """
    xmin, ymin, xmax, ymax = bounds
    sx, symin, sxmax, sy = raster_bounds(ds)
    r = int(ds.GetGeoTransform()[1])
    result = np.zeros((ymax - ymin, xmax - xmin), dtype=bool)
    left, right = max(xmin, sx), min(xmax, sxmax)
    top, bottom = min(ymax, sy), max(ymin, symin)
    if left >= right or bottom >= top:
        return result
    cx0 = (left - sx) // r
    cx1 = (right - sx + r - 1) // r
    cy0 = (sy - top) // r
    cy1 = (sy - bottom + r - 1) // r
    v = read_valid(ds, cx0, cy0, cx1-cx0, cy1-cy0)
    if r != 1:
        v = np.repeat(np.repeat(v, r, axis=0), r, axis=1)
    dx = left - (sx + cx0 * r)
    dy = (sy - cy0 * r) - top
    result[ymax-top:ymax-bottom, left-xmin:right-xmin] = v[
        dy:dy + top-bottom, dx:dx + right-left]
    return result


def coarse_qa(rasters, layer, previous):
    require(previous.get("status") == "PASSED_INTERNAL_MRDEM30_CONSISTENT_MASK_QA_WITH_CAVEATS",
            "MRDEM QA status not accepted")
    checks = previous.get("checks", {})
    require(bool(checks) and all(v is True for v in checks.values()), "MRDEM QA checks failed")
    m, u = rasters["mrdem30_masked"], rasters["mrdem30_uncut"]
    w, h = m.RasterXSize, m.RasterYSize
    fp = mask_grid(layer, m.GetGeoTransform(), w, h)
    mv, uv = read_valid(m, 0, 0, w, h), read_valid(u, 0, 0, w, h)
    ma, ua = m.ReadAsArray(), u.ReadAsArray()
    require(np.array_equal(mv, fp & uv), "MRDEM masked validity does not match AOI + uncut")
    require(np.array_equal(ma[mv].view(np.uint32), ua[mv].view(np.uint32)),
            "MRDEM masked elevations changed relative to uncut")
    count = int(fp.sum())
    require(count == int(mv.sum()) == 1159588, "MRDEM coarse AOI count not reconciled")
    expected_counts = {"raster_pixels_total": w*h, "aoi_pixel_centres": count,
                       "source_valid_inside_aoi": count, "source_missing_inside_aoi": 0,
                       "output_valid_inside_aoi": count, "output_missing_inside_aoi": 0,
                       "output_valid_outside_aoi": 0, "inside_values_differing_from_source": 0}
    require(all(previous["counts"].get(k) == v for k, v in expected_counts.items()),
            "Previous MRDEM QA counts do not reconcile")
    return {**expected_counts, "uncut_valid_native_cells": int(uv.sum()),
            "masked_retained_values_bitwise_equal": True}


def scan(rasters, layer, bounds, strip_rows=STRIP_ROWS, progress=True):
    xmin, ymin, xmax, ymax = bounds
    width, height = xmax-xmin, ymax-ymin
    keys = ("aoi", "hrdem1m_inside", "hrdem2m_inside", "hrdem_overlap_inside",
            "hrdem1m_total", "hrdem2m_total", "hrdem1m_outside_aoi", "hrdem2m_outside_aoi",
            "mrdem30_masked_inside", "mrdem30_uncut_inside", "selected_1m",
            "selected_2m", "selected_30m_masked", "uncovered_masked",
            "selected_30m_uncut", "uncovered_uncut", "extra_uncut_boundary_support",
            "hrdem_gap")
    counts = dict.fromkeys(keys, 0)
    strips = math.ceil(height / strip_rows)
    for i, top in enumerate(range(ymax, ymin, -strip_rows), 1):
        bottom = max(ymin, top-strip_rows)
        bb = (xmin, bottom, xmax, top)
        aoi = mask_grid(layer, (xmin, 1, 0, top, 0, -1), width, top-bottom)
        v1 = native_support(rasters["hrdem1m"], bb)
        v2 = native_support(rasters["hrdem2m"], bb)
        vm = native_support(rasters["mrdem30_masked"], bb)
        vu = native_support(rasters["mrdem30_uncut"], bb)
        require(not np.any(vm & ~vu), "Masked MRDEM has support missing from uncut")
        gap = aoi & ~v1 & ~v2
        def count(a):
            return int(np.count_nonzero(a))
        row = {
            "aoi": count(aoi),
            "hrdem1m_inside": count(aoi & v1),
            "hrdem2m_inside": count(aoi & v2),
            "hrdem_overlap_inside": count(aoi & v1 & v2),
            "hrdem1m_total": count(v1),
            "hrdem2m_total": count(v2),
            "hrdem1m_outside_aoi": count(v1 & ~aoi),
            "hrdem2m_outside_aoi": count(v2 & ~aoi),
            "mrdem30_masked_inside": count(aoi & vm),
            "mrdem30_uncut_inside": count(aoi & vu),
            "selected_1m": count(aoi & v1),
            "selected_2m": count(aoi & ~v1 & v2),
            "selected_30m_masked": count(gap & vm),
            "uncovered_masked": count(gap & ~vm),
            "selected_30m_uncut": count(gap & vu),
            "uncovered_uncut": count(gap & ~vu),
            "extra_uncut_boundary_support": count(gap & vu & ~vm),
            "hrdem_gap": count(gap),
        }
        for k, n in row.items():
            counts[k] += n
        if progress and (i % 10 == 0 or i == strips):
            print(f"COVERAGE_QA_PROGRESS: {i}/{strips}", flush=True)
    return counts


def summarize(c):
    n = c["aoi"]
    require(n > 0, "Empty AOI accounting grid")
    require(c["selected_1m"] + c["selected_2m"] + c["hrdem_gap"] == n,
            "HRDEM priority counts do not balance")
    require(c["selected_2m"] + c["hrdem_overlap_inside"] == c["hrdem2m_inside"],
            "HRDEM overlap counts do not balance")
    require(c["extra_uncut_boundary_support"] == c["uncovered_masked"] - c["uncovered_uncut"],
            "Uncut fallback difference does not balance")
    result = {}
    for name in ("masked", "uncut"):
        parts = {"hrdem_1m": c["selected_1m"], "hrdem_2m": c["selected_2m"],
                 "mrdem_30m": c["selected_30m_" + name],
                 "uncovered": c["uncovered_" + name]}
        require(sum(parts.values()) == n, "Priority categories do not partition AOI")
        result[name] = {k: {"accounting_pixels_1m": v, "area_km2": v / 1e6,
                            "percent_of_aoi": 100*v/n}
                        for k, v in parts.items()}
        result[name]["covered_percent"] = 100*(n-parts["uncovered"])/n
    return result


def main():
    init()
    for key, path in INPUTS.items():
        require(path.is_file(), f"Required input missing: {key}: {path}")
    print("HASHING_INPUTS", flush=True)
    hashes = {k: digest(p) for k, p in INPUTS.items()}
    mem, layer, envelope, area, target = projected_aoi(INPUTS["aoi"])
    rasters = open_rasters(hashes, target)
    previous = json.loads(INPUTS["mrdem30_qa"].read_text(encoding="utf-8"))
    qa30 = coarse_qa(rasters, layer, previous)
    print("MRDEM_MASK_REFERENCE_QA_PASSED: True", flush=True)
    # Include the full DEM extents so expanded support totals can be
    # reconciled against the published native valid pixel inventories.
    xlo, xhi, ylo, yhi = envelope
    boxes = [raster_bounds(ds) for ds in rasters.values()] + [(xlo, ylo, xhi, yhi)]
    bounds = (math.floor(min(b[0] for b in boxes)/1000)*1000,
              math.floor(min(b[1] for b in boxes)/1000)*1000,
              math.ceil(max(b[2] for b in boxes)/1000)*1000,
              math.ceil(max(b[3] for b in boxes)/1000)*1000)
    x0, y0, x1, y1 = bounds
    print("ACCOUNTING_GRID:", json.dumps({"bounds": bounds, "width": x1-x0,
                                          "height": y1-y0, "resolution_m": 1}),
          flush=True)
    counts = scan(rasters, layer, bounds)
    require(counts["hrdem1m_total"] == 464090223, "Published 1 m valid count mismatch")
    require(counts["hrdem2m_total"] == 4*116196821, "Published 2 m valid count mismatch")
    for label in ("hrdem1m", "hrdem2m"):
        require(counts[label+"_inside"] + counts[label+"_outside_aoi"] == counts[label+"_total"],
                f"Unbalanced source support: {label}")
    scenarios = summarize(counts)
    rasters.clear()
    layer = None
    mem = None
    print("RECHECKING_INPUT_HASHES", flush=True)
    require(hashes == {k: digest(p) for k, p in INPUTS.items()}, "Input files changed")
    summary = {
        "method_id": METHOD,
        "script_sha256": digest(Path(__file__)),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "crs": "EPSG:26910",
        "accounting_resolution_m": 1,
        "accounting_bounds": list(bounds),
        "accounting_grid_pixels": (x1-x0)*(y1-y0),
        "aoi_polygon_area_km2": area/1e6,
        "aoi_pixel_centre_area_km2": counts["aoi"]/1e6,
        "counts": counts,
        "scenarios": scenarios,
        "mrdem30_mask_reference_qa": qa30,
        "input_paths": {k: p.relative_to(ROOT).as_posix() for k, p in INPUTS.items()},
        "input_sha256": hashes,
        "input_files_unchanged": True,
        "coverage_qa_completed": True,
        "full_coverage_achieved_masked": counts["uncovered_masked"] == 0,
        "full_coverage_achieved_uncut": counts["uncovered_uncut"] == 0,
        "vertical_datum_compatibility_validated": False,
        "cross_source_elevation_agreement_validated": False,
        "absolute_vertical_accuracy_validated": False,
        "final_integrated_dem_created": False,
        "limitations": [
            "Coverage only: a 1 m accounting pixel represents containing native-cell support, not 1 m elevation information.",
            "MRDEM uncut support is evaluated only inside the final AOI; no extrapolation is permitted.",
            "Coarse-cell support can cross the exact polygon boundary even when the coarse cell centre is outside.",
            "No source elevations are merged and no cross-source elevation adjustment is made.",
            "Source-specific anomaly registers and terrain derivative limitations still apply.",
        ],
    }
    out = Path(tempfile.mkdtemp(prefix="dem_multires_coverage_qa_", dir=DIAG))
    report = out / "multires_coverage_qa.json"
    report.write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print("MULTIRES_COVERAGE_SUMMARY:", json.dumps(summary, indent=2, allow_nan=False))
    print("DIAGNOSTIC_JSON:", report)
    print("MULTIRES_COVERAGE_QA_COMPLETED: True")
    print("EXISTING_INPUT_FILES_MODIFIED: FALSE")
    print("FINAL_INTEGRATED_DEM: FALSE")


if __name__ == "__main__":
    main()
