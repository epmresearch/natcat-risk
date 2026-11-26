"""Cross-resolution elevation agreement diagnostic; never changes DEM inputs.

All differences are lower-resolution elevation minus the area mean of
fully valid, aligned higher-resolution cells over the same footprint.
This is a consistency test, not an independent vertical accuracy test.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
import math
import tempfile

import numpy as np
import step4_2_dem_multires_coverage_qa as c

COVERAGE = c.DIAG / "dem_multires_coverage_qa_xce4r42_" / "multires_coverage_qa.json"
METHOD = "MULTIRES_STRICT_COMMON_SUPPORT_OVERLAP_QA_V1"
THRESHOLDS = (0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10, 20, 50, 100)
DATUM_SOURCES = {
    "hrdem_mosaic": "https://open.canada.ca/data/en/dataset/0fe65119-e96e-4a57-8bfe-9d9245fba06b",
    "mrdem": "https://open.canada.ca/data/en/dataset/18752265-bda3-498c-a4ba-9dfe68cb98da",
}


def aligned_intersection(coarse, *fine):
    """Windows on the coarse native grid, entirely within all raster extents."""
    x0, y0, x1, y1 = c.raster_bounds(coarse)
    r = int(coarse.GetGeoTransform()[1])
    for d in fine:
        fx0, fy0, fx1, fy1 = c.raster_bounds(d)
        fr = int(d.GetGeoTransform()[1])
        c.require(r % fr == 0 and (x0-fx0) % fr == 0 and (y1-fy1) % fr == 0,
                  "Native grids are not integer-aligned")
    boxes = [c.raster_bounds(d) for d in fine]
    left = max([x0]+[b[0] for b in boxes])
    bottom = max([y0]+[b[1] for b in boxes])
    right = min([x1]+[b[2] for b in boxes])
    top = min([y1]+[b[3] for b in boxes])
    col0 = (left-x0+r-1)//r
    col1 = (right-x0)//r
    row0 = (y1-top+r-1)//r
    row1 = (y1-bottom)//r
    c.require(col1 > col0 and row1 > row0, "No common aligned extent")
    return col0, row0, col1, row1


def values(ds, left, top, width, height, resolution):
    """Read native Float32 cells in a window defined on the coarser grid."""
    x, _, _, y = c.raster_bounds(ds)
    r = int(ds.GetGeoTransform()[1])
    c.require(resolution % r == 0 and (left-x) % r == 0 and (y-top) % r == 0,
              "Read window is not aligned")
    factor = resolution//r
    col, row = (left-x)//r, (y-top)//r
    w, h = width*factor, height*factor
    c.require(col >= 0 and row >= 0 and col+w <= ds.RasterXSize
              and row+h <= ds.RasterYSize, "Read would extend outside native raster")
    band = ds.GetRasterBand(1)
    a = band.ReadAsArray(col, row, w, h)
    c.require(a.shape == (h,w) and np.isfinite(a).all(), "Invalid native read")
    v = a != c.NODATA
    if band.GetMaskFlags() != c.gdal.GMF_NODATA:
        c.require(np.array_equal(v, band.GetMaskBand().ReadAsArray(col,row,w,h)!=0),
                  "Mask differs from numeric NoData")
    a = a.reshape(height, factor, width, factor)
    v = v.reshape(height, factor, width, factor)
    count = v.sum(axis=(1,3), dtype=np.uint32)
    mean = a.mean(axis=(1,3), dtype=np.float64)
    return mean, count == factor*factor, count > 0


class Stats:
    def __init__(self, label, resolution):
        self.label = label
        self.resolution = resolution
        self.n = 0
        self.signed = self.absolute = self.squares = 0.0
        self.minimum = math.inf
        self.maximum = -math.inf
        self.over = np.zeros(len(THRESHOLDS), dtype=np.int64)
        self.peak = None

    def add(self, high, low, valid, left, top):
        if not valid.any():
            return
        delta = low[valid] - high[valid]
        ad = np.abs(delta)
        self.n += int(delta.size)
        self.signed += float(delta.sum(dtype=np.float64))
        self.absolute += float(ad.sum(dtype=np.float64))
        self.squares += float(np.dot(delta, delta))
        self.minimum = min(self.minimum, float(delta.min()))
        self.maximum = max(self.maximum, float(delta.max()))
        self.over += np.array([np.count_nonzero(ad > t) for t in THRESHOLDS],
                              dtype=np.int64)
        peak_index = int(np.argmax(ad))
        max_abs = float(ad[peak_index])
        if self.peak is None or max_abs > self.peak["abs_difference_m"]:
            flat = np.flatnonzero(valid)[peak_index]
            i, j = np.unravel_index(flat, valid.shape)
            r = self.resolution
            self.peak = {
                "center_epsg26910": [float(left+(j+0.5)*r), float(top-(i+0.5)*r)],
                "higher_resolution_area_mean_m": float(high[i,j]),
                "lower_resolution_elevation_or_mean_m": float(low[i,j]),
                "signed_difference_m": float(delta[peak_index]),
                "abs_difference_m": max_abs,
            }

    def report(self):
        if self.n == 0:
            return {"comparison": self.label, "strict_common_cells": 0}
        def interval(p):
            # Exact threshold counts give a bounded percentile interval.
            rank = math.ceil(p*self.n)
            lower = 0.0
            for t,n_over in zip(THRESHOLDS, self.over):
                if self.n-int(n_over) >= rank:
                    return [lower, t]
                lower = t
            return [lower, self.peak["abs_difference_m"]]
        return {
            "comparison": self.label, "comparison_grid_resolution_m": self.resolution,
            "strict_common_cells": self.n,
            "strict_common_area_km2": self.n*self.resolution**2/1e6,
            "mean_signed_difference_m": self.signed/self.n,
            "mae_m": self.absolute/self.n,
            "rmse_m": math.sqrt(self.squares/self.n),
            "minimum_signed_difference_m": self.minimum,
            "maximum_signed_difference_m": self.maximum,
            "maximum_absolute_difference": self.peak,
            "absolute_difference_p50_interval_m": interval(0.50),
            "absolute_difference_p95_interval_m": interval(0.95),
            "absolute_difference_p99_interval_m": interval(0.99),
            "counts_over_m": {str(t):int(n) for t,n in zip(THRESHOLDS,self.over)},
            "percent_over_m": {str(t):100*int(n)/self.n for t,n in zip(THRESHOLDS,self.over)},
        }


def compare_2m(one, two, strip_rows=250, progress=True):
    x0, _, _, y1 = c.raster_bounds(two)
    col0,row0,col1,row1 = aligned_intersection(two, one)
    stat = Stats("HRDEM2m minus mean(HRDEM1m 2x2)",2)
    counts = dict.fromkeys(("geometric_candidate_cells", "valid_2m",
                            "valid_2m_with_partial_1m", "valid_2m_with_no_1m"),0)
    total = math.ceil((row1-row0)/strip_rows)
    for i,row in enumerate(range(row0,row1,strip_rows),1):
        h,w = min(strip_rows,row1-row), col1-col0
        left,top = x0+2*col0,y1-2*row
        a,av,any_a = values(one,left,top,w,h,2)
        b,bv,_ = values(two,left,top,w,h,2)
        stat.add(a,b,av&bv,left,top)
        counts["geometric_candidate_cells"] += w*h
        counts["valid_2m"] += int(bv.sum())
        counts["valid_2m_with_partial_1m"] += int((bv&any_a&~av).sum())
        counts["valid_2m_with_no_1m"] += int((bv&~any_a).sum())
        if progress and (i%20==0 or i==total):
            print(f"HRDEM_1M_2M_PROGRESS: {i}/{total}", flush=True)
    c.require(counts["valid_2m"] == stat.n+counts["valid_2m_with_partial_1m"]
              +counts["valid_2m_with_no_1m"],"2 m support does not balance")
    return {"support":counts,"statistics":stat.report()}


def compare_30m(one,two,thirty,strip_rows=12,progress=True):
    x0, _, _, y1 = c.raster_bounds(thirty)
    col0,row0,col1,row1 = aligned_intersection(thirty,one,two)
    labels = [
        "MRDEM30m minus mean(HRDEM1m 30x30)",
        "MRDEM30m minus mean(HRDEM2m 15x15)",
        "MRDEM30m minus mean(HRDEM1m 30x30), triple common",
        "MRDEM30m minus mean(HRDEM2m 15x15), triple common",
        "mean(HRDEM2m 15x15) minus mean(HRDEM1m 30x30), triple common",
    ]
    stats = [Stats(k,30) for k in labels]
    counts = dict.fromkeys(("geometric_candidate_cells","valid_30m",
                            "valid_30m_with_partial_1m","valid_30m_with_no_1m",
                            "valid_30m_with_partial_2m","valid_30m_with_no_2m"),0)
    total = math.ceil((row1-row0)/strip_rows)
    for i,row in enumerate(range(row0,row1,strip_rows),1):
        h,w = min(strip_rows,row1-row),col1-col0
        left,top = x0+30*col0,y1-30*row
        a,av,any_a = values(one,left,top,w,h,30)
        b,bv,any_b = values(two,left,top,w,h,30)
        d,dv,_ = values(thirty,left,top,w,h,30)
        triple=av&bv&dv
        for st,hi,lo,valid in zip(stats,[a,b,a,b,a],[d,d,d,d,b],
                                  [av&dv,bv&dv,triple,triple,triple]):
            st.add(hi,lo,valid,left,top)
        counts["geometric_candidate_cells"] += w*h
        counts["valid_30m"] += int(dv.sum())
        for key,full,any_valid in (("1m",av,any_a),("2m",bv,any_b)):
            counts["valid_30m_with_partial_"+key] += int((dv&any_valid&~full).sum())
            counts["valid_30m_with_no_"+key] += int((dv&~any_valid).sum())
        if progress and (i%20==0 or i==total):
            print(f"HRDEM_MRDEM_30M_PROGRESS: {i}/{total}", flush=True)
    for key,st in zip(("1m","2m"),stats[:2]):
        c.require(counts["valid_30m"] == st.n+counts["valid_30m_with_partial_"+key]
                  +counts["valid_30m_with_no_"+key],"30 m support does not balance")
    c.require(stats[2].n==stats[3].n==stats[4].n,"Triple support mismatch")
    return {"support":counts,"statistics":[st.report() for st in stats]}


def main():
    c.init()
    paths = {**c.INPUTS,"coverage_qa":COVERAGE,
             "coverage_script":Path(c.__file__).resolve()}
    for k,p in paths.items():
        c.require(p.is_file(),f"Required input missing: {k}: {p}")
    print("HASHING_INPUTS",flush=True)
    hashes={k:c.digest(p) for k,p in paths.items()}
    q=json.loads(COVERAGE.read_text(encoding="utf-8"))
    c.require(q["method_id"] == c.METHOD and q["coverage_qa_completed"] is True
              and q["input_files_unchanged"] is True
              and q["full_coverage_achieved_uncut"] is True
              and q["final_integrated_dem_created"] is False,"Coverage QA gate failed")
    c.require(q["script_sha256"]==hashes["coverage_script"],"Coverage script hash mismatch")
    c.require(q["input_sha256"]=={k:hashes[k] for k in c.INPUTS},
              "Input hashes do not match coverage QA")
    c.require(q["counts"]["aoi"]==1043616657 and q["counts"]["selected_1m"]==464090223
              and q["counts"]["selected_2m"]==698372
              and q["counts"]["selected_30m_uncut"]==578828062
              and q["counts"]["uncovered_uncut"]==0,"Unexpected verified coverage inventory")
    mem,layer,envelope,area,target=c.projected_aoi(c.INPUTS["aoi"])
    ds=c.open_rasters(hashes,target)
    print("INPUT_HASH_AND_COVERAGE_GATE_PASSED: True",flush=True)
    p2=compare_2m(ds["hrdem1m"],ds["hrdem2m"])
    p30=compare_30m(ds["hrdem1m"],ds["hrdem2m"],ds["mrdem30_masked"])
    c.require(p2["statistics"]["strict_common_cells"]*4 <= q["counts"]["hrdem_overlap_inside"],
              "2 m strict support exceeds verified overlap")
    c.require(p2["statistics"]["strict_common_cells"]>0
              and all(s["strict_common_cells"]>0 for s in p30["statistics"]),
              "No common support in one or more comparisons")
    ds.clear()
    layer=None
    mem=None
    print("RECHECKING_INPUT_HASHES",flush=True)
    c.require(hashes=={k:c.digest(p) for k,p in paths.items()},"Inputs changed")
    report={
        "method_id":METHOD,"script_sha256":c.digest(Path(__file__)),
        "created_utc":datetime.now(timezone.utc).isoformat(),
        "input_paths":{k:p.relative_to(c.ROOT).as_posix() for k,p in paths.items()},
        "input_sha256":hashes,
        "documented_source_vertical_reference":{
            "hrdem_1m":"CGVD2013","hrdem_2m":"CGVD2013","mrdem_30m":"CGVD2013",
            "evidence_urls":DATUM_SOURCES,
            "product_metadata_compatible":True,
            "local_vertical_crs_or_height_realization_independently_verified":False,
            "vertical_datum_conversion_applied":False,
        },
        "comparison_2m":p2,"comparison_30m":p30,
        "input_files_unchanged":True,"overlap_diagnostic_completed":True,
        "cross_source_elevation_agreement_accepted":False,
        "absolute_vertical_accuracy_validated":False,
        "final_integrated_dem_created":False,
        "limitations":[
            "Lower-resolution values minus strict aligned higher-resolution area means; no interpolation.",
            "30 m comparisons use previously QA-verified AOI-masked MRDEM cells and require full higher-resolution support.",
            "30 m comparison extent is common to both HRDEM raster bounding boxes; pair statistics use separate full-valid masks; triple statistics share exactly the same support.",
            "MRDEM cells may partly intersect the AOI: membership uses the existing native centre mask.",
            "The area mean is a comparison convention, not a reproduction of the provider's resampling algorithm.",
            "Shared LiDAR inputs make agreement non-independent; residuals can reflect resolution, acquisition date and processing.",
            "Threshold counts are descriptive only; no acceptance tolerances, corrections, smoothing, or hazard validation are asserted.",
            "Percentile intervals are derived from exact threshold counts, not exact interpolated quantiles.",
        ],
    }
    out=Path(tempfile.mkdtemp(prefix="dem_multires_overlap_qa_",dir=c.DIAG))
    path=out/"multires_overlap_qa.json"
    path.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("MULTIRES_OVERLAP_SUMMARY:",json.dumps(report,indent=2,allow_nan=False))
    print("DIAGNOSTIC_JSON:",path)
    print("MULTIRES_OVERLAP_QA_COMPLETED: True")
    print("EXISTING_INPUT_FILES_MODIFIED: FALSE")
    print("FINAL_INTEGRATED_DEM: FALSE")


if __name__=="__main__":
    main()
