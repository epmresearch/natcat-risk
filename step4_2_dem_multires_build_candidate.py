"""Build and verify a reference elevation mosaic plus native-source code grid.

This is a diagnostic candidate, not a uniformly resolved 1 m terrain model.
No interpolation, vertical adjustment, smoothing or terrain derivatives.
Existing published products are read-only.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
import math
import shutil
import tempfile

import numpy as np
import step4_2_dem_multires_coverage_qa as c
import step4_2_dem_multires_overlap_qa as o
import step4_2_dem_multires_transition_qa as t

METHOD="MULTIRES_PRIORITY_REFERENCE_CANDIDATE_V1"
TRANSITION=c.DIAG/"dem_multires_transition_qa_weczzfdj"/"multires_transition_qa.json"
OVERLAP_SHA="1fb278da43959fbbb7452daf3e692cef9d3e30f6a56a2be9e3107790e3c009f4"
BOUNDS=t.EXPECTED_BOUNDS
CODES={0:"Outside AOI / NoData",1:"HRDEM 1 m",2:"HRDEM 2 m",30:"MRDEM 30 m uncut"}
DEM_NAME="dem_highway99_multires_reference_candidate.tif"
SOURCE_NAME="dem_highway99_multires_source_candidate.tif"
QA_NAME="mosaic_candidate_qa.json"
REGISTERS={
    "hrdem1m_register":(
        c.ELEV/"hrdem_1m_dtm_highway99_fixed_v1"/"hrdem_1m_dtm_highway99_anomaly_register.gpkg",
        "e2e8840013db83074fc614d3a12d299ee6d199e09465172bf59ffb9e489c76ae"),
    "hrdem2m_register":(
        c.ELEV/"hrdem_2m_dtm_highway99_fixed_v1"/"hrdem_2m_dtm_highway99_anomaly_register.gpkg",
        "e516be189e73f4b8c4b6c48d2905baad857d02fe7dd103fe0ddf09069b8d164d"),
}


def qa_inputs():
    paths={**c.INPUTS,"coverage_qa":o.COVERAGE,"overlap_qa":t.OVERLAP,
           "transition_qa":TRANSITION,"coverage_script":Path(c.__file__).resolve(),
           "overlap_script":Path(o.__file__).resolve(),"transition_script":Path(t.__file__).resolve(),
           **{k:v[0] for k,v in REGISTERS.items()}}
    for k,p in paths.items():
        c.require(p.is_file(),f"Required input missing: {k}: {p}")
    print("HASHING_INPUTS",flush=True)
    hashes={k:c.digest(p) for k,p in paths.items()}
    for key,(_,expected) in REGISTERS.items():
        c.require(hashes[key]==expected,f"Published anomaly register changed: {key}")
    cq=json.loads(paths["coverage_qa"].read_text(encoding="utf-8"))
    oq=json.loads(paths["overlap_qa"].read_text(encoding="utf-8"))
    tq=json.loads(paths["transition_qa"].read_text(encoding="utf-8"))
    c.require(hashes["coverage_qa"]==t.COVERAGE_SHA and hashes["overlap_qa"]==OVERLAP_SHA,
              "Coverage or overlap QA file changed")
    c.require(cq["method_id"]==c.METHOD and cq["coverage_qa_completed"] is True
              and cq["input_files_unchanged"] is True and cq["full_coverage_achieved_uncut"] is True,
              "Coverage QA gate failed")
    c.require(oq["method_id"]==o.METHOD and oq["overlap_diagnostic_completed"] is True
              and oq["input_files_unchanged"] is True,"Overlap diagnostic gate failed")
    c.require(tq["method_id"]==t.METHOD and tq["source_transition_diagnostic_completed"] is True
              and tq["coverage_inventory_reconciled"] is True
              and tq["input_files_unchanged"] is True,"Transition diagnostic gate failed")
    for report,script in ((cq,"coverage_script"),(oq,"overlap_script"),(tq,"transition_script")):
        c.require(report["script_sha256"]==hashes[script],f"Script hash mismatch: {script}")
        c.require(report["final_integrated_dem_created"] is False,"Unexpected historical QA state")
        c.require(all(k in hashes and hashes[k]==v for k,v in report["input_sha256"].items()),
                  f"QA input digest mismatch: {script}")
    c.require(cq["input_sha256"]=={k:hashes[k] for k in c.INPUTS},
              "Incomplete coverage QA lineage")
    c.require(tq["input_sha256"]=={k:hashes[k] for k in
              (*c.INPUTS,"coverage_qa","overlap_qa","coverage_script","overlap_script")},
              "Incomplete transition QA lineage")
    c.require(tuple(cq["accounting_bounds"])==tuple(tq["accounting_bounds"])==BOUNDS
              and tq["aoi_accounting_pixels"]==1043616657
              and tq["uncovered_accounting_pixels"]==0
              and {int(k):v for k,v in tq["selected_source_pixels"].items()}==t.EXPECTED_COUNTS,
              "Priority coverage inventory mismatch")
    c.require(tq["total_transition_edges"]==772181 and
              {k:p["transition_edges"] for k,p in tq["transition_pairs"].items()}==
              {"1_to_2":375077,"1_to_30":23078,"2_to_30":374026},
              "Unexpected transition inventory")
    c.require(tq["source_transitions_accepted_for_terrain_derivatives"] is False and
              tq["cross_source_elevation_agreement_accepted"] is False,
              "Unexpected scientific acceptance flags")
    print("INPUT_QA_CHAIN_AND_REGISTERS_VERIFIED: True",flush=True)
    return paths,hashes,tq


def create(path,bounds,srs,source=False):
    g=c.gdal
    xmin,ymin,xmax,ymax=bounds
    ds=g.GetDriverByName("GTiff").Create(str(path),xmax-xmin,ymax-ymin,1,
              g.GDT_Byte if source else g.GDT_Float32,
              options=["BIGTIFF=YES","TILED=YES","BLOCKXSIZE=256","BLOCKYSIZE=128",
                       "COMPRESS=DEFLATE","PREDICTOR=2" if source else "PREDICTOR=3",
                       "ZLEVEL=6","NUM_THREADS=2"])
    c.require(ds is not None,f"Could not create {path}")
    ds.SetGeoTransform((xmin,1,0,ymax,0,-1))
    ds.SetProjection(srs.ExportToWkt())
    common={
        "METHOD_ID":METHOD,"AREA_OR_POINT":"Area",
        "ACCOUNTING_GRID_SPACING_M":"1",
        "SOURCE_PRIORITY":"HRDEM_1M > HRDEM_2M > MRDEM_30M_UNCUT",
        "SOURCE_ASSIGNMENT":"Containing native cell; AOI pixel-centre mask",
        "VERTICAL_REFERENCE_DOCUMENTED_BY_PROVIDER":"CGVD2013",
        "ABSOLUTE_VERTICAL_ACCURACY_VALIDATED":"NO",
        "CROSS_SOURCE_ELEVATION_AGREEMENT_ACCEPTED":"NO",
        "TERRAIN_DERIVATIVES_FROM_MOSAIC_APPROVED":"NO",
        "PRODUCT_STATUS":"DIAGNOSTIC_CANDIDATE",
        "SOURCE_CODES":"0=NoData; 1=HRDEM1m; 2=HRDEM2m; 30=MRDEM30m",
    }
    ds.SetMetadata(common)
    band=ds.GetRasterBand(1)
    band.SetNoDataValue(0 if source else c.NODATA)
    band.SetDescription("Selected native source resolution code" if source else
                        "Reference terrain elevation from selected native source")
    if not source:
        band.SetUnitType("m")
    return ds


def compose(rasters,aoi,bb):
    """Write rule: assign coarsest first, overwrite with finer valid data."""
    z=np.full(aoi.shape,c.NODATA,np.float32)
    lab=np.zeros(aoi.shape,np.uint8)
    for code in (30,2,1):
        vals=t.native_values(rasters[t.SOURCE_KEYS[code]],bb)
        use=aoi&(vals!=c.NODATA)
        z[use]=vals[use]
        lab[use]=code
    return z,lab


def build(rasters,layer,bounds,dem,src,strip_rows=c.STRIP_ROWS,progress=True):
    xmin,ymin,xmax,ymax=bounds
    w=xmax-xmin
    counts=dict.fromkeys(CODES,0)
    aoi_total=0
    strips=math.ceil((ymax-ymin)/strip_rows)
    for i,top in enumerate(range(ymax,ymin,-strip_rows),1):
        bottom=max(ymin,top-strip_rows)
        bb=xmin,bottom,xmax,top
        aoi=c.mask_grid(layer,(xmin,1,0,top,0,-1),w,top-bottom)
        z,lab=compose(rasters,aoi,bb)
        c.require(np.array_equal(z!=c.NODATA,lab!=0),"Build validity and provenance mismatch")
        c.require(np.array_equal(lab!=0,aoi),"Missing AOI data or outside AOI data")
        c.require(dem.GetRasterBand(1).WriteArray(z,0,ymax-top)==0,"DEM write failed")
        c.require(src.GetRasterBand(1).WriteArray(lab,0,ymax-top)==0,"Source write failed")
        aoi_total+=int(aoi.sum())
        for k in counts:
            counts[k]+=int((lab==k).sum())
        if progress and (i%10==0 or i==strips):
            print(f"MOSAIC_WRITE_PROGRESS: {i}/{strips}",flush=True)
    dem.FlushCache()
    src.FlushCache()
    return counts,aoi_total


def index_native_values(ds,bb):
    """Independent direct native index gather, used for full readback validation."""
    xmin,ymin,xmax,ymax=bb
    sx,_,_,sy=c.raster_bounds(ds)
    r=int(ds.GetGeoTransform()[1])
    # Integer 1 m coordinates plus 0.5 m centres yield these floor indices.
    cols=(np.arange(xmin,xmax,dtype=np.int64)-sx)//r
    rows=(sy-1-np.arange(ymax-1,ymin-1,-1,dtype=np.int64))//r
    # Rows above are equivalent to floor((sy - centre_y) / r).
    okx=(cols>=0)&(cols<ds.RasterXSize)
    oky=(rows>=0)&(rows<ds.RasterYSize)
    out=np.full((ymax-ymin,xmax-xmin),c.NODATA,np.float32)
    if not okx.any() or not oky.any():
        return out
    cx,ry=cols[okx],rows[oky]
    x0,y0=int(cx.min()),int(ry.min())
    w,h=int(cx.max())-x0+1,int(ry.max())-y0+1
    band=ds.GetRasterBand(1)
    a=band.ReadAsArray(x0,y0,w,h)
    c.require(a.shape==(h,w) and a.dtype==np.float32 and np.isfinite(a).all(),
              "Native verification window invalid")
    if band.GetMaskFlags()!=c.gdal.GMF_NODATA:
        c.require(np.array_equal(a!=c.NODATA,band.GetMaskBand().ReadAsArray(x0,y0,w,h)!=0),
                  "Native verification mask mismatch")
    out[np.ix_(oky,okx)]=a[np.ix_(ry-y0,cx-x0)]
    return out


def inspect(ds,bounds,srs,source=False):
    g=c.gdal
    xmin,ymin,xmax,ymax=bounds
    c.require(ds is not None and (ds.RasterXSize,ds.RasterYSize)==(xmax-xmin,ymax-ymin)
              and ds.RasterCount==1 and ds.GetGeoTransform()==(xmin,1,0,ymax,0,-1),
              "Candidate grid mismatch")
    actual=ds.GetSpatialRef()
    actual.SetAxisMappingStrategy(c.osr.OAMS_TRADITIONAL_GIS_ORDER)
    c.require(actual.IsSame(srs),"Candidate CRS mismatch")
    b=ds.GetRasterBand(1)
    c.require(b.DataType==(g.GDT_Byte if source else g.GDT_Float32)
              and b.GetNoDataValue()==(0 if source else c.NODATA),"Candidate band mismatch")
    c.require(ds.GetMetadataItem("METHOD_ID")==METHOD
              and ds.GetMetadataItem("TERRAIN_DERIVATIVES_FROM_MOSAIC_APPROVED")=="NO"
              and ds.GetMetadataItem("SOURCE_CODES")=="0=NoData; 1=HRDEM1m; 2=HRDEM2m; 30=MRDEM30m",
              "Candidate provenance metadata mismatch")
    c.require(b.GetOverviewCount()==0,"Unexpected candidate overviews")


def verify(rasters,layer,bounds,dem,src,strip_rows=c.STRIP_ROWS,progress=True):
    xmin,ymin,xmax,ymax=bounds
    w=xmax-xmin
    counts=dict.fromkeys(CODES,0)
    statistics={k:{"count":0,"minimum_m":math.inf,"maximum_m":-math.inf,
                   "negative_pixels":0} for k in (1,2,30)}
    aoi_total=0
    strips=math.ceil((ymax-ymin)/strip_rows)
    for i,top in enumerate(range(ymax,ymin,-strip_rows),1):
        bottom=max(ymin,top-strip_rows)
        h=top-bottom
        bb=xmin,bottom,xmax,top
        aoi=c.mask_grid(layer,(xmin,1,0,top,0,-1),w,h)
        z=dem.GetRasterBand(1).ReadAsArray(0,ymax-top,w,h)
        lab=src.GetRasterBand(1).ReadAsArray(0,ymax-top,w,h)
        c.require(z.shape==lab.shape==aoi.shape and np.isfinite(z).all(),"Invalid readback")
        c.require(np.array_equal(lab!=0,aoi) and np.array_equal(z!=c.NODATA,aoi),
                  "Candidate mask differs from AOI")
        c.require(np.all(z[~aoi].view(np.uint32)==np.float32(c.NODATA).view(np.uint32)),
                  "Outside AOI elevation not exact NoData")
        remaining=aoi.copy()
        for code in (1,2,30):
            expected=index_native_values(rasters[t.SOURCE_KEYS[code]],bb)
            selected=remaining&(expected!=c.NODATA)
            c.require(np.array_equal(lab==code,selected),
                      f"Source selection mismatch at strip {i}, code {code}")
            actual=z[selected]
            c.require(np.array_equal(actual.view(np.uint32),expected[selected].view(np.uint32)),
                      f"Elevation bitwise mismatch at strip {i}, code {code}")
            n=int(selected.sum())
            st=statistics[code]
            st["count"]+=n
            if n:
                st["minimum_m"]=min(st["minimum_m"],float(actual.min()))
                st["maximum_m"]=max(st["maximum_m"],float(actual.max()))
                st["negative_pixels"]+=int((actual<0).sum())
            remaining &= ~selected
            counts[code]+=n
        c.require(not remaining.any(),"Uncovered AOI in readback")
        counts[0]+=int((~aoi).sum())
        aoi_total+=int(aoi.sum())
        if progress and (i%10==0 or i==strips):
            print(f"MOSAIC_READBACK_PROGRESS: {i}/{strips}",flush=True)
    for st in statistics.values():
        if st["count"]==0:
            st["minimum_m"]=st["maximum_m"]=None
    return counts,aoi_total,statistics


def main():
    c.init()
    paths,hashes,tq=qa_inputs()
    mem,layer,envelope,area,target=c.projected_aoi(c.INPUTS["aoi"])
    rasters=c.open_rasters(hashes,target)
    xmin,ymin,xmax,ymax=BOUNDS
    grid_pixels=(xmax-xmin)*(ymax-ymin)
    free=shutil.disk_usage(c.DIAG).free
    c.require(free>=grid_pixels*5+5*1024**3,"Insufficient disk space for two candidate rasters")
    expected={**t.EXPECTED_COUNTS,0:grid_pixels-sum(t.EXPECTED_COUNTS.values())}
    out=Path(tempfile.mkdtemp(prefix="dem_multires_mosaic_candidate_",dir=c.DIAG))
    print("OUTPUT_DIRECTORY:",out,flush=True)
    print("MOSAIC_PLAN:",json.dumps({"bounds":BOUNDS,"width":xmax-xmin,"height":ymax-ymin,
          "grid_pixels":grid_pixels,"uncompressed_payload_bytes":grid_pixels*5,
          "free_bytes":free}),flush=True)
    dp=out/(DEM_NAME+".partial.tif")
    sp=out/(SOURCE_NAME+".partial.tif")
    dem=create(dp,BOUNDS,target)
    src=create(sp,BOUNDS,target,source=True)
    written,aoi=build(rasters,layer,BOUNDS,dem,src)
    c.require(written==expected and aoi==1043616657,"Build inventory differs from coverage QA")
    dem=None
    src=None
    print("MOSAIC_WRITE_COMPLETED: True",flush=True)
    dem=c.gdal.Open(str(dp),c.gdal.GA_ReadOnly)
    src=c.gdal.Open(str(sp),c.gdal.GA_ReadOnly)
    inspect(dem,BOUNDS,target)
    inspect(src,BOUNDS,target,source=True)
    readback,readback_aoi,stats=verify(rasters,layer,BOUNDS,dem,src)
    c.require(readback==written==expected and readback_aoi==aoi,"Readback inventory mismatch")
    dem=None
    src=None
    rasters.clear()
    layer=None
    mem=None
    print("RECHECKING_INPUT_HASHES",flush=True)
    c.require(hashes=={k:c.digest(p) for k,p in paths.items()},"Existing inputs changed")
    # Promote only within the unique diagnostics directory, never to processed data.
    final_dem,final_src=out/DEM_NAME,out/SOURCE_NAME
    c.require(not final_dem.exists() and not final_src.exists(),"Candidate destination exists")
    dp.rename(final_dem)
    sp.rename(final_src)
    files={k:{"path":p.relative_to(c.ROOT).as_posix(),"sha256":c.digest(p),"bytes":p.stat().st_size}
           for k,p in (("reference_dem",final_dem),("source_provenance",final_src))}
    report={
        "method_id":METHOD,"script_sha256":c.digest(Path(__file__)),
        "created_utc":datetime.now(timezone.utc).isoformat(),
        "crs":"EPSG:26910","accounting_resolution_m":1,"bounds_epsg26910":BOUNDS,
        "width":xmax-xmin,"height":ymax-ymin,"geotransform":[xmin,1,0,ymax,0,-1],
        "total_grid_pixels":grid_pixels,"aoi_pixels":aoi,"valid_pixels":aoi,
        "outside_aoi_nodata_pixels":expected[0],"uncovered_aoi_pixels":0,
        "coverage_percent":100.0,"source_codes":CODES,
        "source_pixel_counts":readback,"source_elevation_statistics":stats,
        "reference_dem_dtype":"Float32","reference_dem_nodata":c.NODATA,
        "source_provenance_dtype":"Byte","source_provenance_nodata":0,
        "input_paths":{k:p.relative_to(c.ROOT).as_posix() for k,p in paths.items()},
        "input_sha256":hashes,"files":files,
        "source_transition_edges":tq["total_transition_edges"],
        "source_transition_qa_completed":True,
        "source_transitions_accepted_for_terrain_derivatives":False,
        "cross_source_elevation_agreement_accepted":False,
        "documented_vertical_reference":"CGVD2013",
        "local_height_realization_independently_verified":False,
        "source_priority_exactly_verified":True,
        "native_selected_elevation_bitwise_equal":True,
        "all_grid_pixels_readback_verified":True,"input_files_unchanged":True,
        "mosaic_candidate_qa_passed":True,
        "final_integrated_dem_published":False,
        "absolute_vertical_accuracy_validated":False,"hazard_model_validated":False,
        "limitations":[
            "Reference elevation only: 1 m grid spacing is not uniform native resolution.",
            "Coarse source values repeated without interpolation; source code raster is required.",
            "Substantial source-transition jumps have been recorded; no correction or smoothing applied.",
            "Do not compute terrain derivatives directly on this mixed-resolution mosaic.",
            "Derivatives require source-native, valid neighbourhood processing and separate QA.",
            "Original source anomaly registers remain applicable; no negative elevations changed.",
            "Complete AOI coverage is not proof of accurate or hazard-ready elevation.",
        ],
    }
    p=out/QA_NAME
    p.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("MULTIRES_MOSAIC_CANDIDATE_SUMMARY:",json.dumps(report,indent=2,allow_nan=False))
    print("REFERENCE_DEM_CANDIDATE:",final_dem)
    print("SOURCE_PROVENANCE_CANDIDATE:",final_src)
    print("DIAGNOSTIC_JSON:",p)
    print("MULTIRES_MOSAIC_CANDIDATE_QA_PASSED: True")
    print("EXISTING_INPUT_FILES_MODIFIED: FALSE")
    print("FINAL_INTEGRATED_DEM_PUBLISHED: FALSE")
    print("TERRAIN_DERIVATIVES_FROM_MOSAIC_APPROVED: FALSE")


if __name__=="__main__":
    main()
