"""Build a readback-verified 2 m mosaic candidate, without resampling.

Run in the natcat-risk project root after the anomaly-register QA.
All outputs are written to a new diagnostics directory.
"""
from pathlib import Path
import json
import shutil
import tempfile

import numpy as np
from osgeo import gdal, osr

import step4_2_hrdem2m_fixed_batch_300 as e
import step4_2_hrdem2m_fixed_mask_batch_300 as m

METHOD = "HRDEM2M_FIXED_MASKED_MOSAIC_V1"
MASK_QA = "hrdem2m_global_masking_qa_opn0z04b/global_masking_qa.json"
REGISTER_QA = "hrdem2m_anomaly_register_erxmut33/anomaly_register_qa.json"
EXPECTED_BOUNDS = (475000, 5463000, 506000, 5557000)
EXPECTED_STATES = {"COMPLETE": 669, "EMPTY": 2}
EXPECTED_VALID = 116196821
EXPECTED_NEGATIVE = 5045476
EXPECTED_FOOTPRINT = 116231518
EXPECTED_INSIDE_NODATA = 34697
EXPECTED_REGISTER_SHA256 = "e516be189e73f4b8c4b6c48d2905baad857d02fe7dd103fe0ddf09069b8d164d"
TILE_M = 1000
RES = 2
SIZE = TILE_M // RES
NODATA = e.NODATA
require = e.require


def inventory():
    return {str(p): e.digest(p) for d in (e.OUT,m.OUT)
            for p in d.iterdir() if p.is_file()}


def layout(bounds):
    xmin,ymin,xmax,ymax = bounds
    require(all(v % TILE_M == 0 for v in bounds) and xmin < xmax and ymin < ymax,
            "Invalid mosaic bounds")
    w,h = (xmax-xmin)//RES,(ymax-ymin)//RES
    return w,h,(xmin,RES,0,ymax,0,-RES)


def tile_array(k, records):
    r = records[k]
    if r["state"] == "EMPTY":
        require(r["masked_valid_pixels"] == 0 and not m.paths(k)[0].is_file(),
                f"Unexpected empty tile: {k}")
        return np.full((SIZE,SIZE),NODATA,dtype=np.float32)
    a,stats = e.inspect_tile(m.paths(k)[0], k)
    require(stats["valid_pixels"] == r["masked_valid_pixels"] and
            stats["nodata_pixels"] == r["masked_nodata_pixels"] and
            stats["min_elevation_m"] == r["min_elevation_m"] and
            stats["max_elevation_m"] == r["max_elevation_m"],
            f"Tile statistics changed: {k}")
    return a


def create_mosaic(path, bounds, records):
    xmin,ymin,xmax,ymax = bounds
    w,h,gt = layout(bounds)
    require(not path.exists(), f"Refusing to overwrite: {path}")
    ds = gdal.GetDriverByName("GTiff").Create(
        str(path),w,h,1,gdal.GDT_Float32,
        options=["TILED=YES","BLOCKXSIZE=256","BLOCKYSIZE=256",
                 "COMPRESS=DEFLATE","PREDICTOR=3","BIGTIFF=YES"])
    require(ds is not None, "Could not create mosaic")
    s = osr.SpatialReference()
    require(s.ImportFromEPSG(26910)==0,"Cannot resolve CRS")
    require(ds.SetGeoTransform(gt)==0 and ds.SetProjection(s.ExportToWkt())==0,
            "Cannot set georeferencing")
    ds.SetMetadataItem("AREA_OR_POINT","Area")
    band = ds.GetRasterBand(1)
    require(band.SetNoDataValue(NODATA)==0,"Cannot set NoData")
    rows = (ymax-ymin)//TILE_M
    written = set()
    try:
        for i,y in enumerate(range(ymax-TILE_M,ymin-1,-TILE_M),1):
            strip = np.full((SIZE,w),NODATA,dtype=np.float32)
            for x in range(xmin,xmax,TILE_M):
                k = (x,y)
                if k in records:
                    col = (x-xmin)//RES
                    strip[:,col:col+SIZE] = tile_array(k,records)
                    written.add(k)
            require(band.WriteArray(strip,0,(i-1)*SIZE)==0,
                    f"Raster write failed at row {i}")
            if i%10==0 or i==rows:
                print("MOSAIC_WRITE_PROGRESS:",f"{i}/{rows}",flush=True)
        require(written==set(records),"Not every candidate tile was written")
        band.FlushCache()
        ds.FlushCache()
    finally:
        band = None
        ds = None


def verify_mosaic(path,bounds,records):
    xmin,ymin,xmax,ymax = bounds
    w,h,gt=layout(bounds)
    ds=gdal.Open(str(path),gdal.GA_ReadOnly)
    require(ds is not None and (ds.RasterXSize,ds.RasterYSize,ds.RasterCount)==(w,h,1),
            "Wrong mosaic dimensions")
    require(ds.GetGeoTransform()==gt and ds.GetSpatialRef().GetAuthorityCode(None)=="26910",
            "Wrong mosaic georeference")
    band=ds.GetRasterBand(1)
    require(band.DataType==gdal.GDT_Float32 and band.GetNoDataValue()==NODATA,
            "Wrong dtype/NoData")
    valid_total=0
    negative_total=0
    min_elev=None
    max_elev=None
    mismatch=0
    invalid_gap=0
    gaps=0
    scanned=0
    checked=set()
    rows=(ymax-ymin)//TILE_M
    for i,y in enumerate(range(ymax-TILE_M,ymin-1,-TILE_M),1):
        # Read all bounding-box pixels, including gaps and empty tiles.
        yoff=(gt[3]-(y+TILE_M))//RES
        strip=band.ReadAsArray(0,yoff,w,SIZE)
        mask=band.GetMaskBand().ReadAsArray(0,yoff,w,SIZE)>0
        require(strip.shape==(SIZE,w) and np.isfinite(strip).all()
                and np.array_equal(mask,e.valid(strip)),"Invalid mosaic mask or nonfinite values")
        scanned+=strip.size
        n=int(mask.sum())
        valid_total+=n
        negative_total+=int((mask&(strip<0)).sum())
        if n:
            zmin,zmax=float(strip[mask].min()),float(strip[mask].max())
            min_elev=zmin if min_elev is None else min(min_elev,zmin)
            max_elev=zmax if max_elev is None else max(max_elev,zmax)
        for x in range(xmin,xmax,TILE_M):
            k=(x,y)
            xoff=(x-gt[0])//RES
            actual=strip[:,xoff:xoff+SIZE]
            if k in records:
                ref=tile_array(k,records)
                mismatch+=int(np.count_nonzero(actual.view(np.uint32)!=ref.view(np.uint32)))
                checked.add(k)
            else:
                gaps+=1
                invalid_gap+=int(np.count_nonzero(actual!=NODATA))
        if i%10==0 or i==rows:
            print("MOSAIC_READBACK_PROGRESS:",f"{i}/{rows}",flush=True)
    band=None
    ds=None
    require(scanned==w*h and checked==set(records),"Incomplete mosaic readback")
    require(gaps==((xmax-xmin)//TILE_M)*rows-len(records),"Gap inventory mismatch")
    require(mismatch==0 and invalid_gap==0,
            f"Readback mismatch={mismatch}, unexpected gap pixels={invalid_gap}")
    return {
        "width":w,"height":h,"bounds_epsg26910":list(bounds),"geotransform":list(gt),
        "crs":"EPSG:26910","dtype":"Float32","nodata":float(NODATA),
        "total_grid_pixels":w*h,"valid_pixels":valid_total,"nodata_pixels":w*h-valid_total,
        "negative_pixels":negative_total,"minimum_elevation_m":min_elev,"maximum_elevation_m":max_elev,
        "candidate_tiles":len(records),"noncandidate_cells":gaps,
        "input_tiffs":sum(r["state"]=="COMPLETE" for r in records.values()),
        "empty_extraction_cells":sum(r["state"]=="EMPTY" for r in records.values()),
        "all_nodata_input_tiffs":sum(r["state"]=="COMPLETE" and r["masked_valid_pixels"]==0 for r in records.values()),
        "input_tile_bitwise_mismatch_pixels":mismatch,
        "unexpected_non_nodata_in_gap_cells":invalid_gap,
        "all_grid_pixels_readback_verified":True,
    }


def main():
    gdal.UseExceptions()
    gdal.SetCacheMax(128*1024*1024)
    tiles,_,ch,mc=m.config()
    mh=e.canonical_hash(mc)
    bounds=tuple(e.load(e.PLAN)["bounds"])
    require(bounds==EXPECTED_BOUNDS,"Production-plan bounds mismatch")
    w,h,gt=layout(bounds)
    require((w,h)==(15500,47000) and len(tiles)==671,"Wrong planned mosaic grid")
    mask_path=e.DIAG/MASK_QA
    register_path=e.DIAG/REGISTER_QA
    register_gpkg=register_path.parent/"anomaly_register.gpkg"
    q=e.load(mask_path)
    reg=e.load(register_path)
    summary=q["summary"]
    require(summary["global_masking_qa_passed"] is True and
            summary["extraction_configuration_sha256"]==ch and
            summary["mask_configuration_sha256"]==mh,"Global masking QA not valid")
    require(reg["register_completed"] is True and reg["all_geometry_validated"] is True
            and reg["extraction_configuration_sha256"]==ch and reg["mask_configuration_sha256"]==mh
            and reg["total_features"]==1110,"Anomaly register QA not valid")
    require(e.digest(register_gpkg)==reg["register_sha256"]==EXPECTED_REGISTER_SHA256,
            "Anomaly register checksum mismatch")
    require(e.digest(mask_path)==reg["input_qa_sha256"]["masking"],"Mask QA hash mismatch")
    provenance={str(e.DIAG/v):reg["input_qa_sha256"][k] for k,v in reg["input_qa_paths"].items()}
    for f,sha in provenance.items():
        require(e.digest(Path(f))==sha,f"Registered QA changed: {f}")
    provenance[str(register_gpkg)]=e.digest(register_gpkg)
    provenance[str(register_path)]=e.digest(register_path)
    mr={tuple(r["tile"]):r for r in q["tiles"]}
    require(len(q["tiles"])==len(mr)==len(tiles)==671 and set(mr)==set(tiles),
            "Candidate inventory mismatch")
    require(all(xmin<=x<xmax and ymin<=y<ymax and x%TILE_M==0 and y%TILE_M==0
                for xmin,ymin,xmax,ymax in [bounds] for x,y in mr),"Candidate outside mosaic")
    states={s:sum(r["state"]==s for r in mr.values()) for s in EXPECTED_STATES}
    require(states==EXPECTED_STATES,"Unexpected tile states")
    require(summary["masked_valid_pixels"]==sum(r["masked_valid_pixels"] for r in mr.values())==EXPECTED_VALID
            and summary["footprint_pixels"]==EXPECTED_FOOTPRINT
            and summary["nodata_inside_footprint"]==EXPECTED_INSIDE_NODATA
            and summary["valid_outside_footprint"]==0,"Masked-pixel totals mismatch")
    before=inventory()
    for i,(k,r) in enumerate(mr.items(),1):
        tif,qa=m.paths(k)
        require(before[str(qa)]==r["mask_qa_sha256"],f"Mask QA hash mismatch: {k}")
        if r["state"]=="COMPLETE":
            require(before[str(tif)]==r["masked_raster_sha256"],f"Raster hash mismatch: {k}")
        else:
            require(not tif.exists(),f"Unexpected empty raster: {k}")
        if i%100==0 or i==len(mr):
            print("INPUT_HASH_PROGRESS:",f"{i}/{len(mr)}",flush=True)
    free=shutil.disk_usage(e.DIAG).free
    raw=w*h*4
    require(free>=2*raw+512*1024*1024,"Insufficient free space for safe mosaic creation")
    out=Path(tempfile.mkdtemp(prefix="hrdem2m_mosaic_candidate_",dir=e.DIAG))
    tmp=out/"hrdem2m_fixed_masked_mosaic_candidate.partial.tif"
    final=out/"hrdem2m_fixed_masked_mosaic_candidate.tif"
    print("OUTPUT_DIRECTORY:",out,flush=True)
    print("MOSAIC_PLAN:",json.dumps({"bounds":list(bounds),"width":w,"height":h,
          "grid_pixels":w*h,"candidate_cells":len(mr),"noncandidate_cells":w*h//(SIZE*SIZE)-len(mr),
          "uncompressed_bytes":raw,"free_bytes":free}),flush=True)
    create_mosaic(tmp,bounds,mr)
    results=verify_mosaic(tmp,bounds,mr)
    require(results["valid_pixels"]==EXPECTED_VALID and results["negative_pixels"]==EXPECTED_NEGATIVE
            and results["minimum_elevation_m"]==summary["minimum_elevation_m"]
            and results["maximum_elevation_m"]==summary["maximum_elevation_m"]
            and results["all_nodata_input_tiffs"]==0,"Mosaic statistics mismatch")
    require(before==inventory(),"Staging changed during mosaic creation")
    require(all(e.digest(Path(p))==sha for p,sha in provenance.items()),"QA inputs changed")
    require(not final.exists(),"Unexpected existing final candidate")
    tmp.rename(final)
    results.update({
        "method_id":METHOD,
        "script_sha256":e.digest(Path(__file__).resolve()),
        "extraction_configuration_sha256":ch,
        "mask_configuration_sha256":mh,
        "input_global_masking_qa":MASK_QA,
        "input_global_masking_qa_sha256":e.digest(mask_path),
        "input_anomaly_register_qa":REGISTER_QA,
        "input_anomaly_register_qa_sha256":e.digest(register_path),
        "register_sha256":EXPECTED_REGISTER_SHA256,
        "input_hashes_verified_before_and_after":True,
        "raster_bytes":final.stat().st_size,
        "mosaic_sha256":e.digest(final),
        "mosaic_candidate_qa_passed":True,
        "absolute_vertical_accuracy_validated":False,
        "negative_source_physical_validity_resolved":False,
        "final_dem_approved":False,
        "final_dem_published":False,
    })
    e.new_json(out/"mosaic_candidate_qa.json",results)
    print("MOSAIC_CANDIDATE_SUMMARY:",json.dumps(results,indent=2))
    print("MOSAIC_CANDIDATE:",final)
    print("DIAGNOSTIC_JSON:",out/"mosaic_candidate_qa.json")
    print("MOSAIC_CANDIDATE_QA_PASSED: True")
    print("FINAL_DEM: FALSE")


if __name__=="__main__":
    main()
