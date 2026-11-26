"""Independent full-grid technical validation of the HRDEM 2 m candidate.

The routine rasterizes the original footprint on the final grid, verifies every
tile's coverage and elevation statistics, and checks the anomaly register
against the candidate. It writes only a new validation report.
"""
from pathlib import Path
import json
import tempfile

import numpy as np
from osgeo import gdal, ogr, osr

import step4_2_hrdem2m_fixed_batch_300 as e
import step4_2_hrdem2m_fixed_mask_batch_300 as m
import step4_2_hrdem2m_build_mosaic_candidate as b
import step4_2_hrdem2m_anomaly_register as ar

METHOD = "HRDEM2M_FINAL_TECHNICAL_VALIDATION_V1"
CANDIDATE_DIR = "hrdem2m_mosaic_candidate_b2mzrrzw"
CANDIDATE_NAME = "hrdem2m_fixed_masked_mosaic_candidate.tif"
EXPECTED_MOSAIC_SHA256 = "5ce4ffa9e13af91a72d19a799ac8dce2dfe1d1b9dea184eeb3908670aa27d6ae"
EXPECTED_REGISTER_SHA256 = b.EXPECTED_REGISTER_SHA256
EXPECTED_EXTRACTION_SHA256 = "6e3e412e3d81794d69ce0ab7257445ae8a966ccbd0c6c9054cfa3e25442da30f"
EXPECTED_MASK_SHA256 = "d7e7d5608bc5a0b3ad0685ffe157cc5029561eed470845262f59d706d154e3c5"
THRESHOLDS = (-1,-5,-10,-20)
require = e.require


def register_inventory(path, expected, mr, negatives):
    ds = ogr.Open(str(path),0)
    require(ds is not None and ds.GetLayerCount()==len(expected),"Invalid register inventory")
    counts = {ly.GetName():ly.GetFeatureCount() for ly in ds}
    require(counts==expected,"Anomaly register layer counts mismatch")
    deep = {}
    for ly in ds:
        name=ly.GetName()
        require(ly.GetSpatialRef().GetAuthorityCode(None)=="26910",f"Invalid register CRS: {name}")
        seen_tiles=set()
        seen_seams=set()
        for feat in ly:
            geom=feat.GetGeometryRef()
            require(geom is not None and not geom.IsEmpty() and geom.IsValid()
                    and geom.GetGeometryType()==ly.GetGeomType(),f"Invalid register geometry: {name}")
            k=(feat.GetField("tile_x"),feat.GetField("tile_y"))
            require(k in mr,f"Unexpected tile in {name}")
            if name=="deep_negative_pixels":
                x,y=geom.GetX(),geom.GetY()
                require((x,y) not in deep and x==int(x) and y==int(y)
                        and k[0]<x<k[0]+1000 and k[1]<y<k[1]+1000
                        and (x-k[0])%2==1 and (y-k[1])%2==1,
                        "Invalid deep-negative pixel coordinate")
                value=feat.GetField("elevation_m")
                require(value<-5 and feat.GetField("source_fidelity_passed")==1
                        and feat.GetField("physical_validated")==0,"Invalid negative point status")
                deep[(int(x),int(y))]=float(value)
            elif name=="tiles_with_internal_nodata":
                require(k not in seen_tiles and geom.Equals(ar.square(k)), "Invalid NoData tile")
                seen_tiles.add(k);r=mr[k]
                require(feat.GetField("footprint_pixels")==r["footprint_pixels"]
                        and feat.GetField("valid_pixels")==r["masked_valid_pixels"]
                        and feat.GetField("nodata_pixels")==r["nodata_inside_footprint"],
                        "NoData register counts mismatch")
            elif name=="negative_elevation_tiles":
                require(k not in seen_tiles and k in negatives and geom.Equals(ar.square(k)),
                        "Invalid negative tile")
                seen_tiles.add(k);r=negatives[k]
                require(feat.GetField("negative_pixels")==r["negative_pixels"]
                        and feat.GetField("minimum_m")==r["minimum_elevation_m"]
                        and feat.GetField("physical_validated")==0,"Negative register statistics mismatch")
                for n in THRESHOLDS:
                    require(feat.GetField(f"below_minus_{abs(n)}")==r["counts_below_m"][str(n)],
                            "Negative threshold count mismatch")
            else:
                direction=feat.GetField("direction");ident=(k,direction)
                require(ident not in seen_seams and geom.Equals(ar.line(k,direction)),
                        "Invalid seam register geometry")
                seen_seams.add(ident)
                require(feat.GetField("physical_validated")==0,"Unexpected physical-validation claim")
                if name=="flagged_seams":
                    require(feat.GetField("bridge_passed")==feat.GetField("native_passed")==1
                            and feat.GetField("strict_crossings")>=20,"Invalid flagged seam status")
        if name=="tiles_with_internal_nodata":
            require(seen_tiles=={k for k,r in mr.items() if r["nodata_inside_footprint"]>0},
                    "NoData tile inventory mismatch")
        if name=="negative_elevation_tiles":
            require(seen_tiles==set(negatives),"Negative tile inventory mismatch")
    feat=geom=ly=None;ds=None
    print("REGISTER_LAYER_COUNTS:",json.dumps(counts),flush=True)
    return counts,deep


def scan(path, footprint, bounds, mr, negatives, deep):
    xmin,ymin,xmax,ymax=bounds
    w,h,gt=b.layout(bounds)
    ds=gdal.Open(str(path),gdal.GA_ReadOnly)
    require(ds is not None and (ds.RasterXSize,ds.RasterYSize,ds.RasterCount)==(w,h,1)
            and ds.GetGeoTransform()==gt and ds.GetSpatialRef().GetAuthorityCode(None)=="26910",
            "Mosaic grid mismatch")
    band=ds.GetRasterBand(1)
    require(band.DataType==gdal.GDT_Float32 and band.GetNoDataValue()==e.NODATA,
            "Mosaic data type or NoData mismatch")
    fpds=ogr.Open(str(footprint),0)
    require(fpds is not None and fpds.GetLayerCount()==1,"Invalid footprint")
    ly=fpds.GetLayer(0)
    require(ly.GetSpatialRef().GetAuthorityCode(None)=="26910","Wrong footprint CRS")
    stripe=gdal.GetDriverByName("MEM").Create("",w,b.SIZE,1,gdal.GDT_Byte)
    stripe.SetProjection(ds.GetProjection())
    sb=stripe.GetRasterBand(1)
    valid_count=fp_count=missing=outside=grid=negative=gap_fp=gap_valid=0
    threshold_counts={str(n):0 for n in THRESHOLDS}
    min_rec=max_rec=None
    checked=set();neg_seen=set();deep_seen={}
    rows=(ymax-ymin)//1000
    for index,y in enumerate(range(ymax-1000,ymin-1,-1000),1):
        yoff=(index-1)*b.SIZE
        a=band.ReadAsArray(0,yoff,w,b.SIZE)
        v=e.valid(a)
        require(np.isfinite(a).all() and np.array_equal(v,band.GetMaskBand().ReadAsArray(0,yoff,w,b.SIZE)>0),
                "Nonfinite pixels or invalid raster mask")
        require(stripe.SetGeoTransform((xmin,2,0,y+1000,0,-2))==0,"Stripe georeference error")
        sb.Fill(0)
        ly.ResetReading()
        require(gdal.RasterizeLayer(stripe,[1],ly,burn_values=[1],options=["ALL_TOUCHED=FALSE"])==0,
                "Footprint rasterization failed")
        fp=sb.ReadAsArray()==1
        grid+=a.size;valid_count+=int(v.sum());fp_count+=int(fp.sum())
        missing+=int((fp&~v).sum());outside+=int((v&~fp).sum())
        for x in range(xmin,xmax,1000):
            k=(x,y);j=(x-xmin)//2
            z=a[:,j:j+b.SIZE];zv=v[:,j:j+b.SIZE];zp=fp[:,j:j+b.SIZE]
            n=int(zv.sum());nf=int(zp.sum());nd=int((zp&~zv).sum())
            if k not in mr:
                gap_fp+=nf;gap_valid+=n
                require(n==nf==0 and np.all(z==e.NODATA),f"Unexpected data in gap {k}")
                continue
            r=mr[k];checked.add(k)
            require(n==r["masked_valid_pixels"] and nf==r["footprint_pixels"]
                    and nd==r["nodata_inside_footprint"] and not (zv&~zp).any()
                    and b.SIZE*b.SIZE-n==r["masked_nodata_pixels"],
                    f"Footprint or coverage statistics mismatch: {k}")
            neg=zv&(z<0);nn=int(neg.sum());negative+=nn
            require(nn==r["negative_pixels"],f"Negative count mismatch: {k}")
            counts={str(t):int((zv&(z<t)).sum()) for t in THRESHOLDS}
            for t,count in counts.items(): threshold_counts[t]+=count
            if n:
                mn=float(z[zv].min());mx=float(z[zv].max())
                require(mn==r["min_elevation_m"] and mx==r["max_elevation_m"],f"Elevation range mismatch: {k}")
                def extreme(value):
                    i,c=np.argwhere(zv&(z==value))[0]
                    return {"value_m":value,"center_epsg26910":[int(x+2*c+1),int(y+999-2*i)],"tile":list(k)}
                if min_rec is None or mn<min_rec["value_m"]: min_rec=extreme(mn)
                if max_rec is None or mx>max_rec["value_m"]: max_rec=extreme(mx)
            else:
                require(r["min_elevation_m"] is None and r["max_elevation_m"] is None,
                        f"Empty tile elevation mismatch: {k}")
            if nn:
                require(k in negatives,"Missing negative tile record")
                rneg=negatives[k];neg_seen.add(k)
                require(nn==rneg["negative_pixels"] and counts==rneg["counts_below_m"]
                        and mn==rneg["minimum_elevation_m"],f"Negative report mismatch: {k}")
                ax,ay=rneg["minimum_center_epsg26910"]
                require(x<ax<x+1000 and y<ay<y+1000 and (ax-x)%2==1 and (ay-y)%2==1
                        and float(z[(y+999-ay)//2,(ax-x-1)//2])==mn,
                        "Incorrect reported minimum location")
            else:
                require(k not in negatives,"Unexpected negative tile record")
            if counts["-5"]:
                rr,cc=np.where(zv&(z<-5))
                for i,c in zip(rr,cc):
                    pos=(int(x+2*c+1),int(y+999-2*i));value=float(z[i,c])
                    require(pos in deep and deep[pos]==value and pos not in deep_seen,
                            f"Deep negative point mismatch: {pos}")
                    deep_seen[pos]=value
        if index%10==0 or index==rows: print("FINAL_READ_PROGRESS:",f"{index}/{rows}",flush=True)
    require(grid==w*h and checked==set(mr) and neg_seen==set(negatives),"Incomplete grid scan")
    require(deep_seen==deep and len(deep_seen)==threshold_counts["-5"],"Incomplete deep negative validation")
    require(outside==gap_fp==gap_valid==0 and fp_count==valid_count+missing,
            "Unexpected footprint coverage")
    sb=stripe=ly=fpds=band=ds=None
    return {
        "grid_pixels":grid,"candidate_cells":len(checked),"valid_pixels":valid_count,
        "nodata_pixels":grid-valid_count,"footprint_pixels":fp_count,
        "nodata_inside_footprint":missing,"valid_outside_footprint":outside,
        "footprint_pixels_in_noncandidate_tiles":gap_fp,
        "valid_pixels_in_noncandidate_tiles":gap_valid,
        "footprint_coverage_percent":100*valid_count/fp_count,
        "minimum":min_rec,"maximum":max_rec,"negative_pixels":negative,
        "negative_tiles":len(neg_seen),"negative_counts_below_m":threshold_counts,
        "deep_negative_points_matching_mosaic":len(deep_seen),
        "all_pixels_scanned":True,
    }


def main():
    gdal.UseExceptions();ogr.UseExceptions();osr.UseExceptions()
    gdal.SetCacheMax(128*1024*1024)
    tiles,cfg,ch,mc=m.config();mh=e.canonical_hash(mc)
    require(ch==EXPECTED_EXTRACTION_SHA256 and mh==EXPECTED_MASK_SHA256,"Configuration mismatch")
    candidate=e.DIAG/CANDIDATE_DIR/CANDIDATE_NAME
    candidate_qa=candidate.parent/"mosaic_candidate_qa.json"
    register_qa=e.DIAG/b.REGISTER_QA
    register=register_qa.parent/"anomaly_register.gpkg"
    mask_qa=e.DIAG/b.MASK_QA
    c=e.load(candidate_qa);reg=e.load(register_qa);q=e.load(mask_qa)
    require(c["mosaic_candidate_qa_passed"] and c["all_grid_pixels_readback_verified"]
            and c["input_tile_bitwise_mismatch_pixels"]==0 and c["unexpected_non_nodata_in_gap_cells"]==0
            and c["extraction_configuration_sha256"]==ch and c["mask_configuration_sha256"]==mh,
            "Candidate QA gate not passed")
    require(e.digest(candidate)==c["mosaic_sha256"]==EXPECTED_MOSAIC_SHA256,"Candidate hash mismatch")
    require(reg["register_completed"] and reg["all_geometry_validated"] and reg["layer_counts"]==ar.EXPECTED
            and reg["extraction_configuration_sha256"]==ch and reg["mask_configuration_sha256"]==mh
            and e.digest(register)==reg["register_sha256"]==c["register_sha256"]==EXPECTED_REGISTER_SHA256,
            "Register QA gate not passed")
    require(e.digest(mask_qa)==c["input_global_masking_qa_sha256"]==reg["input_qa_sha256"]["masking"]
            and e.digest(register_qa)==c["input_anomaly_register_qa_sha256"],"QA input hash mismatch")
    inputs={e.DIAG/v:reg["input_qa_sha256"][k] for k,v in reg["input_qa_paths"].items()}
    inputs.update({candidate:c["mosaic_sha256"],register:reg["register_sha256"],
                   candidate_qa:e.digest(candidate_qa),register_qa:e.digest(register_qa),
                   e.FOOTPRINT:cfg["input_sha256"][e.FOOTPRINT.relative_to(e.ROOT).as_posix()],
                   e.VRT:cfg["input_sha256"][e.VRT.relative_to(e.ROOT).as_posix()],
                   Path(b.__file__).resolve():c["script_sha256"],
                   Path(ar.__file__).resolve():reg["script_sha256"],
                   Path(e.__file__).resolve():cfg["script_sha256"],
                   Path(m.__file__).resolve():mc["mask_script_sha256"]})
    require(all(e.digest(p)==h for p,h in inputs.items()),"Provenance hash mismatch")
    neg=e.load(e.DIAG/reg["input_qa_paths"]["negative"])
    require(neg["summary"]["negative_source_qa_passed"],"Negative-source QA not passed")
    mr={tuple(r["tile"]):r for r in q["tiles"]}
    negatives={tuple(r["tile"]):r for r in neg["negative_tiles"]}
    require(set(mr)==set(tiles) and len(mr)==len(q["tiles"])==671 and len(negatives)==119,
            "QA inventory mismatch")
    counts,deep=register_inventory(register,ar.EXPECTED,mr,negatives)
    bounds=tuple(e.load(e.PLAN)["bounds"])
    require(bounds==b.EXPECTED_BOUNDS and list(bounds)==c["bounds_epsg26910"],"Bounds mismatch")
    summary=scan(candidate,e.FOOTPRINT,bounds,mr,negatives,deep)
    ms=q["summary"];ns=neg["summary"]
    require(summary["grid_pixels"]==728500000 and summary["valid_pixels"]==c["valid_pixels"]==ms["masked_valid_pixels"]==116196821
            and summary["nodata_pixels"]==c["nodata_pixels"]==612303179
            and summary["footprint_pixels"]==ms["footprint_pixels"]==116231518
            and summary["nodata_inside_footprint"]==ms["nodata_inside_footprint"]==34697
            and summary["negative_pixels"]==ns["negative_pixels"]==5045476
            and summary["negative_counts_below_m"]==ns["negative_counts_below_m"]
            and summary["deep_negative_points_matching_mosaic"]==ns["deep_negative_pixels"]==606
            and summary["minimum"]["value_m"]==c["minimum_elevation_m"]
            and summary["maximum"]["value_m"]==c["maximum_elevation_m"],
            "Final totals mismatch")
    require(all(e.digest(p)==h for p,h in inputs.items()),"Input changed during validation")
    summary.update({
        "method_id":METHOD,
        "script_sha256":e.digest(Path(__file__).resolve()),
        "extraction_configuration_sha256":ch,"mask_configuration_sha256":mh,
        "crs":"EPSG:26910","resolution_m":2,
        "candidate_path":candidate.relative_to(e.ROOT).as_posix(),
        "candidate_qa_path":candidate_qa.relative_to(e.ROOT).as_posix(),
        "candidate_qa_sha256":inputs[candidate_qa],
        "register_path":register.relative_to(e.ROOT).as_posix(),
        "register_qa_path":register_qa.relative_to(e.ROOT).as_posix(),
        "register_qa_sha256":inputs[register_qa],
        "register_layer_counts":counts,
        "mosaic_sha256":EXPECTED_MOSAIC_SHA256,"register_sha256":EXPECTED_REGISTER_SHA256,
        "input_sha256":{p.relative_to(e.ROOT).as_posix():h for p,h in inputs.items()},
        "input_hashes_verified_before_and_after":True,
        "final_technical_validation_passed":True,
        "final_dem_published":False,
        "absolute_vertical_accuracy_validated":False,
        "negative_source_physical_validity_resolved":False,
        "hazard_model_validated":False,
    })
    out=Path(tempfile.mkdtemp(prefix="hrdem2m_final_validation_",dir=e.DIAG))
    report=out/"final_validation_qa.json"
    e.new_json(report,summary)
    print("FINAL_TECHNICAL_VALIDATION_SUMMARY:",json.dumps(summary,indent=2))
    print("DIAGNOSTIC_JSON:",report)
    print("FINAL_TECHNICAL_VALIDATION_PASSED: True")
    print("FINAL_DEM_PUBLISHED: FALSE")


if __name__=="__main__":
    main()
