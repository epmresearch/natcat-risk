"""Independent scan of the written multiresolution reference mosaic.

Verifies the candidate bytes, AOI, source inventory and transition statistics.
Writes a diagnostic register of all cardinal source edges with |jump| > 10 m.
The threshold records discrepancies, not hazard or accuracy acceptance.
No published raster is changed or promoted.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
import math
import tempfile

import numpy as np
from osgeo import ogr, osr
import step4_2_dem_multires_build_candidate as b

c=b.c
t=b.t
METHOD="MULTIRES_REFERENCE_FINAL_TECHNICAL_VALIDATION_V1"
CANDIDATE=c.DIAG/"dem_multires_mosaic_candidate_7ed8t26a"
QA=CANDIDATE/b.QA_NAME
EXPECTED_BUILD_SHA="e2ee6dd8261811d7cfb8dfd1a9748f7a8f98e4d72d5e0e4a2df6855d79c65d46"
EXPECTED_TRANSITION_SHA="5a1f067715ee190d9487b9b0d2be50b5d20b47e0152fecc1f6e72063e8dfb1f5"
EXPECTED_FILES={
    "reference_dem":(b.DEM_NAME,"54ab737d96fe940e70cb607ea75899a19830ec4fbce92f1bf6784974d0765063",1054828661),
    "source_provenance":(b.SOURCE_NAME,"775e42da007c5a5d6b23e043610c37d3f9a2644d74ffa35524b1724e435b8393",9956761),
}
EXPECTED_COUNTS={0:2058383343,1:464090223,2:698372,30:578828062}
FLAG_THRESHOLD=10.0
REGISTER_NAME="dem_highway99_multires_transition_register.gpkg"
LAYER="source_transition_jump_gt10m"
EXPECTED_FLAGS={"1_to_2":43,"1_to_30":2754,"2_to_30":48078}


def check_inputs():
    paths,hashes,transition=b.qa_inputs()
    c.require(hashes["transition_qa"]==EXPECTED_TRANSITION_SHA,"Transition QA hash changed")
    paths={**paths,"candidate_qa":QA,"candidate_script":Path(b.__file__).resolve(),
           **{key:CANDIDATE/v[0] for key,v in EXPECTED_FILES.items()}}
    for k,p in paths.items():
        c.require(p.is_file(),f"Missing file: {k}: {p}")
    for k in ("candidate_qa","candidate_script",*EXPECTED_FILES):
        hashes[k]=c.digest(paths[k])
    q=json.loads(QA.read_text(encoding="utf-8"))
    c.require(hashes["candidate_script"]==q["script_sha256"]==EXPECTED_BUILD_SHA,
              "Build script hash mismatch")
    c.require(q["method_id"]==b.METHOD and q["mosaic_candidate_qa_passed"] is True
              and q["all_grid_pixels_readback_verified"] is True
              and q["source_priority_exactly_verified"] is True
              and q["native_selected_elevation_bitwise_equal"] is True
              and q["input_files_unchanged"] is True
              and q["final_integrated_dem_published"] is False,"Candidate QA gate failed")
    c.require(q["source_transitions_accepted_for_terrain_derivatives"] is False
              and q["cross_source_elevation_agreement_accepted"] is False
              and q["absolute_vertical_accuracy_validated"] is False,"Scientific limitations changed")
    c.require(q["input_sha256"]=={k:hashes[k] for k in q["input_paths"]}
              and set(q["input_paths"])==set(paths)-{"candidate_qa","candidate_script",*EXPECTED_FILES},
              "Candidate lineage mismatch")
    c.require(tuple(q["bounds_epsg26910"])==b.BOUNDS
              and {int(k):v for k,v in q["source_pixel_counts"].items()}==EXPECTED_COUNTS
              and q["source_transition_edges"]==772181,
              "Candidate inventory mismatch")
    for key,(name,sha,size) in EXPECTED_FILES.items():
        p=paths[key]
        c.require(hashes[key]==q["files"][key]["sha256"]==sha
                  and p.stat().st_size==q["files"][key]["bytes"]==size
                  and q["files"][key]["path"]==p.relative_to(c.ROOT).as_posix(),
                  f"Candidate raster changed: {key}")
    print("CANDIDATE_AND_UPSTREAM_HASHES_VERIFIED: True",flush=True)
    return paths,hashes,q,transition


def create_register(path,srs):
    ds=ogr.GetDriverByName("GPKG").CreateDataSource(str(path))
    layer=ds.CreateLayer(LAYER,srs,ogr.wkbLineString)
    for name,typ in (("fine_code",ogr.OFTInteger),("coarse_code",ogr.OFTInteger),
                     ("direction",ogr.OFTString),("fine_z_m",ogr.OFTReal),
                     ("coarse_z_m",ogr.OFTReal),("signed_jump_m",ogr.OFTReal),
                     ("abs_jump_m",ogr.OFTReal)):
        c.require(layer.CreateField(ogr.FieldDefn(name,typ))==0,"Could not create register field")
    c.require(ds.StartTransaction()==0,"Could not start register transaction")
    return ds,layer


def add_flag(layer,fi,co,direction,x,y,zf,zc):
    # x,y are the centre of the actual shared 1 m grid edge.
    g=ogr.Geometry(ogr.wkbLineString)
    if direction=="east":
        g.AddPoint_2D(x,y-0.5); g.AddPoint_2D(x,y+0.5)
    else:
        g.AddPoint_2D(x-0.5,y); g.AddPoint_2D(x+0.5,y)
    f=ogr.Feature(layer.GetLayerDefn())
    for k,v in {"fine_code":int(fi),"coarse_code":int(co),"direction":direction,
                "fine_z_m":float(zf),"coarse_z_m":float(zc),
                "signed_jump_m":float(zc-zf),"abs_jump_m":float(abs(zc-zf))}.items():
        f.SetField(k,v)
    f.SetGeometry(g)
    c.require(layer.CreateFeature(f)==0,"Could not write flagged transition")


class EdgeTotals:
    def __init__(self):
        self.n=0
        self.east=self.south=0
        self.sum=self.abs_sum=self.sq_sum=0.0
        self.lo=math.inf;self.hi=-math.inf
        self.thresholds={str(v):0 for v in t.THRESHOLDS}

    def add(self,d,direction):
        if not d.size:
            return
        self.n+=int(d.size)
        if direction=="east": self.east+=int(d.size)
        else: self.south+=int(d.size)
        a=np.abs(d)
        self.sum+=float(d.sum(dtype=np.float64))
        self.abs_sum+=float(a.sum(dtype=np.float64))
        self.sq_sum+=float(np.dot(d,d))
        self.lo=min(self.lo,float(d.min()))
        self.hi=max(self.hi,float(d.max()))
        for value in t.THRESHOLDS:
            self.thresholds[str(value)]+=int(np.count_nonzero(a>value))

    def result(self):
        c.require(self.n>0,"Empty transition pair")
        return {"count":self.n,"direction_counts":{"east":self.east,"south":self.south},
                "mean_signed_m":self.sum/self.n,"mae_m":self.abs_sum/self.n,
                "rmse_m":math.sqrt(self.sq_sum/self.n),
                "minimum_signed_m":self.lo,"maximum_signed_m":self.hi,
                "counts_over_m":self.thresholds}


def edges(z,lab,x0,top,direction,totals,register):
    if direction=="east":
        a,d=lab[:,:-1],lab[:,1:]
        za,zd=z[:,:-1],z[:,1:]
    else:
        a,d=lab[:-1,:],lab[1:,:]
        za,zd=z[:-1,:],z[1:,:]
    rr,cc=np.nonzero((a!=d)&(a!=0)&(d!=0))
    av,dv=a[rr,cc],d[rr,cc]
    for fi,co in t.PAIRS:
        m=((av==fi)&(dv==co))|((av==co)&(dv==fi))
        if not m.any():
            continue
        r,col=rr[m],cc[m]
        forward=av[m]==fi
        a_z=za[r,col].astype(np.float64)
        d_z=zd[r,col].astype(np.float64)
        fine=np.where(forward,a_z,d_z)
        coarse=np.where(forward,d_z,a_z)
        delta=coarse-fine
        totals[(fi,co)].add(delta,direction)
        for j in np.flatnonzero(np.abs(delta)>FLAG_THRESHOLD):
            # Shared edge centre, independent of orientation of source priority.
            x=x0+int(col[j])+ (1.0 if direction=="east" else 0.5)
            y=top-int(r[j])- (0.5 if direction=="east" else 1.0)
            add_flag(register,fi,co,direction,x,y,fine[j],coarse[j])


def scan(dem,src,layer,bounds,register,strip_rows=c.STRIP_ROWS,progress=True):
    xmin,ymin,xmax,ymax=bounds
    w=xmax-xmin
    counts=dict.fromkeys(b.CODES,0)
    stats={k:{"count":0,"minimum_m":math.inf,"maximum_m":-math.inf,"negative_pixels":0}
           for k in (1,2,30)}
    totals={p:EdgeTotals() for p in t.PAIRS}
    prev_z=prev_lab=None
    strips=math.ceil((ymax-ymin)/strip_rows)
    for i,top in enumerate(range(ymax,ymin,-strip_rows),1):
        bottom=max(ymin,top-strip_rows)
        h=top-bottom
        aoi=c.mask_grid(layer,(xmin,1,0,top,0,-1),w,h)
        z=dem.GetRasterBand(1).ReadAsArray(0,ymax-top,w,h)
        lab=src.GetRasterBand(1).ReadAsArray(0,ymax-top,w,h)
        c.require(z.shape==lab.shape==aoi.shape and z.dtype==np.float32 and lab.dtype==np.uint8
                  and np.isfinite(z).all(),"Invalid candidate raster array")
        c.require(np.array_equal(z!=c.NODATA,aoi) and np.array_equal(lab!=0,aoi),
                  "Candidate AOI or data mask differs")
        c.require(np.array_equal(np.isin(lab,[1,2,30]),aoi),"Unknown source code")
        for k in counts:
            m=(lab==k)
            n=int(m.sum())
            counts[k]+=n
            if k and n:
                v=z[m]
                st=stats[k];st["count"]+=n
                st["minimum_m"]=min(st["minimum_m"],float(v.min()))
                st["maximum_m"]=max(st["maximum_m"],float(v.max()))
                st["negative_pixels"]+=int((v<0).sum())
        edges(z,lab,xmin,top,"east",totals,register)
        edges(z,lab,xmin,top,"south",totals,register)
        if prev_lab is not None:
            edges(np.concatenate([prev_z,z[:1]]),np.concatenate([prev_lab,lab[:1]]),
                  xmin,top+1,"south",totals,register)
        prev_z=z[-1:].copy();prev_lab=lab[-1:].copy()
        if progress and (i%10==0 or i==strips):
            print(f"FINAL_READ_PROGRESS: {i}/{strips}",flush=True)
    return counts,stats,{"_to_".join(map(str,p)):s.result() for p,s in totals.items()}


def verify_register(path,srs,expected):
    ds=ogr.Open(str(path),0)
    c.require(ds is not None and ds.GetLayerCount()==1,"Register open or layer count mismatch")
    l=ds.GetLayerByName(LAYER)
    cs=l.GetSpatialRef()
    cs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    c.require(cs.IsSame(srs),"Register CRS mismatch")
    counts=dict.fromkeys(expected,0)
    keys=set()
    for f in l:
        g=f.GetGeometryRef()
        c.require(g is not None and g.IsValid() and not g.IsEmpty() and g.GetPointCount()==2
                  and abs(g.Length()-1.0)<1e-12,"Incorrect transition edge geometry")
        key=f"{f['fine_code']}_to_{f['coarse_code']}"
        c.require(key in counts and f["abs_jump_m"]>FLAG_THRESHOLD
                  and abs(f["signed_jump_m"]-(f["coarse_z_m"]-f["fine_z_m"]))<1e-9
                  and abs(f["abs_jump_m"]-abs(f["signed_jump_m"]))<1e-9,
                  "Incorrect register attributes")
        x0,y0,*_=g.GetPoint(0)
        x1,y1,*_=g.GetPoint(1)
        direction=f["direction"]
        c.require((direction=="east" and x0==x1 and y1-y0==1) or
                  (direction=="south" and y0==y1 and x1-x0==1),"Incorrect edge orientation")
        unique=(x0,y0,x1,y1)
        c.require(unique not in keys,"Duplicate registered edge")
        keys.add(unique);counts[key]+=1
    c.require(counts==expected and l.GetFeatureCount()==sum(expected.values()),"Register counts mismatch")
    l=None;ds=None
    return counts


def compare_previous(counts,stats,actual,q,transition):
    c.require(counts=={int(k):v for k,v in q["source_pixel_counts"].items()},
              "Candidate readback counts differ")
    c.require({str(k):v for k,v in stats.items()}==q["source_elevation_statistics"],
              "Candidate readback statistics differ")
    for key,p in actual.items():
        previous=transition["transition_pairs"][key]
        jump=previous["jump_all"]
        c.require(p["count"]==previous["transition_edges"]==jump["count"]
                  and p["direction_counts"]==previous["direction_counts"]
                  and p["counts_over_m"]==jump["counts_over_m"],
                  f"Transition counts differ: {key}")
        for metric in ("mean_signed_m","mae_m","rmse_m","minimum_signed_m","maximum_signed_m"):
            c.require(math.isclose(p[metric],jump[metric],rel_tol=1e-10,abs_tol=1e-10),
                      f"Transition metric differs: {key}: {metric}")


def main():
    c.init()
    paths,hashes,q,transition=check_inputs()
    mem,layer,envelope,area,target=c.projected_aoi(c.INPUTS["aoi"])
    dem=c.gdal.Open(str(paths["reference_dem"]),c.gdal.GA_ReadOnly)
    src=c.gdal.Open(str(paths["source_provenance"]),c.gdal.GA_ReadOnly)
    b.inspect(dem,b.BOUNDS,target)
    b.inspect(src,b.BOUNDS,target,source=True)
    out=Path(tempfile.mkdtemp(prefix="dem_multires_final_validation_",dir=c.DIAG))
    print("OUTPUT_DIRECTORY:",out,flush=True)
    regpath=out/REGISTER_NAME
    reg,rl=create_register(regpath,target)
    counts,stats,actual=scan(dem,src,layer,b.BOUNDS,rl)
    c.require(reg.CommitTransaction()==0,"Register commit failed")
    rl=None;reg=None;dem=None;src=None
    compare_previous(counts,stats,actual,q,transition)
    c.require(counts==EXPECTED_COUNTS,"Final source counts differ from coverage QA")
    flags=verify_register(regpath,target,EXPECTED_FLAGS)
    print("TRANSITION_REGISTER_VERIFIED:",json.dumps(flags),flush=True)
    layer=None;mem=None
    print("RECHECKING_INPUT_HASHES",flush=True)
    c.require(hashes=={k:c.digest(p) for k,p in paths.items()},"Inputs changed during final validation")
    report={
        "method_id":METHOD,"script_sha256":c.digest(Path(__file__)),
        "created_utc":datetime.now(timezone.utc).isoformat(),
        "crs":"EPSG:26910","accounting_resolution_m":1,"bounds_epsg26910":b.BOUNDS,
        "total_grid_pixels":sum(counts.values()),"aoi_pixels":sum(counts[k] for k in (1,2,30)),
        "source_pixel_counts":counts,"source_elevation_statistics":stats,
        "uncovered_aoi_pixels":0,"valid_outside_aoi_pixels":0,"coverage_percent":100.0,
        "all_grid_pixels_scanned":True,"source_counts_and_statistics_reconciled":True,
        "transition_statistics_reproduced_from_candidate":True,
        "total_transition_edges":sum(p["count"] for p in actual.values()),
        "transition_pairs":actual,"register":{
            "path":regpath.relative_to(c.ROOT).as_posix(),"sha256":c.digest(regpath),
            "bytes":regpath.stat().st_size,"layer":LAYER,"threshold_abs_jump_m_strict_gt":FLAG_THRESHOLD,
            "counts":flags,"features":sum(flags.values()),"geometry_and_attributes_verified":True,
            "threshold_is_diagnostic_not_hazard_or_accuracy_criterion":True,
        },
        "input_paths":{k:p.relative_to(c.ROOT).as_posix() for k,p in paths.items()},
        "input_sha256":hashes,"input_files_unchanged":True,
        "final_technical_validation_passed":True,
        "scientific_use":"Reference elevation with mandatory source provenance only",
        "cross_source_elevation_agreement_accepted":False,
        "source_transitions_accepted_for_terrain_derivatives":False,
        "absolute_vertical_accuracy_validated":False,"hazard_model_validated":False,
        "final_integrated_dem_published":False,
        "limitations":q["limitations"]+[
            "The cross-source register contains only edges with absolute jump strictly greater than 10 m. Other edges are not certified safe.",
            "Flagged jumps may reflect terrain variation, source resolution and processing; they are not proven errors.",
            "Source-specific anomaly registers and transition diagnostics must accompany the reference DEM.",
        ],
    }
    p=out/"final_validation_qa.json"
    p.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("FINAL_TECHNICAL_VALIDATION_SUMMARY:",json.dumps(report,indent=2,allow_nan=False))
    print("TRANSITION_REGISTER:",regpath)
    print("DIAGNOSTIC_JSON:",p)
    print("FINAL_TECHNICAL_VALIDATION_PASSED: True")
    print("EXISTING_INPUT_FILES_MODIFIED: FALSE")
    print("FINAL_INTEGRATED_DEM_PUBLISHED: FALSE")


if __name__=="__main__":
    main()
