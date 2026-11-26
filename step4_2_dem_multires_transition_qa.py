"""Read-only source-transition diagnostic for the proposed priority DEM.

Evaluates all cardinal, 1 m accounting-grid adjacencies between distinct
selected source codes (1, 2, 30). Native values are assigned by containing
cell, without interpolation; no mosaic, smoothing or slope is produced.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
import math
import tempfile

import numpy as np
import step4_2_dem_multires_coverage_qa as c
import step4_2_dem_multires_overlap_qa as o

METHOD = "MULTIRES_PRIORITY_SOURCE_TRANSITION_QA_V1"
OVERLAP = c.DIAG / "dem_multires_overlap_qa_cp_t07sl" / "multires_overlap_qa.json"
COVERAGE = o.COVERAGE
COVERAGE_SHA = "60b726149550cc3f50e01cad1f0e146228afcf0b40b20b931a48a423f4e7bf18"
SOURCE_KEYS = {1:"hrdem1m", 2:"hrdem2m", 30:"mrdem30_uncut"}
PAIRS = ((1,2),(1,30),(2,30))
THRESHOLDS = (0.1,0.5,1,2,5,10,20,50,100)
EXPECTED_COUNTS = {1:464090223,2:698372,30:578828062}
EXPECTED_BOUNDS = (475000,5463000,508000,5557000)


def native_values(ds, bounds):
    """Map each 1 m accounting centre to containing native value, not interpolation."""
    xmin,ymin,xmax,ymax = bounds
    sx,symin,sxmax,sy = c.raster_bounds(ds)
    r = int(ds.GetGeoTransform()[1])
    out = np.full((ymax-ymin,xmax-xmin),c.NODATA,dtype=np.float32)
    left,right = max(xmin,sx),min(xmax,sxmax)
    top,bottom = min(ymax,sy),max(ymin,symin)
    if left>=right or bottom>=top:
        return out
    cx0 = (left-sx)//r
    cx1 = (right-sx+r-1)//r
    cy0 = (sy-top)//r
    cy1 = (sy-bottom+r-1)//r
    band = ds.GetRasterBand(1)
    a = band.ReadAsArray(cx0,cy0,cx1-cx0,cy1-cy0)
    c.require(a.shape==(cy1-cy0,cx1-cx0) and a.dtype==np.float32 and np.isfinite(a).all(),
              "Invalid native raster window")
    if band.GetMaskFlags()!=c.gdal.GMF_NODATA:
        m=band.GetMaskBand().ReadAsArray(cx0,cy0,cx1-cx0,cy1-cy0)!=0
        c.require(np.array_equal(m,a!=c.NODATA),"Native masks disagree")
    if r!=1:
        a=np.repeat(np.repeat(a,r,axis=0),r,axis=1)
    dx,dy=left-(sx+cx0*r),(sy-cy0*r)-top
    out[ymax-top:ymax-bottom,left-xmin:right-xmin]=a[dy:dy+top-bottom,dx:dx+right-left]
    return out


def choose(values,aoi):
    label=np.zeros(aoi.shape,np.uint8)
    for code in (30,2,1):
        label[aoi & (values[code]!=c.NODATA)]=code
    return label


class Metrics:
    def __init__(self):
        self.n=0
        self.total=self.abs_total=self.sq_total=0.0
        self.lo=math.inf
        self.hi=-math.inf
        self.over=np.zeros(len(THRESHOLDS),np.int64)
        self.peak=None

    def add(self,d,sample):
        if d.size==0:
            return
        c.require(np.isfinite(d).all(),"Nonfinite differences")
        ad=np.abs(d)
        self.n+=int(d.size)
        self.total+=float(d.sum(dtype=np.float64))
        self.abs_total+=float(ad.sum(dtype=np.float64))
        self.sq_total+=float(np.dot(d,d))
        self.lo=min(self.lo,float(d.min()))
        self.hi=max(self.hi,float(d.max()))
        self.over+=np.array([np.count_nonzero(ad>t) for t in THRESHOLDS],dtype=np.int64)
        k=int(np.argmax(ad))
        if self.peak is None or float(ad[k])>self.peak["absolute_difference_m"]:
            self.peak={**sample(k),"signed_difference_m":float(d[k]),"absolute_difference_m":float(ad[k])}

    def report(self):
        if not self.n:
            return {"count":0}
        def interval(p):
            rank=math.ceil(p*self.n)
            low=0
            for t,n in zip(THRESHOLDS,self.over):
                if self.n-int(n)>=rank:
                    return [low,t]
                low=t
            return [low,self.peak["absolute_difference_m"]]
        return {
            "count":self.n,"mean_signed_m":self.total/self.n,
            "mae_m":self.abs_total/self.n,"rmse_m":math.sqrt(self.sq_total/self.n),
            "minimum_signed_m":self.lo,"maximum_signed_m":self.hi,
            "abs_p50_interval_m":interval(.50),"abs_p95_interval_m":interval(.95),
            "abs_p99_interval_m":interval(.99),
            "counts_over_m":{str(t):int(n) for t,n in zip(THRESHOLDS,self.over)},
            "percent_over_m":{str(t):100*int(n)/self.n for t,n in zip(THRESHOLDS,self.over)},
            "maximum_absolute":self.peak,
        }


class Pair:
    def __init__(self):
        self.edges=0
        self.directions={"east":0,"south":0}
        self.coarse_context=0
        self.same_native_coarse_cell=0
        self.jump=Metrics()
        self.jump_with_context=Metrics()
        self.coarse_own_step=Metrics()
        self.source_offset_at_fine=Metrics()

    def report(self):
        c.require(self.edges==sum(self.directions.values())==self.jump.n,"Transition edge totals mismatch")
        c.require(self.coarse_context==self.jump_with_context.n==self.coarse_own_step.n
                  ==self.source_offset_at_fine.n,"Context metric count mismatch")
        return {
            "transition_edges":self.edges,
            "accounting_grid_edge_length_km":self.edges/1000,
            "direction_counts":self.directions,
            "coarse_context_available_edges":self.coarse_context,
            "coarse_context_missing_edges":self.edges-self.coarse_context,
            "same_native_coarse_cell_edges":self.same_native_coarse_cell,
            "jump_all":self.jump.report(),
            "jump_on_coarse_context_support":self.jump_with_context.report(),
            "coarse_native_step_on_coarse_context_support":self.coarse_own_step.report(),
            "cross_source_offset_at_fine_endpoint":self.source_offset_at_fine.report(),
        }


def process_edges(label,values,xmin,top,direction,states,grids):
    """Count each east/south undirected edge once, oriented fine -> coarse."""
    if direction=="east":
        a,b=label[:,:-1],label[:,1:]
        dr,dc=0,1
    else:
        a,b=label[:-1,:],label[1:,:]
        dr,dc=1,0
    rr,cc=np.nonzero((a!=b)&(a!=0)&(b!=0))
    av,bv=a[rr,cc],b[rr,cc]
    c.require(np.isin(av,[1,2,30]).all() and np.isin(bv,[1,2,30]).all(),
              "Unexpected source code")
    for fid,cid in PAIRS:
        hits=((av==fid)&(bv==cid))|((av==cid)&(bv==fid))
        if not hits.any():
            continue
        r,col=rr[hits],cc[hits]
        forward=av[hits]==fid
        fr=r+(~forward)*dr
        fc=col+(~forward)*dc
        cr=r+forward*dr
        cc2=col+forward*dc
        fine=values[fid][fr,fc].astype(np.float64)
        coarse=values[cid][cr,cc2].astype(np.float64)
        cf=values[cid][fr,fc].astype(np.float64)
        c.require(np.all(fine!=c.NODATA) and np.all(coarse!=c.NODATA),
                  "Selected transition endpoint is NoData")
        c.require(np.all(values[fid][cr,cc2]==c.NODATA),"Fine source priority violated")
        jump=coarse-fine
        valid=cf!=c.NODATA
        n=int(jump.size)
        st=states[(fid,cid)]
        st.edges+=n
        st.directions[direction]+=n
        st.coarse_context+=int(valid.sum())
        gx,_,_,gy=grids[cid]
        fine_x=xmin+fc+0.5
        fine_y=top-fr-0.5
        coarse_x=xmin+cc2+0.5
        coarse_y=top-cr-0.5
        same=((np.floor((fine_x-gx)/cid)==np.floor((coarse_x-gx)/cid))
              &(np.floor((gy-fine_y)/cid)==np.floor((gy-coarse_y)/cid)))
        c.require(np.all(~same|valid),"One native cell has inconsistent validity")
        st.same_native_coarse_cell+=int(same.sum())

        def sample(k):
            return {
                "direction":direction,
                "fine_center_epsg26910":[float(fine_x[k]),float(fine_y[k])],
                "coarse_center_epsg26910":[float(coarse_x[k]),float(coarse_y[k])],
                "fine_source_resolution_m":fid,
                "coarse_source_resolution_m":cid,
                "fine_elevation_m":float(fine[k]),
                "coarse_elevation_m":float(coarse[k]),
                "coarse_at_fine_endpoint_m":float(cf[k]) if valid[k] else None,
                "edge_center_epsg26910":[float((fine_x[k]+coarse_x[k])/2),
                                         float((fine_y[k]+coarse_y[k])/2)],
            }
        st.jump.add(jump,sample)
        vi=np.flatnonzero(valid)
        def context_sample(k):
            return sample(int(vi[k]))
        own=coarse[valid]-cf[valid]
        offset=cf[valid]-fine[valid]
        c.require(np.allclose(jump[valid],own+offset,atol=1e-9,rtol=0),
                  "Transition decomposition does not reconcile")
        st.jump_with_context.add(jump[valid],context_sample)
        st.coarse_own_step.add(own,context_sample)
        st.source_offset_at_fine.add(offset,context_sample)


def scan(rasters,layer,bounds,strip_rows=c.STRIP_ROWS,progress=True):
    xmin,ymin,xmax,ymax=map(int,bounds)
    width,height=xmax-xmin,ymax-ymin
    grids={k:c.raster_bounds(rasters[v]) for k,v in SOURCE_KEYS.items()}
    states={p:Pair() for p in PAIRS}
    counts={1:0,2:0,30:0}
    aoi_count=uncovered=0
    prev_lab=prev_values=None
    strips=math.ceil(height/strip_rows)
    for i,top in enumerate(range(ymax,ymin,-strip_rows),1):
        bottom=max(ymin,top-strip_rows)
        bb=(xmin,bottom,xmax,top)
        aoi=c.mask_grid(layer,(xmin,1,0,top,0,-1),width,top-bottom)
        values={k:native_values(rasters[v],bb) for k,v in SOURCE_KEYS.items()}
        label=choose(values,aoi)
        aoi_count+=int(aoi.sum())
        uncovered+=int((aoi&(label==0)).sum())
        c.require(not np.any(~aoi&(label!=0)),"Selection outside AOI")
        for k in counts:
            counts[k]+=int((label==k).sum())
        process_edges(label,values,xmin,top,"east",states,grids)
        process_edges(label,values,xmin,top,"south",states,grids)
        if prev_lab is not None:
            # Include the one set of north-south edges across strip breaks.
            seam_lab=np.concatenate((prev_lab,label[:1]),axis=0)
            seam_vals={k:np.concatenate((prev_values[k],values[k][:1]),axis=0)
                       for k in SOURCE_KEYS}
            process_edges(seam_lab,seam_vals,xmin,top+1,"south",states,grids)
        prev_lab=label[-1:].copy()
        prev_values={k:v[-1:].copy() for k,v in values.items()}
        if progress and (i%10==0 or i==strips):
            print(f"TRANSITION_QA_PROGRESS: {i}/{strips}",flush=True)
    c.require(sum(counts.values())+uncovered==aoi_count,"Selection inventory mismatch")
    return counts,aoi_count,uncovered,{"_to_".join(map(str,p)):st.report() for p,st in states.items()}


def main():
    c.init()
    paths={**c.INPUTS,"coverage_qa":COVERAGE,"overlap_qa":OVERLAP,
           "coverage_script":Path(c.__file__).resolve(),"overlap_script":Path(o.__file__).resolve()}
    for k,p in paths.items():
        c.require(p.is_file(),f"Required input missing: {k}: {p}")
    print("HASHING_INPUTS",flush=True)
    hashes={k:c.digest(p) for k,p in paths.items()}
    cq=json.loads(COVERAGE.read_text(encoding="utf-8"))
    oq=json.loads(OVERLAP.read_text(encoding="utf-8"))
    c.require(hashes["coverage_qa"]==COVERAGE_SHA,"Unexpected coverage QA bytes")
    c.require(cq["coverage_qa_completed"] is True and cq["input_files_unchanged"] is True
              and cq["full_coverage_achieved_uncut"] is True
              and tuple(cq["accounting_bounds"])==EXPECTED_BOUNDS
              and cq["script_sha256"]==hashes["coverage_script"],"Coverage gate failed")
    c.require(oq["method_id"]==o.METHOD and oq["overlap_diagnostic_completed"] is True
              and oq["input_files_unchanged"] is True
              and oq["script_sha256"]==hashes["overlap_script"],"Overlap gate failed")
    c.require(oq["input_sha256"]=={k:hashes[k] for k in (*c.INPUTS,"coverage_qa","coverage_script")}
              and cq["input_sha256"]=={k:hashes[k] for k in c.INPUTS},
              "Input hashes do not match QA chain")
    c.require(oq["comparison_2m"]["statistics"]["strict_common_cells"]==115963824
              and oq["comparison_30m"]["statistics"][0]["strict_common_cells"]==508727,
              "Unexpected overlap inventory")
    mem,layer,envelope,area,target=c.projected_aoi(c.INPUTS["aoi"])
    rasters=c.open_rasters(hashes,target)
    print("INPUT_QA_CHAIN_VERIFIED: True",flush=True)
    counts,aoi,missing,pairs=scan(rasters,layer,EXPECTED_BOUNDS)
    c.require(counts==EXPECTED_COUNTS and aoi==1043616657 and missing==0,
              "Priority source inventory differs from coverage QA")
    c.require(all(p["transition_edges"]>0 for p in pairs.values()),"Unexpected missing transition pair")
    rasters.clear()
    layer=None
    mem=None
    print("RECHECKING_INPUT_HASHES",flush=True)
    c.require(hashes=={k:c.digest(p) for k,p in paths.items()},"Input files changed")
    report={
        "method_id":METHOD,"script_sha256":c.digest(Path(__file__)),
        "created_utc":datetime.now(timezone.utc).isoformat(),
        "crs":"EPSG:26910","accounting_grid_resolution_m":1,
        "accounting_bounds":EXPECTED_BOUNDS,
        "input_sha256":hashes,"input_paths":{k:p.relative_to(c.ROOT).as_posix() for k,p in paths.items()},
        "selection_rule":"AOI centre, valid HRDEM1m > HRDEM2m > MRDEM30m uncut; containing cell",
        "aoi_accounting_pixels":aoi,"uncovered_accounting_pixels":missing,
        "selected_source_pixels":counts,
        "total_transition_edges":sum(p["transition_edges"] for p in pairs.values()),
        "transition_pairs":pairs,
        "metrics_definition":{
            "jump":"coarse selected elevation minus fine selected elevation at adjacent accounting-centre endpoints",
            "coarse_native_step":"coarse source at coarse endpoint minus that same source at fine endpoint",
            "source_offset":"coarse source at fine endpoint minus fine selected elevation at fine endpoint",
            "decomposition":"jump = coarse_native_step + source_offset on coarse-context-supported edges",
            "edge_length":"sum of cardinal 1 m shared cell edges; not vector perimeter length or independent samples",
        },
        "coverage_inventory_reconciled":True,"input_files_unchanged":True,
        "source_transition_diagnostic_completed":True,
        "cross_source_elevation_agreement_accepted":False,
        "source_transitions_accepted_for_terrain_derivatives":False,
        "absolute_vertical_accuracy_validated":False,
        "final_integrated_dem_created":False,
        "limitations":[
            "Adjacent centres are 1 m apart, but this is an accounting grid, not uniformly resolved 1 m terrain.",
            "The cross-source jump includes natural local relief, coarse-cell placement and source differences; it is not itself a measurement error.",
            "Coarse-context offset uses nearest containing native cell, not an interpolated matched point or block mean.",
            "Samples are spatially dependent, especially multiple edges referring to the same 30 m native cell.",
            "Source priority precludes the fine source at the coarse endpoint; the counterfactual fine gradient is unavailable.",
            "Diagonal edges are not counted. No correction, smoothing, slope, or hazard model is produced.",
            "Thresholds are descriptive, not pass/fail or independently justified hazard tolerances.",
        ],
    }
    out=Path(tempfile.mkdtemp(prefix="dem_multires_transition_qa_",dir=c.DIAG))
    p=out/"multires_transition_qa.json"
    p.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("MULTIRES_TRANSITION_SUMMARY:",json.dumps(report,indent=2,allow_nan=False))
    print("DIAGNOSTIC_JSON:",p)
    print("SOURCE_TRANSITION_QA_COMPLETED: True")
    print("EXISTING_INPUT_FILES_MODIFIED: FALSE")
    print("FINAL_INTEGRATED_DEM: FALSE")


if __name__=="__main__":
    main()
