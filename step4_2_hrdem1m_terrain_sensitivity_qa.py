"""Read-only terrain derivative sensitivity audit on masked HRDEM, with real tile halos.
Outputs JSON only. No interpolation, warp, edge extrapolation, or raster modification.
"""
from pathlib import Path
from functools import lru_cache
import json
import tempfile
import numpy as np
import step4_2_hrdem1m_fixed_batch_300 as e
import step4_2_hrdem1m_global_masking_qa as q
import step4_2_hrdem1m_targeted_source_qa as t

PRIOR = 'hrdem1m_targeted_source_qa_eobrsvvs/targeted_source_qa.json'
N = 1000
SLOPE_THRESHOLDS = (15, 30, 45, 60)
DIFF_THRESHOLDS = (1, 5, 10, 20)
ASPECT_MIN_SLOPE = 5.0


def derivatives(a):
    # Require the complete 3x3 neighborhood. Array rows increase southward.
    a = np.asarray(a, dtype=np.float64)
    v = np.isfinite(a) & (a != q.NODATA)
    h,w = a.shape
    ok = np.ones((h-2,w-2), dtype=bool)
    for r in range(3):
        for c in range(3):
            ok &= v[r:r+h-2,c:c+w-2]
    nw,n,ne = a[:-2,:-2],a[:-2,1:-1],a[:-2,2:]
    west,z,east = a[1:-1,:-2],a[1:-1,1:-1],a[1:-1,2:]
    sw,s,se = a[2:,:-2],a[2:,1:-1],a[2:,2:]
    dx = ((ne-nw)+2*(east-west)+(se-sw))/8.0
    dy = ((nw-sw)+2*(n-s)+(ne-se))/8.0
    slope = np.degrees(np.arctan(np.hypot(dx,dy)))
    aspect = np.degrees(np.arctan2(-dx,-dy)) % 360.0
    # Five-point elevation Laplacian, m/m^2, NOT plan/profile curvature.
    lap = (east-z)+(west-z)+(n-z)+(s-z)
    return ok,slope,aspect,lap


def metric(a,b,thresholds,quantiles=False):
    d = b-a; ad = np.abs(d); n = int(d.size)
    return dict(n=n, signed_sum=float(d.sum()), abs_sum=float(ad.sum()),
                squared_sum=float(np.dot(d,d)), max_abs=float(ad.max()) if n else None,
                counts_over={str(x):int((ad>x).sum()) for x in thresholds},
                p95_abs=float(np.percentile(ad,95)) if n and quantiles else None)


def finalize(r):
    n=r['n']
    return dict(**r,mean_signed=r['signed_sum']/n if n else None,
                mae=r['abs_sum']/n if n else None,rmse=(r['squared_sum']/n)**.5 if n else None)


def compare(a,b,selection=None,quantiles=False):
    va,sa,aa,la=derivatives(a); vb,sb,ab,lb=derivatives(b)
    q.require(np.array_equal(va,vb),'Derivative support mismatch')
    centers=(a[1:-1,1:-1]!=q.NODATA)&np.isfinite(a[1:-1,1:-1])
    if selection is not None:
        va &= selection; centers &= selection
    av,bv=sa[va],sb[va]
    asp=va&(sa>=ASPECT_MIN_SLOPE)&(sb>=ASPECT_MIN_SLOPE)
    circular=(ab[asp]-aa[asp]+180)%360-180
    # Signed aspect difference is the shortest angular displacement.
    return dict(valid_centers=int(centers.sum()), strict_3x3_centers=int(va.sum()),
        slope_deg=finalize(metric(av,bv,DIFF_THRESHOLDS,quantiles)),
        aspect_deg=finalize(metric(np.zeros(circular.size),circular,(5,15,30,90),quantiles)),
        laplacian_per_m=finalize(metric(la[va],lb[va],(.1,.5,1,5),quantiles)),
        slope_thresholds={str(k):dict(original_ge=int((av>=k).sum()),fixed_ge=int((bv>=k).sum()),
            upward=int(((av<k)&(bv>=k)).sum()),downward=int(((av>=k)&(bv<k)).sum())) for k in SLOPE_THRESHOLDS})


def aggregate(rows):
    out={k:sum(r[k] for r in rows) for k in ('valid_centers','strict_3x3_centers')}
    for name in ('slope_deg','aspect_deg','laplacian_per_m'):
        items=[r[name] for r in rows]; maxima=[v['max_abs'] for v in items if v['max_abs'] is not None]
        total={k:sum(v[k] for v in items) for k in ('n','signed_sum','abs_sum','squared_sum')}
        total.update(max_abs=max(maxima,default=None),p95_abs=None,
            counts_over={k:sum(v['counts_over'][k] for v in items) for k in items[0]['counts_over']})
        out[name]=finalize(total)
    out['slope_thresholds']={str(k):{n:sum(r['slope_thresholds'][str(k)][n] for r in rows)
        for n in ('original_ge','fixed_ge','upward','downward')} for k in SLOPE_THRESHOLDS}
    return out


def patch(read,bbox,n=N):
    # Assemble by geographic position, including real adjacent tile rows and corners.
    xmin,ymin,xmax,ymax=bbox
    out=[np.full((ymax-ymin,xmax-xmin),q.NODATA,np.float32) for _ in range(2)]
    for y in range(ymin//n*n,ymax,n):
        for x in range(xmin//n*n,xmax,n):
            pair=read((x,y))
            if pair is None: continue
            l,b,r,t=max(x,xmin),max(y,ymin),min(x+n,xmax),min(y+n,ymax)
            dest=(slice(ymax-t,ymax-b),slice(l-xmin,r-xmin))
            orig=(slice(y+n-t,y+n-b),slice(l-x,r-x))
            for a,src in zip(out,pair): a[dest]=src[orig]
    return out


def main():
    e.legacy.gdal.UseExceptions(); e.legacy.ogr.UseExceptions()
    root=e.legacy.project_root(); coverage,_,srcdir,maskdir=e.paths(root)
    diag=srcdir.parent; cfg=e.config(root)
    pp=diag/PRIOR; prior=e.read_json(pp)
    q.require(prior['configuration']==cfg and prior['script_sha256']==e.sha(Path(t.__file__))
        and prior['summary']['all_sites_numerically_consistent'] is True,'Targeted source QA mismatch')
    lp=root/prior['prior_local_report']; q.require(e.sha(lp)==prior['prior_local_report_sha256'],'Local report hash mismatch')
    local=e.read_json(lp); sp=root/local['prior_seam_report']
    q.require(e.sha(sp)==local['prior_seam_report_sha256'],'Seam report hash mismatch')
    seam=e.read_json(sp); mp=root/seam['prior_masking_qa']
    q.require(e.sha(mp)==seam['prior_masking_qa_sha256'],'Mask report hash mismatch')
    mask=e.read_json(mp)
    q.require(mask['configuration']==cfg and mask['summary']['masking_qa_passed'] is True,'Mask QA mismatch')
    lookup={(x,y):g for x,y,g in e.legacy.candidate_tiles(coverage)}
    records={tuple(r['tile']):r for r in mask['tiles']}
    q.require(len(lookup)==656 and set(lookup)==set(records),'Candidate inventory mismatch')
    olddir=diag/'hrdem1m_masked_staging'
    stems={e.stem(*k) for k,r in records.items() if r['state']=='COMPLETE'}
    for folder in (olddir,maskdir):
        inv=q.inventory(folder,stems,masked=True)
        q.require(not inv['missing_entries'] and not inv['unexpected_entries'],f'Inventory mismatch {folder}')
    verified={}
    @lru_cache(maxsize=48)
    def read(k):
        if k not in records or records[k]['state']=='EMPTY': return None
        x,y=k; base=e.stem(x,y); fp=maskdir/(base+'_masked.tif'); op=olddir/(base+'_masked.tif')
        side=e.read_json(maskdir/(base+'_masked.qa.json')); prov=side['provenance']
        q.require(prov['extraction_configuration']==cfg and prov['mask_method_id']==mask['mask_method_id']
            and all(prov.get(n)==v for n,v in mask['mask_script_hashes'].items()),f'Provenance mismatch {k}')
        if k not in verified:
            fixed_hash=e.sha(fp); q.require(side['output_sha256']==fixed_hash,f'Fixed hash mismatch {k}')
            verified[k]=dict(tile=list(k),fixed_sha256=fixed_hash,original_current_sha256=e.sha(op))
        a,_=q.read_grid(op,x,y); b,_=q.read_grid(fp,x,y)
        va=a!=q.NODATA; vb=b!=q.NODATA
        q.require(np.array_equal(va,vb) and int(vb.sum())==records[k]['masked_valid_pixels'],f'Validity mismatch {k}')
        return a,b
    rows=[]
    for i,k in enumerate(sorted(lookup,key=lambda p:(-p[1],p[0])),1):
        x,y=k
        a,b=patch(read,(x-1,y-1,x+N+1,y+N+1))
        r=compare(a,b); r.update(tile=list(k),source_group=lookup[k]); rows.append(r)
        if i%25==0 or i==len(lookup): print(f'TERRAIN_PROGRESS: {i}/{len(lookup)}',flush=True)
    summary=aggregate(rows)
    q.require(summary['valid_centers']==mask['summary']['masked_valid_pixels'] and len(verified)==655,'Global center count mismatch')
    summary['valid_centers_excluded_by_3x3']=summary['valid_centers']-summary['strict_3x3_centers']
    print('===== GLOBAL TERRAIN SENSITIVITY ====='); print(json.dumps(summary,indent=2),flush=True)
    sites=[]
    for site in local['sites']:
        x,y=site['tile']; d=site['direction']; bbox=t.seam_box(x,y,d)
        l,b,r,u=bbox; a,f=patch(read,(l-1,b-1,r+1,u+1))
        shape=(u-b,r-l); near=np.zeros(shape,bool)
        # Four center columns/rows straddle the seam: distances 0.5 and 1.5 m.
        if d=='east': near[:,shape[1]//2-2:shape[1]//2+2]=True
        else: near[shape[0]//2-2:shape[0]//2+2,:]=True
        rr=dict(kind=site['kind'],tile=[x,y],direction=d,bbox=list(bbox),
            center_band=compare(a,f,near,True),full_strip=compare(a,f,quantiles=True))
        sites.append(rr)
        print('SEAM_TERRAIN:',site['kind'],(x,y,d),'strict_centers:',rr['center_band']['strict_3x3_centers'],
            'slope_p95_abs_deg:',rr['center_band']['slope_deg']['p95_abs'],
            'slope_max_abs_deg:',rr['center_band']['slope_deg']['max_abs'],flush=True)
    hot=[]
    for site in prior['sites']:
        if 'hotspot' not in site: continue
        l,b,r,u=site['bbox']; a,f=patch(read,(l-1,b-1,r+1,u+1))
        rr=dict(label=site['label'],bbox=site['bbox'],metrics=compare(a,f,quantiles=True))
        hot.append(rr); print('HOTSPOT_TERRAIN:',json.dumps(rr),flush=True)
    q.require(len(sites)==82 and len(hot)==2,'Local site count mismatch')
    # Site strips can overlap; report site-level distributions, not unique-pixel totals.
    ss={}
    for kind in ('flagged','control'):
        group=[r['center_band']['slope_deg'] for r in sites if r['kind']==kind]
        p95=[r['p95_abs'] for r in group if r['n']>=20]
        ss[kind]=dict(sites=len(group),sites_with_ge20_strict_centers=len(p95),
            median_site_slope_p95_abs_deg=float(np.median(p95)) if p95 else None,
            sites_with_p95_slope_change_over_5deg=sum(v>5 for v in p95))
    print('===== LOCAL TERRAIN SUMMARY ====='); print(json.dumps(ss,indent=2),flush=True)
    report=Path(tempfile.mkdtemp(prefix='hrdem1m_terrain_sensitivity_',dir=diag))/'terrain_sensitivity_qa.json'
    e.save_json(report,dict(configuration=cfg,prior_targeted_source_qa=str(pp.relative_to(root)),
        prior_targeted_source_qa_sha256=e.sha(pp),script_sha256=e.sha(Path(__file__)),
        summary=summary,site_summary=ss,tiles=rows,seam_sites=sites,hotspots=hot,verified_hashes=list(verified.values()),
        methods=dict(slope='Horn 3x3, 1m spacing, degrees',aspect='downslope clockwise from grid north, shortest circular difference, both slopes >=5 degrees',
            laplacian='five-point elevation Laplacian, inverse meters, strict 3x3 support',slope_thresholds_deg=list(SLOPE_THRESHOLDS)),
        limitations=['NoData and outer edges excluded; actual neighboring pixels are used across tile boundaries.',
            'Slope thresholds are illustrative sensitivity bins, not landslide thresholds. Laplacian is not plan/profile curvature.',
            'Site strips may overlap and are purposive. No inferential test or independent vertical truth.',
            'Original files have no historical hashes; current hashes recorded. Differences do not establish physical accuracy.',
            'This is not factor-of-safety or hazard validation.'], final_dem_approved=False))
    print('DIAGNOSTIC_JSON:',report); print('TERRAIN_SENSITIVITY_QA_COMPLETED: True'); print('FINAL_DEM: FALSE')

if __name__=='__main__': main()
