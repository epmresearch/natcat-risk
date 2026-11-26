"""Local parallel-line controls for fixed HRDEM seams. No warp or raster writes."""
from pathlib import Path
from collections import Counter
from functools import lru_cache
import json
import tempfile
import numpy as np
import step4_2_hrdem1m_global_fixed_seam_qa as s
import step4_2_hrdem1m_fixed_batch_300 as e
import step4_2_hrdem1m_global_masking_qa as q
import step4_2_hrdem1m_fixed_mask_batch_300 as m

PRIOR='hrdem1m_global_fixed_seam_qa_9rfgcq9o/global_fixed_seam_qa.json'
OFFSETS=(-50,-20,-10,-5,-2,-1,0,1,2,5,10,20,50)
W=52

def strip(a,b,direction):
    # Rows are along the seam; columns increase eastward or northward.
    if direction=='east':
        return np.concatenate((a[:,-W:],b[:,:W]),axis=1)
    if direction=='north':
        return np.concatenate((a[:W,:][::-1,:].T,b[-W:,:][::-1,:].T),axis=1)
    raise ValueError(direction)

def scan(z):
    vals={}; masks={}
    for o in OFFSETS:
        v=z[:,W+o-2:W+o+2]
        masks[o]=(np.isfinite(v)&(v!=s.NODATA)).all(axis=1)
        vals[o]=s.residual(v)
    central=masks[0]
    joint=np.logical_and.reduce(list(masks.values()))
    rows=[]
    for o in OFFSETS:
        both=central&masks[o]
        rows.append(dict(offset_m=o,all_valid=s.describe(vals[o][masks[o]]),
                         paired_center=s.describe(vals[0][both]),
                         paired_offset=s.describe(vals[o][both]),
                         joint_support=s.describe(vals[o][joint])))
    # The full common support allows like-for-like comparisons of all lines.
    j=int(joint.sum()); center=s.describe(vals[0][joint]); distant=[r for r in rows if abs(r['offset_m'])>=5]
    peak=max((r['joint_support']['p95_abs_m'] for r in distant),default=None) if j else None
    ratio=(center['p95_abs_m']/peak) if j>=20 and peak is not None and peak>0 else None
    return dict(center=s.describe(vals[0][central]),joint_crossings=j,
                joint_fraction_of_center=j/int(central.sum()) if central.any() else None,
                joint_center_p95_m=center['p95_abs_m'],joint_max_offset_p95_m=peak,
                joint_center_to_max_offset_p95_ratio=ratio,offsets=rows)

def key(r): return (*r['tile'],r['direction'])

def controls(flags, pool):
    # Unique nearest controls; same orientation/group and enough central crossings.
    used=set(); result=[]
    for r in sorted(flags,key=key):
        x,y,d=key(r); minimum=min(500,r['common_crossings'])
        eligible=[c for c in pool if key(c) not in used and c['direction']==d
                  and c['source_groups']==r['source_groups'] and c['common_crossings']>=minimum]
        if not eligible:
            result.append((r,None)); continue
        c=min(eligible,key=lambda c:((c['tile'][0]-x)**2+(c['tile'][1]-y)**2,key(c)))
        used.add(key(c)); result.append((r,c))
    return result

def main():
    m.legacy.gdal.UseExceptions(); m.legacy.ogr.UseExceptions()
    root=e.legacy.project_root(); coverage,_,_,folder=e.paths(root); diag=folder.parent
    prior_path=diag/PRIOR; prior=e.read_json(prior_path); cfg=e.config(root)
    q.require(prior['configuration']==cfg and prior['summary']['seam_audit_completed'] is True,'Prior seam QA/configuration mismatch')
    q.require(prior['audit_script_sha256']==e.sha(Path(s.__file__)),'Seam QA script changed')
    mask_path=root/prior['prior_masking_qa']; mask=e.read_json(mask_path)
    q.require(e.sha(mask_path)==prior['prior_masking_qa_sha256'] and mask['summary']['masking_qa_passed'] is True,'Masking report mismatch')
    plan=e.legacy.candidate_tiles(coverage); lookup={(x,y):g for x,y,g in plan}
    records={tuple(r['tile']):r for r in mask['tiles']}
    q.require(len(lookup)==656 and set(records)==set(lookup),'Candidate inventory mismatch')
    rows=prior['seams']; flags=[r for r in rows if r['fixed']['flagged']]
    pool=[r for r in rows if r['fixed']['assessable'] and not r['fixed']['flagged'] and not r['old']['flagged']]
    q.require(len(rows)==1127 and len(flags)==41 and len(pool)==1034,'Unexpected seam counts')
    pairs=controls(flags,pool); missing=[list(key(r)) for r,c in pairs if c is None]
    verified=set()
    @lru_cache(maxsize=8)
    def read(k):
        x,y=k; base=e.stem(x,y); fp=folder/(base+'_masked.tif'); side=e.read_json(folder/(base+'_masked.qa.json'))
        prov=side['provenance']
        q.require(prov['extraction_configuration']==cfg and prov['mask_method_id']==m.MASK_ID and all(prov.get(n)==v for n,v in mask['mask_script_hashes'].items()),f'Provenance mismatch {k}')
        if k not in verified:
            q.require(side['output_sha256']==e.sha(fp),f'Hash mismatch {k}'); verified.add(k)
        a,_=q.read_grid(fp,x,y)
        q.require(int((a!=s.NODATA).sum())==records[k]['masked_valid_pixels'],f'Valid count mismatch {k}')
        return a
    results=[]
    for r,c in pairs:
        for kind,row in [('flagged',r),('control',c)]:
            if row is None: continue
            k=tuple(row['tile']); other=tuple(row['neighbor']); d=row['direction']
            q.require(other==((k[0]+1000,k[1]) if d=='east' else (k[0],k[1]+1000)),f'Neighbor mismatch {k}')
            result=scan(strip(read(k),read(other),d))
            central=result['center']; old=row['fixed']
            q.require(central['crossings']==old['crossings'] and central['flagged']==old['flagged'] and all(abs(central[n]-old[n])<=1e-9 for n in ('p50_abs_m','p95_abs_m','max_abs_m','signed_mean_m')),f'Central residual mismatch {k} {d}')
            result.update(kind=kind,tile=list(k),neighbor=list(other),direction=d,source_groups=row['source_groups'],historical_flag=row['historical_flag'],paired_flag=list(key(r)))
            results.append(result)
            print('SITE:',kind,key(row),'center_p95:',central['p95_abs_m'],'joint_crossings:',result['joint_crossings'],'ratio:',result['joint_center_to_max_offset_p95_ratio'],flush=True)
    summary={}
    for kind in ('flagged','control'):
        rr=[r for r in results if r['kind']==kind]; ratios=[r['joint_center_to_max_offset_p95_ratio'] for r in rr if r['joint_center_to_max_offset_p95_ratio'] is not None and r['joint_fraction_of_center']>=.8]
        summary[kind]=dict(sites=len(rr),joint_support_ge20_and_ge80pct=sum(r['joint_crossings']>=20 and r['joint_fraction_of_center']>=.8 for r in rr),
                           comparable_nonzero_offset_ratios=len(ratios),center_p95_above_all_distant_offsets=sum(v>1 for v in ratios),center_p95_at_least_twice_all_distant_offsets=sum(v>=2 for v in ratios),median_ratio=float(np.median(ratios)) if ratios else None)
    summary.update(missing_controls=missing,unique_tiles_hash_verified=len(verified),central_metrics_reproduced=True,final_dem_approved=False)
    report=Path(tempfile.mkdtemp(prefix='hrdem1m_seam_local_context_',dir=diag))/'seam_local_context_qa.json'
    e.save_json(report,dict(configuration=cfg,prior_seam_report=str(prior_path.relative_to(root)),prior_seam_report_sha256=e.sha(prior_path),script_sha256=e.sha(Path(__file__)),offsets_m=list(OFFSETS),summary=summary,sites=results,limitations=['Local controls are not independent terrain truth. Nearest controls are selected greedily without replacement, not matched for terrain slope or land cover.','Ratios compare center P95 with the maximum P95 at offsets of at least 5 m using identical crossing rows. Summary ratios require >=20 joint rows and >=80% central support; thresholds are descriptive only.','A localized signal does not establish source artifact or natural origin. No raster is modified.']))
    print('===== LOCAL CONTEXT SUMMARY ====='); print(json.dumps(summary,indent=2)); print('DIAGNOSTIC_JSON:',report); print('LOCAL_CONTEXT_QA_COMPLETED: True'); print('FINAL_DEM: FALSE')

if __name__=='__main__': main()
