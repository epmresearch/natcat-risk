"""Compare all adjacent HRDEM tile seams using the historical four-pixel residual.
Read local rasters only; write one new diagnostic JSON. Run serially.
"""
from __future__ import annotations
from collections import Counter
from pathlib import Path
import json
import tempfile
import numpy as np
import step4_2_hrdem1m_fixed_batch_300 as e
import step4_2_hrdem1m_fixed_mask_batch_300 as m
import step4_2_hrdem1m_global_masking_qa as q

PRIOR = 'hrdem1m_global_masking_qa_t3koqz9w/global_masking_qa.json'
N = 1000
NODATA = -32767.0
HISTORICAL = {
    (485000,5490000,'east'), (484000,5488000,'east'),
    (483000,5503000,'east'), (485000,5489000,'north'),
    (485000,5474000,'north'), (483000,5467000,'north'),
    (484000,5505000,'east'), (483000,5470000,'east'),
    (485000,5472000,'east'), (485000,5475000,'north'),
    (484000,5487000,'east'), (489000,5503000,'east'),
    (490000,5536000,'north'), (486000,5506000,'east'),
}


def edges(a):
    # Copies avoid retaining whole tile arrays in memory.
    return {'west': a[:,:2].copy(), 'east': a[:,-2:].copy(),
            'north': a[:2,:].copy(), 'south': a[-2:,:].copy()}


def stencil(a, b, direction):
    # Each row contains four elevations in west-to-east or south-to-north order.
    if direction == 'east':
        return np.concatenate((a['east'], b['west']), axis=1)
    if direction == 'north':
        return np.concatenate((a['north'][::-1,:].T, b['south'][::-1,:].T), axis=1)
    raise ValueError(direction)


def residual(z):
    z = z.astype(np.float64, copy=False)
    return (z[:,2]-z[:,1])-0.5*((z[:,1]-z[:,0])+(z[:,3]-z[:,2]))


def describe(r):
    n = int(r.size)
    if not n:
        return dict(crossings=0,p50_abs_m=None,p95_abs_m=None,max_abs_m=None,
                    signed_mean_m=None,flagged=False,assessable=False)
    a = np.abs(r)
    p50,p95 = map(float,np.percentile(a,[50,95]))
    mean = float(r.mean())
    return dict(crossings=n,p50_abs_m=p50,p95_abs_m=p95,max_abs_m=float(a.max()),
                signed_mean_m=mean,assessable=n>=20,
                flagged=n>=20 and (p95>0.5 or abs(mean)>0.2))


def compare(old, fixed):
    ov = np.isfinite(old) & (old!=NODATA)
    fv = np.isfinite(fixed) & (fixed!=NODATA)
    q.require(np.array_equal(ov,fv),'Seam stencil validity mismatch')
    ok = ov.all(axis=1)
    a,b = residual(old[ok]),residual(fixed[ok])
    d = b-a
    return dict(old=describe(a),fixed=describe(b),
                common_crossings=int(ok.sum()),
                max_abs_residual_change_m=float(np.abs(d).max()) if d.size else None)


def main():
    m.legacy.gdal.UseExceptions(); m.legacy.ogr.UseExceptions()
    root=e.legacy.project_root()
    coverage,_,_,fixed_dir=e.paths(root)
    diag=fixed_dir.parent; old_dir=diag/'hrdem1m_masked_staging'
    prior_path=diag/PRIOR; prior=e.read_json(prior_path)
    cfg=e.config(root)
    q.require(prior['configuration']==cfg and prior['summary']['masking_qa_passed'] is True,
              'Prior masking QA or extraction configuration mismatch')
    q.require(prior['mask_method_id']==m.MASK_ID,'Mask method mismatch')
    mask_hashes=dict(legacy_mask_sha256=e.sha(Path(m.legacy.__file__)),
                     fixed_mask_sha256=e.sha(Path(m.__file__)))
    q.require(prior['mask_script_hashes']==mask_hashes,'Masking script hashes changed')
    plan=e.legacy.candidate_tiles(coverage); lookup={(x,y):g for x,y,g in plan}
    recs={tuple(r['tile']):r for r in prior['tiles']}
    q.require(len(plan)==len(lookup)==len(recs)==656 and set(lookup)==set(recs),'Candidate inventory mismatch')
    empty={k for k,r in recs.items() if r['state']=='EMPTY'}
    q.require(empty=={(497000,5549000)},'Unexpected empty inventory')
    complete=set(lookup)-empty; stems={e.stem(*k) for k in complete}
    for folder in (old_dir,fixed_dir):
        inv=q.inventory(folder,stems,masked=True)
        q.require(not inv['missing_entries'] and not inv['unexpected_entries'],f'Inventory mismatch: {folder}')
        print('MASKED_INVENTORY:',folder.name,json.dumps(inv),flush=True)
    old_edges={}; fixed_edges={}; zero_retained=[]
    for i,(x,y,g) in enumerate(plan,1):
        k=(x,y); r=recs[k]
        q.require(r['source_group']==g,f'Source group mismatch: {k}')
        if k in empty:
            a=np.full((N,N),NODATA,dtype=np.float32)
            old_edges[k]=edges(a); fixed_edges[k]=edges(a)
        else:
            base=e.stem(x,y); fn=base+'_masked.tif'
            fp=fixed_dir/fn; side=e.read_json(fixed_dir/(base+'_masked.qa.json'))
            prov=side['provenance']
            q.require(prov['extraction_configuration']==cfg and prov['mask_method_id']==m.MASK_ID
                      and all(prov.get(n)==v for n,v in prior['mask_script_hashes'].items())
                      and side['output_sha256']==e.sha(fp),f'Fixed masked integrity mismatch: {k}')
            a,_=q.read_grid(old_dir/fn,x,y); b,_=q.read_grid(fp,x,y)
            va=a!=NODATA; vb=b!=NODATA
            q.require(np.array_equal(va,vb) and int(vb.sum())==r['masked_valid_pixels']
                      and side['masked_valid_pixels']==r['masked_valid_pixels'],f'Validity/count mismatch: {k}')
            if not vb.any():
                zero_retained.append(dict(tile=list(k),footprint_pixels=r['footprint_pixels'],
                    source_valid_pixels=r['source_valid_pixels'],nodata_inside_footprint=r['nodata_inside_footprint']))
            old_edges[k]=edges(a); fixed_edges[k]=edges(b)
        if i%50==0 or i==len(plan): print(f'EDGE_READ_PROGRESS: {i}/{len(plan)}',flush=True)
    rows=[]; seen=set()
    for x,y in lookup:
        k=(x,y)
        for direction,other in (('east',(x+N,y)),('north',(x,y+N))):
            if other not in lookup: continue
            key=(x,y,direction); seen.add(key)
            row=compare(stencil(old_edges[k],old_edges[other],direction),
                        stencil(fixed_edges[k],fixed_edges[other],direction))
            rows.append(dict(tile=[x,y],neighbor=list(other),direction=direction,
                source_groups=[lookup[k],lookup[other]],both_actual_tiffs=k not in empty and other not in empty,
                historical_flag=key in HISTORICAL,**row))
    q.require(HISTORICAL<=seen,'Historical flagged pair absent')
    def key(r): return (*r['tile'],r['direction'])
    old_flags={key(r) for r in rows if r['old']['flagged']}
    new_flags={key(r) for r in rows if r['fixed']['flagged']}
    transitions=Counter(('flagged' if r['old']['flagged'] else 'not_flagged')+'_to_'+
                        ('flagged' if r['fixed']['flagged'] else 'not_flagged')
                        for r in rows if r['old']['assessable'])
    summary=dict(candidate_pairs=len(rows),actual_tiff_pairs=sum(r['both_actual_tiffs'] for r in rows),
        east_pairs=sum(r['direction']=='east' for r in rows),north_pairs=sum(r['direction']=='north' for r in rows),
        assessable_pairs=sum(r['old']['assessable'] for r in rows),
        zero_crossing_pairs=sum(r['common_crossings']==0 for r in rows),
        low_crossing_pairs_1_to_19=sum(0<r['common_crossings']<20 for r in rows),
        common_valid_crossings=sum(r['common_crossings'] for r in rows),
        old_flagged=len(old_flags),fixed_flagged=len(new_flags),
        persistent_flags=len(old_flags & new_flags),cleared_flags=len(old_flags-new_flags),
        newly_flagged=len(new_flags-old_flags),historical_flags_reproduced=old_flags==HISTORICAL,
        historical_flags_not_reproduced=sorted(HISTORICAL-old_flags),
        extra_old_flags=sorted(old_flags-HISTORICAL),transitions=dict(transitions),
        old_vs_fixed_validity_mismatch=0,zero_retained_source_tiles=zero_retained,
        seam_audit_completed=True,final_dem_approved=False)
    print('===== GLOBAL SEAM SUMMARY ====='); print(json.dumps(summary,indent=2))
    print('===== HISTORICAL 14 SEAMS =====')
    for r in rows:
        if r['historical_flag']: print(json.dumps(r))
    print('===== NEWLY FLAGGED SEAMS =====')
    for r in rows:
        if key(r) in new_flags-old_flags: print(json.dumps(r))
    print('===== TOP 10 FIXED SEAMS BY P95 (ASSESSABLE) =====')
    top=sorted((r for r in rows if r['fixed']['assessable']),key=lambda r:r['fixed']['p95_abs_m'],reverse=True)[:10]
    for r in top: print(json.dumps(r))
    report=Path(tempfile.mkdtemp(prefix='hrdem1m_global_fixed_seam_qa_',dir=diag))/'global_fixed_seam_qa.json'
    e.save_json(report,dict(configuration=cfg,prior_masking_qa=str(prior_path.relative_to(root)),
        prior_masking_qa_sha256=e.sha(prior_path),audit_script_sha256=e.sha(Path(__file__)),
        criteria=dict(minimum_crossings=20,p95_abs_m_strictly_above=0.5,abs_signed_mean_m_strictly_above=0.2),
        summary=summary,seams=rows,top_fixed_p95=top,
        limitations=['Four valid pixels required; only adjacent candidate pairs are tested.',
                     'Flags indicate local residuals, not verified artifacts or absolute elevation errors.',
                     'A cleared flag is not proof of accuracy; a new flag is not proof of a processing seam.',
                     'Original rasters have no historical hashes. No source data or final mosaic is modified.']))
    print('DIAGNOSTIC_JSON:',report)
    print('GLOBAL_SEAM_AUDIT_COMPLETED: True'); print('FINAL_DEM: FALSE')

if __name__=='__main__':
    main()
