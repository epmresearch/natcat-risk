"""Audit all original and fixed-scale HRDEM staging tiles without changing rasters.
Run serially from the project root. Only a new diagnostic JSON is written.
"""
from __future__ import annotations

from collections import Counter
import json
from pathlib import Path
import tempfile

import numpy as np
import step4_2_hrdem1m_fixed_batch_300 as e

NODATA = -32767.0
EXPECTED_EMPTY = {(497000, 5549000)}
THRESHOLDS = (0.001, 0.01, 0.1, 0.5, 1.0, 5.0, 10.0)


def inventory(folder, expected):
    tifs = {p.stem for p in folder.glob('*_1km.tif')}
    qa = {p.name[:-8] for p in folder.glob('*_1km.qa.json')}
    allowed = {k + suffix for k in expected for suffix in ('.tif', '.qa.json')}
    extra = sorted(p.name for p in folder.iterdir() if p.name not in allowed)
    return dict(tif_count=len(tifs), qa_count=len(qa),
                missing_qa=sorted(expected-qa), unexpected_entries=extra,
                tif_without_qa=sorted(tifs-qa))


def load_array(path):
    ds = e.legacy.gdal.Open(str(path), e.legacy.gdal.GA_ReadOnly)
    z = ds.GetRasterBand(1).ReadAsArray()
    ds = None
    if z.shape != (1000, 1000):
        raise RuntimeError(f'Unexpected array dimensions: {path}')
    return z


def valid(z):
    return np.isfinite(z) & (z != NODATA)


def main():
    e.legacy.gdal.UseExceptions()
    e.legacy.ogr.UseExceptions()
    root = e.legacy.project_root()
    coverage, vrt, fixed, _ = e.paths(root)
    old = root / 'results/diagnostics/hrdem1m_batch_staging'
    plan = e.legacy.candidate_tiles(coverage)
    expected = {e.stem(x,y) for x,y,g in plan}
    if len(plan) != 656 or len(expected) != 656:
        raise RuntimeError('Candidate inventory changed')
    cfg = e.config(root)
    inventories = {name: inventory(folder,expected) for name,folder in [('original',old),('fixed',fixed)]}
    print('===== GLOBAL INVENTORY =====', flush=True)
    for name, inv in inventories.items():
        print(name.upper(), json.dumps(inv), flush=True)
        if inv['missing_qa'] or inv['unexpected_entries'] or inv['tif_without_qa']:
            raise RuntimeError(f'Inventory mismatch: {name}; investigate before masking')
    records = []
    old_states, new_states = Counter(), Counter()
    old_empty, new_empty = set(), set()
    old_valid_total = new_valid_total = common_total = 0
    old_only_total = new_only_total = nonfinite_total = 0
    signed_sum = abs_sum = squared_sum = 0.0
    max_abs = 0.0
    threshold_counts = {str(t): 0 for t in THRESHOLDS}
    print('===== FULL TILE QA =====', flush=True)
    for i,(x,y,g) in enumerate(plan,1):
        stem = e.stem(x,y)
        state, new_qa = e.state(root,x,y,g,cfg)
        if state not in ('COMPLETE','EMPTY'):
            raise RuntimeError(f'Unprocessed candidate: {stem}')
        new_states[state] += 1
        old_qa = e.read_json(old/(stem+'.qa.json'))
        old_path, new_path = old/(stem+'.tif'), fixed/(stem+'.tif')
        if (old_qa.get('tile_x'),old_qa.get('tile_y'),old_qa.get('source_group')) != (x,y,g):
            raise RuntimeError(f'Original QA identity mismatch: {stem}')
        if old_qa.get('status') == e.EMPTY:
            if old_path.exists() or old_qa.get('valid_pixels') != 0 or old_qa.get('nodata_pixels') != 1000000:
                raise RuntimeError(f'Invalid original empty marker: {stem}')
            old_state = 'EMPTY'
            old_empty.add((x,y))
            a = np.full((1000,1000), NODATA, dtype=np.float32)
        else:
            if not old_path.exists():
                raise RuntimeError(f'Original TIFF missing: {stem}')
            qa = e.legacy.inspect_tiff(old_path,x,y)
            if qa['valid_pixels'] == 0 or any(old_qa.get(k) != v for k,v in qa.items()):
                raise RuntimeError(f'Original structural QA mismatch: {stem}')
            old_state = 'COMPLETE'
            a = load_array(old_path)
        old_states[old_state] += 1
        if state == 'EMPTY':
            new_empty.add((x,y))
            b = np.full((1000,1000), NODATA, dtype=np.float32)
        else:
            b = load_array(new_path)
        va,vb = valid(a),valid(b)
        both = va & vb
        na,nb,n = int(va.sum()),int(vb.sum()),int(both.sum())
        oa,ob = int((va & ~vb).sum()),int((vb & ~va).sum())
        nf_a,nf_b = int((~np.isfinite(a)).sum()),int((~np.isfinite(b)).sum())
        if na != old_qa['valid_pixels'] or nb != new_qa['valid_pixels']:
            raise RuntimeError(f'Valid count mismatch: {stem}')
        d = b[both].astype(np.float64)-a[both].astype(np.float64)
        ad = np.abs(d)
        counts = {str(t):int((ad>t).sum()) for t in THRESHOLDS}
        s,sa,ss = float(d.sum()),float(ad.sum()),float(np.dot(d,d))
        p = [float(v) for v in np.percentile(ad,[50,95,100])] if n else [None,None,None]
        record = dict(tile=[x,y],source_group=g,original_state=old_state,fixed_state=state,
                      original_valid=na,fixed_valid=nb,common_valid=n,
                      original_only_valid=oa,fixed_only_valid=ob,validity_mismatch=oa+ob,
                      original_nonfinite=nf_a,fixed_nonfinite=nf_b,
                      abs_diff_p50_p95_max_m=p,mean_signed_fixed_minus_original_m=s/n if n else None,
                      mae_m=sa/n if n else None,rmse_m=(ss/n)**0.5 if n else None,
                      counts_over_m=counts)
        records.append(record)
        old_valid_total+=na; new_valid_total+=nb; common_total+=n
        old_only_total+=oa; new_only_total+=ob; nonfinite_total+=nf_a+nf_b
        signed_sum+=s; abs_sum+=sa; squared_sum+=ss
        if n: max_abs=max(max_abs,p[2])
        for k,v in counts.items(): threshold_counts[k]+=v
        if i%25==0 or i==len(plan): print(f'QA_PROGRESS: {i}/{len(plan)}',flush=True)
        if oa or ob or nf_a or nf_b or old_state!=state:
            print('REVIEW_TILE:',json.dumps(record),flush=True)
    grid_pixels=len(plan)*1000000
    integrity_ok = (new_states['COMPLETE']+new_states['EMPTY']==len(plan)
                    and inventories['fixed']['tif_count']==new_states['COMPLETE']
                    and inventories['original']['tif_count']==old_states['COMPLETE'])
    same_empty = old_empty==new_empty==EXPECTED_EMPTY
    same_validity = old_only_total==new_only_total==0
    review_gate = integrity_ok and same_empty and same_validity and nonfinite_total==0 and common_total>0
    summary = dict(candidates=len(plan),original_states=dict(old_states),fixed_states=dict(new_states),
                   original_empty_tiles=sorted(old_empty),fixed_empty_tiles=sorted(new_empty),
                   total_grid_pixels=grid_pixels,original_valid_pixels=old_valid_total,
                   fixed_valid_pixels=new_valid_total,original_nodata_or_invalid_pixels=grid_pixels-old_valid_total,
                   fixed_nodata_or_invalid_pixels=grid_pixels-new_valid_total,common_valid_pixels=common_total,
                   original_only_valid=old_only_total,fixed_only_valid=new_only_total,
                   validity_mismatch=old_only_total+new_only_total,nonfinite_pixels_across_both_versions=nonfinite_total,
                   tiles_with_validity_mismatch=sum(r['validity_mismatch']>0 for r in records),
                   tiles_with_elevation_difference_over_0_5m=sum(r['counts_over_m']['0.5']>0 for r in records),
                   mean_signed_fixed_minus_original_m=signed_sum/common_total if common_total else None,
                   mae_m=abs_sum/common_total if common_total else None,
                   rmse_m=(squared_sum/common_total)**0.5 if common_total else None,
                   maximum_abs_difference_m=max_abs if common_total else None,
                   counts_over_m=threshold_counts,inventory_and_fixed_provenance_passed=integrity_ok,
                   empty_inventory_unchanged=same_empty,validity_masks_equal=same_validity,
                   masking_review_gate_passed=review_gate)
    top = sorted((r for r in records if r['common_valid']),key=lambda r:r['abs_diff_p50_p95_max_m'][2],reverse=True)[:10]
    print('===== GLOBAL SUMMARY =====',flush=True)
    print(json.dumps(summary,indent=2),flush=True)
    print('===== TOP 10 TILES BY MAXIMUM ABSOLUTE DIFFERENCE =====',flush=True)
    for r in top: print(json.dumps(r),flush=True)
    # Write only a new report; never modify the original or fixed staging folders.
    report_dir=Path(tempfile.mkdtemp(prefix='hrdem1m_global_extraction_qa_',dir=fixed.parent))
    report=report_dir/'global_extraction_qa.json'
    payload=dict(configuration=cfg,inventories=inventories,summary=summary,tiles=records,
                 top_10_by_maximum_difference=top,
                 limitations=['Differences are fixed minus original, not vertical accuracy errors.',
                              'Per-tile percentiles are not pooled global percentiles.',
                              'Empty markers represent all-NoData grids, not independently reread source rasters.',
                              'The VRT hash does not authenticate remote COG contents.',
                              'Original rasters have no recorded historical SHA-256; fixed hashes are checked.',
                              'This audit does not establish absence of source artifacts or validate the final DEM.'])
    with report.open('x',encoding='utf-8') as f: json.dump(payload,f,indent=2,allow_nan=False); f.write('\n')
    print('DIAGNOSTIC_JSON:',report,flush=True)
    print('GLOBAL_EXTRACTION_QA_PASSED:',review_gate,flush=True)
    print('FINAL_DEM: FALSE',flush=True)
    if not review_gate: raise SystemExit(2)


if __name__=='__main__':
    main()
