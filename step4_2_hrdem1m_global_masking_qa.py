"""Read-only audit of the full HRDEM masking stage; write a new JSON report only.
Run serially in the original project root, with the production scripts unchanged.
"""
from __future__ import annotations

from collections import Counter
import json
from pathlib import Path
import tempfile

import numpy as np
import step4_2_hrdem1m_fixed_batch_300 as e
import step4_2_hrdem1m_fixed_mask_batch_300 as m

L = m.legacy
N = 1000
NODATA = L.NODATA
EMPTY_TILE = (497000, 5549000)
THRESHOLDS = (0.001, 0.01, 0.1, 0.5, 1.0, 5.0, 10.0)
PRIOR = 'hrdem1m_global_extraction_qa_v0n9af73/global_extraction_qa.json'


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def inventory(folder, stems, masked=False):
    suffixes = ('_masked.tif', '_masked.qa.json') if masked else ('.tif', '.qa.json')
    expected = {s + suffix for s in stems for suffix in suffixes}
    if not masked:
        expected.add(e.stem(*EMPTY_TILE) + '.qa.json')
    actual = {p.name for p in folder.iterdir() if p.is_file()}
    extra = sorted({p.name for p in folder.iterdir()} - expected)
    return dict(tif_count=sum(n.endswith('.tif') for n in actual),
                qa_count=sum(n.endswith('.qa.json') for n in actual),
                missing_entries=sorted(expected-actual), unexpected_entries=extra)


def read_grid(path, x, y, coverage=None):
    ds = L.gdal.Open(str(path), L.gdal.GA_ReadOnly)
    band = L.validate_grid(ds, x, y, str(path))
    a = band.ReadAsArray()
    require(a is not None and a.shape == (N,N), f'Invalid array: {path}')
    require(np.isfinite(a).all(), f'Nonfinite raster values: {path}')
    inside = L.rasterized_footprint(coverage, ds, x, y) if coverage is not None else None
    band = None
    ds = None
    return a, inside


def empty_footprint(coverage, x, y):
    # Construct the missing tile grid in memory; no source raster is re-extracted.
    ds = L.gdal.GetDriverByName('MEM').Create('', N, N, 1, L.gdal.GDT_Byte)
    srs = L.osr.SpatialReference()
    srs.ImportFromEPSG(26910)
    ds.SetProjection(srs.ExportToWkt())
    ds.SetGeoTransform((x, 1, 0, y+N, 0, -1))
    inside = L.rasterized_footprint(coverage, ds, x, y)
    ds = None
    return inside


def mask_stats(valid, inside):
    kept = valid & inside
    return dict(source_valid_pixels=int(valid.sum()), footprint_pixels=int(inside.sum()),
                valid_inside_footprint=int(kept.sum()), nodata_inside_footprint=int((inside & ~valid).sum()),
                valid_outside_footprint_removed=int((valid & ~inside).sum()),
                masked_valid_pixels=int(kept.sum()), masked_nodata_pixels=N*N-int(kept.sum()))


def audit_mask(source, masked, inside, x, y, label):
    valid = source != NODATA
    kept = valid & inside
    actual = masked != NODATA
    require(np.array_equal(actual, kept), f'{label}: validity differs from source AND footprint at {(x,y)}')
    require(np.array_equal(masked[kept], source[kept]), f'{label}: retained elevations changed at {(x,y)}')
    return mask_stats(valid, inside)


def differences(old, fixed, x, y):
    common = (old != NODATA) & (fixed != NODATA)
    d = fixed[common].astype(np.float64)-old[common].astype(np.float64)
    ad = np.abs(d)
    n = int(d.size)
    s, a, ss = float(d.sum()), float(ad.sum()), float(np.dot(d,d))
    peak = float(ad.max()) if n else None
    coord = None
    if n:
        # Locate the first maximum in row-major order on the common valid grid.
        idx = int(np.flatnonzero(common)[int(ad.argmax())])
        row, col = divmod(idx,N)
        coord = [x+col+0.5,y+N-row-0.5]
    counts = {str(t): int((ad>t).sum()) for t in THRESHOLDS}
    rec = dict(common_valid=n, mean_signed_fixed_minus_original_m=s/n if n else None,
               mae_m=a/n if n else None, rmse_m=(ss/n)**0.5 if n else None,
               abs_diff_p50_p95_max_m=[float(v) for v in np.percentile(ad,[50,95,100])] if n else [None]*3,
               first_max_abs_difference_pixel_center_epsg26910=coord, counts_over_m=counts)
    return rec, (s,a,ss)


def main():
    L.gdal.UseExceptions(); L.ogr.UseExceptions()
    root = e.legacy.project_root()
    coverage, _, fixed_src, fixed_mask = e.paths(root)
    diag = fixed_src.parent
    old_src, old_mask = diag/'hrdem1m_batch_staging', diag/'hrdem1m_masked_staging'
    cfg = e.config(root)
    prior_path = diag/PRIOR
    prior = e.read_json(prior_path)
    require(prior['configuration']==cfg and prior['summary']['masking_review_gate_passed'] is True,
            'Prior global extraction QA is not applicable to current configuration')
    plan = e.legacy.candidate_tiles(coverage)
    keys = {(x,y) for x,y,g in plan}
    require(len(plan)==len(keys)==656 and EMPTY_TILE in keys, 'Candidate inventory changed')
    stems = {e.stem(x,y) for x,y in keys if (x,y)!=EMPTY_TILE}
    inventories = {name:inventory(folder,stems,masked) for name,folder,masked in (
        ('original_source',old_src,False),('fixed_source',fixed_src,False),
        ('original_masked',old_mask,True),('fixed_masked',fixed_mask,True))}
    print('===== GLOBAL MASKING INVENTORY =====',flush=True)
    for name, inv in inventories.items():
        print(name.upper(),json.dumps(inv),flush=True)
        require(not inv['missing_entries'] and not inv['unexpected_entries'], f'Inventory mismatch: {name}')
    mask_hashes = dict(legacy_mask_sha256=e.sha(Path(L.__file__)),fixed_mask_sha256=e.sha(Path(m.__file__)))
    prior_tiles = {tuple(r['tile']):r for r in prior['tiles']}
    require(set(prior_tiles)==keys,'Prior QA candidate mismatch')
    totals = Counter(); group_totals = {}; counts = Counter(); states = Counter()
    signed = absolute = squared = 0.0
    records = []
    print('===== FULL MASKED TILE QA =====',flush=True)
    for i,(x,y,g) in enumerate(plan,1):
        base=e.stem(x,y)
        state, source_qa=e.state(root,x,y,g,cfg)
        states[state]+=1
        require(state==('EMPTY' if (x,y)==EMPTY_TILE else 'COMPLETE'),f'Unexpected extraction state: {base}')
        old_q=e.read_json(old_src/(base+'.qa.json'))
        require((old_q['tile_x'],old_q['tile_y'],old_q['source_group'])==(x,y,g),f'Original source QA identity: {base}')
        require(prior_tiles[x,y]['fixed_state']==state and prior_tiles[x,y]['original_state']==state, f'Prior state changed: {base}')
        if state=='EMPTY':
            require(old_q['status']==e.EMPTY and old_q['valid_pixels']==0 and old_q['nodata_pixels']==N*N,f'Original empty QA mismatch: {base}')
            inside=empty_footprint(coverage,x,y)
            st=mask_stats(np.zeros((N,N),dtype=bool),inside)
            rec=dict(tile=[x,y],source_group=g,state=state,**st,common_valid=0,
                     validity_mismatch=0,retained_values_exact=True, counts_over_m={str(t):0 for t in THRESHOLDS},
                     abs_diff_p50_p95_max_m=[None]*3,first_max_abs_difference_pixel_center_epsg26910=None)
            print('EMPTY_TILE_FOOTPRINT:',json.dumps(rec),flush=True)
        else:
            fsrc,inside=read_grid(fixed_src/(base+'.tif'),x,y,coverage)
            osrc,_=read_grid(old_src/(base+'.tif'),x,y)
            require(np.array_equal(fsrc!=NODATA,osrc!=NODATA),f'Source validity mismatch: {base}')
            fpath=fixed_mask/(base+'_masked.tif'); opath=old_mask/(base+'_masked.tif')
            fa,_=read_grid(fpath,x,y); oa,_=read_grid(opath,x,y)
            st=audit_mask(fsrc,fa,inside,x,y,'FIXED')
            ost=audit_mask(osrc,oa,inside,x,y,'ORIGINAL')
            require(st==ost,f'Original and fixed masking counts differ: {base}')
            require(st['source_valid_pixels']==source_qa['valid_pixels']==old_q['valid_pixels']==prior_tiles[x,y]['fixed_valid']==prior_tiles[x,y]['original_valid'],f'Source valid count mismatch: {base}')
            fq=e.read_json(fixed_mask/(base+'_masked.qa.json'))
            oq=e.read_json(old_mask/(base+'_masked.qa.json'))
            expected_prov=dict(mask_method_id=m.MASK_ID,extraction_configuration=cfg,
                               source_sha256=source_qa['output_sha256'],source_qa_sha256=e.sha(fixed_src/(base+'.qa.json')),
                               footprint_sha256=cfg['footprint_sha256'],**mask_hashes)
            require(fq.get('provenance')==expected_prov and fq.get('output_sha256')==e.sha(fpath),f'Fixed mask provenance or hash mismatch: {base}')
            for q,src,label in ((fq,fixed_src,'FIXED'),(oq,old_src,'ORIGINAL')):
                require((q.get('tile_x'),q.get('tile_y'))==(x,y) and q.get('status')=='MASKED_STAGING_QA_PASSED_WITH_CAVEATS' and all(q.get(k)==v for k,v in st.items()) and q.get('original_values_unchanged_for_retained_pixels') is True and q.get('source_staging_tif')==str((src/(base+'.tif')).relative_to(root)) and q.get('mask_source')==str(coverage.relative_to(root)),f'{label} mask sidecar mismatch: {base}')
            require(np.array_equal(oa!=NODATA,fa!=NODATA),f'Masked validity mismatch: {base}')
            diff, sums=differences(oa,fa,x,y)
            signed+=sums[0]; absolute+=sums[1]; squared+=sums[2]
            rec=dict(tile=[x,y],source_group=g,state=state,**st,**diff,validity_mismatch=0,retained_values_exact=True)
        records.append(rec); totals.update(st); counts.update(rec['counts_over_m'])
        group_totals.setdefault(g,Counter()).update(st)
        if i%25==0 or i==len(plan): print(f'QA_PROGRESS: {i}/{len(plan)}',flush=True)
    require(states=={'COMPLETE':655,'EMPTY':1},'Unexpected final state counts')
    source=totals['source_valid_pixels']; retained=totals['masked_valid_pixels']; footprint=totals['footprint_pixels']; missing=totals['nodata_inside_footprint']; removed=totals['valid_outside_footprint_removed']
    require(source==prior['summary']['fixed_valid_pixels']==502332071,'Global source valid count changed')
    require(source==retained+removed and footprint==retained+missing and retained+totals['masked_nodata_pixels']==len(plan)*N*N and retained==sum(r['common_valid'] for r in records) and footprint>0 and retained>0,'Global pixel accounting mismatch')
    top=sorted((r for r in records if r['common_valid']),key=lambda r:r['abs_diff_p50_p95_max_m'][2],reverse=True)[:10]
    empty=next(r for r in records if r['state']=='EMPTY')
    summary=dict(candidates=len(plan),states=dict(states),actual_masked_tiffs=655,grid_pixels_including_empty=len(plan)*N*N,
                 **dict(totals),empty_tile_footprint_pixels=empty['footprint_pixels'],empty_tile_nodata_inside_footprint=empty['nodata_inside_footprint'],
                 footprint_area_m2_on_1m_grid=footprint,valid_area_m2_on_1m_grid=retained,missing_area_m2_on_1m_grid=missing,
                 footprint_coverage_percent=100*retained/footprint,missing_inside_footprint_percent=100*missing/footprint,
                 actual_masked_tiffs_nodata_pixels=totals['masked_nodata_pixels']-N*N,
                 tiles_with_nodata_inside_footprint=sum(r['nodata_inside_footprint']>0 for r in records),
                 nonempty_source_tiles_with_zero_retained=sum(r['state']=='COMPLETE' and r['masked_valid_pixels']==0 for r in records),
                 original_vs_fixed_masked_validity_mismatch=0,retained_values_modified=0,valid_pixels_outside_footprint_after_masking=0,
                 fixed_mask_provenance_and_hashes_passed=True,all_mask_arrays_exact=True,
                 mean_signed_fixed_minus_original_m=signed/retained,mae_m=absolute/retained,rmse_m=(squared/retained)**0.5,
                 maximum_abs_difference_m=top[0]['abs_diff_p50_p95_max_m'][2],counts_over_m=dict(counts),
                 tiles_with_elevation_difference_over_0_5m=sum(r['counts_over_m']['0.5']>0 for r in records),masking_qa_passed=True)
    print('===== GLOBAL MASKING SUMMARY =====',flush=True); print(json.dumps(summary,indent=2),flush=True)
    print('===== MASKING TOTALS BY SOURCE GROUP =====',flush=True); print(json.dumps(group_totals,indent=2),flush=True)
    print('===== TOP 10 MASKED TILES BY MAXIMUM ABSOLUTE DIFFERENCE =====',flush=True)
    for r in top: print(json.dumps(r),flush=True)
    report_dir=Path(tempfile.mkdtemp(prefix='hrdem1m_global_masking_qa_',dir=diag))
    report=report_dir/'global_masking_qa.json'
    payload=dict(configuration=cfg,mask_method_id=m.MASK_ID,mask_script_hashes=mask_hashes,
                 prior_extraction_qa=str(prior_path.relative_to(root)),prior_extraction_qa_sha256=e.sha(prior_path),
                 audit_script_sha256=e.sha(Path(__file__)),inventories=inventories,summary=summary,
                 source_group_totals=group_totals,tiles=records,top_10_by_maximum_difference=top,
                 limitations=['Footprint coverage is on the candidate grid using the existing pixel-center rasterization rule, not exact vector area or coverage of the entire AOI.',
                              'The empty extraction is represented by an all-NoData grid; its source was not fetched again.',
                              'The reported masked_nodata_pixels includes the virtual empty grid; actual TIFF NoData is separately reported.',
                              'Original files lack historical SHA-256; old mask values are verified against current old source arrays.',
                              'Original versus fixed differences are not absolute vertical accuracy errors.',
                              'This checks masking, not source artifacts, seam quality or final DEM validity.'])
    e.save_json(report,payload)
    print('DIAGNOSTIC_JSON:',report,flush=True)
    print('GLOBAL_MASKING_QA_PASSED: True',flush=True); print('FINAL_DEM: FALSE',flush=True)


if __name__=='__main__':
    main()
