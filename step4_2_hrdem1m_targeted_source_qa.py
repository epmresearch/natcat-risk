"""Targeted source-fidelity audit; reads local TIFFs and official HRDEM COG windows.
No GDAL warp, raster writes, or changes to production files.
"""
from pathlib import Path
from functools import lru_cache
import json
import tempfile
import numpy as np
from osgeo import gdal
import step4_2_hrdem1m_fixed_batch_300 as e
import step4_2_hrdem1m_global_masking_qa as q
import step4_2_hrdem1m_global_fixed_seam_qa as s
import step4_2_hrdem1m_seam_local_context_qa as local
import step4_2_hrdem1m_multisite_fixed_qa as multi
import step4_2_hrdem1m_cross_mosaic_qa as cross

PRIOR = 'hrdem1m_seam_local_context_p_79r4to/seam_local_context_qa.json'
SEAMS = [
    (483000,5503000,'east'),   # Strong localization, historical
    (484000,5505000,'east'),   # Largest individual residual
    (490000,5536000,'north'),  # Strong localization, historical
    (487000,5507000,'north'),  # Newly flagged, localized
    (484000,5488000,'east'),   # Severe, incomplete common support
]
CONTROLS = [(482000,5503000,'east'),(489000,5536000,'north')]
HOTSPOTS = [(490453.5,5503923.5),(484884.5,5505072.5)]
TOL = 0.001
W = 52

def valid(a):
    return np.isfinite(a) & (a != q.NODATA)

def seam_box(x,y,d):
    if d == 'east': return (x+1000-W,y,x+1000+W,y+1000)
    if d == 'north': return (x,y+1000-W,x+1000,y+1000+W)
    raise ValueError(d)

def center_stencil(a,d):
    if d=='east': return a[:,W-2:W+2]
    if d=='north': return a[W-2:W+2,:][::-1,:].T
    raise ValueError(d)

def compare_reference(a, ref, support):
    use = support & valid(a)
    delta = np.abs(a[use].astype(np.float64)-ref[use])
    return dict(strict_support=int(support.sum()), common_valid=int(use.sum()),
                strict_missing=int((support&~valid(a)).sum()),
                max_abs_m=float(delta.max()) if delta.size else None,
                mae_m=float(delta.mean()) if delta.size else None,
                over_0_001m=int((delta>TOL).sum()))

def source_fidelity_ok(rec):
    return rec['common_valid']>=1000 and rec['strict_missing']==0 and rec['over_0_001m']==0

def main():
    gdal.UseExceptions(); e.legacy.ogr.UseExceptions()
    root=e.legacy.project_root(); coverage,vrt,srcdir,maskdir=e.paths(root)
    diag=srcdir.parent; cfg=e.config(root)
    prior_path=diag/PRIOR; prior=e.read_json(prior_path)
    q.require(prior['configuration']==cfg and prior['script_sha256']==e.sha(Path(local.__file__)),'Local QA configuration mismatch')
    q.require(prior['summary']['central_metrics_reproduced'] is True,'Prior local QA incomplete')
    seam_path=root/prior['prior_seam_report']
    q.require(e.sha(seam_path)==prior['prior_seam_report_sha256'],'Seam report hash mismatch')
    seam_report=e.read_json(seam_path)
    mask_path=root/seam_report['prior_masking_qa']
    q.require(e.sha(mask_path)==seam_report['prior_masking_qa_sha256'],'Mask report hash mismatch')
    mask_report=e.read_json(mask_path)
    q.require(mask_report['configuration']==cfg and mask_report['summary']['masking_qa_passed'] is True,'Mask QA mismatch')
    lookup={(x,y):g for x,y,g in e.legacy.candidate_tiles(coverage)}
    sites={(*r['tile'],r['direction']):r for r in prior['sites']}
    tiles={tuple(r['tile']):r for r in mask_report['tiles']}
    q.require(len(lookup)==656 and set(lookup)==set(tiles),'Inventory mismatch')
    for k in SEAMS+CONTROLS:
        q.require(k in sites and sites[k]['kind']==('flagged' if k in SEAMS else 'control'),f'Unexpected site {k}')
    for key,value in {'GDAL_DISABLE_READDIR_ON_OPEN':'YES','GDAL_HTTP_MAX_RETRY':'3',
                      'GDAL_HTTP_RETRY_DELAY':'2','GDAL_HTTP_TIMEOUT':'45'}.items():
        gdal.SetConfigOption(key,value)
    src=gdal.Open(str(vrt),gdal.GA_ReadOnly)
    q.require((src.RasterXSize,src.RasterYSize)==(500000,1000000) and
              np.allclose(src.GetGeoTransform(),(-2000000,1,0,1000000,0,-1),rtol=0,atol=1e-8)
              and src.GetSpatialRef().GetAuthorityCode(None)=='3979'
              and src.GetRasterBand(1).GetNoDataValue()==q.NODATA,'Native VRT geometry mismatch')
    @lru_cache(maxsize=8)
    def read(k):
        x,y=k; stem=e.stem(x,y)
        state,rec=e.state(root,x,y,lookup[k],cfg)
        q.require(state=='COMPLETE',f'Unexpected tile state {k}')
        u,_=q.read_grid(srcdir/(stem+'.tif'),x,y)
        mp=maskdir/(stem+'_masked.tif'); side=e.read_json(maskdir/(stem+'_masked.qa.json'))
        prov=side['provenance']
        q.require(prov['extraction_configuration']==cfg and all(prov.get(n)==v for n,v in mask_report['mask_script_hashes'].items())
                  and side['output_sha256']==e.sha(mp),'Masked provenance/hash mismatch')
        m,_=q.read_grid(mp,x,y)
        o,_=q.read_grid(diag/'hrdem1m_masked_staging'/(stem+'_masked.tif'),x,y)
        q.require(np.array_equal(valid(o),valid(m)) and int(valid(m).sum())==tiles[k]['masked_valid_pixels']
                  and np.all(valid(u)[valid(m)]) and np.array_equal(u[valid(m)],m[valid(m)]),'Local mask mismatch')
        return u,m,o
    def patch(bbox):
        xmin,ymin,xmax,ymax=bbox
        out=[np.full((ymax-ymin,xmax-xmin),q.NODATA,dtype=np.float32) for _ in range(3)]
        covered=np.zeros(out[0].shape,dtype=bool)
        for y in range(ymin//1000*1000,ymax,1000):
            for x in range(xmin//1000*1000,xmax,1000):
                l,b,r,t=max(x,xmin),max(y,ymin),min(x+1000,xmax),min(y+1000,ymax)
                if l>=r or b>=t: continue
                a=read((x,y))
                dest=(slice(ymax-t,ymax-b),slice(l-xmin,r-xmin))
                orig=(slice(y+1000-t,y+1000-b),slice(l-x,r-x))
                for dest_a,src_a in zip(out,a): dest_a[dest]=src_a[orig]
                covered[dest]=True
        q.require(covered.all(),'Patch not fully covered by production tiles')
        return out
    results=[]
    def audit(label,bbox,d=None,point=None):
        print('AUDITING_SOURCE:',label,bbox,flush=True)
        u,m,o=patch(bbox)
        ref,ok,window,region=multi.native_reference(src,bbox)
        direct=cross.direct_source_audit(src,window)
        direct_ok=bool(direct and all(r['validity_mismatch']==0 and r['exact_equal'] for r in direct))
        unmasked=compare_reference(u,ref,ok)
        masked=compare_reference(m,ref,ok&valid(m))
        old=compare_reference(o,ref,ok&valid(m))
        rec=dict(label=label,bbox=list(bbox),source_window=window,source_region=region,
                 direct_source_audit=direct,direct_source_equal=direct_ok,
                 fixed_unmasked_vs_native=unmasked,fixed_masked_vs_native=masked,
                 original_masked_vs_native=old,
                 masked_valid_pixels=int(valid(m).sum()),
                 masked_valid_without_strict_reference=int((valid(m)&~ok).sum()),
                 numerical_fidelity_passed=bool(direct_ok and source_fidelity_ok(unmasked)))
        if d:
            a,b,c=(center_stencil(v,d) for v in (m,o,ref))
            support=center_stencil(ok,d).all(axis=1)&valid(a).all(axis=1)&valid(b).all(axis=1)
            rr=[s.residual(v[support]) for v in (a,b,c)]
            central=s.describe(s.residual(a[valid(a).all(axis=1)]))
            rec.update(central_fixed_all=central,
                       strict_common_crossings=int(support.sum()),
                       strict_fraction_of_masked_crossings=float(support.sum()/central['crossings']) if central['crossings'] else None,
                       fixed_residual_common=s.describe(rr[0]),old_residual_common=s.describe(rr[1]),
                       native_residual_common=s.describe(rr[2]),
                       max_abs_fixed_native_residual_difference_m=float(np.max(np.abs(rr[0]-rr[2]))) if support.any() else None)
            rec['residual_support_sufficient']=bool(support.sum()>=20 and rec['strict_fraction_of_masked_crossings']>=.8)
        if point:
            x,y=point; row=int(bbox[3]-y); col=int(x-bbox[0])
            q.require(valid(m)[row,col] and valid(o)[row,col],'Hotspot no longer valid')
            rec['hotspot']={'center_epsg26910':list(point),'original_m':float(o[row,col]),'fixed_m':float(m[row,col]),
                            'native_m':float(ref[row,col]) if ok[row,col] else None,
                            'strict_reference_valid':bool(ok[row,col]),
                            'fixed_minus_original_m':float(m[row,col])-float(o[row,col]),
                            'fixed_minus_native_m':float(m[row,col])-float(ref[row,col]) if ok[row,col] else None}
        results.append(rec)
        print('SITE_RESULT:',json.dumps(rec,allow_nan=False),flush=True)
        return rec
    for k in SEAMS+CONTROLS:
        x,y,d=k; rec=audit(f'{sites[k]["kind"]}:{x},{y},{d}',seam_box(*k),d)
        prior_center=sites[k]['center']
        q.require(rec['central_fixed_all']['crossings']==prior_center['crossings'] and
                  abs(rec['central_fixed_all']['p95_abs_m']-prior_center['p95_abs_m'])<=1e-9,'Central metric changed')
    for x,y in HOTSPOTS:
        ix,iy=int(x),int(y)
        rec=audit(f'hotspot:{x},{y}',(ix-32,iy-32,ix+33,iy+33),point=(x,y))
        expected=tiles[(ix//1000*1000,iy//1000*1000)]
        q.require(expected['first_max_abs_difference_pixel_center_epsg26910']==[x,y] and
                  abs(abs(rec['hotspot']['fixed_minus_original_m'])-expected['abs_diff_p50_p95_max_m'][2])<=1e-9,'Hotspot metric changed')
    summary=dict(sites=len(results),direct_source_equal=sum(r['direct_source_equal'] for r in results),
                 numerical_fidelity_passed=sum(r['numerical_fidelity_passed'] for r in results),
                 strict_reference_pixels=sum(r['fixed_unmasked_vs_native']['strict_support'] for r in results),
                 strict_reference_missing=sum(r['fixed_unmasked_vs_native']['strict_missing'] for r in results),
                 fixed_unmasked_over_0_001m=sum(r['fixed_unmasked_vs_native']['over_0_001m'] for r in results),
                 max_fixed_unmasked_native_difference_m=max((r['fixed_unmasked_vs_native']['max_abs_m'] for r in results if r['fixed_unmasked_vs_native']['max_abs_m'] is not None),default=None),
                 seams_with_sufficient_strict_support=sum(r.get('residual_support_sufficient',False) for r in results),
                 hotspots_with_strict_reference=sum(r.get('hotspot',{}).get('strict_reference_valid',False) for r in results),
                 all_sites_numerically_consistent=all(r['numerical_fidelity_passed'] for r in results),
                 final_dem_approved=False)
    out=Path(tempfile.mkdtemp(prefix='hrdem1m_targeted_source_qa_',dir=diag))/'targeted_source_qa.json'
    e.save_json(out,dict(configuration=cfg,prior_local_report=str(prior_path.relative_to(root)),prior_local_report_sha256=e.sha(prior_path),
                         script_sha256=e.sha(Path(__file__)),helper_hashes={p.__name__:e.sha(Path(p.__file__)) for p in (multi,cross,q,s,local)},
                         summary=summary,sites=results,
                         limitations=['Targeted purposive sample, not exhaustive.','Strict four-valid-neighbor native reference from the same source; not independent terrain truth.','Direct source audit verifies VRT pixel values for sampled windows only. No source version guarantee.','Source fidelity does not distinguish natural terrain from internal source artifacts.']))
    print('===== TARGETED SOURCE SUMMARY ====='); print(json.dumps(summary,indent=2))
    print('DIAGNOSTIC_JSON:',out)
    print('TARGETED_SOURCE_QA_COMPLETED: True')
    print('FINAL_DEM: FALSE')
    return 0 if summary['all_sites_numerically_consistent'] else 2

if __name__=='__main__':
    raise SystemExit(main())
