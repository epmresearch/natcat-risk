"""Targeted negative-elevation source audit. Read-only inputs; diagnostic JSON only."""
from pathlib import Path
import json
import tempfile
import numpy as np
from osgeo import gdal
import step4_2_hrdem1m_fixed_batch_300 as e
import step4_2_hrdem1m_global_masking_qa as q
import step4_2_hrdem1m_multisite_fixed_qa as multi
import step4_2_hrdem1m_cross_mosaic_qa as cross
import step4_2_hrdem1m_final_validation as f

PRIOR = 'hrdem1m_final_validation_ebk_g0f9/final_validation_qa.json'
DEEP = {(489000, 5506000): 1806, (479000, 5468000): 626}
CONTROL = (477000, 5473000)
TOL = 0.001


def valid(a):
    return np.isfinite(a) & (a != q.NODATA)


def deep_bounds(a, x, y):
    rr, cc = np.where(valid(a) & (a < -5))
    if not rr.size:
        raise RuntimeError('No deep-negative pixels')
    # A two-pixel context margin, clipped to the existing tile.
    l, r = max(0, int(cc.min())-2), min(1000, int(cc.max())+3)
    t, b = max(0, int(rr.min())-2), min(1000, int(rr.max())+3)
    return (x+l, y+1000-b, x+r, y+1000-t)


def compare(a, ref, strict, select):
    eligible = select & strict
    d = np.abs(a[eligible].astype(np.float64)-ref[eligible])
    return dict(selected=int(select.sum()), strict=int(eligible.sum()),
                without_strict=int((select & ~strict).sum()),
                max_abs_m=float(d.max()) if d.size else None,
                over_0_001m=int((d > TOL).sum()))


def main():
    gdal.UseExceptions(); e.legacy.ogr.UseExceptions()
    root=e.legacy.project_root(); coverage,vrt,srcdir,maskdir=e.paths(root)
    diag=srcdir.parent; cfg=e.config(root)
    p=diag/PRIOR; ph=e.sha(p); prev=e.read_json(p); ps=prev['summary']
    q.require(prev['configuration']==cfg and prev['script_sha256']==e.sha(Path(f.__file__)), 'Final QA configuration mismatch')
    q.require(ps['final_technical_validation_passed'] and not ps['final_dem_published'], 'Unexpected final QA state')
    q.require(ps['negative_pixels']==20054860 and ps['negative_counts_below_m']=={'-1':18408762,'-5':2432,'-10':1685,'-20':484}, 'Negative inventory changed')
    cp=root/prev['candidate_report']; q.require(e.sha(cp)==prev['candidate_report_sha256'],'Candidate report changed')
    candidate=e.read_json(cp); cs=candidate['summary']; mp=root/prev['mosaic']; gp=root/prev['register']
    q.require(e.sha(mp)==ps['mosaic_sha256']==cs['mosaic_sha256']==f.MOSAIC_SHA and e.sha(gp)==prev['register_sha256']==f.REGISTER_SHA,'Mosaic/register changed')
    hashes={tuple(r['tile']):r['sha256'] for r in candidate['input_tile_sha256']}
    rows={tuple(r['tile']):r for r in prev['negative_tile_details']}
    q.require({k:r['below_minus_5'] for k,r in rows.items() if r['below_minus_5']}==DEEP,'Deep negative tile inventory mismatch')
    lookup={(x,y):g for x,y,g in e.legacy.candidate_tiles(coverage)}
    for key,value in {'GDAL_DISABLE_READDIR_ON_OPEN':'YES','GDAL_HTTP_MAX_RETRY':'3','GDAL_HTTP_RETRY_DELAY':'2','GDAL_HTTP_TIMEOUT':'45'}.items():
        gdal.SetConfigOption(key,value)
    src=gdal.Open(str(vrt),gdal.GA_ReadOnly)
    q.require((src.RasterXSize,src.RasterYSize)==(500000,1000000) and src.GetGeoTransform()==(-2000000,1,0,1000000,0,-1) and src.GetSpatialRef().GetAuthorityCode(None)=='3979' and src.GetRasterBand(1).GetNoDataValue()==q.NODATA,'Native VRT mismatch')
    ds=gdal.Open(str(mp),gdal.GA_ReadOnly); band=ds.GetRasterBand(1)
    results=[]
    for k in [*DEEP, CONTROL]:
        x,y=k; stem=e.stem(*k); r=rows[k]
        state,_=e.state(root,x,y,lookup[k],cfg); q.require(state=='COMPLETE','Extraction state changed')
        path=maskdir/(stem+'_masked.tif'); q.require(e.sha(path)==hashes[k], 'Masked tile changed')
        m,_=q.read_grid(path,x,y); u,_=q.read_grid(srcdir/(stem+'.tif'),x,y)
        mo=band.ReadAsArray(x-475000,5557000-y-1000,1000,1000)
        q.require(np.array_equal(m.view(np.uint32),mo.view(np.uint32)) and np.array_equal(m[valid(m)].view(np.uint32),u[valid(m)].view(np.uint32)),'Local transfer mismatch')
        q.require(int((valid(m)&(m < -5)).sum())==r['below_minus_5'],'Deep count changed')
        px,py=r['minimum']['center_epsg26910']; col=int(px-x); row=int(y+1000-py)
        q.require(float(m[row,col])==r['minimum']['value_m'],'Minimum changed')
        if k in DEEP:
            bbox=deep_bounds(m,x,y)
        else:
            # Shallow-negative comparison around the reported tile minimum.
            bbox=(max(x,int(px)-32),max(y,int(py)-32),min(x+1000,int(px)+33),min(y+1000,int(py)+33))
        l,b,rr,t=bbox; sl=(slice(y+1000-t,y+1000-b),slice(l-x,rr-x))
        a=m[sl]; un=u[sl]
        print('AUDITING_NEGATIVE_SOURCE:',k,bbox,flush=True)
        ref,ok,window,region=multi.native_reference(src,bbox)
        direct=cross.direct_source_audit(src,window)
        direct_ok=bool(direct and all(v['exact_equal'] and v['validity_mismatch']==0 for v in direct))
        selected=valid(a)&(a < (-5 if k in DEEP else 0))
        target=compare(a,ref,ok,selected)
        context=compare(un,ref,ok,valid(un))
        strict_missing=int((ok&~valid(un)).sum())
        ri,ci=int(t-py),int(px-l)
        minimum=dict(center_epsg26910=[px,py],fixed_m=float(a[ri,ci]),strict=bool(ok[ri,ci]),native_reference_m=float(ref[ri,ci]) if ok[ri,ci] else None)
        minimum['fixed_minus_reference_m']=minimum['fixed_m']-minimum['native_reference_m'] if minimum['strict'] else None
        q.require(target['selected']==DEEP.get(k,target['selected']) and target['selected']>0,'Selected negative pixels missing')
        rec=dict(tile=list(k),kind='deep' if k in DEEP else 'shallow_control',bbox=list(bbox),source_window=window,source_region=region,direct_source_audit=direct,direct_source_equal=direct_ok,negative_comparison=target,unmasked_context_comparison=context,strict_missing_from_unmasked=strict_missing,minimum=minimum)
        rec['numerical_fidelity_passed']=bool(direct_ok and target['strict']>0 and target['over_0_001m']==0 and context['over_0_001m']==0 and strict_missing==0 and minimum['strict'] and abs(minimum['fixed_minus_reference_m'])<=TOL)
        results.append(rec); print('SITE_RESULT:',json.dumps(rec,allow_nan=False),flush=True)
    # Record what the source itself declares; do not infer a vertical datum from a horizontal EPSG code.
    metadata=[]
    for url in sorted(cross.halo.EXPECTED_VRT):
        cog=gdal.Open(url,gdal.GA_ReadOnly); srs=cog.GetSpatialRef(); cb=cog.GetRasterBand(1)
        metadata.append(dict(url=url,crs_wkt=srs.ExportToWkt(),band_unit=cb.GetUnitType(),dataset_metadata=cog.GetMetadata(),band_metadata=cb.GetMetadata()))
        cb=None; cog=None
    band=None; ds=None; src=None
    deep=[r for r in results if r['kind']=='deep']
    summary=dict(sites=len(results),sites_numerically_consistent=sum(r['numerical_fidelity_passed'] for r in results),deep_negative_pixels=sum(r['negative_comparison']['selected'] for r in deep),deep_negative_strict=sum(r['negative_comparison']['strict'] for r in deep),deep_negative_without_strict=sum(r['negative_comparison']['without_strict'] for r in deep),deep_negative_over_0_001m=sum(r['negative_comparison']['over_0_001m'] for r in deep),all_sites_numerically_consistent=all(r['numerical_fidelity_passed'] for r in results),absolute_vertical_accuracy_validated=False,negative_source_physical_validity_resolved=False,final_dem_published=False)
    q.require(summary['deep_negative_pixels']==2432,'Deep inventory incomplete')
    q.require(e.sha(p)==ph and e.sha(mp)==f.MOSAIC_SHA and e.sha(gp)==f.REGISTER_SHA and e.sha(vrt)==cfg['vrt_sha256'],'Inputs changed during audit')
    folder=Path(tempfile.mkdtemp(prefix='hrdem1m_negative_source_qa_',dir=diag)); report=folder/'negative_source_qa.json'
    e.save_json(report,dict(configuration=cfg,script_sha256=e.sha(Path(__file__)),prior_final_report=str(p.relative_to(root)),prior_final_report_sha256=ph,summary=summary,sites=results,source_metadata=metadata,limitations=['Strict reference uses four valid native neighbors; partial stencils remain unverified.','Native source fidelity is not independent physical validation.','Shallow control is a targeted patch, not a representative sample.','No elevation modification or publication.']))
    print('===== NEGATIVE SOURCE SUMMARY ====='); print(json.dumps(summary,indent=2))
    print('SOURCE_METADATA:',json.dumps(metadata,allow_nan=False))
    print('DIAGNOSTIC_JSON:',report); print('NEGATIVE_SOURCE_QA_COMPLETED: True'); print('FINAL_DEM_PUBLISHED: FALSE')
    q.require(summary['all_sites_numerically_consistent'],'Numerical fidelity review required; see report')


if __name__=='__main__': main()
