"""Diagnostic cross-mosaic QA; never modifies existing staging rasters.
Requires the existing multisite, bridge QA and halo64 pilot modules.
"""
from __future__ import annotations
import json
import os
import tempfile
from pathlib import Path
import numpy as np
from osgeo import gdal, osr
import step4_2_hrdem1m_multisite_fixed_qa as multi
import step4_2_hrdem1m_flagged_seams_bridge_qa as qa
import step4_2_hrdem1m_halo64_pilot as halo

CANDIDATES = [(481000,5480000),(481000,5481000),(482000,5481000),
              (483000,5481000),(483000,5482000),(484000,5482000)]
BOUNDARY = 500000  # Native row boundary and native northing coincide here.
HALF = 64
TOL = 0.001
MIN_REF = 1000
MIN_CROSS = 20
MIN_SIDE = 100


def transform(a, b):
    s, t = osr.SpatialReference(), osr.SpatialReference()
    s.ImportFromEPSG(a); t.ImportFromEPSG(b)
    for crs in (s,t): crs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return osr.CoordinateTransformation(s,t)


def screen():
    # Rank local patches by existing valid coverage, without touching remote sources.
    forward, backward = transform(26910,3979), transform(3979,26910)
    choices = []
    for x,y in CANDIDATES:
        old = halo.inspect(halo.STAGING/f'hrdem1m_x{x}_y{y}_utm10_1km.tif',x,y,1000,1000)
        m = halo.valid(old)
        corners = np.asarray(forward.TransformPoints([(x+dx,y+dy) for dx in (0,1000) for dy in (0,1000)]))
        xx = np.arange(np.floor(corners[:,0].min()),np.ceil(corners[:,0].max())+1,10.)
        points = np.asarray(backward.TransformPoints([(float(v),float(BOUNDARY)) for v in xx]))
        ranked = []
        for cx,cy in np.rint(points[:,:2]).astype(int):
            if not (x+HALF<=cx<=x+1000-HALF and y+HALF<=cy<=y+1000-HALF): continue
            col,row = cx-HALF-x,y+1000-(cy+HALF)
            score = int(m[row:row+2*HALF,col:col+2*HALF].sum())
            if score: ranked.append((score,int(cx),int(cy)))
        selected = []
        for score,cx,cy in sorted(ranked, key=lambda q:(-q[0],q[1],q[2])):
            if all((cx-u)**2+(cy-v)**2>=64**2 for _,u,v in selected):
                selected.append((score,cx,cy))
            if len(selected)==3: break
        print('LOCAL_CANDIDATE:',(x,y),'PATCHES:',selected,flush=True)
        choices.extend((score,x,y,cx,cy) for score,cx,cy in selected)
    return sorted(choices,key=lambda q:(-q[0],*q[1:]))


def reference(src,bbox):
    ref,ok,window,region = multi.native_reference(src,bbox)
    xmin,ymin,xmax,ymax = bbox
    rr,cc = np.indices(ref.shape)
    xy = np.asarray(transform(26910,3979).TransformPoints(list(zip(
        (xmin+cc.ravel()+.5).tolist(),(ymax-rr.ravel()-.5).tolist()))))
    if not np.isfinite(xy[:,:2]).all(): raise RuntimeError('Invalid coordinates')
    inv = gdal.InvGeoTransform(src.GetGeoTransform())
    iy = np.floor(inv[3]+inv[4]*xy[:,0]+inv[5]*xy[:,1]-.5).astype(np.int64).reshape(ref.shape)
    cross = iy==BOUNDARY-1
    north = iy+1<BOUNDARY
    south = iy>=BOUNDARY
    counts = {'strict_reference_valid':int(ok.sum()),'strict_cross_source_stencils':int((ok&cross).sum()),
              'strict_2_4_only':int((ok&north).sum()),'strict_2_3_only':int((ok&south).sum())}
    return ref,ok,cross,north,south,window,region,counts


def direct_source_audit(src,window):
    x,y,w,h = window
    vb = src.GetRasterBand(1)
    result = []
    for region,lo,hi in [('2_4',0,BOUNDARY),('2_3',BOUNDARY,1000000)]:
        a,b = max(y,lo),min(y+h,hi)
        if b<=a: continue
        url = next(u for u in halo.EXPECTED_VRT if f'{region}-mosaic' in u)
        ds = gdal.Open(url,gdal.GA_ReadOnly)
        expected = (-2000000.,1.,0.,1000000.-lo,0.,-1.)
        if ((ds.RasterXSize,ds.RasterYSize)!=(500000,500000) or
            not np.allclose(ds.GetGeoTransform(),expected,rtol=0,atol=1e-8) or
            not ds.GetSpatialRef().IsSame(src.GetSpatialRef())):
            raise RuntimeError(f'Unexpected source geometry: {region}')
        sb = ds.GetRasterBand(1)
        if sb.GetNoDataValue()!=qa.NODATA: raise RuntimeError('Unexpected source NoData')
        win_v=(x,a,w,b-a); win_s=(x,a-lo,w,b-a)
        v=vb.ReadAsArray(*win_v); s=sb.ReadAsArray(*win_s)
        vv=qa.valid(v)&(vb.GetMaskBand().ReadAsArray(*win_v)>0)
        sv=qa.valid(s)&(sb.GetMaskBand().ReadAsArray(*win_s)>0)
        m=vv&sv
        record={'region':region,'common_valid':int(m.sum()),'validity_mismatch':int((vv!=sv).sum()),
                'exact_equal':bool(m.any() and np.array_equal(v[m],s[m]))}
        result.append(record)
    return result


def main():
    gdal.UseExceptions(); qa.validate_vrt()
    choices=screen()
    opts={'GDAL_DISABLE_READDIR_ON_OPEN':'YES','GDAL_HTTP_MAX_RETRY':'3',
          'GDAL_HTTP_RETRY_DELAY':'2','GDAL_HTTP_TIMEOUT':'45'}
    env=os.environ.copy(); env.update(opts)
    for k,v in opts.items(): gdal.SetConfigOption(k,v)
    src=gdal.Open(str(qa.VRT),gdal.GA_ReadOnly)
    if ((src.RasterXSize,src.RasterYSize)!=(500000,1000000) or
        not np.allclose(src.GetGeoTransform(),(-2000000.,1.,0.,1000000.,0.,-1.),rtol=0,atol=1e-8) or
        src.GetSpatialRef().GetAuthorityCode(None)!='3979' or src.GetRasterBand(1).GetNoDataValue()!=qa.NODATA):
        raise RuntimeError('Unexpected native VRT grid')
    out=Path(tempfile.mkdtemp(prefix='hrdem1m_cross_mosaic_',dir=qa.DIAG))
    path=out/'cross_mosaic_qa.json'
    report={'gdal_version':gdal.VersionInfo('--version'),'tolerance_m':TOL,
            'min_reference':MIN_REF,'min_cross_source_stencils':MIN_CROSS,'min_each_side':MIN_SIDE,
            'candidate_screen':'Up to three locally high-coverage 128m patches per candidate; not exhaustive',
            'screened_patches':[],'commands':[],'passed':False}
    print('OUTPUT_DIRECTORY:',out,flush=True)
    for score,x,y,cx,cy in choices:
        bbox=(cx-HALF,cy-HALF,cx+HALF,cy+HALF)
        ref,ok,cross,north,south,window,region,counts=reference(src,bbox)
        r={'tile':[x,y],'bbox':list(bbox),'source_window':window,'source_region':region,**counts}
        report['screened_patches'].append(r)
        print('NATIVE_SUPPORT:',json.dumps(r),flush=True)
        qa.write_json_atomic(path,report)
        if (region=='crosses_VRT_boundary' and counts['strict_reference_valid']>=MIN_REF and
            counts['strict_cross_source_stencils']>=MIN_CROSS and counts['strict_2_4_only']>=MIN_SIDE and
            counts['strict_2_3_only']>=MIN_SIDE): break
    else:
        print('CROSS_MOSAIC_QA_PASSED: FALSE (no qualifying sampled patch)')
        qa.write_json_atomic(path,report); print('DIAGNOSTIC_JSON:',path); return 2
    audit=direct_source_audit(src,window)
    print('DIRECT_SOURCE_AUDIT:',json.dumps(audit),flush=True)
    report['selected_patch']=r; report['direct_source_audit']=audit
    qa.write_json_atomic(path,report)
    if len(audit)!=2 or not all(q['validity_mismatch']==0 and q['exact_equal'] for q in audit):
        print('CROSS_MOSAIC_QA_PASSED: FALSE (source audit failed)'); return 2
    wide=(x-64,y-64,x+1064,y+1064)
    a=multi.warp(out/'wide_fixed.tif',wide,env,report['commands'])
    b=multi.warp(out/'narrow_fixed.tif',bbox,env,report['commands'])
    col,row=bbox[0]-wide[0],wide[3]-bbox[3]
    a=a[row:row+b.shape[0],col:col+b.shape[1]]
    if a.shape!=b.shape: raise RuntimeError('Wide crop mismatch')
    extent=multi.stats(a,b)
    va,vb=qa.valid(a),qa.valid(b)
    comparisons={}
    for name,support in [('all_strict',ok),('cross_source',ok&cross),('2_4_only',ok&north),('2_3_only',ok&south)]:
        m=support&va&vb
        comparisons[name]={'strict_count':int(support.sum()),'missing_wide':int((support&~va).sum()),
                           'missing_narrow':int((support&~vb).sum()),
                           'wide_vs_native':multi.stats(a,ref,m,m),'narrow_vs_native':multi.stats(b,ref,m,m)}
    passed=(extent['common_valid']>=MIN_REF and extent['validity_mismatch']==0 and extent['over_0_001m']==0 and
            all(q['strict_count']>0 and q['missing_wide']==q['missing_narrow']==0 and
                q['wide_vs_native']['over_0_001m']==q['narrow_vs_native']['over_0_001m']==0 for q in comparisons.values()))
    report.update({'wide_bbox':list(wide),'wide_vs_narrow':extent,'reference_comparisons':comparisons,'passed':bool(passed),
                   'note':'Numerical source fidelity only, not independent vertical accuracy.'})
    qa.write_json_atomic(path,report)
    print('WIDE_VS_NARROW:',json.dumps(extent))
    for name,q in comparisons.items(): print('REFERENCE_'+name.upper()+':',json.dumps(q))
    print('CROSS_MOSAIC_QA_PASSED:',bool(passed)); print('DIAGNOSTIC_JSON:',path); print('ORIGINAL_TILES_MODIFIED: FALSE')
    return 0 if passed else 2

if __name__=='__main__':
    raise SystemExit(main())
