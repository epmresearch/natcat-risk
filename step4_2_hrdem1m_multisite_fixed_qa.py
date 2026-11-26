"""Diagnostic multi-site fixed-scale QA. Run from natcat-risk root.
No existing staging, masked, or diagnostic rasters are overwritten.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from osgeo import gdal, osr
import step4_2_hrdem1m_flagged_seams_bridge_qa as qa

# Three previously flagged sites, excluding the extensively tested first site.
FLAGGED = [(489000, 5503000, 'east'), (485000, 5489000, 'north'),
           (483000, 5467000, 'north')]
CONTROL_TARGETS = [(489000, 5503000, 'east'), (483000, 5467000, 'north')]
TOL = 0.001  # Numerical agreement threshold, not vertical accuracy.
MIN_CONTROL = 500
MIN_REFERENCE = 1000
NODATA = qa.NODATA


def stats(a, b, va=None, vb=None):
    va = qa.valid(a) if va is None else va
    vb = qa.valid(b) if vb is None else vb
    m = va & vb
    delta = np.abs(a[m].astype(np.float64) - b[m].astype(np.float64))
    return {'common_valid': int(m.sum()), 'validity_mismatch': int((va != vb).sum()),
            'first_only_valid': int((va & ~vb).sum()),
            'second_only_valid': int((vb & ~va).sum()),
            'abs_diff_p50_p95_max_m': np.percentile(delta, [50, 95, 100]).tolist() if delta.size else None,
            'over_0_001m': int((delta > TOL).sum()),
            'over_0_5m': int((delta > 0.5).sum())}


def select_controls():
    pat = re.compile(r'hrdem1m_x(\d+)_y(\d+)_utm10_1km_masked\.tif$')
    tiles = {tuple(map(int, m.groups())) for p in qa.MASK.glob('*.tif')
             if (m := pat.fullmatch(p.name))}
    known = set(qa.FLAGS)
    controls = []
    for tx, ty, direction in CONTROL_TARGETS:
        dx, dy = (1000, 0) if direction == 'east' else (0, 1000)
        candidates = [(x, y, direction) for x, y in tiles
                      if (x + dx, y + dy) in tiles and (x, y, direction) not in known]
        candidates.sort(key=lambda t: ((t[0]-tx)**2+(t[1]-ty)**2, t[0], t[1]))
        for site in candidates:
            old = qa.combined_old(*site)
            metric = qa.distribution(qa.seam_residual(old, direction, qa.valid(old)))
            if (metric['crossings'] >= MIN_CONTROL and metric['p95_abs_m'] <= 0.5
                    and abs(metric['signed_mean_m']) <= 0.2):
                controls.append(site)
                print('CONTROL_SELECTED:', site, metric, flush=True)
                break
        else:
            raise RuntimeError(f'No eligible control found for {direction}')
    return controls


def native_reference(src, bbox):
    xmin, ymin, xmax, ymax = bbox
    h, w = ymax-ymin, xmax-xmin
    rr, cc = np.indices((h, w))
    s = src.GetSpatialRef()
    t = osr.SpatialReference()
    t.ImportFromEPSG(26910)
    s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    t.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    ct = osr.CoordinateTransformation(t, s)
    xy = np.asarray(ct.TransformPoints(list(zip((xmin+cc.ravel()+0.5).tolist(),
                                               (ymax-rr.ravel()-0.5).tolist()))))
    if not np.isfinite(xy[:, :2]).all():
        raise RuntimeError('Nonfinite transformed coordinates')
    inv = gdal.InvGeoTransform(src.GetGeoTransform())
    px = (inv[0]+inv[1]*xy[:, 0]+inv[2]*xy[:, 1]-0.5).reshape(h, w)
    py = (inv[3]+inv[4]*xy[:, 0]+inv[5]*xy[:, 1]-0.5).reshape(h, w)
    ix, iy = np.floor(px).astype(np.int64), np.floor(py).astype(np.int64)
    u, v = px-ix, py-iy
    x0, y0 = int(ix.min()), int(iy.min())
    ww, hh = int(ix.max()-x0+2), int(iy.max()-y0+2)
    if not (x0 >= 0 and y0 >= 0 and x0+ww <= src.RasterXSize and y0+hh <= src.RasterYSize):
        raise RuntimeError('Reference window outside VRT')
    band = src.GetRasterBand(1)
    z = band.ReadAsArray(x0, y0, ww, hh).astype(np.float64)
    sv = qa.valid(z) & (band.GetMaskBand().ReadAsArray(x0, y0, ww, hh) > 0)
    j, i = ix-x0, iy-y0
    ok = sv[i,j] & sv[i,j+1] & sv[i+1,j] & sv[i+1,j+1]
    ref = (1-u)*(1-v)*z[i,j]+u*(1-v)*z[i,j+1]+(1-u)*v*z[i+1,j]+u*v*z[i+1,j+1]
    ref = np.where(ok, ref, NODATA)
    region = '2_4' if y0+hh <= 500000 else ('2_3' if y0 >= 500000 else 'crosses_VRT_boundary')
    return ref, ok, [x0, y0, ww, hh], region


def warp(path, bbox, env, commands):
    cmd = ['gdalwarp', '-q', '-of', 'GTiff', '-t_srs', 'EPSG:26910',
           '-te', *map(str,bbox), '-tr','1','1','-tap','-r','bilinear',
           '-et','0','-ovr','NONE','-srcnodata','-32767','-dstnodata','-32767',
           '-wo','XSCALE=1','-wo','YSCALE=1','-co','TILED=YES','-co','COMPRESS=DEFLATE',
           '-co','PREDICTOR=3','-co','BIGTIFF=IF_SAFER',str(qa.VRT),str(path)]
    commands.append(cmd)
    print('WARPING:', path.name, flush=True)
    subprocess.run(cmd, env=env, check=True)
    return qa.raster(path, *bbox)


def main():
    gdal.UseExceptions()
    qa.validate_vrt()
    opts = {'GDAL_DISABLE_READDIR_ON_OPEN':'YES','GDAL_HTTP_MAX_RETRY':'3',
            'GDAL_HTTP_RETRY_DELAY':'2','GDAL_HTTP_TIMEOUT':'45'}
    env = os.environ.copy(); env.update(opts)
    for k,v in opts.items(): gdal.SetConfigOption(k,v)
    src = gdal.Open(str(qa.VRT), gdal.GA_ReadOnly)
    if (src.RasterXSize,src.RasterYSize) != (500000,1000000) or not np.allclose(src.GetGeoTransform(),(-2000000,1,0,1000000,0,-1),rtol=0,atol=1e-8):
        raise RuntimeError('Unexpected native VRT grid')
    s = src.GetSpatialRef(); s.AutoIdentifyEPSG()
    if s.GetAuthorityCode(None) != '3979' or src.GetRasterBand(1).GetNoDataValue() != NODATA:
        raise RuntimeError('Unexpected VRT CRS or NoData')
    controls = select_controls()
    out = Path(tempfile.mkdtemp(prefix='hrdem1m_multisite_fixed_',dir=qa.DIAG))
    print('OUTPUT_DIRECTORY:',out,flush=True)
    report = {'gdal_version':gdal.VersionInfo('--version'), 'gdalwarp_version':subprocess.check_output(['gdalwarp','--version'],text=True).strip(),
              'tolerance_m':TOL,'min_reference_pixels':MIN_REFERENCE,'reference':'strict four-valid-neighbor bilinear from same HRDEM VRT, not independent terrain truth',
              'control_selection':'nearest eligible old unflagged join to each predetermined target; >=500 crossings',
              'commands':[],'sites':[],'complete':False}
    path = out/'multisite_fixed_qa.json'
    for site in FLAGGED+controls:
        x,y,direction=site
        kind = 'flagged' if site in FLAGGED else 'control'
        key=f'{kind}_x{x}_y{y}_{direction}'
        print('SITE:',key,flush=True)
        narrow=qa.box_for(*site)
        wide=(x-64,y-64,x+1064,y+1064)
        ww=warp(out/(key+'_wide.tif'),wide,env,report['commands'])
        n=warp(out/(key+'_narrow.tif'),narrow,env,report['commands'])
        c0=narrow[0]-wide[0]; r0=wide[3]-narrow[3]
        w=ww[r0:r0+n.shape[0],c0:c0+n.shape[1]]
        if w.shape != n.shape: raise RuntimeError('Common strip outside wide raster')
        ref,ok,window,region=native_reference(src,narrow)
        old=qa.combined_old(*site)
        m=qa.valid(w)&qa.valid(n)&ok
        extent=stats(w,n)
        native_w=stats(w,ref,m,m); native_n=stats(n,ref,m,m)
        compare=qa.valid(old)&m
        seam={'old':qa.distribution(qa.seam_residual(old,direction,compare)),
              'fixed_wide':qa.distribution(qa.seam_residual(w,direction,compare)),
              'fixed_narrow':qa.distribution(qa.seam_residual(n,direction,compare)),
              'native_reference':qa.distribution(qa.seam_residual(ref,direction,compare))}
        record={'site':list(site),'kind':kind,'source_window':window,'source_region':region,
                'wide_vs_narrow':extent,'strict_reference_valid':int(ok.sum()),
                'strict_reference_missing_in_wide':int((ok & ~qa.valid(w)).sum()),
                'strict_reference_missing_in_narrow':int((ok & ~qa.valid(n)).sum()),
                'wide_vs_native_common':native_w,'narrow_vs_native_common':native_n,
                'seam_on_old_and_reference_common_support':seam,
                'extent_consistent':extent['common_valid']>=MIN_REFERENCE and extent['validity_mismatch']==0 and extent['over_0_001m']==0,
                'native_consistent':int(m.sum())>=MIN_REFERENCE and not (ok & ~qa.valid(w)).any() and not (ok & ~qa.valid(n)).any() and native_w['over_0_001m']==0 and native_n['over_0_001m']==0}
        report['sites'].append(record)
        qa.write_json_atomic(path,report)
        print('SITE_RESULT:',json.dumps(record),flush=True)
    report['complete']=True
    report['all_extent_consistent']=all(r['extent_consistent'] for r in report['sites'])
    report['all_native_consistent']=all(r['native_consistent'] for r in report['sites'])
    qa.write_json_atomic(path,report)
    print('ALL_SITES_COMPLETED:',len(report['sites']))
    print('ALL_EXTENT_CONSISTENT:',report['all_extent_consistent'])
    print('ALL_NATIVE_CONSISTENT:',report['all_native_consistent'])
    print('DIAGNOSTIC_JSON:',path)
    print('ORIGINAL_TILES_MODIFIED: FALSE')

if __name__=='__main__':
    main()
