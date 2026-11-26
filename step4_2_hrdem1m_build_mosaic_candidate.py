"""Build a lossless HRDEM mosaic candidate, not an approved final DEM.
Run serially in the project root. No warp, interpolation or input edits.
"""
from pathlib import Path
import json
import os
import shutil
import tempfile
import numpy as np
import step4_2_hrdem1m_fixed_batch_300 as e
import step4_2_hrdem1m_global_masking_qa as q
import step4_2_hrdem1m_anomaly_register as ar

N = 1000
REGISTER = 'hrdem1m_anomaly_register_q_06k0gb/anomaly_register_qa.json'
REGISTER_SHA = 'e2e8840013db83074fc614d3a12d299ee6d199e09465172bf59ffb9e489c76ae'
EXPECTED = 464090223
OPTIONS = ['TILED=YES', 'BLOCKXSIZE=256', 'BLOCKYSIZE=256',
           'COMPRESS=DEFLATE', 'PREDICTOR=3', 'BIGTIFF=YES']


def bounds(keys):
    return min(x for x,y in keys), min(y for x,y in keys), max(x for x,y in keys)+N, max(y for x,y in keys)+N


def strip(y, xmin, xmax, complete, read):
    # Every gap, including the empty extraction, is explicitly written as NoData.
    a = np.full((N, xmax-xmin), q.NODATA, dtype=np.float32)
    for x in range(xmin,xmax,N):
        if (x,y) in complete:
            a[:,x-xmin:x-xmin+N] = read((x,y))
    return a


def verify_strip(actual, y, xmin, xmax, complete, read):
    # Independent tile-window checks, including all noncandidate cells.
    q.require(actual.shape==(N,xmax-xmin) and actual.dtype==np.float32, 'Readback shape/type mismatch')
    valid=0; gaps=0; lo=None; hi=None
    for x in range(xmin,xmax,N):
        k=(x,y); block=actual[:,x-xmin:x-xmin+N]
        q.require(np.isfinite(block).all(), f'Nonfinite mosaic cell: {k}')
        if k in complete:
            expected=read(k)
            q.require(np.array_equal(block.view(np.uint32),expected.view(np.uint32)), f'Bitwise mosaic mismatch: {k}')
            v=block[block!=q.NODATA]; valid+=int(v.size)
            if v.size:
                low,high=float(v.min()),float(v.max())
                lo=low if lo is None else min(lo,low); hi=high if hi is None else max(hi,high)
        else:
            q.require(np.all(block==q.NODATA), f'Unexpected data in gap/empty cell: {k}')
            gaps+=1
    return valid,gaps,lo,hi


def main():
    gdal,ogr,osr=q.L.gdal,q.L.ogr,q.L.osr
    gdal.UseExceptions(); ogr.UseExceptions()
    root=e.legacy.project_root(); coverage,_,src,maskdir=e.paths(root); diag=src.parent
    cfg=e.config(root); rp=diag/REGISTER; register=e.read_json(rp); rs=register['summary']
    q.require(register['configuration']==cfg and register['script_sha256']==e.sha(Path(ar.__file__)), 'Register configuration/script mismatch')
    gpkg=rp.parent/'anomaly_register.gpkg'
    q.require(rs['register_completed'] and not rs['final_dem_approved'] and rs['register_sha256']==REGISTER_SHA==e.sha(gpkg), 'Register hash/status mismatch')
    chain=register['report_chain_sha256']
    for rel,h in chain.items():
        q.require(e.sha(root/rel)==h, f'Prior report changed: {rel}')
    def report_named(name):
        matches=[root/p for p in chain if Path(p).name==name]
        q.require(len(matches)==1, f'Expected exactly one prior report: {name}')
        r=e.read_json(matches[0]); q.require(r['configuration']==cfg, f'Prior config mismatch: {name}'); return r
    mask=report_named('global_masking_qa.json'); terrain=report_named('terrain_sensitivity_qa.json')
    q.require(mask['summary']['masking_qa_passed'] and mask['summary']['masked_valid_pixels']==EXPECTED, 'Masking QA not accepted')
    q.require(terrain['summary']['valid_centers']==EXPECTED and terrain['summary']['strict_3x3_centers']==463636335, 'Terrain QA mismatch')
    q.require(rs['layer_counts']==dict(slope_change_gt10deg=7250,aspect_change_gt90deg=233,elevation_change_gt10m=45,flagged_seams=41,insufficient_seams=52,tiles_with_internal_nodata=265), 'Register layer counts mismatch')
    q.require(register['decision']['retain_elevations_unchanged'] and not register['decision']['exclude_flagged_pixels_automatically'], 'Anomaly handling decision mismatch')
    lookup={(x,y):g for x,y,g in e.legacy.candidate_tiles(coverage)}
    rows={tuple(r['tile']):r for r in mask['tiles']}
    hashes={tuple(r['tile']):r['fixed_sha256'] for r in terrain['verified_hashes']}
    complete={k for k,r in rows.items() if r['state']=='COMPLETE'}
    empty={k for k,r in rows.items() if r['state']=='EMPTY'}
    q.require(len(lookup)==656 and set(lookup)==set(rows)==complete|empty and len(complete)==655 and set(hashes)==complete and empty=={(497000,5549000)}, 'Candidate inventory mismatch')
    q.require({k for k in complete if rows[k]['masked_valid_pixels']==0}=={(484000,5482000)}, 'All-NoData TIFF inventory mismatch')
    inv=q.inventory(maskdir,{e.stem(*k) for k in complete},masked=True)
    q.require(not inv['missing_entries'] and not inv['unexpected_entries'], 'Unexpected mask inventory: '+str(inv))
    def path(k): return maskdir/(e.stem(*k)+'_masked.tif')
    def check_hash(k):
        p=path(k); side=p.with_suffix('.qa.json'); r=e.read_json(side)
        q.require(r['output_sha256']==hashes[k]==e.sha(p) and r['provenance']['extraction_configuration']==cfg and r['masked_valid_pixels']==rows[k]['masked_valid_pixels'], f'Input hash/provenance mismatch: {k}')
    def read(k):
        a,_=q.read_grid(path(k),*k)
        q.require(a.dtype==np.float32 and int((a!=q.NODATA).sum())==rows[k]['masked_valid_pixels'],f'Input pixel count/type mismatch: {k}')
        return a
    for i,k in enumerate(sorted(complete),1):
        check_hash(k)
        if i%100==0 or i==len(complete): print(f'INPUT_HASH_PROGRESS: {i}/{len(complete)}',flush=True)
    xmin,ymin,xmax,ymax=bounds(lookup); w=xmax-xmin; h=ymax-ymin
    q.require(w%N==0 and h%N==0 and all((x-xmin)%N==0 and (y-ymin)%N==0 for x,y in lookup), 'Off-grid candidate')
    cells=(w//N)*(h//N); total=w*h; raw=total*4
    free=shutil.disk_usage(diag).free
    q.require(free>raw+2*1024**3, f'Insufficient conservative free-space reserve: {free} bytes')
    gt=(xmin,1,0,ymax,0,-1); crs=osr.SpatialReference(); crs.ImportFromEPSG(26910)
    print('MOSAIC_PLAN:',json.dumps(dict(bounds=[xmin,ymin,xmax,ymax],width=w,height=h,grid_pixels=total,candidate_cells=len(lookup),noncandidate_cells=cells-len(lookup),uncompressed_bytes=raw,free_bytes=free)),flush=True)
    folder=Path(tempfile.mkdtemp(prefix='hrdem1m_mosaic_candidate_',dir=diag)); tmp=folder/'hrdem1m_fixed_masked_mosaic.partial.tif'; out=folder/'hrdem1m_fixed_masked_mosaic_candidate.tif'
    ds=gdal.GetDriverByName('GTiff').Create(str(tmp),w,h,1,gdal.GDT_Float32,options=OPTIONS)
    q.require(ds is not None,'Cannot create candidate TIFF')
    ds.SetGeoTransform(gt); ds.SetProjection(crs.ExportToWkt()); band=ds.GetRasterBand(1); band.SetNoDataValue(q.NODATA)
    for i,y in enumerate(range(ymax-N,ymin-1,-N),1):
        a=strip(y,xmin,xmax,complete,read)
        q.require(band.WriteArray(a,0,ymax-y-N)==0, f'Write failed at y={y}')
        if i%10==0 or i==h//N: print(f'MOSAIC_WRITE_PROGRESS: {i}/{h//N}',flush=True)
    band.FlushCache(); band=None; ds.FlushCache(); ds=None
    # Reopen the physical GeoTIFF and read every pixel, without resampling.
    ds=gdal.Open(str(tmp),gdal.GA_ReadOnly)
    q.require(ds is not None and (ds.RasterXSize,ds.RasterYSize,ds.RasterCount)==(w,h,1) and ds.GetGeoTransform()==gt and ds.GetSpatialRef().IsSame(crs),'Mosaic grid/CRS mismatch')
    band=ds.GetRasterBand(1); q.require(band.DataType==gdal.GDT_Float32 and band.GetNoDataValue()==q.NODATA,'Mosaic dtype/NoData mismatch')
    count=gaps=0; minimum=maximum=None
    for i,y in enumerate(range(ymax-N,ymin-1,-N),1):
        a=band.ReadAsArray(0,ymax-y-N,w,N)
        q.require(a is not None,'Mosaic read failed')
        n,g,lo,hi=verify_strip(a,y,xmin,xmax,complete,read); count+=n; gaps+=g
        if lo is not None:
            minimum=lo if minimum is None else min(minimum,lo); maximum=hi if maximum is None else max(maximum,hi)
        if i%10==0 or i==h//N: print(f'MOSAIC_READBACK_PROGRESS: {i}/{h//N}',flush=True)
    band=None; ds=None
    q.require(count==EXPECTED and gaps==cells-len(complete),'Mosaic valid/gap accounting mismatch')
    # Confirm the source inputs stayed unchanged throughout construction/readback.
    for k in sorted(complete): check_hash(k)
    q.require(e.sha(gpkg)==REGISTER_SHA and all(e.sha(root/p)==v for p,v in chain.items()),'QA inputs changed during construction')
    q.require(not any(Path(str(tmp)+s).exists() for s in ('.aux.xml','.msk','.ovr')),'Unexpected TIFF sidecar; investigate')
    os.replace(tmp,out)
    summary=dict(candidate_tiles=len(lookup),input_tiffs=len(complete),noncandidate_cells=cells-len(lookup),empty_extraction_cells=len(empty),all_nodata_input_tiffs=1,
        bounds_epsg26910=[xmin,ymin,xmax,ymax],width=w,height=h,geotransform=list(gt),crs='EPSG:26910',dtype='Float32',nodata=q.NODATA,
        total_grid_pixels=total,valid_pixels=count,nodata_pixels=total-count,minimum_elevation_m=minimum,maximum_elevation_m=maximum,
        input_tile_bitwise_mismatch_pixels=0,unexpected_valid_in_gap_cells=0,all_grid_pixels_readback_verified=True,input_hashes_verified_before_and_after=True,
        raster_bytes=out.stat().st_size,mosaic_sha256=e.sha(out),mosaic_candidate_qa_passed=True,final_dem_approved=False)
    report=folder/'mosaic_candidate_qa.json'
    e.save_json(report,dict(configuration=cfg,script_sha256=e.sha(Path(__file__)),register_report=str(rp.relative_to(root)),register_report_sha256=e.sha(rp),register_sha256=REGISTER_SHA,
        report_chain_sha256=chain,input_tile_sha256=[dict(tile=list(k),sha256=hashes[k]) for k in sorted(complete)],creation_options=OPTIONS,summary=summary,
        mosaic=str(out.relative_to(root)),limitations=['Candidate bounding rectangle includes all candidate cells and background NoData; it is not the exact vector AOI.','No new resampling, filling or smoothing; no claim of source vertical accuracy.','Anomaly register remains separate; physical landslide model sensitivity and final publication remain pending.']))
    print('===== MOSAIC CANDIDATE SUMMARY ====='); print(json.dumps(summary,indent=2)); print('MOSAIC_CANDIDATE:',out); print('DIAGNOSTIC_JSON:',report); print('MOSAIC_CANDIDATE_QA_PASSED: True'); print('FINAL_DEM: FALSE')


if __name__=='__main__': main()
