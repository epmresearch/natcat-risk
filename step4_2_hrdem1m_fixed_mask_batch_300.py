"""Mask only provenance-verified fixed-scale HRDEM tiles. Default: read-only plan."""
from __future__ import annotations
import argparse
import json
import os
import sys
from pathlib import Path
import step4_2_hrdem1m_fixed_batch_300 as e
import step4_2_hrdem1m_mask_batch_300 as legacy

MASK_ID='HRDEM1M_FIXED_PIXEL_CENTER_MASK_V1'

def mask_one(root,x,y,group,cfg):
    coverage,_,srcdir,out=e.paths(root); base=e.stem(x,y)
    status,record=e.state(root,x,y,group,cfg)
    if status=='PENDING': raise RuntimeError(f'Extraction pending: {base}')
    tif=out/(base+'_masked.tif'); q=out/(base+'_masked.qa.json'); tmp=out/(base+'_masked.partial.tif')
    if tmp.exists() or q.with_name(q.name+'.partial').exists(): raise RuntimeError(f'Incomplete mask output: {base}')
    if status=='EMPTY':
        if tif.exists() or q.exists(): raise RuntimeError(f'Mask output for empty extraction: {base}')
        print('SKIPPED_EMPTY_EXTRACTION:',base); return 'EMPTY'
    source=srcdir/(base+'.tif')
    provenance={'mask_method_id':MASK_ID,'extraction_configuration':cfg,
                'source_sha256':record['output_sha256'],'source_qa_sha256':e.sha(srcdir/(base+'.qa.json')),
                'footprint_sha256':cfg['footprint_sha256'],
                'legacy_mask_sha256':e.sha(Path(legacy.__file__)), 'fixed_mask_sha256':e.sha(Path(__file__))}
    expected,stats,projection,geo=legacy.inspect_source_and_mask(source,coverage,x,y)
    if tif.exists()!=q.exists(): raise RuntimeError(f'Incomplete masked TIFF/QA pair: {base}')
    if tif.exists():
        previous=e.read_json(q)
        if (previous.get('provenance')!=provenance or previous.get('tile_x')!=x or previous.get('tile_y')!=y
            or previous.get('status')!='MASKED_STAGING_QA_PASSED_WITH_CAVEATS' or previous.get('output_sha256')!=e.sha(tif)
            or any(previous.get(k)!=v for k,v in stats.items())):
            raise RuntimeError(f'Masked provenance or stats mismatch: {base}')
        legacy.verify_masked_file(tif,expected,x,y)
        print('SKIPPED_VERIFIED_MASK:',base); return 'SKIPPED'
    out.mkdir(parents=True,exist_ok=True)
    ds=None
    try:
        ds=legacy.gdal.GetDriverByName('GTiff').Create(str(tmp),1000,1000,1,legacy.gdal.GDT_Float32,options=['TILED=YES','COMPRESS=DEFLATE','PREDICTOR=3'])
        if ds is None: raise RuntimeError(f'Cannot create {tmp}')
        ds.SetGeoTransform(geo); ds.SetProjection(projection); ds.GetRasterBand(1).SetNoDataValue(legacy.NODATA)
        ds.GetRasterBand(1).WriteArray(expected); ds.FlushCache(); ds=None
        legacy.verify_masked_file(tmp,expected,x,y)
        digest=e.sha(tmp)
        if tif.exists() or q.exists(): raise RuntimeError('Unexpected existing output; run serially')
        os.replace(tmp,tif)
        e.save_json(q,{'status':'MASKED_STAGING_QA_PASSED_WITH_CAVEATS','tile_x':x,'tile_y':y,
                       'source_staging_tif':str(source.relative_to(root)),'mask_source':str(coverage.relative_to(root)),
                       'mask_rule':'pixel-center rasterization of AOI-clipped official HRDEM footprint',
                       'crs':'EPSG:26910','pixel_size_m':1,'provenance':provenance,'output_sha256':digest,
                       'original_values_unchanged_for_retained_pixels':True,**stats,
                       'caveats':['Missing heights preserved as NoData','Staging only; not vertical accuracy or final DEM']})
        print('FIXED_MASK_QA_PASSED:',base,json.dumps(stats),flush=True); return 'SUCCESS'
    finally:
        ds=None
        if tmp.exists(): tmp.unlink()

def main():
    p=argparse.ArgumentParser(description=__doc__); mode=p.add_mutually_exclusive_group()
    mode.add_argument('--pilot',action='store_true'); mode.add_argument('--tile',nargs=2,type=int); mode.add_argument('--batch',type=int); a=p.parse_args()
    if a.batch is not None and not 1<=a.batch<=300: p.error('--batch must be 1..300')
    legacy.gdal.UseExceptions(); legacy.ogr.UseExceptions(); root=e.legacy.project_root()
    coverage,_,srcdir,out=e.paths(root); plan=e.legacy.candidate_tiles(coverage); lookup={(x,y):g for x,y,g in plan}
    if len(plan)!=656 or len(lookup)!=656 or any(k not in lookup for k in e.PILOT): raise RuntimeError('Candidate inventory changed')
    cfg=e.config(root)
    print('MASK_METHOD:',MASK_ID); print('INPUT_DIRECTORY:',srcdir); print('OUTPUT_DIRECTORY:',out)
    print('FIXED_STAGING_TIFS:',len(list(srcdir.glob('*_1km.tif')))); print('FIXED_MASKED_TIFS:',len(list(out.glob('*_masked.tif'))))
    if not a.pilot and a.tile is None and a.batch is None:
        print('PLAN_ONLY: no network raster reads or output writes'); return
    if a.pilot: chosen=[(x,y,lookup[x,y]) for x,y in e.PILOT]
    elif a.tile is not None:
        x,y=a.tile; chosen=[(x,y,lookup[x,y])]
    else:
        chosen=[]
        for x,y,g in plan:
            status,_=e.state(root,x,y,g,cfg)
            base=e.stem(x,y); tif=out/(base+'_masked.tif'); q=out/(base+'_masked.qa.json')
            if (out/(base+'_masked.partial.tif')).exists() or q.with_name(q.name+'.partial').exists():
                raise RuntimeError(f'Incomplete mask output: {base}')
            if status in ('PENDING','EMPTY'):
                if tif.exists() or q.exists(): raise RuntimeError(f'Mask without nonempty extraction: {base}')
                continue
            if tif.exists() or q.exists(): mask_one(root,x,y,g,cfg); continue
            chosen.append((x,y,g))
            if len(chosen)==a.batch: break
    counts={'SUCCESS':0,'EMPTY':0,'SKIPPED':0}
    for x,y,g in chosen: counts[mask_one(root,x,y,g,cfg)]+=1
    print('FIXED_MASK_FINISHED:',counts); print('FINAL_DEM: FALSE')

if __name__=='__main__':
    try: main()
    except (Exception,KeyboardInterrupt) as ex:
        print('ERROR:',ex,file=sys.stderr); sys.exit(1)
