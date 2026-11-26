"""Fixed-scale HRDEM extraction. Default: local, read-only plan. Run serially."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import step4_2_hrdem1m_batch_vrt_300 as legacy

METHOD = 'HRDEM1M_FIXED_BILINEAR_V1'
PILOT = [(485000,5490000),(486000,5490000),(485000,5489000),
         (488000,5503000),(489000,5503000),(490000,5503000),
         (481000,5481000),(497000,5549000)]
WARP = ['-of','GTiff','-t_srs','EPSG:26910','-tr','1','1','-tap',
        '-r','bilinear','-et','0','-ovr','NONE','-wo','XSCALE=1',
        '-wo','YSCALE=1','-ot','Float32','-srcnodata','-32767',
        '-dstnodata','-32767','-co','TILED=YES','-co','COMPRESS=DEFLATE',
        '-co','PREDICTOR=3','-co','BIGTIFF=IF_SAFER']
EMPTY = 'EMPTY_RASTER_NO_VALID_PIXELS'
SUCCESS = 'STRUCTURAL_QA_PASSED_WITH_CAVEATS'

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):
            h.update(b)
    return h.hexdigest()

def paths(root):
    d=root/'results/diagnostics'
    return (root/'data/processed/elevation/hrdem_1m_dtm_coverage_highway99.gpkg',
            d/'step4_2_hrdem1m_boundary_virtual.vrt',
            d/'hrdem1m_fixed_batch_staging',d/'hrdem1m_fixed_masked_staging')

def stem(x,y):
    return f'hrdem1m_x{x}_y{y}_utm10_1km'

def config(root):
    coverage,vrt,_,_=paths(root)
    legacy.check_vrt_file(vrt)
    # A local WAL could make a main-file hash insufficient for provenance.
    for suffix in ('-wal','-journal'):
        if Path(str(coverage)+suffix).exists():
            raise RuntimeError('Footprint has a database sidecar; close writers first')
    cli=subprocess.check_output(['gdalwarp','--version'],text=True).strip()
    version=legacy.gdal.VersionInfo('RELEASE_NAME')
    if not cli.startswith('GDAL '+version+','):
        raise RuntimeError(f'Python/CLI GDAL version mismatch: {version} / {cli}')
    return {'method_id':METHOD,'warp_options':WARP,'gdal_python_version':version,
            'gdal_cli_version':cli,'vrt_sha256':sha(vrt),'footprint_sha256':sha(coverage),
            'legacy_extraction_sha256':sha(Path(legacy.__file__)),
            'fixed_extraction_sha256':sha(Path(__file__))}

def command(vrt,partial,x,y):
    return ['gdalwarp',*WARP,'-te',str(x),str(y),str(x+1000),str(y+1000),str(vrt),str(partial)]

def read_json(p):
    return json.loads(p.read_text(encoding='utf-8'))

def save_json(p,record):
    tmp=p.with_name(p.name+'.partial')
    if p.exists() or tmp.exists():
        raise RuntimeError(f'Existing QA or partial; investigate: {p}')
    with tmp.open('x',encoding='utf-8') as f:
        json.dump(record,f,indent=2,allow_nan=False); f.write('\n')
    os.replace(tmp,p)

def state(root,x,y,group,cfg):
    _,vrt,out,_=paths(root)
    base=stem(x,y); tif=out/(base+'.tif'); q=out/(base+'.qa.json')
    partial=out/(base+'.partial.tif')
    if partial.exists() or q.with_name(q.name+'.partial').exists():
        raise RuntimeError(f'Incomplete output; investigate before retry: {base}')
    if not q.exists():
        if tif.exists(): raise RuntimeError(f'TIFF without QA: {tif}')
        return 'PENDING',None
    r=read_json(q)
    if (r.get('configuration')!=cfg or r.get('tile_x')!=x or r.get('tile_y')!=y
        or r.get('source_group')!=group or r.get('warp_command')!=command(vrt,partial,x,y)):
        raise RuntimeError(f'QA provenance mismatch: {q}')
    if r.get('status')==EMPTY:
        if tif.exists() or r.get('valid_pixels')!=0 or r.get('nodata_pixels')!=1000000 or r.get('output_sha256') is not None:
            raise RuntimeError(f'Invalid empty marker: {q}')
        return 'EMPTY',r
    if r.get('status')!=SUCCESS or not tif.exists() or r.get('output_sha256')!=sha(tif):
        raise RuntimeError(f'Invalid TIFF/QA pair: {q}')
    qa=legacy.inspect_tiff(tif,x,y)
    if qa['valid_pixels']==0 or any(r.get(k)!=v for k,v in qa.items()):
        raise RuntimeError(f'QA stats mismatch: {q}')
    return 'COMPLETE',r

def process(root,x,y,group,cfg):
    status,_=state(root,x,y,group,cfg)
    if status!='PENDING':
        print('SKIPPED_VERIFIED:',stem(x,y),status,flush=True); return 'SKIPPED'
    _,vrt,out,_=paths(root); out.mkdir(parents=True,exist_ok=True)
    base=stem(x,y); tif=out/(base+'.tif'); q=out/(base+'.qa.json'); tmp=out/(base+'.partial.tif')
    cmd=command(vrt,tmp,x,y)
    env=os.environ.copy(); env.update({'GDAL_DISABLE_READDIR_ON_OPEN':'YES','GDAL_HTTP_MAX_RETRY':'3','GDAL_HTTP_RETRY_DELAY':'2','GDAL_HTTP_TIMEOUT':'45'})
    print('EXTRACTING_FIXED:',base,group,flush=True)
    try:
        subprocess.run(cmd,check=True,env=env)
        qa=legacy.inspect_tiff(tmp,x,y)
        record={'tile_x':x,'tile_y':y,'source_group':group,'configuration':cfg,
                'warp_command':cmd,'source_mode':'SINGLE_TWO_SOURCE_VRT',**qa}
        if qa['valid_pixels']==0:
            save_json(q,{**record,'status':EMPTY,'output_sha256':None})
            print('EMPTY_FIXED:',base,flush=True); return 'EMPTY'
        digest=sha(tmp)
        if tif.exists() or q.exists(): raise RuntimeError('Unexpected existing output; run serially')
        os.replace(tmp,tif)
        save_json(q,{**record,'status':SUCCESS,'output_sha256':digest})
        print('FIXED_STRUCTURAL_QA_PASSED:',base,'valid_pixels:',qa['valid_pixels'],flush=True)
        return 'SUCCESS'
    finally:
        if tmp.exists(): tmp.unlink()

def main():
    p=argparse.ArgumentParser(description=__doc__); mode=p.add_mutually_exclusive_group()
    mode.add_argument('--pilot',action='store_true'); mode.add_argument('--tile',nargs=2,type=int)
    mode.add_argument('--batch',type=int); a=p.parse_args()
    if a.batch is not None and not 1<=a.batch<=300: p.error('--batch must be 1..300')
    legacy.gdal.UseExceptions(); legacy.ogr.UseExceptions(); root=legacy.project_root()
    coverage,vrt,out,_=paths(root); plan=legacy.candidate_tiles(coverage); lookup={(x,y):g for x,y,g in plan}
    if len(plan)!=656 or len(lookup)!=656 or any(k not in lookup for k in PILOT): raise RuntimeError('Candidate inventory changed')
    cfg=config(root)
    print('METHOD:',METHOD); print('CANDIDATES:',len(plan)); print('OUTPUT_DIRECTORY:',out)
    print('CONFIGURATION:',json.dumps(cfg)); print('PILOT_TILES:',[(x,y,lookup[x,y]) for x,y in PILOT])
    print('EXISTING_FIXED_TIFS:',len(list(out.glob('*_1km.tif')))); print('EXISTING_FIXED_QA:',len(list(out.glob('*_1km.qa.json'))))
    if not a.pilot and a.tile is None and a.batch is None:
        print('PLAN_ONLY: no network raster reads or output writes'); return
    if a.pilot: chosen=[(x,y,lookup[x,y]) for x,y in PILOT]
    elif a.tile is not None:
        x,y=a.tile; chosen=[(x,y,lookup[x,y])]
    else:
        chosen=[]
        for x,y,g in plan:
            if state(root,x,y,g,cfg)[0]=='PENDING': chosen.append((x,y,g))
            if len(chosen)==a.batch: break
    counts={'SUCCESS':0,'EMPTY':0,'SKIPPED':0}
    for x,y,g in chosen: counts[process(root,x,y,g,cfg)]+=1
    print('FIXED_EXTRACTION_FINISHED:',counts); print('FINAL_DEM: FALSE')

if __name__=='__main__':
    try: main()
    except (Exception,KeyboardInterrupt) as e:
        print('ERROR:',e,file=sys.stderr); sys.exit(1)
