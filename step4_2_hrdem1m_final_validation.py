"""Validate the immutable HRDEM mosaic candidate before publication.
Read-only raster/vector audit; write a new diagnostic JSON only. No remote reads.
Negative elevations are reported, not corrected or classified as errors.
"""
from pathlib import Path
import json
import tempfile
import numpy as np
import step4_2_hrdem1m_fixed_batch_300 as e
import step4_2_hrdem1m_global_masking_qa as q
import step4_2_hrdem1m_build_mosaic_candidate as b

PRIOR = 'hrdem1m_mosaic_candidate_hjcps_iq/mosaic_candidate_qa.json'
MOSAIC_SHA = 'eedd95b99c4c0a56067a1c15dcce413cd8ef64b56d68b7145ec159a94643ffb6'
REGISTER_SHA = 'e2e8840013db83074fc614d3a12d299ee6d199e09465172bf59ffb9e489c76ae'
N = 1000


def block_stats(a, inside, x, y, nodata):
    if a.shape != inside.shape or not np.isfinite(a).all():
        raise RuntimeError(f'Invalid mosaic block: {(x,y)}')
    valid = a != nodata
    negative = valid & (a < 0)
    n = int(valid.sum())
    def extreme(low):
        if not n:
            return None
        v = np.where(valid, a, np.inf if low else -np.inf)
        i = int(v.argmin() if low else v.argmax())
        r,c = divmod(i,a.shape[1])
        return dict(value_m=float(a[r,c]), center_epsg26910=[x+c+.5,y+a.shape[0]-r-.5],tile=[x,y])
    return dict(valid=n, footprint=int(inside.sum()), missing=int((inside & ~valid).sum()),
                outside_valid=int((valid & ~inside).sum()), negative=int(negative.sum()),
                below_minus_1=int((valid & (a < -1)).sum()),
                below_minus_5=int((valid & (a < -5)).sum()),
                below_minus_10=int((valid & (a < -10)).sum()),
                below_minus_20=int((valid & (a < -20)).sum()),
                minimum=extreme(True), maximum=extreme(False)), negative


def main():
    gdal,ogr,osr=q.L.gdal,q.L.ogr,q.L.osr
    gdal.UseExceptions(); ogr.UseExceptions()
    gdal.SetConfigOption('GDAL_PAM_ENABLED','NO')
    root=e.legacy.project_root(); coverage,_,src,maskdir=e.paths(root); diag=src.parent
    cfg=e.config(root); p=diag/PRIOR; digest=e.sha(p); prior=e.read_json(p); s=prior['summary']
    q.require(prior['configuration']==cfg and prior['script_sha256']==e.sha(Path(b.__file__)), 'Candidate configuration/script mismatch')
    q.require(s['mosaic_candidate_qa_passed'] and s['all_grid_pixels_readback_verified'] and s['input_hashes_verified_before_and_after'] and not s['final_dem_approved'], 'Prior mosaic QA not accepted')
    q.require(s['input_tile_bitwise_mismatch_pixels']==s['unexpected_valid_in_gap_cells']==0,'Prior mismatch')
    path=root/prior['mosaic']; q.require(path.parent==p.parent and e.sha(path)==s['mosaic_sha256']==MOSAIC_SHA and path.stat().st_size==s['raster_bytes'],'Mosaic identity mismatch')
    rp=root/prior['register_report']; q.require(e.sha(rp)==prior['register_report_sha256'],'Register report changed')
    register=e.read_json(rp); rs=register['summary']; gpkg=rp.parent/'anomaly_register.gpkg'
    q.require(register['configuration']==cfg and rs['register_completed'] and rs['register_sha256']==prior['register_sha256']==REGISTER_SHA==e.sha(gpkg),'Register mismatch')
    chain=prior['report_chain_sha256']; q.require(chain==register['report_chain_sha256'],'Report chain mismatch')
    for rel,h in chain.items(): q.require(e.sha(root/rel)==h,f'Prior report changed: {rel}')
    matches=[root/rel for rel in chain if Path(rel).name=='global_masking_qa.json']
    q.require(len(matches)==1,'Missing/ambiguous masking QA'); mask=e.read_json(matches[0]); ms=mask['summary']
    q.require(mask['configuration']==cfg and ms['masking_qa_passed'],'Masking report mismatch')
    lookup={(x,y):g for x,y,g in e.legacy.candidate_tiles(coverage)}
    rows={tuple(r['tile']):r for r in mask['tiles']}; complete={k for k,r in rows.items() if r['state']=='COMPLETE'}
    hashes={tuple(r['tile']):r['sha256'] for r in prior['input_tile_sha256']}
    q.require(len(lookup)==656 and set(lookup)==set(rows) and len(complete)==655 and set(hashes)==complete,'Candidate inventory mismatch')
    q.require(set(rows)-complete=={(497000,5549000)},'Empty candidate mismatch')
    # Inspect actual register layers, independently of the previous reported counts.
    crs=osr.SpatialReference(); crs.ImportFromEPSG(26910)
    crs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    expected={'slope_change_gt10deg':(7250,ogr.wkbPoint),'aspect_change_gt90deg':(233,ogr.wkbPoint),'elevation_change_gt10m':(45,ogr.wkbPoint),'flagged_seams':(41,ogr.wkbLineString),'insufficient_seams':(52,ogr.wkbLineString),'tiles_with_internal_nodata':(265,ogr.wkbPolygon)}
    for suffix in ('-wal','-journal'): q.require(not Path(str(gpkg)+suffix).exists(),'Register database sidecar present')
    vec=ogr.Open(str(gpkg),0); q.require(vec is not None and vec.GetLayerCount()==len(expected),'Register layer inventory mismatch')
    layer_counts={}
    for name,(count,geomtype) in expected.items():
        ly=vec.GetLayerByName(name)
        q.require(ly is not None and ly.GetSpatialRef() is not None and ly.GetSpatialRef().IsSame(crs) and ly.GetGeomType()==geomtype,f'Register layer geometry/CRS mismatch: {name}')
        q.require(ly.GetFeatureCount()==count,f'Register feature count mismatch: {name}'); actual=0
        for f in ly:
            geom=f.GetGeometryRef(); q.require(geom is not None and not geom.IsEmpty() and geom.GetGeometryType()==geomtype,f'Invalid register feature: {name}'); actual+=1
        q.require(actual==count,f'Register iteration count mismatch: {name}'); layer_counts[name]=actual
        f=None; geom=None; ly=None
    vec=None
    q.require(layer_counts==rs['layer_counts'],'Register count/report mismatch')
    print('REGISTER_LAYER_COUNTS:',json.dumps(layer_counts),flush=True)
    xmin,ymin,xmax,ymax=b.bounds(lookup); w=xmax-xmin; h=ymax-ymin; gt=(xmin,1,0,ymax,0,-1)
    q.require((xmin,ymin,xmax,ymax)==(475000,5464000,506000,5557000) and list(gt)==s['geotransform'],'Unexpected extent')
    for suffix in ('.aux.xml','.msk','.ovr'): q.require(not Path(str(path)+suffix).exists(),f'Unexpected mosaic sidecar: {suffix}')
    ds=gdal.Open(str(path),gdal.GA_ReadOnly)
    q.require(ds is not None and ds.GetDriver().ShortName=='GTiff' and (ds.RasterXSize,ds.RasterYSize,ds.RasterCount)==(w,h,1) and ds.GetGeoTransform()==gt and ds.GetSpatialRef().IsSame(crs),'Mosaic format/grid/CRS mismatch')
    band=ds.GetRasterBand(1); q.require(band.DataType==gdal.GDT_Float32 and band.GetNoDataValue()==q.NODATA and band.GetMaskFlags()==gdal.GMF_NODATA,'Mosaic dtype/NoData/mask semantics mismatch')
    q.require(band.GetBlockSize()==[256,256] or tuple(band.GetBlockSize())==(256,256),'Unexpected block size')
    image_structure=ds.GetMetadata('IMAGE_STRUCTURE'); q.require(image_structure.get('COMPRESSION')=='DEFLATE','Unexpected compression')
    with path.open('rb') as fh: magic=fh.read(4)
    q.require(magic in (b'II+\x00',b'MM\x00+'),'Not a BigTIFF')
    total=dict(valid=0,footprint=0,missing=0,outside_valid=0,negative=0,below_minus_1=0,below_minus_5=0,below_minus_10=0,below_minus_20=0)
    minimum=maximum=None; negative_tiles=[]; negative_sources=0; cells=0
    for i,y in enumerate(range(ymax-N,ymin-1,-N),1):
        strip=band.ReadAsArray(0,ymax-y-N,w,N)
        q.require(strip is not None and strip.shape==(N,w) and strip.dtype==np.float32 and np.isfinite(strip).all(),f'Invalid strip: {y}')
        for x in range(xmin,xmax,N):
            k=(x,y); a=strip[:,x-xmin:x-xmin+N]; cells+=1
            if k not in lookup:
                q.require(np.all(a==q.NODATA),f'Unexpected noncandidate data: {k}'); continue
            inside=q.empty_footprint(coverage,x,y); d,neg=block_stats(a,inside,x,y,q.NODATA); r=rows[k]
            q.require(d['outside_valid']==0 and d['valid']==r['masked_valid_pixels'] and d['footprint']==r['footprint_pixels'] and d['missing']==r['nodata_inside_footprint'],f'Footprint or count mismatch: {k}')
            for key in total: total[key]+=d[key]
            if d['minimum'] is not None:
                if minimum is None or d['minimum']['value_m']<minimum['value_m']: minimum=d['minimum']
                if maximum is None or d['maximum']['value_m']>maximum['value_m']: maximum=d['maximum']
            if d['negative']:
                # Trace negatives to the hashed local fixed-scale extraction, not independent ground truth.
                status,_=e.state(root,x,y,lookup[k],cfg); q.require(status=='COMPLETE',f'Invalid extraction state: {k}')
                fp=maskdir/(e.stem(x,y)+'_masked.tif'); q.require(e.sha(fp)==hashes[k],f'Masked input changed: {k}')
                unmasked,_=q.read_grid(src/(e.stem(x,y)+'.tif'),x,y); masked,_=q.read_grid(fp,x,y)
                q.require(np.array_equal(a[neg].view(np.uint32),unmasked[neg].view(np.uint32)) and np.array_equal(a[neg].view(np.uint32),masked[neg].view(np.uint32)),f'Negative elevation transfer mismatch: {k}')
                negative_sources+=d['negative']; negative_tiles.append(dict(tile=list(k),source_group=lookup[k],**{key:d[key] for key in ('negative','below_minus_1','below_minus_5','below_minus_10','below_minus_20','minimum')}))
        if i%10==0 or i==h//N: print(f'FINAL_READ_PROGRESS: {i}/{h//N}',flush=True)
    band=None; ds=None
    q.require(cells==2883 and w*h==2883000000 and total['valid']==ms['masked_valid_pixels']==464090223 and total['footprint']==ms['footprint_pixels']==464219963 and total['missing']==ms['nodata_inside_footprint']==129740 and total['outside_valid']==0,'Global spatial accounting mismatch')
    q.require(minimum['value_m']==s['minimum_elevation_m'] and maximum['value_m']==s['maximum_elevation_m'] and negative_sources==total['negative'],'Elevation accounting mismatch')
    q.require(e.sha(path)==MOSAIC_SHA and e.sha(p)==digest and e.sha(gpkg)==REGISTER_SHA and e.sha(rp)==prior['register_report_sha256'] and e.sha(coverage)==cfg['footprint_sha256'] and all(e.sha(root/rel)==v for rel,v in chain.items()),'QA inputs changed during validation')
    negative_tiles.sort(key=lambda d:d['minimum']['value_m'])
    summary=dict(grid_pixels=w*h,candidate_cells=len(lookup),valid_pixels=total['valid'],nodata_pixels=w*h-total['valid'],footprint_pixels=total['footprint'],nodata_inside_footprint=total['missing'],valid_outside_footprint=total['outside_valid'],footprint_coverage_percent=100*total['valid']/total['footprint'],minimum=minimum,maximum=maximum,negative_pixels=total['negative'],negative_tiles=len(negative_tiles),negative_counts_below_m={str(-k):total[f'below_minus_{k}'] for k in (1,5,10,20)},negative_pixels_matching_local_fixed_extraction_and_mask=negative_sources,register_layer_counts=layer_counts,mosaic_sha256=MOSAIC_SHA,all_pixels_scanned=True,final_technical_validation_passed=True,final_dem_published=False,absolute_vertical_accuracy_validated=False,negative_source_physical_validity_resolved=False)
    folder=Path(tempfile.mkdtemp(prefix='hrdem1m_final_validation_',dir=diag)); report=folder/'final_validation_qa.json'
    e.save_json(report,dict(configuration=cfg,script_sha256=e.sha(Path(__file__)),candidate_report=str(p.relative_to(root)),candidate_report_sha256=digest,mosaic=str(path.relative_to(root)),register=str(gpkg.relative_to(root)),register_sha256=REGISTER_SHA,summary=summary,negative_tile_details=negative_tiles,image_structure=image_structure,limitations=['Technical validation only. Negative elevations were traced to local fixed extraction, not independently validated against native COG or ground measurements.','Source anomalies remain unresolved; anomaly register is diagnostic, not an exclusion mask.','No publication, interpolation or correction performed.']))
    print('===== FINAL TECHNICAL VALIDATION SUMMARY ====='); print(json.dumps(summary,indent=2))
    print('===== TOP 10 NEGATIVE TILES BY MINIMUM =====')
    for r in negative_tiles[:10]: print(json.dumps(r))
    print('DIAGNOSTIC_JSON:',report); print('FINAL_TECHNICAL_VALIDATION_PASSED: True'); print('FINAL_DEM_PUBLISHED: FALSE')


if __name__=='__main__': main()
