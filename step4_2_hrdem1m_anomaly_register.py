"""Spatial diagnostic register for the selected fixed-scale HRDEM. No raster writes.
The output layers are review indicators, not an exclusion mask or a terrain correction.
"""
from pathlib import Path
from functools import lru_cache
from collections import Counter
import json
import tempfile
import os
import numpy as np
import step4_2_hrdem1m_fixed_batch_300 as e
import step4_2_hrdem1m_global_masking_qa as q
import step4_2_hrdem1m_terrain_sensitivity_qa as t

PRIOR = 'hrdem1m_terrain_sensitivity_u4w4q9l1/terrain_sensitivity_qa.json'
N = 1000


def key(r):
    return (*r['tile'], r['direction'])


def point_xy(x, y, row, col):
    return x+int(col)+0.5, y+N-int(row)-0.5


def seam_endpoints(x, y, direction):
    if direction == 'east':
        return [(x+N,y),(x+N,y+N)]
    if direction == 'north':
        return [(x,y+N),(x+N,y+N)]
    raise ValueError(direction)


def localization(r):
    ratio = r['joint_center_to_max_offset_p95_ratio']
    if r['joint_crossings'] < 20 or r['joint_fraction_of_center'] < .8 or ratio is None:
        return 'insufficient_support'
    return 'strong' if ratio >= 2 else 'moderate' if ratio > 1 else 'no_dominant_peak'


def main():
    gdal,ogr,osr=q.L.gdal,q.L.ogr,q.L.osr
    gdal.UseExceptions(); ogr.UseExceptions()
    root=e.legacy.project_root(); coverage,_,src,maskdir=e.paths(root)
    diag=src.parent; cfg=e.config(root); pp=diag/PRIOR
    terrain=e.read_json(pp)
    q.require(terrain['configuration']==cfg and terrain['script_sha256']==e.sha(Path(t.__file__)), 'Terrain report/config mismatch')
    chain={str(pp.relative_to(root)):e.sha(pp)}
    def previous(parent,pathkey,hashkey):
        p=root/parent[pathkey]; digest=e.sha(p)
        q.require(digest==parent[hashkey],f'Prior report hash mismatch: {p}')
        chain[str(p.relative_to(root))]=digest
        r=e.read_json(p); q.require(r['configuration']==cfg,f'Prior config mismatch: {p}')
        return r
    native=previous(terrain,'prior_targeted_source_qa','prior_targeted_source_qa_sha256')
    local=previous(native,'prior_local_report','prior_local_report_sha256')
    seam=previous(local,'prior_seam_report','prior_seam_report_sha256')
    mask=previous(seam,'prior_masking_qa','prior_masking_qa_sha256')
    q.require(native['summary']['all_sites_numerically_consistent'] and local['summary']['central_metrics_reproduced'] and seam['summary']['seam_audit_completed'] and mask['summary']['masking_qa_passed'],'Prior QA incomplete')
    previous(mask,'prior_extraction_qa','prior_extraction_qa_sha256')
    lookup={(x,y):group for x,y,group in e.legacy.candidate_tiles(coverage)}
    mr={tuple(r['tile']):r for r in mask['tiles']}; tr={tuple(r['tile']):r for r in terrain['tiles']}
    hashes={tuple(r['tile']):r for r in terrain['verified_hashes']}
    complete={k for k,r in mr.items() if r['state']=='COMPLETE'}
    q.require(len(lookup)==656 and set(lookup)==set(mr)==set(tr) and len(complete)==655 and set(hashes)==complete,'Tile inventory mismatch')
    olddir=diag/'hrdem1m_masked_staging'
    stems={e.stem(*k) for k in complete}
    for folder in (olddir,maskdir):
        inv=q.inventory(folder,stems,masked=True)
        q.require(not inv['missing_entries'] and not inv['unexpected_entries'],f'Inventory mismatch: {folder}')
    def paths(k):
        name=e.stem(*k)+'_masked.tif'
        return olddir/name,maskdir/name
    # Verify every file against the terrain audit before trusting its per-tile counts.
    for i,k in enumerate(sorted(complete),1):
        op,fp=paths(k); h=hashes[k]
        q.require(e.sha(op)==h['original_current_sha256'] and e.sha(fp)==h['fixed_sha256'],f'Raster hash changed: {k}')
        if i%100==0 or i==len(complete): print(f'HASH_PROGRESS: {i}/{len(complete)}',flush=True)
    @lru_cache(maxsize=24)
    def read(k):
        if k not in complete: return None
        op,fp=paths(k); a,_=q.read_grid(op,*k); b,_=q.read_grid(fp,*k)
        q.require(np.array_equal(a!=q.NODATA,b!=q.NODATA) and int((b!=q.NODATA).sum())==mr[k]['masked_valid_pixels'],f'Validity mismatch: {k}')
        return a,b
    ss=terrain['summary']; ms=mask['summary']
    expected={'slope_gt10':ss['slope_deg']['counts_over']['10'], 'slope_gt20':ss['slope_deg']['counts_over']['20'], 'aspect_gt90':ss['aspect_deg']['counts_over']['90'], 'elevation_gt10':ms['counts_over_m']['10.0']}
    q.require(expected==dict(slope_gt10=7250,slope_gt20=1154,aspect_gt90=233,elevation_gt10=45),'Unexpected previous extreme counts')
    # Only compute derivatives for tiles known to contain one of the selected indicators.
    selected=[k for k in sorted(lookup) if tr[k]['slope_deg']['counts_over']['10'] or tr[k]['aspect_deg']['counts_over']['90'] or mr[k]['counts_over_m']['10.0']]
    folder=Path(tempfile.mkdtemp(prefix='hrdem1m_anomaly_register_',dir=diag))
    temp=folder/'anomaly_register.partial.gpkg'; final=folder/'anomaly_register.gpkg'
    driver=ogr.GetDriverByName('GPKG'); q.require(driver is not None,'GPKG driver unavailable')
    ds=driver.CreateDataSource(str(temp)); q.require(ds is not None,'Cannot create GPKG')
    crs=osr.SpatialReference(); crs.ImportFromEPSG(26910)
    crs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    I,R,S=ogr.OFTInteger,ogr.OFTReal,ogr.OFTString
    layers={}
    def layer(name,geom,fields):
        ly=ds.CreateLayer(name,crs,geom); q.require(ly is not None,f'Layer creation failed: {name}')
        for n,typ in fields:
            q.require(ly.CreateField(ogr.FieldDefn(n,typ))==0,f'Field creation failed: {name}.{n}')
        layers[name]=ly
    xyfields=[('tile_x',I),('tile_y',I),('source_group',S)]
    points=xyfields+[('original',R),('fixed',R),('signed_delta',R),('abs_delta',R)]
    layer('slope_change_gt10deg',ogr.wkbPoint,points+[('gt20deg',I)])
    layer('aspect_change_gt90deg',ogr.wkbPoint,points)
    layer('elevation_change_gt10m',ogr.wkbPoint,points)
    seamfields=xyfields+[('direction',S),('neighbor_x',I),('neighbor_y',I),('neighbor_group',S),('crossings',I),('old_p95_m',R),('fixed_p95_m',R),('fixed_mean_m',R),('fixed_max_m',R)]
    layer('flagged_seams',ogr.wkbLineString,seamfields+[('historical',I),('local_class',S),('local_ratio',R),('joint_fraction',R),('slope_p95_deg',R),('slope_max_deg',R)])
    layer('insufficient_seams',ogr.wkbLineString,seamfields)
    layer('tiles_with_internal_nodata',ogr.wkbPolygon,xyfields+[('state',S),('footprint_px',I),('missing_px',I),('missing_pct',R)])
    counts=Counter()
    def write(name,geom,fields):
        ly=layers[name]; f=ogr.Feature(ly.GetLayerDefn())
        for k,v in fields.items():
            if v is not None: f.SetField(k,v)
        f.SetGeometry(geom); q.require(ly.CreateFeature(f)==0,f'Feature write failed: {name}')
        counts[name]+=1
    def base(k): return dict(tile_x=k[0],tile_y=k[1],source_group=lookup[k])
    def point(name,k,row,col,a,b,d,extra=None):
        geom=ogr.Geometry(ogr.wkbPoint); geom.AddPoint_2D(*point_xy(*k,row,col))
        fields=dict(**base(k),original=float(a),fixed=float(b),signed_delta=float(d),abs_delta=abs(float(d)))
        if extra: fields.update(extra)
        write(name,geom,fields)
    localrows={key(r):r for r in local['sites'] if r['kind']=='flagged'}
    terrainrows={key(r):r for r in terrain['seam_sites'] if r['kind']=='flagged'}
    flagkeys={key(r) for r in seam['seams'] if r['fixed']['flagged']}
    q.require(len(flagkeys)==41 and flagkeys==set(localrows)==set(terrainrows),'Flagged seam inventory mismatch')
    classes=Counter()
    for r in seam['seams']:
        fixed=r['fixed']; old=r['old']; k=tuple(r['tile']); direction=r['direction']
        if fixed['assessable'] and not fixed['flagged']: continue
        geom=ogr.Geometry(ogr.wkbLineString)
        for xy in seam_endpoints(*k,direction): geom.AddPoint_2D(*xy)
        fields=dict(**base(k),direction=direction,neighbor_x=r['neighbor'][0],neighbor_y=r['neighbor'][1],neighbor_group=r['source_groups'][1],crossings=r['common_crossings'],old_p95_m=old['p95_abs_m'],fixed_p95_m=fixed['p95_abs_m'],fixed_mean_m=fixed['signed_mean_m'],fixed_max_m=fixed['max_abs_m'])
        if fixed['flagged']:
            lc=localrows[key(r)]; tc=terrainrows[key(r)]['center_band']['slope_deg']; cl=localization(lc); classes[cl]+=1
            fields.update(historical=int(r['historical_flag']),local_class=cl,local_ratio=lc['joint_center_to_max_offset_p95_ratio'],joint_fraction=lc['joint_fraction_of_center'],slope_p95_deg=tc['p95_abs'],slope_max_deg=tc['max_abs'])
            write('flagged_seams',geom,fields)
        else: write('insufficient_seams',geom,fields)
    for k,r in mr.items():
        if not r['nodata_inside_footprint']: continue
        x,y=k; ring=ogr.Geometry(ogr.wkbLinearRing)
        for xy in [(x,y),(x+N,y),(x+N,y+N),(x,y+N),(x,y)]: ring.AddPoint_2D(*xy)
        geom=ogr.Geometry(ogr.wkbPolygon); geom.AddGeometry(ring)
        write('tiles_with_internal_nodata',geom,dict(**base(k),state=r['state'],footprint_px=r['footprint_pixels'],missing_px=r['nodata_inside_footprint'],missing_pct=100*r['nodata_inside_footprint']/r['footprint_pixels']))
    maxima={}; total=Counter()
    for i,k in enumerate(selected,1):
        x,y=k; a,b=t.patch(read,(x-1,y-1,x+N+1,y+N+1))
        va,sa,aa,_=t.derivatives(a); vb,sb,ab,_=t.derivatives(b)
        q.require(np.array_equal(va,vb) and int(va.sum())==tr[k]['strict_3x3_centers'],f'Derivative support mismatch: {k}')
        ca=a[1:-1,1:-1]; cb=b[1:-1,1:-1]
        v=(ca!=q.NODATA)&(cb!=q.NODATA)
        sd=sb-sa; ad=(ab-aa+180)%360-180; ed=cb.astype(np.float64)-ca.astype(np.float64)
        asp=va&(sa>=t.ASPECT_MIN_SLOPE)&(sb>=t.ASPECT_MIN_SLOPE)
        tests=[('slope_change_gt10deg',va&(np.abs(sd)>10),sa,sb,sd,'slope_gt10'),('aspect_change_gt90deg',asp&(np.abs(ad)>90),aa,ab,ad,'aspect_gt90'),('elevation_change_gt10m',v&(np.abs(ed)>10),ca,cb,ed,'elevation_gt10')]
        tilecounts={'slope_gt20':int((va&(np.abs(sd)>20)).sum())}
        for name,sel,old,fixed,delta,tag in tests:
            rr,cc=np.nonzero(sel); tilecounts[tag]=len(rr)
            for row,col in zip(rr,cc):
                d=float(delta[row,col]); extra={'gt20deg':int(abs(d)>20)} if tag=='slope_gt10' else None
                point(name,k,row,col,old[row,col],fixed[row,col],d,extra)
                if tag not in maxima or abs(d)>maxima[tag]['abs_delta']:
                    maxima[tag]=dict(tile=list(k),center_epsg26910=list(point_xy(*k,row,col)),original=float(old[row,col]),fixed=float(fixed[row,col]),signed_delta=d,abs_delta=abs(d))
        expected_tile=dict(slope_gt10=tr[k]['slope_deg']['counts_over']['10'],slope_gt20=tr[k]['slope_deg']['counts_over']['20'],aspect_gt90=tr[k]['aspect_deg']['counts_over']['90'],elevation_gt10=mr[k]['counts_over_m']['10.0'])
        q.require(tilecounts==expected_tile,f'Per-tile extreme count mismatch: {k}: {tilecounts} != {expected_tile}')
        total.update(tilecounts)
        if i%10==0 or i==len(selected): print(f'EXTREME_TILE_PROGRESS: {i}/{len(selected)}',flush=True)
    q.require(dict(total)==expected,'Global extreme count mismatch')
    q.require(dict(classes)==dict(strong=17,moderate=8,no_dominant_peak=7,insufficient_support=9),'Localization count mismatch')
    expected_layers=dict(slope_change_gt10deg=7250,aspect_change_gt90deg=233,elevation_change_gt10m=45,flagged_seams=41,insufficient_seams=52,tiles_with_internal_nodata=265)
    q.require(dict(counts)==expected_layers,'Feature count mismatch')
    for tag,reference in [('slope_gt10',ss['slope_deg']['max_abs']),('aspect_gt90',ss['aspect_deg']['max_abs']),('elevation_gt10',ms['maximum_abs_difference_m'])]:
        q.require(abs(maxima[tag]['abs_delta']-reference)<1e-9,f'Maximum mismatch: {tag}')
    # Close all layer handles before reopening and renaming the database (Windows).
    layers.clear(); ds=None
    check=ogr.Open(str(temp),0); q.require(check is not None,'Cannot reopen GPKG')
    q.require(check.GetLayerCount()==len(expected_layers),'Layer inventory mismatch')
    for name,n in expected_layers.items():
        ly=check.GetLayerByName(name)
        q.require(ly is not None and ly.GetFeatureCount()==n and ly.GetSpatialRef().IsSame(crs),'Reopened GPKG layer mismatch: '+name)
        ly=None
    check=None
    q.require(not any(Path(str(temp)+suffix).exists() for suffix in ('-wal','-journal','-shm')),'Database sidecar remains')
    os.replace(temp,final)
    summary=dict(raster_pairs_hash_verified=len(complete),selected_extreme_tiles=len(selected),extreme_counts=dict(total),layer_counts=dict(counts),localization_classes=dict(classes),maximum_locations=maxima,register_sha256=e.sha(final),register_completed=True,final_dem_approved=False)
    report=folder/'anomaly_register_qa.json'
    e.save_json(report,dict(configuration=cfg,report_chain_sha256=chain,script_sha256=e.sha(Path(__file__)),summary=summary,decision=dict(preferred_resampling='fixed_scale_bilinear',retain_elevations_unchanged=True,exclude_flagged_pixels_automatically=False,source_artifact_origin='unresolved',physical_model_validation='pending'),limitations=['Point layers indicate sensitivity between processing methods, not measured error. Thresholds are descriptive and omit smaller differences.','Seam lines depict full nominal tile edges, not confirmed artifact lengths or exact crossing locations.','NoData polygons depict full tiles containing internal NoData, not the actual missing-pixel geometry.','Layers may overlap; feature counts cannot be added as unique affected pixels.','No raster edited; no exclusion, smoothing, filling or final DEM approval.']))
    print('===== ANOMALY REGISTER SUMMARY ====='); print(json.dumps(summary,indent=2))
    print('REGISTER_GPKG:',final); print('DIAGNOSTIC_JSON:',report)
    print('ANOMALY_REGISTER_COMPLETED: True'); print('FINAL_DEM: FALSE')

if __name__=='__main__':
    main()
