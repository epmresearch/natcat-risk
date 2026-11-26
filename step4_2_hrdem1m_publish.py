"""Publish the validated HRDEM research input as a versioned, immutable copy.
Default is a read-only plan; --publish explicitly enables copying. No GDAL writes.
Run serially from the project root. Existing releases are never replaced.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import tempfile
from datetime import datetime, timezone

MOSAIC_SHA = 'eedd95b99c4c0a56067a1c15dcce413cd8ef64b56d68b7145ec159a94643ffb6'
REGISTER_SHA = 'e2e8840013db83074fc614d3a12d299ee6d199e09465172bf59ffb9e489c76ae'
NEGATIVE_REPORT = 'results/diagnostics/hrdem1m_negative_source_qa__cuibs6i/negative_source_qa.json'
FINAL_REPORT = 'results/diagnostics/hrdem1m_final_validation_ebk_g0f9/final_validation_qa.json'
MOSAIC_REPORT = 'results/diagnostics/hrdem1m_mosaic_candidate_hjcps_iq/mosaic_candidate_qa.json'
REGISTER_REPORT = 'results/diagnostics/hrdem1m_anomaly_register_q_06k0gb/anomaly_register_qa.json'
RELEASE = 'hrdem_1m_dtm_highway99_fixed_v1'
DEM = RELEASE + '.tif'
REGISTER = 'hrdem_1m_dtm_highway99_anomaly_register.gpkg'
PRODUCT_URL = 'https://open.canada.ca/data/en/dataset/0fe65119-e96e-4a57-8bfe-9d9245fba06b'
LAYERS = {'slope_change_gt10deg':7250,'aspect_change_gt90deg':233,'elevation_change_gt10m':45,'flagged_seams':41,'insufficient_seams':52,'tiles_with_internal_nodata':265}
REPORT_SCRIPTS = {
    'global_extraction_qa.json':'step4_2_hrdem1m_global_extraction_qa.py',
    'global_masking_qa.json':'step4_2_hrdem1m_global_masking_qa.py',
    'global_fixed_seam_qa.json':'step4_2_hrdem1m_global_fixed_seam_qa.py',
    'seam_local_context_qa.json':'step4_2_hrdem1m_seam_local_context_qa.py',
    'targeted_source_qa.json':'step4_2_hrdem1m_targeted_source_qa.py',
    'terrain_sensitivity_qa.json':'step4_2_hrdem1m_terrain_sensitivity_qa.py',
    'anomaly_register_qa.json':'step4_2_hrdem1m_anomaly_register.py',
    'mosaic_candidate_qa.json':'step4_2_hrdem1m_build_mosaic_candidate.py',
    'final_validation_qa.json':'step4_2_hrdem1m_final_validation.py',
    'negative_source_qa.json':'step4_2_hrdem1m_negative_source_qa.py',
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def local(root, rel):
    p = (root / rel).resolve()
    require(p.is_relative_to(root) and p.is_file(), f'Missing or external input: {p}')
    return p


def check_no_sidecars(path, suffixes):
    for suffix in suffixes:
        require(not Path(str(path)+suffix).exists(), f'Unexpected sidecar: {path}{suffix}')


def copy_verified(source, dest, expected):
    require(not dest.exists(), f'Destination exists: {dest}')
    dest.parent.mkdir(parents=True, exist_ok=True)
    require(sha(source)==expected, f'Source changed before copying: {source}')
    with source.open('rb') as inp, dest.open('xb') as out:
        shutil.copyfileobj(inp, out, length=1024*1024)
    require(sha(dest)==expected and sha(source)==expected, f'Copy hash mismatch: {dest}')


def write_new(path, content):
    with path.open('x', encoding='utf-8', newline='\n') as f:
        f.write(content)


def validate(root, cfg):
    negp=local(root,NEGATIVE_REPORT); fp=local(root,FINAL_REPORT)
    cp=local(root,MOSAIC_REPORT); rp=local(root,REGISTER_REPORT)
    neg,final,candidate,reg=map(read,(negp,fp,cp,rp))
    ns,fs,cs,rs=[r['summary'] for r in (neg,final,candidate,reg)]
    for r in (neg,final,candidate,reg):
        require(r['configuration']==cfg,'Configuration mismatch')
    require(local(root,neg['prior_final_report'])==fp and sha(fp)==neg['prior_final_report_sha256'],'Final report link mismatch')
    require(local(root,final['candidate_report'])==cp and sha(cp)==final['candidate_report_sha256'],'Mosaic report link mismatch')
    require(local(root,candidate['register_report'])==rp and sha(rp)==candidate['register_report_sha256'],'Register report link mismatch')
    require(ns['all_sites_numerically_consistent'] and ns['sites']==ns['sites_numerically_consistent']==3 and ns['deep_negative_pixels']==ns['deep_negative_strict']==2432 and ns['deep_negative_without_strict']==ns['deep_negative_over_0_001m']==0,'Negative QA not accepted')
    require(fs['final_technical_validation_passed'] and fs['all_pixels_scanned'] and fs['valid_pixels']==464090223 and fs['footprint_pixels']==464219963 and fs['nodata_inside_footprint']==129740 and fs['valid_outside_footprint']==0 and fs['grid_pixels']==2883000000,'Final validation not accepted')
    require(fs['negative_pixels']==fs['negative_pixels_matching_local_fixed_extraction_and_mask']==20054860 and fs['negative_counts_below_m']=={'-1':18408762,'-5':2432,'-10':1685,'-20':484},'Negative inventory mismatch')
    require(cs['mosaic_candidate_qa_passed'] and cs['all_grid_pixels_readback_verified'] and cs['input_hashes_verified_before_and_after'] and cs['input_tile_bitwise_mismatch_pixels']==cs['unexpected_valid_in_gap_cells']==0,'Mosaic QA not accepted')
    require(cs['crs']=='EPSG:26910' and cs['dtype']=='Float32' and cs['nodata']==-32767 and cs['geotransform']==[475000,1,0,5557000,0,-1] and (cs['width'],cs['height'])==(31000,93000),'Mosaic grid mismatch')
    require(rs['register_completed'] and rs['layer_counts']==fs['register_layer_counts']==LAYERS,'Register inventory mismatch')
    require(reg['decision']['retain_elevations_unchanged'] and not reg['decision']['exclude_flagged_pixels_automatically'],'Anomaly handling mismatch')
    require(not fs['absolute_vertical_accuracy_validated'] and not ns['negative_source_physical_validity_resolved'],'Unexpected scientific approval status')
    dem=local(root,final['mosaic']); gpkg=local(root,final['register'])
    require(dem==local(root,candidate['mosaic']) and sha(dem)==fs['mosaic_sha256']==cs['mosaic_sha256']==MOSAIC_SHA and dem.stat().st_size==cs['raster_bytes'],'Mosaic identity mismatch')
    require(gpkg==rp.parent/'anomaly_register.gpkg' and sha(gpkg)==rs['register_sha256']==final['register_sha256']==candidate['register_sha256']==REGISTER_SHA,'Register identity mismatch')
    check_no_sidecars(dem,('.aux.xml','.msk','.ovr'))
    check_no_sidecars(gpkg,('-wal','-journal','-shm'))
    chain=candidate['report_chain_sha256']; require(chain==reg['report_chain_sha256'],'QA chain mismatch')
    reports={local(root,p):h for p,h in chain.items()}
    reports.update({p:sha(p) for p in (negp,fp,cp,rp)})
    require(len(reports)==10 and {p.name for p in reports}==set(REPORT_SCRIPTS),'QA report inventory mismatch')
    for p,h in reports.items():
        require(sha(p)==h,f'QA report hash mismatch: {p}')
        r=read(p)
        require(r['configuration']==cfg,f'QA configuration mismatch: {p}')
        # The original extraction audit did not record its own script hash.
        historical_script_hash=r.get('script_sha256',r.get('audit_script_sha256'))
        if p.name=='global_extraction_qa.json':
            require(r['summary']['masking_review_gate_passed'],'Extraction QA not accepted')
        else:
            require(historical_script_hash is not None and sha(local(root,REPORT_SCRIPTS[p.name]))==historical_script_hash,f'QA script changed: {p.name}')
    return dem,gpkg,reports,final,candidate,neg


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--publish',action='store_true',help='Create the versioned research release')
    args=parser.parse_args()
    import step4_2_hrdem1m_fixed_batch_300 as e
    root=e.legacy.project_root().resolve(); cfg=e.config(root)
    coverage,vrt,_,_=e.paths(root)
    dest=coverage.parent/RELEASE
    require(not dest.exists() and not dest.is_symlink(),f'Release already exists; no overwrite: {dest}')
    require(not list(dest.parent.glob('.'+RELEASE+'.partial_*')),'Previous partial release exists; inspect it first')
    dem,gpkg,reports,final,candidate,neg=validate(root,cfg)
    entries={DEM:(dem,MOSAIC_SHA),REGISTER:(gpkg,REGISTER_SHA)}
    for path,digest in reports.items():
        entries['qa/'+path.name]=(path,digest)
    entries['provenance/source_vrt_snapshot.vrt']=(vrt,cfg['vrt_sha256'])
    payload=sum(p.stat().st_size for p,_ in entries.values())
    require(shutil.disk_usage(dest.parent).free>payload+100*1024**2,'Insufficient free disk space')
    print('PUBLICATION_PLAN:',json.dumps({'destination':str(dest),'files_to_copy':len(entries),'payload_bytes':payload,'publish_requested':args.publish}),flush=True)
    if not args.publish:
        print('FINAL_DEM_PUBLISHED: FALSE'); return
    stage=Path(tempfile.mkdtemp(prefix='.'+RELEASE+'.partial_',dir=dest.parent))
    inventory={}
    for name,(source,digest) in entries.items():
        copy_verified(source,stage/name,digest)
        inventory[name]={'sha256':digest,'bytes':(stage/name).stat().st_size,'original_project_path':source.relative_to(root).as_posix()}
        print('COPY_VERIFIED:',name,flush=True)
    cs=candidate['summary']; fs=final['summary']
    manifest={
      'release_id':RELEASE,'created_utc':datetime.now(timezone.utc).isoformat(),
      'status':'technically_validated_research_input','configuration':cfg,
      'source':{'publisher':'Natural Resources Canada','product':'High Resolution Digital Elevation Model Mosaic (HRDEM Mosaic) - CanElevation Series, DTM','product_documentation':PRODUCT_URL,'licence':'Open Government Licence - Canada','native_metadata':neg['source_metadata'], 'remote_assets_version_pinned':False},
      'horizontal_crs':'EPSG:26910','native_horizontal_crs':'EPSG:3979','vertical_datum':'CGVD2013 (official product documentation; not encoded in audited source WKT)','height_type':'orthometric','elevation_unit':'m (product convention; source band unit empty)','vertical_transformation_applied':False,
      'raster':{k:cs[k] for k in ('bounds_epsg26910','width','height','geotransform','dtype','nodata','total_grid_pixels','valid_pixels','nodata_pixels','minimum_elevation_m','maximum_elevation_m')},
      'footprint':{'project_path':coverage.relative_to(root).as_posix(),'sha256':cfg['footprint_sha256'],'rasterized_pixels':fs['footprint_pixels'],'missing_pixels':fs['nodata_inside_footprint'],'coverage_percent':fs['footprint_coverage_percent']},
      'negative_elevations':{'pixels':fs['negative_pixels'],'counts_below_m':fs['negative_counts_below_m'],'deep_pixels_native_reference_verified':2432,'deep_review_sites':[r for r in neg['sites'] if r['kind']=='deep']},
      'anomaly_register_layer_counts':LAYERS,'file_inventory':inventory,
      'publication_script_sha256':sha(Path(__file__)),
      'qa_script_hashes':{name:{'path':script,'current_sha256':sha(local(root,script)), 'historical_hash_available':name!='global_extraction_qa.json'} for name,script in REPORT_SCRIPTS.items()},
      'limits':['Technical validation and numerical source fidelity are not independent vertical accuracy validation.','The physical origin of source anomalies and negative elevations remains unresolved.','Spatial register is diagnostic, not an exclusion mask or measured-error map; seam lines and NoData tile polygons are nominal geometries.','Physical landslide model and risk sensitivity remain unvalidated.','Coverage is relative to the rasterized HRDEM footprint, not necessarily the complete AOI.','Native-source audit is targeted; remote source assets were not version-pinned.','Historical QA JSON paths refer to the original project paths; copied evidence is indexed here.','The global extraction QA report did not record its own script hash; its report hash is verified through the downstream chain.'],
      'elevation_values_modified_during_publication':False,'absolute_vertical_accuracy_validated':False,'negative_source_physical_validity_resolved':False,'hazard_model_validated':False,
      'publication_qa_passed':True,
    }
    readme=f'''# HRDEM 1 m DTM — Highway 99 — Fixed Scale v1

Technically validated research input. Not an independently validated terrain survey or landslide hazard model.

- DEM: `{DEM}`
- Spatial diagnostic register: `{REGISTER}`
- CRS: NAD83 / UTM zone 10N (EPSG:26910); 1 m grid; Float32; NoData -32767.
- Vertical reference: CGVD2013 orthometric heights, according to official product documentation. The audited source WKT is horizontal-only and the source band unit field is empty. No vertical transformation was applied.
- Source: Natural Resources Canada, HRDEM Mosaic, CanElevation Series, DTM. Open Government Licence - Canada. {PRODUCT_URL}
- Fixed-scale point bilinear resampling (XSCALE=YSCALE=1, exact transformer, no source overviews), pixel-center footprint masking and exact pixel mosaicking.
- Valid pixels: 464,090,223. Rasterized footprint pixels: 464,219,963. Missing inside footprint: 129,740 (99.972052% coverage). This is not whole-AOI coverage.
- All elevation values, including negative values, retained. No smoothing, gap filling, seam replacement or automatic anomaly exclusion.
- All 2,432 pixels below -5 m agree with the strict native bilinear reference within 0.001 m. This does not establish physical validity.
- Diagnostic register: 41 flagged seams, 52 insufficient-data seams, 265 tiles with internal NoData, and extreme differences between processing methods. Do not treat these features as a terrain-error or exclusion mask.
- Absolute vertical accuracy, source artifact origin and physical landslide model sensitivity remain unresolved.

`publication_manifest.json` records hashes and lineage. `qa/` contains the ten unchanged QA reports. The VRT snapshot is provenance evidence, not the final raster or a pinned copy of the remote sources. Original diagnostics and staging files remain untouched.
'''
    write_new(stage/'README.md',readme)
    manifest['readme_sha256']=sha(stage/'README.md')
    write_new(stage/'publication_manifest.json',json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    # Recheck all copied payloads and inputs before the directory is exposed.
    for name,(source,digest) in entries.items():
        require(sha(stage/name)==sha(source)==digest,f'Final hash mismatch: {name}')
    require(sha(coverage)==cfg['footprint_sha256'],'Footprint changed during publication')
    check_no_sidecars(dem,('.aux.xml','.msk','.ovr')); check_no_sidecars(gpkg,('-wal','-journal','-shm'))
    require(not dest.exists() and not dest.is_symlink(),'Destination appeared during publication')
    stage.rename(dest)  # Same-filesystem rename; run serially, Windows refuses existing destination.
    for name,(_,digest) in entries.items():
        require(sha(dest/name)==digest,f'Published hash mismatch: {name}')
    require(read(dest/'publication_manifest.json')==manifest and sha(dest/'README.md')==manifest['readme_sha256'],'Published metadata mismatch')
    print('PUBLISHED_DEM:',dest/DEM)
    print('PUBLISHED_REGISTER:',dest/REGISTER)
    print('PUBLICATION_MANIFEST:',dest/'publication_manifest.json')
    print('PUBLISHED_MOSAIC_SHA256:',sha(dest/DEM))
    print('PUBLISHED_REGISTER_SHA256:',sha(dest/REGISTER))
    print('PUBLICATION_QA_PASSED: True')
    print('FINAL_DEM_PUBLISHED: TRUE')
    print('ABSOLUTE_VERTICAL_ACCURACY_VALIDATED: FALSE')
    print('HAZARD_MODEL_VALIDATED: FALSE')


if __name__=='__main__':
    main()
