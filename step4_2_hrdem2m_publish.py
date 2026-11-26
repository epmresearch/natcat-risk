"""Publish the validated HRDEM 2 m research input as an immutable v1 package.

Without --publish this script verifies inputs and prints a plan only.
With --publish it copies into a temporary sibling directory, verifies all
copies by SHA-256, and renames the directory only after the checks succeed.
No existing publication or staging raster is overwritten.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import tempfile

import step4_2_hrdem2m_fixed_batch_300 as e
import step4_2_hrdem2m_fixed_mask_batch_300 as m
import step4_2_hrdem2m_anomaly_register as ar
import step4_2_hrdem2m_build_mosaic_candidate as b
import step4_2_hrdem2m_final_validation as v

METHOD = "HRDEM2M_RESEARCH_INPUT_PUBLICATION_V1"
FINAL_QA = "hrdem2m_final_validation_ij5c2ve2/final_validation_qa.json"
EXTRACTION_QA = "hrdem2m_global_extraction_qa_uvxf80ty/global_extraction_qa.json"
DEST_NAME = "hrdem_2m_dtm_highway99_fixed_v1"
DEM_NAME = DEST_NAME + ".tif"
REGISTER_NAME = "hrdem_2m_dtm_highway99_anomaly_register.gpkg"
ONE_M = "data/processed/elevation/hrdem_1m_dtm_highway99_fixed_v1"
PROTECTED_1M = {
    ONE_M+"/hrdem_1m_dtm_highway99_fixed_v1.tif":
        "eedd95b99c4c0a56067a1c15dcce413cd8ef64b56d68b7145ec159a94643ffb6",
    ONE_M+"/hrdem_1m_dtm_highway99_anomaly_register.gpkg":
        "e2e8840013db83074fc614d3a12d299ee6d199e09465172bf59ffb9e489c76ae",
}
require = e.require


def relative(p):
    return p.relative_to(e.ROOT).as_posix()


def verify_hashes(hashes):
    for name,sha in hashes.items():
        p=e.ROOT/name
        require(p.is_file() and e.digest(p)==sha,f"Input hash mismatch: {name}")


def build_plan():
    e.gdal.UseExceptions()
    tiles,cfg,ch,mc=m.config()
    mh=e.canonical_hash(mc)
    require(ch==v.EXPECTED_EXTRACTION_SHA256 and mh==v.EXPECTED_MASK_SHA256,
            "Extraction/mask configuration mismatch")
    final_path=e.DIAG/FINAL_QA
    final=e.load(final_path)
    require(final["final_technical_validation_passed"] is True
            and final["all_pixels_scanned"] is True
            and final["input_hashes_verified_before_and_after"] is True
            and final["final_dem_published"] is False
            and final["extraction_configuration_sha256"]==ch
            and final["mask_configuration_sha256"]==mh
            and final["mosaic_sha256"]==v.EXPECTED_MOSAIC_SHA256
            and final["register_sha256"]==v.EXPECTED_REGISTER_SHA256,
            "Final validation gate not passed")
    require(final["grid_pixels"]==728500000 and final["candidate_cells"]==671
            and final["valid_pixels"]==116196821 and final["nodata_inside_footprint"]==34697
            and final["footprint_pixels"]==116231518 and final["valid_outside_footprint"]==0
            and final["deep_negative_points_matching_mosaic"]==606
            and final["register_layer_counts"]==ar.EXPECTED
            and final["absolute_vertical_accuracy_validated"] is False
            and final["negative_source_physical_validity_resolved"] is False
            and final["hazard_model_validated"] is False,"Unexpected final statistics or claims")
    inputs=dict(final["input_sha256"])
    inputs[relative(Path(v.__file__).resolve())]=final["script_sha256"]
    inputs[relative(final_path)]=e.digest(final_path)
    verify_hashes(inputs)
    candidate=e.ROOT/final["candidate_path"]
    register=e.ROOT/final["register_path"]
    candidate_qa=e.ROOT/final["candidate_qa_path"]
    register_qa=e.ROOT/final["register_qa_path"]
    reg=e.load(register_qa)
    cand=e.load(candidate_qa)
    require(cand["mosaic_candidate_qa_passed"] and reg["register_completed"]
            and cand["mosaic_sha256"]==final["mosaic_sha256"]
            and reg["register_sha256"]==final["register_sha256"],
            "Candidate or register gate not passed")
    verify_hashes(PROTECTED_1M)
    print("EXISTING_1M_HASHES_VERIFIED: True",flush=True)
    require(e.load(e.OUT/"configuration.json")==cfg and e.load(m.OUT/"configuration.json")==mc,
            "Stored configuration mismatch")
    ex_path=e.DIAG/EXTRACTION_QA
    ex=e.load(ex_path)
    mask=e.load(e.DIAG/b.MASK_QA)
    require(ex["summary"]["global_extraction_qa_passed"]
            and ex["summary"]["configuration_sha256"]==ch
            and ex["summary"]["valid_pixels"]==125764085
            and ex["summary"]["states"]==dict(b.EXPECTED_STATES, PENDING=0),
            "Extraction QA gate not passed")
    er={tuple(r["tile"]):r for r in ex["tiles"]}
    mr={tuple(r["tile"]):r for r in mask["tiles"]}
    require(set(er)==set(mr)==set(tiles) and len(tiles)==671
            and all(er[k]["state"]==mr[k]["state"]
                    and er[k]["valid_pixels"]==mr[k]["source_valid_pixels"] for k in er),
            "Extraction/mask inventory mismatch")
    # Verify the complete staging inventory before publishing; the final
    # package itself retains the compact aggregate QA reports, not the tiles.
    for index,k in enumerate(tiles,1):
        et,eq=e.paths(k);mt,mq=m.paths(k)
        require(e.digest(eq)==er[k]["qa_sha256"] and e.digest(mq)==mr[k]["mask_qa_sha256"],
                f"Staging QA mismatch: {k}")
        if er[k]["state"]=="COMPLETE":
            require(e.digest(et)==er[k]["raster_sha256"]
                    and e.digest(mt)==mr[k]["masked_raster_sha256"],f"Staging raster mismatch: {k}")
        else:
            require(not et.exists() and not mt.exists(),f"Unexpected empty-tile raster: {k}")
        if index%100==0 or index==len(tiles):
            print("STAGING_HASH_PROGRESS:",f"{index}/{len(tiles)}",flush=True)

    files={
        DEM_NAME:candidate,
        REGISTER_NAME:register,
        "qa/global_extraction_qa.json":ex_path,
        "qa/global_masking_qa.json":e.DIAG/b.MASK_QA,
        "qa/global_seam_qa.json":e.DIAG/ar.QA_PATHS["screening"],
        "qa/targeted_seam_bridge_qa.json":e.DIAG/ar.QA_PATHS["bridge"],
        "qa/targeted_native_source_qa.json":e.DIAG/ar.QA_PATHS["native"],
        "qa/negative_source_qa.json":e.DIAG/ar.QA_PATHS["negative"],
        "qa/anomaly_register_qa.json":register_qa,
        "qa/mosaic_candidate_qa.json":candidate_qa,
        "qa/final_validation_qa.json":final_path,
        "qa/scale_comparison_qa.json":e.TEST/"scale_comparison_qa.json",
        "qa/native_bilinear_reference_qa.json":e.TEST/"native_bilinear_reference_qa.json",
        "qa/pilot_extraction_qa.json":e.OUT/"pilot_extraction_qa.json",
        "qa/pilot_seam_edge_qa.json":m.SEAM_QA,
        "qa/pilot_mask_qa.json":m.OUT/"pilot_mask_qa.json",
        "provenance/source_vrt_snapshot.vrt":e.VRT,
        "provenance/footprint_snapshot.gpkg":e.FOOTPRINT,
        "provenance/production_plan_qa.json":e.PLAN,
        "provenance/warp_commands.json":e.TEST/"warp_commands.json",
        "provenance/extraction_configuration.json":e.OUT/"configuration.json",
        "provenance/mask_configuration.json":m.OUT/"configuration.json",
    }
    for name in e.SITES:
        files[f"provenance/pilot_reference/{name}_narrow_fixed.tif"]=e.TEST/(name+"_narrow_fixed.tif")
    expected={**inputs,**cfg["input_sha256"]}
    expected[relative(e.OUT/"pilot_extraction_qa.json")]=mc["extraction_pilot_qa_sha256"]
    expected[relative(m.SEAM_QA)]=mc["seam_edge_qa_sha256"]
    inventory={}
    for name,src in files.items():
        require(src.is_file(),f"Missing publication input: {src}")
        rel=relative(src)
        digest=e.digest(src)
        if rel in expected:
            require(digest==expected[rel],f"Publication input changed: {rel}")
        inventory[name]={"original_project_path":rel,
                         "sha256":digest,"bytes":src.stat().st_size}
    scripts={name:{"path":name,"sha256":sha} for name,sha in inputs.items() if name.endswith(".py")}
    scripts[Path(__file__).name]={"path":relative(Path(__file__).resolve()),
                                 "sha256":e.digest(Path(__file__).resolve())}
    manifest={
        "publication_method":METHOD,
        "version":"fixed_v1",
        "release_status":"technically_validated_research_input",
        "final_dem_published":True,
        "absolute_vertical_accuracy_validated":False,
        "negative_source_physical_validity_resolved":False,
        "hazard_model_validated":False,
        "integrated_multiresolution_dem":False,
        "source_products":[
            "https://canelevation-dem.s3.ca-central-1.amazonaws.com/hrdem-mosaic-2m/2_3-mosaic-2m-dtm.tif",
            "https://canelevation-dem.s3.ca-central-1.amazonaws.com/hrdem-mosaic-2m/2_4-mosaic-2m-dtm.tif"],
        "extraction_configuration_sha256":ch,
        "mask_configuration_sha256":mh,
        "configuration":cfg,
        "mask_configuration":mc,
        "footprint":{"project_path":relative(e.FOOTPRINT),"sha256":e.digest(e.FOOTPRINT)},
        "mosaic_sha256":final["mosaic_sha256"],
        "register_sha256":final["register_sha256"],
        "technical_validation":{k:final[k] for k in (
            "grid_pixels","valid_pixels","nodata_pixels","footprint_pixels",
            "nodata_inside_footprint","valid_outside_footprint","footprint_coverage_percent",
            "minimum","maximum","negative_pixels","negative_tiles","negative_counts_below_m",
            "deep_negative_points_matching_mosaic","register_layer_counts")},
        "grid":{"crs":"EPSG:26910","resolution_m":2,"dtype":"Float32","nodata":-32767,
                "width":cand["width"],"height":cand["height"],
                "bounds_epsg26910":cand["bounds_epsg26910"],"geotransform":cand["geotransform"]},
        "qa_script_hashes":scripts,
        "publication_script_sha256":e.digest(Path(__file__).resolve()),
        "protected_1m_sha256":PROTECTED_1M,
        "file_inventory":inventory,
    }
    return files,manifest


def readme(manifest):
    t=manifest["technical_validation"]
    return f"""# Highway 99 HRDEM 2 m — fixed v1

## Status

Technically validated research input for terrain processing in the natcat-risk
project. This is a separate 2 m dataset, not the integrated 1 m / 2 m / 30 m DEM
and not a validated hazard model. No published HRDEM 1 m files were changed.

## Processing

- Source: NRCan HRDEM 2 m DTM regional mosaics 2_3 and 2_4.
- Native source CRS: EPSG:3979. Output horizontal CRS: EPSG:26910.
- Grid: 2 m, 15,500 columns by 47,000 rows; Float32, NoData -32767.
- Fixed-scale bilinear interpolation: XSCALE=1, YSCALE=1, error threshold 0,
  source overviews disabled, aligned target grid.
- Pixel-center footprint mask; 671 candidate cells, 669 extracted TIFFs and
  two empty extraction markers.
- Local mosaic assembly used direct pixel placement, without further resampling.
- Full-grid bitwise readback and independent footprint validation passed.

## Validation

- Valid cells: {t['valid_pixels']:,}; footprint coverage:
  {t['footprint_coverage_percent']:.8f}%.
- NoData inside footprint: {t['nodata_inside_footprint']:,} cells.
- Valid outside footprint: {t['valid_outside_footprint']}.
- Elevation range: {t['minimum']['value_m']:.6f} to {t['maximum']['value_m']:.6f} m.
- Negative cells: {t['negative_pixels']:,}; 606 cells below -5 m.
- All 20 flagged seams reproduced numerically in continuous-source bridge
  tests and independent four-neighbor bilinear reference tests.
- The 606 deep negative cells were numerically consistent with the native
  source. They were not removed, filled or smoothed.
- The anomaly register contains 20 flagged seams, 92 insufficient-support
  boundaries, 273 tiles containing footprint NoData, 119 negative-elevation
  tiles, and 606 deep-negative pixel-center points.
- Seam lines and tile polygons represent reporting units, not the precise
  spatial extent of every affected pixel.

## Limitations

Absolute vertical accuracy, physical validity of negative source elevations,
and hazard model accuracy have NOT been validated. Numerical source fidelity
does not establish that source elevations are true ground elevations.
Coverage gaps remain NoData. The dataset is intended for research with
documented limitations and should not be treated as a hazard prediction.

## Files and reproducibility

The GeoTIFF and GeoPackage are byte-identical copies of the validated candidate
and register. The QA reports retain their original contents; their pre-publication
status flags describe the stage at which each report was generated.
The publication manifest records original relative paths, SHA-256 checksums,
configuration and the final publication status. Large staging tiles and
intermediate seam rasters remain outside this release. The provenance directory
contains the input VRT, footprint, processing configurations, production plan
and small pilot references. Original source URLs are recorded in the manifest;
remote data can change independently of this local release.
"""


def verify_package(directory, inventory, readme_sha):
    for name,r in inventory.items():
        path=directory/name
        require(path.is_file() and path.stat().st_size==r["bytes"]
                and e.digest(path)==r["sha256"],f"Published payload hash mismatch: {name}")
    require(e.digest(directory/"README.md")==readme_sha,"README hash mismatch")


def publish_package(destination,files,manifest,protect=None):
    require(not destination.exists(),f"Refusing to overwrite {destination}")
    parent=destination.parent
    require(parent.is_dir(),"Publication directory parent does not exist")
    size=sum(r["bytes"] for r in manifest["file_inventory"].values())
    require(shutil.disk_usage(parent).free>=size+64*1024*1024,"Insufficient free space")
    tmp=Path(tempfile.mkdtemp(prefix=destination.name+"_partial_",dir=parent))
    print("PUBLICATION_STAGING:",tmp,flush=True)
    for name,src in files.items():
        target=tmp/name
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(src,target)
        require(e.digest(target)==manifest["file_inventory"][name]["sha256"],
                f"Copy hash mismatch: {name}")
        print("COPY_VERIFIED:",name,flush=True)
    (tmp/"README.md").write_text(readme(manifest),encoding="utf-8")
    manifest["readme_sha256"]=e.digest(tmp/"README.md")
    manifest["published_at_utc"]=datetime.now(timezone.utc).isoformat()
    manifest["publication_qa_passed"]=True
    manifest["existing_1m_unchanged"]=True
    e.new_json(tmp/"publication_manifest.json",manifest)
    verify_package(tmp,manifest["file_inventory"],manifest["readme_sha256"])
    for name,src in files.items():
        require(e.digest(src)==manifest["file_inventory"][name]["sha256"],
                f"Source changed during publication: {name}")
    if protect:
        protect()
    require(not destination.exists(),"Destination appeared during publication")
    tmp.rename(destination)
    verify_package(destination,manifest["file_inventory"],manifest["readme_sha256"])
    require(e.load(destination/"publication_manifest.json")==manifest,"Manifest readback mismatch")
    return destination


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--publish",action="store_true",help="Create the validated v1 publication")
    args=parser.parse_args()
    destination=e.ROOT/"data/processed/elevation"/DEST_NAME
    require(not destination.exists(),f"Refusing to overwrite {destination}")
    files,manifest=build_plan()
    plan={"destination":str(destination),"files_to_copy":len(files),
          "payload_bytes":sum(x["bytes"] for x in manifest["file_inventory"].values()),
          "publish_requested":args.publish}
    print("PUBLICATION_PLAN:",json.dumps(plan),flush=True)
    if not args.publish:
        print("PLAN_ONLY: no publication files written")
        return
    publish_package(destination,files,manifest,lambda:verify_hashes(PROTECTED_1M))
    print("PUBLISHED_DEM:",destination/DEM_NAME)
    print("PUBLISHED_REGISTER:",destination/REGISTER_NAME)
    print("PUBLICATION_MANIFEST:",destination/"publication_manifest.json")
    print("PUBLISHED_MOSAIC_SHA256:",manifest["mosaic_sha256"])
    print("PUBLISHED_REGISTER_SHA256:",manifest["register_sha256"])
    print("EXISTING_1M_UNCHANGED: TRUE")
    print("PUBLICATION_QA_PASSED: True")
    print("FINAL_DEM_PUBLISHED: TRUE")
    print("ABSOLUTE_VERTICAL_ACCURACY_VALIDATED: FALSE")
    print("HAZARD_MODEL_VALIDATED: FALSE")


if __name__=="__main__":
    main()
