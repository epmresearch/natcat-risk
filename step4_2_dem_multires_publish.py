"""Publish the validated multiresolution reference DEM v1, with provenance.

Plan only unless --publish is passed. All payloads are byte-identical to
the validated inputs. The package is reference elevation, not a hazard-ready
terrain model or a uniformly resolved 1 m DEM. No existing product is replaced.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import json
import shutil
import tempfile

import step4_2_dem_multires_final_validation as f

c=f.c
b=f.b
METHOD="MULTIRES_REFERENCE_PUBLICATION_V1"
DEST_NAME="dem_highway99_multires_reference_v1"
DEM_NAME="dem_highway99_multires_reference_v1.tif"
SOURCE_NAME="dem_highway99_multires_source_provenance_v1.tif"
FINAL_DIR=c.DIAG/"dem_multires_final_validation_lpz6xn6d"
FINAL_QA=FINAL_DIR/"final_validation_qa.json"
REGISTER=FINAL_DIR/f.REGISTER_NAME
EXPECTED_FINAL_SCRIPT="1c245f57a1f4296c95068fbe213f80fc780437849555a6b7e83e9768c968bbd8"
EXPECTED_REGISTER="fb7e72e3c423208a18f95331500e878a053c380ca16e6a5cadcd7f6d561316f5"
SOURCE_MANIFESTS={
    "hrdem1m_publication_manifest":c.ELEV/"hrdem_1m_dtm_highway99_fixed_v1"/"publication_manifest.json",
    "hrdem2m_publication_manifest":c.ELEV/"hrdem_2m_dtm_highway99_fixed_v1"/"publication_manifest.json",
}


def build_plan():
    paths,hashes,candidate,transition=f.check_inputs()
    original_keys=set(paths)
    paths.update({
        "final_validation_qa":FINAL_QA,"final_validation_script":Path(f.__file__).resolve(),
        "transition_register":REGISTER,"publication_script":Path(__file__).resolve(),
        **SOURCE_MANIFESTS,
    })
    for k in set(paths)-original_keys:
        c.require(paths[k].is_file(),f"Required file missing: {k}: {paths[k]}")
        hashes[k]=c.digest(paths[k])
    q=json.loads(FINAL_QA.read_text(encoding="utf-8"))
    c.require(q["method_id"]==f.METHOD and q["final_technical_validation_passed"] is True
              and q["input_files_unchanged"] is True and q["all_grid_pixels_scanned"] is True
              and q["source_counts_and_statistics_reconciled"] is True
              and q["transition_statistics_reproduced_from_candidate"] is True
              and q["final_integrated_dem_published"] is False,"Final technical QA gate failed")
    c.require(q["script_sha256"]==hashes["final_validation_script"]==EXPECTED_FINAL_SCRIPT,
              "Final technical QA script hash mismatch")
    c.require(q["input_sha256"]=={k:hashes[k] for k in original_keys}
              and q["input_paths"]=={k:paths[k].relative_to(c.ROOT).as_posix() for k in original_keys},
              "Final technical QA lineage mismatch")
    c.require(q["source_transitions_accepted_for_terrain_derivatives"] is False
              and q["cross_source_elevation_agreement_accepted"] is False
              and q["absolute_vertical_accuracy_validated"] is False
              and q["hazard_model_validated"] is False,
              "Scientific-use restrictions must stay unchanged")
    c.require(tuple(q["bounds_epsg26910"])==b.BOUNDS and q["total_grid_pixels"]==3102000000
              and q["aoi_pixels"]==1043616657 and q["uncovered_aoi_pixels"]==0
              and q["valid_outside_aoi_pixels"]==0 and q["coverage_percent"]==100.0
              and q["total_transition_edges"]==772181,"Final coverage/transition inventory mismatch")
    counts={int(k):v for k,v in q["source_pixel_counts"].items()}
    stats={int(k):v for k,v in q["source_elevation_statistics"].items()}
    c.require(counts==f.EXPECTED_COUNTS,"Final source counts changed")
    f.compare_previous(counts,stats,q["transition_pairs"],candidate,transition)
    r=q["register"]
    c.require(r["sha256"]==hashes["transition_register"]==EXPECTED_REGISTER
              and r["bytes"]==REGISTER.stat().st_size==10416128
              and r["path"]==REGISTER.relative_to(c.ROOT).as_posix()
              and r["counts"]==f.EXPECTED_FLAGS and r["features"]==50875
              and r["threshold_abs_jump_m_strict_gt"]==10.0
              and r["geometry_and_attributes_verified"] is True
              and r["threshold_is_diagnostic_not_hazard_or_accuracy_criterion"] is True,
              "Transition register QA gate failed")
    mem,layer,envelope,area,target=c.projected_aoi(c.INPUTS["aoi"])
    c.require(f.verify_register(REGISTER,target,f.EXPECTED_FLAGS)==f.EXPECTED_FLAGS,
              "Register read-only review failed")
    layer=None;mem=None
    print("FINAL_VALIDATION_AND_REGISTER_GATE_PASSED: True",flush=True)

    # Preserve the complete diagnostic chain and key local provenance.
    # Larger published HRDEM rasters are referenced, not duplicated.
    layout={
        DEM_NAME:"reference_dem",SOURCE_NAME:"source_provenance",f.REGISTER_NAME:"transition_register",
        "qa/mrdem30_consistent_mask_qa.json":"mrdem30_qa",
        "qa/multires_coverage_qa.json":"coverage_qa",
        "qa/multires_overlap_qa.json":"overlap_qa",
        "qa/multires_transition_qa.json":"transition_qa",
        "qa/mosaic_candidate_qa.json":"candidate_qa",
        "qa/final_validation_qa.json":"final_validation_qa",
        "provenance/highway99_terrain_search_5000m.geojson":"aoi",
        "provenance/mrdem30_exact_bbox.tif":"mrdem30_uncut",
        "provenance/hrdem1m_publication_manifest.json":"hrdem1m_publication_manifest",
        "provenance/hrdem2m_publication_manifest.json":"hrdem2m_publication_manifest",
        "source_anomaly_registers/hrdem_1m_dtm_highway99_anomaly_register.gpkg":"hrdem1m_register",
        "source_anomaly_registers/hrdem_2m_dtm_highway99_anomaly_register.gpkg":"hrdem2m_register",
    }
    for key in ("coverage_script","overlap_script","transition_script","candidate_script",
                "final_validation_script","publication_script"):
        layout["provenance/scripts/"+paths[key].name]=key
    files={name:paths[key] for name,key in layout.items()}
    inventory={name:{"original_project_path":paths[key].relative_to(c.ROOT).as_posix(),
                     "sha256":hashes[key],"bytes":paths[key].stat().st_size}
               for name,key in layout.items()}
    manifest={
        "method_id":METHOD,"release_id":DEST_NAME,"publication_script_sha256":hashes["publication_script"],
        "product_role":"Technically validated multiresolution reference elevation with mandatory source provenance",
        "crs":"EPSG:26910","accounting_grid_spacing_m":1,
        "bounds_epsg26910":list(b.BOUNDS),"geotransform":[475000,1,0,5557000,0,-1],
        "width":33000,"height":94000,"total_grid_pixels":3102000000,
        "aoi_pixels":1043616657,"coverage_percent":100.0,
        "reference_dem_dtype":"Float32","reference_dem_nodata":c.NODATA,
        "source_provenance_dtype":"Byte","source_provenance_nodata":0,
        "source_codes":{str(k):v for k,v in b.CODES.items()},
        "source_pixel_counts":q["source_pixel_counts"],
        "source_elevation_statistics":q["source_elevation_statistics"],
        "source_assignment":"AOI pixel centre, valid HRDEM 1m > HRDEM 2m > MRDEM 30m uncut; containing native cell",
        "mosaic_sha256":hashes["reference_dem"],"source_provenance_sha256":hashes["source_provenance"],
        "register_sha256":hashes["transition_register"],"transition_edges":772181,
        "transition_register_features":50875,"transition_register_counts":r["counts"],
        "transition_register_threshold_m_strict_gt":10.0,
        "documented_vertical_reference":"CGVD2013",
        "source_vertical_reference_evidence_urls":b.o.DATUM_SOURCES,
        "local_height_realization_independently_verified":False,
        "final_technical_validation_passed":True,
        "cross_source_elevation_agreement_accepted":False,
        "source_transitions_accepted_for_terrain_derivatives":False,
        "terrain_derivatives_from_mosaic_approved":False,
        "absolute_vertical_accuracy_validated":False,
        "hazard_model_validated":False,
        "limitations":q["limitations"]+[
            "Candidate GeoTIFF metadata is retained byte-for-byte; publication_manifest.json defines the local v1 release status.",
            "Original HRDEM rasters remain in their separately published directories; their hashes are recorded here.",
        ],
        "input_paths":{k:p.relative_to(c.ROOT).as_posix() for k,p in paths.items()},
        "input_sha256":hashes,"file_inventory":inventory,
        "final_integrated_dem_published":False,"publication_qa_passed":False,
    }
    return paths,files,manifest


def readme(m):
    return """# Highway 99 multiresolution reference DEM — v1

## Intended use

Technically validated reference terrain elevation, with mandatory
source-provenance raster and source-transition/anomaly records.
This is not a uniformly resolved 1 m terrain model or a validated
landslide hazard product.

## Main products

- `dem_highway99_multires_reference_v1.tif`: Float32 elevation, NoData -32767.
- `dem_highway99_multires_source_provenance_v1.tif`: Byte source codes:
  0 outside AOI/NoData; 1 HRDEM 1 m; 2 HRDEM 2 m; 30 MRDEM 30 m.
- `dem_highway99_multires_transition_register.gpkg`: 50,875 cardinal
  source edges with absolute elevation jump strictly greater than 10 m.
  This is a descriptive diagnostic threshold, not an accuracy or hazard
  criterion. Other edges are not certified safe.
- `source_anomaly_registers/`: original HRDEM source anomaly registers.
- `qa/`: unchanged diagnostic reports; historical status flags refer to
  their original processing stages.
- `provenance/`: AOI, uncut MRDEM input, source publication manifests and
  processing scripts. The published HRDEM input rasters are referenced by
  path and SHA-256 in the release manifest rather than duplicated.

## Grid and coverage

NAD83 / UTM zone 10N (EPSG:26910), 33000 x 94000 cells, 1 m accounting
grid spacing, bounds (475000, 5463000, 508000, 5557000).

AOI pixel-centres: 1,043,616,657; covered 100%; no valid pixels outside AOI.
Selected source counts:

| Source | Accounting pixels | Area (km2) | Percent of AOI |
|---|---:|---:|---:|
| HRDEM 1 m | 464,090,223 | 464.090223 | 44.469415 |
| HRDEM 2 m | 698,372 | 0.698372 | 0.066918 |
| MRDEM 30 m | 578,828,062 | 578.828062 | 55.463666 |

Native cell values are copied exactly in priority order. No interpolation,
blending, smoothing, vertical correction, or negative-elevation removal
is applied at the integration stage. Source metadata documents CGVD2013,
but the local vertical realization and absolute vertical accuracy were
not independently verified.

## Critical limitations

Do not compute slope, aspect, curvature or other terrain derivatives
directly on this mixed-resolution mosaic. Use source-native neighbourhoods
with separate derivative coverage and consistency QA. Significant elevation
jumps occur at source changes and coarse native-cell boundaries. Cross-source
agreement, absolute vertical accuracy and hazard model validation remain
unaccepted or unvalidated. Negative elevations retain their source values;
the anomaly registers and QA reports must remain associated with this product.

The two release GeoTIFFs retain their validated candidate metadata, including
the historical product-status label, to preserve exact file SHA-256 checksums.
The publication manifest is authoritative for the local release status.
There has been no modification to the original published HRDEM products.
"""


def verify_package(directory,manifest):
    inventory=manifest["file_inventory"]
    for name,r in inventory.items():
        p=directory/name
        c.require(p.is_file() and p.stat().st_size==r["bytes"] and c.digest(p)==r["sha256"],
                  f"Package byte verification failed: {name}")
    c.require(c.digest(directory/"README.md")==manifest["readme_sha256"],"README hash mismatch")
    c.require(json.loads((directory/"publication_manifest.json").read_text(encoding="utf-8"))==manifest,
              "Publication manifest mismatch")
    expected=set(inventory)|{"README.md","publication_manifest.json"}
    observed={p.relative_to(directory).as_posix() for p in directory.rglob("*") if p.is_file()}
    c.require(expected==observed,"Unexpected publication payload inventory")


def publish_package(destination,paths,files,manifest):
    c.require(not destination.exists(),f"Refusing to replace existing release: {destination}")
    c.require(destination.parent.is_dir(),"Destination parent missing")
    size=sum(v["bytes"] for v in manifest["file_inventory"].values())
    c.require(shutil.disk_usage(destination.parent).free>=size+64*1024**2,"Insufficient disk space")
    tmp=Path(tempfile.mkdtemp(prefix=destination.name+"_partial_",dir=destination.parent))
    print("PUBLICATION_STAGING:",tmp,flush=True)
    for name,src in files.items():
        p=tmp/name
        p.parent.mkdir(parents=True,exist_ok=True)
        c.require(not p.exists(),f"Unexpected staging destination: {p}")
        shutil.copyfile(src,p)
        row=manifest["file_inventory"][name]
        c.require(p.stat().st_size==row["bytes"] and c.digest(p)==row["sha256"],
                  f"Copy hash mismatch: {name}")
        print("COPY_VERIFIED:",name,flush=True)
    manifest=dict(manifest)
    (tmp/"README.md").write_text(readme(manifest),encoding="utf-8")
    manifest["readme_sha256"]=c.digest(tmp/"README.md")
    # Protect every checked input, not only copied payloads.
    print("RECHECKING_ALL_EXISTING_INPUT_HASHES",flush=True)
    c.require(manifest["input_sha256"]=={k:c.digest(v) for k,v in paths.items()},
              "Input changed during publication")
    manifest.update({"published_at_utc":datetime.now(timezone.utc).isoformat(),
                     "input_files_unchanged":True,"final_integrated_dem_published":True,
                     "publication_qa_passed":True})
    (tmp/"publication_manifest.json").write_text(json.dumps(manifest,indent=2,allow_nan=False)+"\n",
                                               encoding="utf-8")
    verify_package(tmp,manifest)
    c.require(not destination.exists(),"Publication destination appeared during staging")
    tmp.rename(destination)
    verify_package(destination,manifest)
    print("EXISTING_INPUT_HASHES_UNCHANGED: True",flush=True)
    return manifest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--publish",action="store_true",help="Create local validated v1 release")
    args=parser.parse_args()
    c.init()
    destination=c.ELEV/DEST_NAME
    c.require(not destination.exists(),f"Refusing to replace existing release: {destination}")
    paths,files,manifest=build_plan()
    print("PUBLICATION_PLAN:",json.dumps({
        "destination":str(destination),"files_to_copy":len(files),
        "payload_bytes":sum(v["bytes"] for v in manifest["file_inventory"].values()),
        "publish_requested":args.publish}),flush=True)
    if not args.publish:
        print("PLAN_ONLY: no publication files written")
        return
    m=publish_package(destination,paths,files,manifest)
    print("PUBLISHED_REFERENCE_DEM:",destination/DEM_NAME)
    print("PUBLISHED_SOURCE_PROVENANCE:",destination/SOURCE_NAME)
    print("PUBLISHED_TRANSITION_REGISTER:",destination/f.REGISTER_NAME)
    print("PUBLICATION_MANIFEST:",destination/"publication_manifest.json")
    print("PUBLISHED_MOSAIC_SHA256:",m["mosaic_sha256"])
    print("PUBLISHED_PROVENANCE_SHA256:",m["source_provenance_sha256"])
    print("PUBLISHED_REGISTER_SHA256:",m["register_sha256"])
    print("PUBLICATION_QA_PASSED: True")
    print("FINAL_INTEGRATED_DEM_PUBLISHED: TRUE")
    print("TERRAIN_DERIVATIVES_FROM_MOSAIC_APPROVED: FALSE")
    print("ABSOLUTE_VERTICAL_ACCURACY_VALIDATED: FALSE")
    print("HAZARD_MODEL_VALIDATED: FALSE")


if __name__=="__main__":
    main()
