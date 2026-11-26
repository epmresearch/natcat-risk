"""Pixel-center footprint masking for the verified HRDEM 2 m extraction.

Default: read-only plan. Run --pilot before --batch.
Empty extraction candidates get a QA marker; completed extractions retain a
masked TIFF even if all its pixels are NoData. No interpolation is performed.
Never run concurrent instances or edit the extraction script after staging.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
from uuid import uuid4

import numpy as np
from osgeo import gdal, ogr
import step4_2_hrdem2m_fixed_batch_300 as e

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "diagnostics" / "hrdem2m_fixed_masked_staging"
SEAM_QA = ROOT / "results" / "diagnostics" / "hrdem2m_pilot_seam_edge_gn4y2wf_" / "pilot_seam_edge_qa.json"
METHOD = "HRDEM2M_FIXED_PIXEL_CENTER_MASK_V1"
NODATA = e.NODATA
SIZE = e.SIZE
TOTAL = SIZE * SIZE


def require(condition, message):
    e.require(condition, message)


def config():
    tiles, extraction = e.configuration()
    ch = e.canonical_hash(extraction)
    require(e.load(e.OUT / "configuration.json") == extraction, "Extraction configuration mismatch")
    pilot = e.load(e.OUT / "pilot_extraction_qa.json")
    seam = e.load(SEAM_QA)
    require(pilot.get("configuration_sha256") == ch and pilot.get("production_pilot_passed") is True,
            "Production pilot gate not passed")
    require(seam.get("configuration_sha256") == ch and seam.get("pilot_seam_edge_qa_passed") is True and
            seam.get("staging_unchanged") is True, "Pilot seam/edge gate not passed")
    cfg = {
        "mask_method": METHOD,
        "extraction_configuration_sha256": ch,
        "footprint_sha256": e.digest(e.FOOTPRINT),
        "extraction_pilot_qa_sha256": e.digest(e.OUT / "pilot_extraction_qa.json"),
        "seam_edge_qa_sha256": e.digest(SEAM_QA),
        "mask_script_sha256": e.digest(Path(__file__).resolve()),
        "gdal_version": gdal.VersionInfo("RELEASE_NAME"),
        "mask_rule": "GDAL_PIXEL_CENTER_FOOTPRINT_RASTERIZATION_ALL_TOUCHED_FALSE",
        "resolution_m": 2,
        "crs": "EPSG:26910",
        "nodata": NODATA,
        "absolute_vertical_accuracy_validated": False,
    }
    return tiles, extraction, ch, cfg


def paths(tile):
    stem = e.paths(tile)[0].stem + "_masked"
    return OUT / (stem + ".tif"), OUT / (stem + ".qa.json")


def mask_for_tile(layer, tile):
    x, y = tile
    ds = gdal.GetDriverByName("MEM").Create("", SIZE, SIZE, 1, gdal.GDT_Byte)
    ds.SetGeoTransform((x, 2, 0, y+e.TILE_M, 0, -2))
    ds.SetProjection(layer.GetSpatialRef().ExportToWkt())
    ds.GetRasterBand(1).Fill(0)
    layer.ResetReading()
    require(gdal.RasterizeLayer(ds, [1], layer, burn_values=[1], options=["ALL_TOUCHED=FALSE"]) == 0,
            f"Footprint rasterization failed: {tile}")
    raw = ds.ReadAsArray()
    require(raw.shape == (SIZE, SIZE) and np.isin(raw, [0,1]).all(), f"Invalid mask: {tile}")
    return raw == 1


def expected_values(layer, tile, source_state):
    footprint = mask_for_tile(layer, tile)
    if source_state == "COMPLETE":
        a, _ = e.inspect_tile(e.paths(tile)[0], tile)
    else:
        require(source_state == "EMPTY", f"Extraction is not ready: {tile}")
        a = np.full((SIZE, SIZE), NODATA, dtype=np.float32)
    valid = e.valid(a)
    keep = valid & footprint
    result = np.where(keep, a, NODATA).astype(np.float32)
    n = int(keep.sum())
    stats = {
        "source_valid_pixels": int(valid.sum()),
        "footprint_pixels": int(footprint.sum()),
        "valid_inside_footprint": n,
        "nodata_inside_footprint": int((footprint & ~valid).sum()),
        "valid_outside_footprint_removed": int((valid & ~footprint).sum()),
        "masked_valid_pixels": n,
        "masked_nodata_pixels": TOTAL-n,
        "valid_outside_footprint": int((e.valid(result) & ~footprint).sum()),
        "retained_values_bitwise_equal": bool(np.array_equal(result[keep].view(np.uint32), a[keep].view(np.uint32))),
    }
    require(stats["source_valid_pixels"] == n + stats["valid_outside_footprint_removed"] and
            stats["footprint_pixels"] == n + stats["nodata_inside_footprint"] and
            stats["valid_outside_footprint"] == 0 and stats["retained_values_bitwise_equal"], f"Mask invariants failed: {tile}")
    return result, stats


def state(tile, source_state, ch, mh):
    tif, qa = paths(tile)
    if not qa.exists():
        require(not tif.exists(), f"Orphan masked TIFF: {tif}")
        return "PENDING"
    require(source_state in ("COMPLETE", "EMPTY"), f"Masked result without completed extraction: {tile}")
    rec = e.load(qa)
    require(rec.get("tile") == list(tile) and rec.get("mask_method") == METHOD and
            rec.get("mask_configuration_sha256") == mh and rec.get("extraction_configuration_sha256") == ch and
            rec.get("source_state") == source_state, f"Mask provenance mismatch: {qa}")
    source_tif, source_qa = e.paths(tile)
    require(rec.get("source_qa_sha256") == e.digest(source_qa), f"Source QA hash changed: {source_qa}")
    source_hash = e.load(source_qa).get("raster_sha256")
    require(rec.get("source_raster_sha256") == source_hash, f"Source raster provenance changed: {tile}")
    require(rec.get("source_valid_pixels") == e.load(source_qa).get("valid_pixels"), f"Source pixel count mismatch: {tile}")
    n = rec.get("masked_valid_pixels")
    require(isinstance(n, int) and 0 <= n <= TOTAL and rec.get("masked_nodata_pixels") == TOTAL-n and
            rec.get("valid_inside_footprint") == n and rec.get("valid_outside_footprint") == 0 and
            rec.get("retained_values_bitwise_equal") is True and
            rec.get("source_valid_pixels") == n + rec.get("valid_outside_footprint_removed", -1) and
            rec.get("footprint_pixels") == n + rec.get("nodata_inside_footprint", -1),
            f"Mask QA accounting mismatch: {qa}")
    if source_state == "EMPTY":
        require(rec.get("status") == "EMPTY_EXTRACTION" and n == 0 and not tif.exists() and
                rec.get("raster_sha256") is None, f"Invalid empty mask state: {qa}")
        return "EMPTY"
    require(rec.get("status") == "COMPLETE" and tif.is_file() and
            e.digest(tif) == rec.get("raster_sha256"), f"Masked raster hash or state mismatch: {tif}")
    return "COMPLETE"


def verify(path, expected, tile):
    actual, _ = e.inspect_tile(path, tile)
    require(np.array_equal(actual.view(np.uint32), expected.view(np.uint32)), f"Masked readback mismatch: {path}")


def mask_one(layer, tile, source_state, ch, mh):
    require(source_state in ("COMPLETE", "EMPTY"), f"Extraction incomplete: {tile}")
    source_tif, source_qa = e.paths(tile)
    require(e.state(tile, ch) == source_state, f"Extraction state changed: {tile}")
    source_qa_sha = e.digest(source_qa)
    source_sha = e.load(source_qa)["raster_sha256"]
    tif, qa = paths(tile)
    existing = state(tile, source_state, ch, mh)
    expected, stats = expected_values(layer, tile, source_state)
    require(stats["source_valid_pixels"] == e.load(source_qa)["valid_pixels"], f"Extraction count mismatch: {tile}")
    if existing != "PENDING":
        rec = e.load(qa)
        require(all(rec.get(k) == v for k,v in stats.items()), f"Existing mask QA differs: {qa}")
        if existing == "COMPLETE":
            verify(tif, expected, tile)
        print("SKIPPED_VERIFIED_MASK:", tile, existing, flush=True)
        return "SKIPPED"
    tmp = tif.with_name(tif.stem + "." + uuid4().hex + ".partial.tif")
    ds = None
    try:
        raster_hash = None
        if source_state == "COMPLETE":
            x, y = tile
            ds = gdal.GetDriverByName("GTiff").Create(str(tmp), SIZE, SIZE, 1, gdal.GDT_Float32,
                                                     options=["TILED=YES","COMPRESS=DEFLATE","PREDICTOR=3"])
            ds.SetGeoTransform((x,2,0,y+e.TILE_M,0,-2))
            ds.SetProjection(layer.GetSpatialRef().ExportToWkt())
            band = ds.GetRasterBand(1)
            band.SetNoDataValue(NODATA)
            band.WriteArray(expected)
            band.FlushCache()
            band = None
            ds.FlushCache()
            ds = None
            verify(tmp, expected, tile)
            raster_hash = e.digest(tmp)
        require(e.digest(source_qa) == source_qa_sha and
                (source_state == "EMPTY" or e.digest(source_tif) == source_sha), f"Input changed during masking: {tile}")
        require(not tif.exists() and not qa.exists(), "Output appeared during processing; stop concurrent writers")
        if source_state == "COMPLETE":
            os.replace(tmp, tif)
        rec = {
            "tile": list(tile),
            "mask_method": METHOD,
            "mask_configuration_sha256": mh,
            "extraction_configuration_sha256": ch,
            "source_state": source_state,
            "source_qa_sha256": source_qa_sha,
            "source_raster_sha256": source_sha,
            "raster_sha256": raster_hash,
            "status": "COMPLETE" if source_state == "COMPLETE" else "EMPTY_EXTRACTION",
            "qa_level": "EXACT_MASK_READBACK_NOT_VERTICAL_ACCURACY",
            **stats,
        }
        e.new_json(qa, rec)
        print("FIXED_MASK_QA_PASSED:" if source_state == "COMPLETE" else "EMPTY_EXTRACTION_ACCOUNTED:",
              tile, json.dumps(stats), flush=True)
        return "SUCCESS" if source_state == "COMPLETE" else "EMPTY"
    finally:
        ds = None
        tmp.unlink(missing_ok=True)


def inventory_hashes(directory):
    return {p.name:e.digest(p) for p in directory.iterdir() if p.is_file()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--pilot", action="store_true")
    group.add_argument("--batch", action="store_true")
    parser.add_argument("--limit", type=int, default=300)
    args = parser.parse_args()
    require(1 <= args.limit <= 300, "--limit must be within 1..300")
    gdal.UseExceptions()
    ogr.UseExceptions()
    tiles, extraction, ch, cfg = config()
    mh = e.canonical_hash(cfg)
    pilots = e.pilot_tiles(tiles)
    source_states = {k:e.state(k,ch) for k in tiles}
    print("MASK_METHOD:", METHOD, flush=True)
    print("MASK_CONFIGURATION_SHA256:", mh, flush=True)
    print("INPUT_DIRECTORY:", e.OUT, flush=True)
    print("OUTPUT_DIRECTORY:", OUT, flush=True)
    print("EXTRACTION_STATE_COUNTS:", {s:list(source_states.values()).count(s) for s in ("COMPLETE","EMPTY","PENDING")}, flush=True)
    print("PILOT_TILES:", pilots, flush=True)
    conf = OUT / "configuration.json"
    if OUT.exists():
        require(conf.is_file() and e.load(conf) == cfg, "Mask staging configuration mismatch")
        require(not list(OUT.glob("*partial*")), "Incomplete masked files require review")
        known = {p.name for k in tiles for p in paths(k)}
        known.update(("configuration.json","pilot_mask_qa.json"))
        require(all(p.name in known for p in OUT.iterdir()), "Unexpected files in mask staging")
    elif args.pilot or args.batch:
        require(args.pilot, "Run --pilot first")
        OUT.mkdir(parents=True, exist_ok=False)
        e.new_json(conf, cfg)
    masks = {k:state(k,source_states[k],ch,mh) for k in tiles} if OUT.exists() else {k:"PENDING" for k in tiles}
    print("MASK_STATE_COUNTS:", {s:list(masks.values()).count(s) for s in ("COMPLETE","EMPTY","PENDING")}, flush=True)
    if not (args.pilot or args.batch):
        print("PLAN_ONLY: no network raster reads or output writes")
        return
    if args.pilot:
        selected = pilots
        require(all(source_states[k] in ("COMPLETE","EMPTY") for k in selected), "Extraction pilot is incomplete")
        before = inventory_hashes(e.OUT)
    else:
        require(all(v in ("COMPLETE","EMPTY") for v in source_states.values()), "Complete extraction before full masking")
        review = e.load(OUT/"pilot_mask_qa.json")
        require(review.get("mask_configuration_sha256") == mh and review.get("pilot_mask_qa_passed") is True,
                "Mask pilot gate not passed")
        selected = [k for k in tiles if masks[k] == "PENDING"][:args.limit]
    source = ogr.Open(str(e.FOOTPRINT), 0)
    require(source is not None, "Cannot open footprint")
    layer = source.GetLayer(0)
    require(layer.GetSpatialRef() is not None and layer.GetSpatialRef().GetAuthorityCode(None) == "26910",
            "Wrong footprint CRS")
    results = [mask_one(layer,k,source_states[k],ch,mh) for k in selected]
    print("FIXED_MASK_FINISHED:", {s:results.count(s) for s in ("SUCCESS","EMPTY","SKIPPED")}, flush=True)
    if args.pilot:
        unchanged = before == inventory_hashes(e.OUT)
        reports = [e.load(paths(k)[1]) for k in pilots]
        passed = unchanged and all(r["retained_values_bitwise_equal"] and r["valid_outside_footprint"] == 0 for r in reports)
        summary = {
            "mask_configuration_sha256": mh,
            "extraction_configuration_sha256": ch,
            "pilot_tile_count": len(pilots),
            "source_staging_unchanged": unchanged,
            "total_footprint_pixels": sum(r["footprint_pixels"] for r in reports),
            "total_masked_valid_pixels": sum(r["masked_valid_pixels"] for r in reports),
            "total_nodata_inside_footprint": sum(r["nodata_inside_footprint"] for r in reports),
            "total_removed_valid_outside_footprint": sum(r["valid_outside_footprint_removed"] for r in reports),
            "total_valid_outside_footprint": sum(r["valid_outside_footprint"] for r in reports),
            "all_retained_values_bitwise_equal": all(r["retained_values_bitwise_equal"] for r in reports),
            "pilot_mask_qa_passed": passed,
            "final_dem_published": False,
        }
        report = OUT / "pilot_mask_qa.json"
        if report.exists():
            require(e.load(report) == summary, "Existing pilot mask summary differs")
        else:
            e.new_json(report, summary)
        print("PILOT_MASK_SUMMARY:", json.dumps(summary), flush=True)
        print("PILOT_MASK_QA_PASSED:", passed, flush=True)
        require(passed, "Mask pilot QA failed")
    print("FINAL_DEM: FALSE")
    print("EXISTING_1M_DEM_MODIFIED: FALSE")


if __name__ == "__main__":
    try:
        main()
    except (Exception, KeyboardInterrupt) as exc:
        print("ERROR:", exc, file=sys.stderr)
        sys.exit(1)
