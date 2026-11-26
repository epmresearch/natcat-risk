"""HRDEM 2 m fixed-bilinear extraction, with a resumable, isolated pilot.

Run without flags for a read-only plan. Run --pilot before --batch.
The script does not mask tiles, change legacy rasters, or publish a DEM.
Do not run multiple instances against the same staging directory.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
from uuid import uuid4

import numpy as np
from osgeo import gdal

ROOT = Path(__file__).resolve().parent
DIAG = ROOT / "results" / "diagnostics"
TEST = DIAG / "hrdem2m_scale_control_1zla68yg"
PLAN = TEST / "production_plan_qa.json"
FOOTPRINT = ROOT / "data" / "processed" / "elevation" / "hrdem_2m_dtm_coverage_highway99_corrected.gpkg"
VRT = DIAG / "step4_2_hrdem2m_boundary_virtual.vrt"
OUT = DIAG / "hrdem2m_fixed_batch_staging"
METHOD = "HRDEM2M_FIXED_BILINEAR_V1"
NODATA = -32767.0
SIZE = 500
TILE_M = 1000
SITES = {
    "pilot_2_3": (484264, 5472056),
    "pilot_2_4": (492112, 5517584),
    "pilot_cross_boundary": (482894, 5481648),
}
OPTIONS = [
    "-of", "GTiff", "-t_srs", "EPSG:26910",
    "-tr", "2", "2", "-tap", "-r", "bilinear", "-et", "0",
    "-ovr", "NONE", "-ot", "Float32",
    "-srcnodata", "-32767", "-dstnodata", "-32767",
    "-wo", "XSCALE=1", "-wo", "YSCALE=1",
    "-co", "TILED=YES", "-co", "COMPRESS=DEFLATE",
    "-co", "PREDICTOR=3", "-co", "BIGTIFF=IF_SAFER",
]
HTTP = {
    "GDAL_DISABLE_READDIR_ON_OPEN": "YES",
    "GDAL_HTTP_MAX_RETRY": "3",
    "GDAL_HTTP_RETRY_DELAY": "2",
    "GDAL_HTTP_TIMEOUT": "45",
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    h = sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def canonical_hash(data):
    return sha256(json.dumps(data, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def new_json(path, data):
    require(not path.exists(), f"Refusing to replace existing JSON: {path}")
    tmp = path.with_name(path.name + "." + uuid4().hex + ".partial")
    try:
        tmp.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        require(not path.exists(), f"Unexpected existing output: {path}")
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def site_tile(x, y):
    return x // TILE_M * TILE_M, y // TILE_M * TILE_M


def reference_tiles(center):
    cx, cy = center
    xlo, ylo = site_tile(cx-32, cy-128)
    xhi, yhi = site_tile(cx+31, cy+127)
    return [(x, y) for y in range(ylo, yhi+1, TILE_M)
            for x in range(xlo, xhi+1, TILE_M)]


def pilot_tiles(tiles):
    pool = set(tiles)
    anchors = [site_tile(*xy) for xy in SITES.values()]
    result = list(dict.fromkeys(k for xy in SITES.values() for k in reference_tiles(xy)))
    require(set(result) <= pool, "A validated reference window extends outside the production plan")
    for x, y in anchors:
        neighbors = [(x+TILE_M, y), (x-TILE_M, y), (x, y+TILE_M), (x, y-TILE_M)]
        if any(k in result for k in neighbors):
            continue
        n = next((k for k in neighbors if k in pool), None)
        require(n is not None, f"No adjacent candidate tile for {(x,y)}")
        result.append(n)
    for edge in (min(tiles, key=lambda k: (k[1], k[0])), max(tiles, key=lambda k: (k[1], k[0]))):
        if edge not in result:
            result.append(edge)
    return result


def configuration():
    plan = load(PLAN)
    require(plan.get("method_id") == METHOD, "Unexpected plan method")
    require(plan.get("crs") == "EPSG:26910" and plan.get("resolution_m") == 2, "Wrong plan grid")
    require(plan.get("tile_size_m") == TILE_M and plan.get("tile_dimensions") == [SIZE, SIZE],
            "Wrong plan tile dimensions")
    require(plan.get("candidate_rule") == "ALL_TOUCHED_1KM_FOOTPRINT_GRID", "Wrong candidate rule")
    tiles = [tuple(k) for k in plan["candidate_tiles"]]
    require(len(tiles) == plan["candidate_count"] == len(set(tiles)) == 671, "Unexpected candidate inventory")
    require(all(len(k) == 2 and all(isinstance(v, int) and v % TILE_M == 0 for v in k) for k in tiles),
            "Invalid tile coordinates")
    xmin, ymin, xmax, ymax = plan["bounds"]
    require(plan["bounds"] == [475000, 5463000, 506000, 5557000], "Unexpected published planning bounds")
    require(all(xmin <= x < xmax and ymin <= y < ymax for x, y in tiles), "A tile lies outside the plan")
    require(digest(FOOTPRINT) == plan["footprint_sha256"], "Footprint has changed since planning")
    require(digest(VRT) == plan["vrt_sha256"], "VRT has changed since planning")
    extent_qa = TEST / "scale_comparison_qa.json"
    native_qa = TEST / "native_bilinear_reference_qa.json"
    require(load(extent_qa).get("fixed_extent_consistency_passed") is True, "Extent QA gate not passed")
    require(load(native_qa).get("native_bilinear_reference_passed") is True, "Native reference QA gate not passed")
    files = [PLAN, FOOTPRINT, VRT, extent_qa, native_qa, TEST / "warp_commands.json"]
    files += [TEST / (name + "_narrow_fixed.tif") for name in SITES]
    cfg = {
        "method_id": METHOD,
        "input_sha256": {f.relative_to(ROOT).as_posix(): digest(f) for f in files},
        "script_sha256": digest(Path(__file__).resolve()),
        "gdal_python_version": gdal.VersionInfo("RELEASE_NAME"),
        "gdal_cli_version": subprocess.check_output(["gdalwarp", "--version"], text=True).strip(),
        "warp_options": OPTIONS,
        "http_config": HTTP,
        "pilot_tiles": [list(k) for k in pilot_tiles(tiles)],
        "candidate_count": len(tiles),
        "absolute_vertical_accuracy_validated": False,
    }
    return tiles, cfg


def paths(tile):
    x, y = tile
    stem = f"hrdem2m_x{x}_y{y}_utm10_1km"
    return OUT / (stem + ".tif"), OUT / (stem + ".qa.json")


def valid(a):
    return np.isfinite(a) & (a != NODATA)


def inspect_tile(path, tile):
    x, y = tile
    ds = gdal.Open(str(path), gdal.GA_ReadOnly)
    require(ds is not None, f"Cannot open raster: {path}")
    require((ds.RasterXSize, ds.RasterYSize, ds.RasterCount) == (SIZE, SIZE, 1), f"Wrong dimensions: {path}")
    require(ds.GetGeoTransform() == (x, 2, 0, y+TILE_M, 0, -2), f"Wrong geotransform: {path}")
    require(ds.GetSpatialRef() is not None and ds.GetSpatialRef().GetAuthorityCode(None) == "26910",
            f"Wrong CRS: {path}")
    band = ds.GetRasterBand(1)
    require(band.DataType == gdal.GDT_Float32 and band.GetNoDataValue() == NODATA, f"Wrong dtype/NoData: {path}")
    a = band.ReadAsArray()
    require(np.isfinite(a).all(), f"Nonfinite pixels: {path}")
    m = valid(a)
    require(np.array_equal(m, band.GetMaskBand().ReadAsArray() > 0), f"Validity mask mismatch: {path}")
    n = int(m.sum())
    stats = {
        "total_pixels": SIZE * SIZE,
        "valid_pixels": n,
        "nodata_pixels": SIZE * SIZE - n,
        "min_elevation_m": float(a[m].min()) if n else None,
        "max_elevation_m": float(a[m].max()) if n else None,
    }
    band = None
    ds = None
    return a, stats


def state(tile, cfg_sha):
    tif, qa = paths(tile)
    if not qa.exists():
        require(not tif.exists(), f"Orphan raster without QA. Review manually: {tif}")
        return "PENDING"
    rec = load(qa)
    require(rec.get("tile") == list(tile) and rec.get("method_id") == METHOD and
            rec.get("configuration_sha256") == cfg_sha, f"QA provenance mismatch: {qa}")
    n = rec.get("valid_pixels")
    require(isinstance(n, int) and 0 <= n <= SIZE*SIZE and
            rec.get("nodata_pixels") == SIZE*SIZE-n and rec.get("total_pixels") == SIZE*SIZE,
            f"Invalid QA counts: {qa}")
    if rec.get("status") == "EMPTY":
        require(n == 0 and not tif.exists() and rec.get("raster_sha256") is None, f"Invalid empty state: {qa}")
        return "EMPTY"
    require(rec.get("status") == "COMPLETE" and n > 0 and tif.is_file(),
            f"Missing or inconsistent complete state: {qa}")
    require(digest(tif) == rec.get("raster_sha256"), f"Raster hash mismatch: {tif}")
    return "COMPLETE"


def process(tile, cfg_sha):
    existing = state(tile, cfg_sha)
    if existing != "PENDING":
        print("SKIPPED_VERIFIED:", tile, existing, flush=True)
        return "SKIPPED"
    x, y = tile
    tif, qa = paths(tile)
    tmp = tif.with_name(tif.stem + "." + uuid4().hex + ".partial.tif")
    cmd = ["gdalwarp"] + OPTIONS + ["-te", str(x), str(y), str(x+TILE_M), str(y+TILE_M), str(VRT), str(tmp)]
    try:
        print("EXTRACTING_FIXED:", tif.stem, flush=True)
        subprocess.run(cmd, env={**os.environ, **HTTP}, check=True)
        _, stats = inspect_tile(tmp, tile)
        n = stats["valid_pixels"]
        rec = {
            "method_id": METHOD,
            "configuration_sha256": cfg_sha,
            "tile": list(tile),
            "status": "COMPLETE" if n else "EMPTY",
            "qa_level": "STRUCTURAL_ONLY_NOT_VERTICAL_ACCURACY",
            **stats,
            "raster_sha256": digest(tmp) if n else None,
        }
        require(not tif.exists() and not qa.exists(), "Output appeared during processing; stop concurrent writers")
        if n:
            os.replace(tmp, tif)
        new_json(qa, rec)
        print("FIXED_STRUCTURAL_QA_PASSED:" if n else "EMPTY_FIXED:", tif.stem, "valid_pixels:", n, flush=True)
        return "SUCCESS" if n else "EMPTY"
    finally:
        tmp.unlink(missing_ok=True)


def review_pilot(cfg_sha, pilots):
    records = []
    for name, center in SITES.items():
        tiles = reference_tiles(center)
        ref_ds = gdal.Open(str(TEST / (name + "_narrow_fixed.tif")), gdal.GA_ReadOnly)
        cx, cy = center
        require((ref_ds.RasterXSize, ref_ds.RasterYSize, ref_ds.RasterCount) == (32,128,1)
                and ref_ds.GetGeoTransform() == (cx-32,2,0,cy+128,0,-2)
                and ref_ds.GetRasterBand(1).GetNoDataValue() == NODATA
                and ref_ds.GetRasterBand(1).DataType == gdal.GDT_Float32
                and ref_ds.GetSpatialRef().GetAuthorityCode(None) == "26910", f"Wrong reference raster: {name}")
        b = ref_ds.ReadAsArray()
        a = np.full((128, 32), NODATA, dtype=np.float32)
        filled = np.zeros(a.shape, dtype=bool)
        for tile in tiles:
            require(state(tile, cfg_sha) == "COMPLETE", f"Reference pilot unexpectedly empty: {tile}")
            data, _ = inspect_tile(paths(tile)[0], tile)
            tx, ty = tile
            xmin, xmax = max(tx, cx-32), min(tx+TILE_M, cx+32)
            ymin, ymax = max(ty, cy-128), min(ty+TILE_M, cy+128)
            width, height = (xmax-xmin)//2, (ymax-ymin)//2
            sx, sy = (xmin-tx)//2, (ty+TILE_M-ymax)//2
            dx, dy = (xmin-(cx-32))//2, (cy+128-ymax)//2
            dest = np.s_[dy:dy+height, dx:dx+width]
            require(not filled[dest].any(), f"Overlapping reference tile windows: {name}")
            a[dest] = data[sy:sy+height, sx:sx+width]
            filled[dest] = True
        require(filled.all(), f"Reference window has missing tile coverage: {name}")
        va, vb = valid(a), valid(b)
        m = va & vb
        delta = np.abs(a[m].astype(np.float64)-b[m].astype(np.float64))
        rec = {
            "site": name,
            "tiles": [list(k) for k in tiles],
            "common_valid": int(m.sum()),
            "validity_mismatch": int((va != vb).sum()),
            "max_abs_diff_m": float(delta.max()) if delta.size else None,
            "over_0_001m": int((delta > 0.001).sum()),
            "all_raw_values_exact": bool(np.array_equal(a, b)),
        }
        rec["passed"] = (rec["common_valid"] >= 1000 and rec["validity_mismatch"] == 0 and
                         rec["over_0_001m"] == 0 and bool(np.isfinite(b).all()))
        records.append(rec)
        print("PILOT_REFERENCE_COMPARISON:", json.dumps(rec), flush=True)
        ref_ds = None
    tile_states = [{"tile": list(k), "state": state(k, cfg_sha),
                    "valid_pixels": load(paths(k)[1])["valid_pixels"]} for k in pilots]
    passed = all(r["passed"] for r in records) and all(v["state"] in ("COMPLETE", "EMPTY") for v in tile_states)
    summary = {"configuration_sha256": cfg_sha, "tile_states": tile_states,
               "reference_comparisons": records, "production_pilot_passed": passed,
               "masking_validated": False, "seam_qa_completed": False,
               "absolute_vertical_accuracy_validated": False, "final_dem_published": False}
    p = OUT / "pilot_extraction_qa.json"
    if p.exists():
        require(load(p) == summary, f"Pilot review differs from existing report: {p}")
    else:
        new_json(p, summary)
    print("PILOT_QA_JSON:", p, flush=True)
    print("PRODUCTION_PILOT_PASSED:", passed, flush=True)
    require(passed, "Pilot reference QA failed; do not start batch extraction")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    group = ap.add_mutually_exclusive_group()
    group.add_argument("--pilot", action="store_true", help="Extract and verify pilot tiles")
    group.add_argument("--batch", action="store_true", help="Extract pending tiles after pilot QA")
    ap.add_argument("--limit", type=int, default=300, help="Maximum pending tiles per batch")
    args = ap.parse_args()
    require(args.limit > 0, "--limit must be positive")
    gdal.UseExceptions()
    tiles, cfg = configuration()
    cfg_sha = canonical_hash(cfg)
    pilots = pilot_tiles(tiles)
    print("METHOD:", METHOD, "\nCANDIDATES:", len(tiles), "\nOUTPUT_DIRECTORY:", OUT, flush=True)
    print("CONFIGURATION_SHA256:", cfg_sha, flush=True)
    print("PILOT_TILES:", pilots, flush=True)
    config_file = OUT / "configuration.json"
    if OUT.exists():
        require(config_file.is_file(), f"Staging exists without configuration: {OUT}")
        require(load(config_file) == cfg, "Staging configuration differs; do not mix runs")
        partials = list(OUT.glob("*partial*"))
        require(not partials, f"Incomplete prior run requires review: {partials}")
        known = {p.name for k in tiles for p in paths(k)}
        known.update(("configuration.json", "pilot_extraction_qa.json"))
        unexpected = [p.name for p in OUT.iterdir() if p.name not in known]
        require(not unexpected, f"Unexpected staging entries: {unexpected}")
    elif args.pilot or args.batch:
        require(args.pilot, "Run --pilot first")
        OUT.mkdir(parents=True, exist_ok=False)
        new_json(config_file, cfg)
    if not (args.pilot or args.batch):
        states = [state(k, cfg_sha) for k in tiles] if OUT.exists() else ["PENDING"] * len(tiles)
        print("STATE_COUNTS:", {s: states.count(s) for s in ("COMPLETE","EMPTY","PENDING")})
        print("PLAN_ONLY: no network raster reads or output writes")
        return
    if args.pilot:
        selected = pilots
    else:
        review = load(OUT / "pilot_extraction_qa.json")
        require(review.get("configuration_sha256") == cfg_sha and
                review.get("production_pilot_passed") is True, "Pilot QA gate not passed")
        review_pilot(cfg_sha, pilots)
        selected = [k for k in tiles if state(k, cfg_sha) == "PENDING"][:args.limit]
    results = [process(k, cfg_sha) for k in selected]
    print("FIXED_EXTRACTION_FINISHED:", {s: results.count(s) for s in ("SUCCESS","EMPTY","SKIPPED")}, flush=True)
    if args.pilot:
        review_pilot(cfg_sha, pilots)
    print("FINAL_DEM: FALSE")
    print("EXISTING_1M_DEM_MODIFIED: FALSE")


if __name__ == "__main__":
    try:
        main()
    except (Exception, KeyboardInterrupt) as exc:
        print("ERROR:", exc, file=sys.stderr)
        sys.exit(1)
