"""Create the HRDEM 2 m spatial anomaly register from completed QA reports.

Run from the natcat-risk project root, after all QA stages have passed.
This script reads existing inputs and creates a new diagnostics directory.
It does not edit the extraction, masking, or HRDEM 1 m products.
"""
from pathlib import Path
import json
import math
import tempfile

import numpy as np
from osgeo import gdal, ogr, osr

import step4_2_hrdem2m_fixed_batch_300 as e
import step4_2_hrdem2m_fixed_mask_batch_300 as m

METHOD = "HRDEM2M_ANOMALY_REGISTER_V1"
QA_PATHS = {
    "masking": "hrdem2m_global_masking_qa_opn0z04b/global_masking_qa.json",
    "screening": "hrdem2m_global_seam_qa_lp3mi24_/global_seam_qa.json",
    "bridge": "hrdem2m_targeted_seam_bridge_mjlpoefq/targeted_seam_bridge_qa.json",
    "native": "hrdem2m_targeted_native_source_x08sit7j/targeted_native_source_qa.json",
    "negative": "hrdem2m_negative_source_qa_7f428dzp/negative_source_qa.json",
}
GATES = {
    "masking": "global_masking_qa_passed",
    "screening": "global_seam_screening_completed",
    "bridge": "targeted_seam_bridge_qa_passed",
    "native": "targeted_native_source_qa_passed",
    "negative": "negative_source_qa_passed",
}
EXPECTED = {
    "flagged_seams": 20,
    "insufficient_seams": 92,
    "tiles_with_internal_nodata": 273,
    "negative_elevation_tiles": 119,
    "deep_negative_pixels": 606,
}
FLAGGED = {"localized_high_residual", "high_residual_nonlocalized"}
require = e.require


def inventory():
    return {str(p): e.digest(p) for d in (e.OUT, m.OUT)
            for p in d.iterdir() if p.is_file()}


def seam_key(r):
    return tuple(r["tile"]), r["direction"]


def spatial_ref():
    s = osr.SpatialReference()
    require(s.ImportFromEPSG(26910) == 0, "Cannot create EPSG:26910")
    s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return s


def line(k, d):
    x, y = k
    require(d in ("east", "north"), f"Invalid seam direction: {d}")
    pts = [(x+1000, y), (x+1000, y+1000)] if d == "east" else [
        (x, y+1000), (x+1000, y+1000)]
    g = ogr.Geometry(ogr.wkbLineString)
    for xx, yy in pts:
        g.AddPoint_2D(xx, yy)
    return g


def square(k):
    x, y = k
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for xx, yy in ((x,y), (x+1000,y), (x+1000,y+1000), (x,y+1000), (x,y)):
        ring.AddPoint_2D(xx, yy)
    g = ogr.Geometry(ogr.wkbPolygon)
    g.AddGeometry(ring)
    return g


def put(layer, geom, attrs):
    f = ogr.Feature(layer.GetLayerDefn())
    require(f.SetGeometry(geom) == 0, "Cannot assign geometry")
    for name, value in attrs.items():
        if value is not None:
            f.SetField(name, value)
    require(layer.CreateFeature(f) == 0, f"Cannot write {layer.GetName()}")


def create_layers(ds):
    s = spatial_ref()
    I, R, T = ogr.OFTInteger, ogr.OFTReal, ogr.OFTString
    base = [("tile_x", I), ("tile_y", I)]
    sf = base + [
        ("neighbor_x", I), ("neighbor_y", I), ("direction", T),
        ("classification", T), ("valid_crossings", I), ("context_crossings", I),
        ("seam_p95_m", R), ("matched_p95_m", R), ("baseline_p95_m", R),
        ("localization_ratio", R), ("bridge_passed", I), ("native_passed", I),
        ("strict_crossings", I), ("max_ref_diff_m", R),
        ("physical_validated", I), ("status", T),
    ]
    defs = {
        "flagged_seams": (ogr.wkbLineString, sf),
        "insufficient_seams": (ogr.wkbLineString, sf),
        "tiles_with_internal_nodata": (ogr.wkbPolygon, base + [
            ("footprint_pixels", I), ("valid_pixels", I), ("nodata_pixels", I),
            ("nodata_percent", R)]),
        "negative_elevation_tiles": (ogr.wkbPolygon, base + [
            ("negative_pixels", I), ("below_minus_1", I), ("below_minus_5", I),
            ("below_minus_10", I), ("below_minus_20", I), ("minimum_m", R),
            ("minimum_x", I), ("minimum_y", I), ("physical_validated", I)]),
        "deep_negative_pixels": (ogr.wkbPoint, base + [
            ("elevation_m", R), ("source_fidelity_passed", I),
            ("physical_validated", I)]),
    }
    layers = {}
    for name, (geom, fields) in defs.items():
        lyr = ds.CreateLayer(name, srs=s, geom_type=geom)
        require(lyr is not None, f"Cannot create layer {name}")
        for field, kind in fields:
            require(lyr.CreateField(ogr.FieldDefn(field, kind)) == 0,
                    f"Cannot create field {field}")
        layers[name] = lyr
    return layers


def write_register(gpkg, mr, flagged, insufficient, br, nr, negatives, deep):
    require(not gpkg.exists(), f"Refusing to overwrite {gpkg}")
    ds = ogr.GetDriverByName("GPKG").CreateDataSource(str(gpkg))
    require(ds is not None, "Cannot create GeoPackage")
    layers = create_layers(ds)
    require(ds.StartTransaction() == 0, "Cannot start transaction")
    try:
        for name, rows in (("flagged_seams", flagged), ("insufficient_seams", insufficient)):
            for r in rows:
                k, n, d = tuple(r["tile"]), r["neighbor"], r["direction"]
                b, a = br.get(seam_key(r)), nr.get(seam_key(r))
                put(layers[name], line(k, d), {
                    "tile_x": k[0], "tile_y": k[1], "neighbor_x": n[0], "neighbor_y": n[1],
                    "direction": d, "classification": r["classification"],
                    "valid_crossings": r["valid_crossings"],
                    "context_crossings": r["matched_context_crossings"],
                    "seam_p95_m": r["seam_p95_m"],
                    "matched_p95_m": r["matched_seam_p95_m"],
                    "baseline_p95_m": r["local_baseline_p95_m"],
                    "localization_ratio": r["localization_ratio"],
                    "bridge_passed": int(b["passed"]) if b else None,
                    "native_passed": int(a["passed"]) if a else None,
                    "strict_crossings": a["strict_seam_crossings"] if a else None,
                    "max_ref_diff_m": a["maximum_abs_reference_difference_m"] if a else None,
                    "physical_validated": 0,
                    "status": "numerically_source_consistent" if a else "insufficient_support",
                })
        for k, r in mr.items():
            nd = r["nodata_inside_footprint"]
            if nd > 0:
                put(layers["tiles_with_internal_nodata"], square(k), {
                    "tile_x": k[0], "tile_y": k[1], "footprint_pixels": r["footprint_pixels"],
                    "valid_pixels": r["masked_valid_pixels"], "nodata_pixels": nd,
                    "nodata_percent": 100 * nd / r["footprint_pixels"],
                })
        for r in negatives:
            k, c = tuple(r["tile"]), r["counts_below_m"]
            x, y = r["minimum_center_epsg26910"]
            put(layers["negative_elevation_tiles"], square(k), {
                "tile_x": k[0], "tile_y": k[1], "negative_pixels": r["negative_pixels"],
                "below_minus_1": c["-1"], "below_minus_5": c["-5"],
                "below_minus_10": c["-10"], "below_minus_20": c["-20"],
                "minimum_m": r["minimum_elevation_m"], "minimum_x": x,
                "minimum_y": y, "physical_validated": 0,
            })
        for k, r in deep.items():
            a, _ = e.inspect_tile(m.paths(k)[0], k)
            rr, cc = np.where(e.valid(a) & (a < -5))
            require(len(rr) == r["selected_pixels"], f"Deep negative count mismatch: {k}")
            for i, j in zip(rr, cc):
                g = ogr.Geometry(ogr.wkbPoint)
                g.AddPoint_2D(int(k[0]+2*j+1), int(k[1]+999-2*i))
                put(layers["deep_negative_pixels"], g, {
                    "tile_x": k[0], "tile_y": k[1], "elevation_m": float(a[i,j]),
                    "source_fidelity_passed": 1, "physical_validated": 0,
                })
            print("DEEP_NEGATIVE_TILE_REGISTERED:", k, len(rr), flush=True)
        require(ds.CommitTransaction() == 0, "Cannot commit register")
    except BaseException:
        ds.RollbackTransaction()
        raise
    finally:
        layers.clear()
        ds = None


def validate_register(gpkg, expected, mr, flagged, insufficient, negatives, deep):
    ds = ogr.Open(str(gpkg), 0)
    require(ds is not None and ds.GetLayerCount() == len(expected), "Unexpected layer count")
    counts = {l.GetName(): l.GetFeatureCount() for l in ds}
    require(counts == expected, f"Unexpected feature counts: {counts}")
    seam_ids = {
        "flagged_seams": {seam_key(r) for r in flagged},
        "insufficient_seams": {seam_key(r) for r in insufficient},
    }
    expected_tiles = {
        "tiles_with_internal_nodata": {k for k, r in mr.items() if r["nodata_inside_footprint"] > 0},
        "negative_elevation_tiles": {tuple(r["tile"]) for r in negatives},
    }
    pixel_counts = {k: 0 for k in deep}
    point_ids = set()
    for l in ds:
        name = l.GetName()
        require(l.GetSpatialRef().GetAuthorityCode(None) == "26910", f"Wrong CRS: {name}")
        seen = set()
        for f in l:
            g = f.GetGeometryRef()
            require(g is not None and not g.IsEmpty() and g.IsValid()
                    and g.GetGeometryType() == l.GetGeomType(), f"Invalid geometry: {name}")
            k = (f.GetField("tile_x"), f.GetField("tile_y"))
            if name in seam_ids:
                d = f.GetField("direction")
                ident = (k, d)
                require(ident not in seen and g.Equals(line(k, d)), "Incorrect seam geometry")
                seen.add(ident)
                if name == "flagged_seams":
                    require(f.GetField("bridge_passed") == f.GetField("native_passed") == 1
                            and f.GetField("physical_validated") == 0, "Incorrect seam QA status")
            elif name in expected_tiles:
                require(k not in seen and g.Equals(square(k)), "Incorrect tile geometry")
                seen.add(k)
            else:
                x, y = g.GetX(), g.GetY()
                require(k in deep and k[0] < x < k[0]+1000 and k[1] < y < k[1]+1000
                        and (x-k[0]) % 2 == 1 and (y-k[1]) % 2 == 1
                        and math.isfinite(f.GetField("elevation_m"))
                        and f.GetField("elevation_m") < -5
                        and f.GetField("source_fidelity_passed") == 1
                        and f.GetField("physical_validated") == 0, "Incorrect deep negative point")
                require((x,y) not in point_ids, "Duplicate pixel coordinate")
                point_ids.add((x,y))
                pixel_counts[k] += 1
        if name in seam_ids:
            require(seen == seam_ids[name], f"Unexpected seams in {name}")
        elif name in expected_tiles:
            require(seen == expected_tiles[name], f"Unexpected tiles in {name}")
        print("LAYER_VERIFIED:", name, l.GetFeatureCount(), flush=True)
    require(pixel_counts == {k: r["selected_pixels"] for k,r in deep.items()},
            "Deep negative point counts do not match source audit")
    f = g = l = None
    ds = None
    return counts


def main():
    gdal.UseExceptions()
    ogr.UseExceptions()
    osr.UseExceptions()
    tiles, _, ch, mc = m.config()
    mh = e.canonical_hash(mc)
    paths = {k: e.DIAG / v for k,v in QA_PATHS.items()}
    input_hashes = {k: e.digest(p) for k,p in paths.items()}
    q = {k: e.load(p) for k,p in paths.items()}
    for k, flag in GATES.items():
        require(q[k]["summary"].get(flag) is True and
                q[k]["summary"].get("mask_configuration_sha256") == mh,
                f"QA gate or mask configuration mismatch: {k}")
    for k in ("masking", "bridge", "native", "negative"):
        require(q[k]["summary"]["extraction_configuration_sha256"] == ch,
                f"Extraction configuration mismatch: {k}")
    before = inventory()
    mr = {tuple(r["tile"]): r for r in q["masking"]["tiles"]}
    require(set(mr) == set(tiles) and len(mr) == len(tiles) == 671, "Tile inventory mismatch")
    for k,r in mr.items():
        tif, qa = m.paths(k)
        require(before[str(qa)] == r["mask_qa_sha256"], f"Mask QA hash mismatch: {k}")
        require(r["state"] in ("COMPLETE", "EMPTY"), f"Invalid mask state: {k}")
        if r["state"] == "COMPLETE":
            require(before[str(tif)] == r["masked_raster_sha256"], f"Masked raster hash mismatch: {k}")
        else:
            require(not tif.exists(), f"Unexpected empty-tile raster: {k}")
    br = {seam_key(r): r for r in q["bridge"]["seams"]}
    nr = {seam_key(r): r for r in q["native"]["seams"]}
    flagged = [r for r in q["screening"]["seams"] if r["classification"] in FLAGGED]
    insufficient = [r for r in q["screening"]["seams"] if r["classification"].startswith("insufficient")]
    require(len(flagged) == 20 and len(insufficient) == 92 and
            {seam_key(r) for r in flagged} == set(br) == set(nr) and
            all(r["passed"] for r in list(br.values()) + list(nr.values())),
            "Seam audit inventory mismatch")
    require(all(r["strict_seam_crossings"] >= 20 for r in nr.values()), "Insufficient native seam support")
    negatives = q["negative"]["negative_tiles"]
    deep = {tuple(r["tile"]): r for r in q["negative"]["source_audits"] if r["kind"] == "deep_negative"}
    require(len(deep) == 2 and sum(r["selected_pixels"] for r in deep.values()) == 606 and
            all(r["numerical_fidelity_passed"] and r["without_strict_support"] == 0
                and r["strict_pixels"] == r["selected_pixels"] for r in deep.values()),
            "Deep negative source audit mismatch")
    require(len(negatives) == 119 and sum(r["negative_pixels"] for r in negatives) == 5045476,
            "Negative inventory mismatch")
    require(sum(r["nodata_inside_footprint"] > 0 for r in mr.values()) == 273, "NoData tile count mismatch")
    print("ALL_INPUT_QA_AND_MASK_HASHES_VERIFIED:", len(mr), flush=True)
    out = Path(tempfile.mkdtemp(prefix="hrdem2m_anomaly_register_", dir=e.DIAG))
    gpkg = out / "anomaly_register.gpkg"
    print("OUTPUT_DIRECTORY:", out, flush=True)
    write_register(gpkg, mr, flagged, insufficient, br, nr, negatives, deep)
    counts = validate_register(gpkg, EXPECTED, mr, flagged, insufficient, negatives, deep)
    require(before == inventory(), "Staging changed during register creation")
    require(input_hashes == {k: e.digest(p) for k,p in paths.items()}, "QA reports changed")
    summary = {
        "method_id": METHOD,
        "script_sha256": e.digest(Path(__file__).resolve()),
        "extraction_configuration_sha256": ch,
        "mask_configuration_sha256": mh,
        "crs": "EPSG:26910",
        "layer_counts": counts,
        "total_features": sum(counts.values()),
        "input_qa_paths": QA_PATHS,
        "input_qa_sha256": input_hashes,
        "register_sha256": e.digest(gpkg),
        "staging_unchanged": True,
        "all_geometry_validated": True,
        "register_completed": True,
        "absolute_vertical_accuracy_validated": False,
        "negative_source_physical_validity_resolved": False,
        "final_dem_published": False,
    }
    e.new_json(out / "anomaly_register_qa.json", summary)
    print("ANOMALY_REGISTER_SUMMARY:", json.dumps(summary, indent=2))
    print("REGISTER_GPKG:", gpkg)
    print("DIAGNOSTIC_JSON:", out / "anomaly_register_qa.json")
    print("ANOMALY_REGISTER_COMPLETED: True")
    print("FINAL_DEM: FALSE")


if __name__ == "__main__":
    main()
