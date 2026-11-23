"""Independent, offline vector coverage QA for Highway 99 DEM source discovery.

This validates *coverage polygons*, not the DEM raster pixel values or vertical accuracy.
No downloads or edits to input files are made.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import geopandas as gpd

# Rounded values recorded during manual QGIS checks (km2)
EXPECTED = {
    "aoi_km2": 1043.616776,
    "hrdem_1m_km2": 464.219678,
    "hrdem_2m_km2": 462.717825,
    "combined_km2": 926.937503,
}
# 0.01 km2 = 10,000 m2. Meaningful discrepancies should not be hidden by rounding.
AREA_TOLERANCE_KM2 = 0.01

SOURCES = {
    "aoi": Path("data/processed/study_area/highway99_terrain_search_5000m.geojson"),
    "hrdem_1m": Path("data/processed/elevation/hrdem_1m_dtm_coverage_highway99.gpkg"),
    "hrdem_2m": Path("data/processed/elevation/hrdem_2m_dtm_coverage_highway99.gpkg"),
    "combined": Path("data/processed/elevation/hrdem_combined_dtm_coverage_highway99.gpkg"),
}


def area_km2(geometry: Any) -> float:
    return geometry.area / 1_000_000.0


def coverage_geometry(path: Path, label: str, target_epsg: int = 26910) -> tuple[Any, dict[str, Any], gpd.GeoDataFrame]:
    gdf = gpd.read_file(path)
    if gdf.empty:
        raise ValueError(f"{label}: no features")
    if gdf.crs is None:
        raise ValueError(f"{label}: missing CRS")
    if gdf.geometry.isna().any() or gdf.geometry.is_empty.any():
        raise ValueError(f"{label}: missing or empty geometry")
    if not gdf.geometry.is_valid.all():
        bad = int((~gdf.geometry.is_valid).sum())
        raise ValueError(f"{label}: {bad} invalid geometries")
    original_epsg = gdf.crs.to_epsg()
    projected = gdf.to_crs(epsg=target_epsg)
    geom = projected.geometry.union_all()
    if geom.is_empty or not geom.is_valid:
        raise ValueError(f"{label}: union invalid/empty")
    info = {
        "relative_path": str(SOURCES[label]).replace('\\', '/'),
        "source_epsg": original_epsg,
        "features": int(len(gdf)),
        "geometry_type": geom.geom_type,
        "area_km2": round(area_km2(geom), 9),
        "columns": [str(c) for c in gdf.columns if c != 'geometry'],
    }
    return geom, info, gdf


def geometry_metrics(aoi: Any, one: Any, two: Any, saved_combined: Any) -> dict[str, float]:
    """Compute overlap and symmetric differences independently of QGIS attributes."""
    independent_union = one.union(two)
    common = one.intersection(two)
    return {
        "aoi_km2": area_km2(aoi),
        "hrdem_1m_km2": area_km2(one),
        "hrdem_2m_km2": area_km2(two),
        "independent_union_km2": area_km2(independent_union),
        "saved_combined_km2": area_km2(saved_combined),
        "overlap_km2": area_km2(common),
        "only_1m_km2": area_km2(one.difference(two)),
        "only_2m_km2": area_km2(two.difference(one)),
        "uncovered_km2": area_km2(aoi.difference(independent_union)),
        "pct_1m": 100.0 * one.area / aoi.area,
        "pct_2m": 100.0 * two.area / aoi.area,
        "pct_combined": 100.0 * independent_union.area / aoi.area,
        "independent_union_vs_saved_symmetric_difference_km2": area_km2(independent_union.symmetric_difference(saved_combined)),
        "coverage_outside_aoi_km2": area_km2(independent_union.difference(aoi)),
        "saved_combined_outside_aoi_km2": area_km2(saved_combined.difference(aoi)),
        "union_minus_sum_identity_error_km2": abs(area_km2(independent_union) - (area_km2(one) + area_km2(two) - area_km2(common))),
    }


def check_attr(gdf: gpd.GeoDataFrame, field: str, expected: float, tolerance: float) -> tuple[bool | None, str]:
    if field not in gdf.columns:
        return None, f"Field '{field}' missing (geometry still verified independently)"
    if len(gdf) != 1:
        return None, f"Field '{field}' has {len(gdf)} records: not a single summary value"
    try:
        val = float(gdf.iloc[0][field])
    except (TypeError, ValueError):
        return False, f"Field '{field}' not numeric"
    if not math.isfinite(val):
        return False, f"Field '{field}' is not finite"
    return abs(val - expected) <= tolerance, f"Field '{field}'={val}, geometry-based={expected:.6f}"


def validate_project(root: Path) -> dict[str, Any]:
    errors: list[str] = []
    warnings: list[str] = []
    info: dict[str, Any] = {}
    for label, relative in SOURCES.items():
        if not (root / relative).is_file():
            errors.append(f"Missing input: {relative.as_posix()}")
    qgz = root / "qgis/highway99_dem_coverage.qgz"
    if qgz.is_file():
        import zipfile
        if not zipfile.is_zipfile(qgz):
            errors.append("QGIS .qgz exists but is not a valid ZIP/QGZ archive")
    else:
        warnings.append("Saved QGIS project not found at qgis/highway99_dem_coverage.qgz")
    if errors:
        return {"status": "FAILED_MISSING_FILES", "errors": errors, "warnings": warnings, "inputs": {}, "metrics": {}}

    geoms: dict[str, Any] = {}
    frames: dict[str, gpd.GeoDataFrame] = {}
    try:
        for label, relative in SOURCES.items():
            geom, metadata, gdf = coverage_geometry(root / relative, label)
            geoms[label] = geom
            frames[label] = gdf
            info[label] = metadata
    except Exception as exc:
        errors.append(str(exc))
        return {"status": "FAILED_GEOMETRY_READ", "errors": errors, "warnings": warnings, "inputs": info, "metrics": {}}

    for label in ("hrdem_1m", "hrdem_2m", "combined"):
        if info[label]["source_epsg"] != 26910:
            errors.append(f"{label}: saved layer CRS is EPSG:{info[label]['source_epsg']} not EPSG:26910")
        if info[label]["features"] != 1:
            errors.append(f"{label}: expected one dissolved feature, got {info[label]['features']}")

    try:
        m = geometry_metrics(geoms["aoi"], geoms["hrdem_1m"], geoms["hrdem_2m"], geoms["combined"])
    except Exception as exc:
        errors.append(f"Overlay calculation failed: {exc}")
        return {"status": "FAILED_SPATIAL_OVERLAY", "errors": errors, "warnings": warnings, "inputs": info, "metrics": {}}

    for metric, expected in EXPECTED.items():
        actual = m["independent_union_km2"] if metric == "combined_km2" else m[metric]
        if abs(actual - expected) > AREA_TOLERANCE_KM2:
            errors.append(f"{metric}: computed {actual:.6f}, manual QGIS value {expected:.6f} (km2)")

    if m["independent_union_vs_saved_symmetric_difference_km2"] > 0.001:
        errors.append("Saved combined geometry does not match independent union (symmetric difference > 0.001 km2)")
    if m["coverage_outside_aoi_km2"] > 0.001 or m["saved_combined_outside_aoi_km2"] > 0.001:
        errors.append("Coverage contains unexpected area outside the 5km study AOI (>0.001 km2)")
    if m["union_minus_sum_identity_error_km2"] > 0.000001:
        errors.append("Area identity (union=1m+2m-overlap) not satisfied")
    if m["independent_union_km2"] > m["aoi_km2"] + 0.001:
        errors.append("Unique HRDEM coverage exceeds AOI area")

    attribute_checks = {
        "hrdem_1m_pct": check_attr(frames["hrdem_1m"], "coverage_pct", m["pct_1m"], 0.02),
        "hrdem_2m_pct": check_attr(frames["hrdem_2m"], "coverage_pct", m["pct_2m"], 0.02),
        "hrdem_2m_area": check_attr(frames["hrdem_2m"], "coverage_km2", m["hrdem_2m_km2"], 0.001),
        "combined_area": check_attr(frames["combined"], "combined_km2", m["independent_union_km2"], 0.001),
        "combined_pct": check_attr(frames["combined"], "combined_pct", m["pct_combined"], 0.02),
    }
    for key, (passed, detail) in attribute_checks.items():
        if passed is None:
            warnings.append(key + ': ' + detail)
        elif not passed:
            errors.append(key + ': ' + detail)

    if "coverage_km2" not in frames["hrdem_1m"].columns:
        warnings.append("1m summary has no exact 'coverage_km2' field (likely a naming typo in QGIS); numeric geometry is validated separately")
    if abs(m["overlap_km2"]) < 0.000001:
        warnings.append("1m/2m polygons show zero overlap at reported precision; visually similar mosaics may be adjacent, but this applies to vector coverage polygons only")

    return {
        "status": "PASSED_INTERNAL_VECTOR_COVERAGE_QA_WITH_CAVEATS" if not errors else "FAILED_INTERNAL_VECTOR_COVERAGE_QA",
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": info,
        "metrics": {key: round(val, 9) for key, val in m.items()},
        "attribute_checks": {key: {"passed": passed, "detail": detail} for key, (passed, detail) in attribute_checks.items()},
        "errors": errors,
        "warnings": warnings,
        "scientific_scope": "STAC-derived DTM vector coverage only; not a DEM raster NoData, vertical accuracy, landslide hazard, or official roadway validation",
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    default_root = Path(__file__).resolve().parents[2]
    p.add_argument("--root", type=Path, default=default_root, help="Path to natcat-risk repository")
    p.add_argument("--out", type=Path, default=None, help="Output JSON path")
    args = p.parse_args()
    root = args.root.resolve()
    report = validate_project(root)
    report["project_root_name"] = root.name
    target = args.out or root / "results/diagnostics/step4_1_dem_coverage_final_qa.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print("STEP 4.1 — INDEPENDENT DEM VECTOR COVERAGE QA")
    print("Status:", report["status"])
    if report.get("metrics"):
        m = report["metrics"]
        for name in ["aoi_km2", "hrdem_1m_km2", "hrdem_2m_km2", "overlap_km2", "independent_union_km2", "pct_combined", "uncovered_km2", "independent_union_vs_saved_symmetric_difference_km2", "coverage_outside_aoi_km2"]:
            print(f"  {name}: {m[name]:.6f}")
    for msg in report["errors"]:
        print("ERROR:", msg)
    for msg in report["warnings"]:
        print("NOTE:", msg)
    print("Report:", target)
    print("Scope: polygon coverage, NOT raster pixel QA or model validation")
    return 0 if not report["errors"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
