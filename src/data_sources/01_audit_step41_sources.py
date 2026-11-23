"""Offline Step 4.1 data-source registry and study area AOI audit.

This script does not download, authenticate to, or validate spatial coverage of
any third-party dataset. It only validates the catalog and local study geometry.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlparse

import geopandas as gpd

ROOT = Path(__file__).resolve().parents[2]
EXPECTED_FIELDS = [
    "source_id", "domain", "dataset_name", "publisher", "priority", "intended_use",
    "native_resolution_or_scale", "stated_geographic_coverage",
    "highway99_coverage_verified", "access", "licence_status",
    "licence_or_rights", "dataset_url", "verification_date_utc", "limitations_or_action",
]
PRIORITIES = {"P0", "P1", "P2", "REFERENCE"}


def read_and_validate_catalog(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != EXPECTED_FIELDS:
            raise ValueError("Unexpected catalog columns or column order")
        records = list(reader)
    if not records:
        raise ValueError("Source registry is empty")
    ids: set[str] = set()
    for index, record in enumerate(records, start=2):
        if any(not (record.get(field) or "").strip() for field in EXPECTED_FIELDS):
            raise ValueError(f"Missing required field on CSV row {index}")
        identifier = record["source_id"]
        if identifier in ids:
            raise ValueError(f"Duplicate source_id: {identifier}")
        ids.add(identifier)
        if record["priority"] not in PRIORITIES:
            raise ValueError(f"Invalid priority: {record['priority']}")
        if record["highway99_coverage_verified"] not in {"NO", "YES"}:
            raise ValueError("Invalid coverage flag")
        parsed = urlparse(record["dataset_url"])
        if parsed.scheme != "https" or not parsed.netloc:
            raise ValueError(f"Invalid source URL: {record['dataset_url']}")
    return records


def audit_aoi(aoi_file: Path, expected_quality_file: Path | None = None) -> dict:
    frame = gpd.read_file(aoi_file)
    if frame.empty or frame.crs is None:
        raise ValueError("AOI is empty or its CRS is missing")
    if not frame.geometry.is_valid.all() or frame.geometry.is_empty.any():
        raise ValueError("AOI contains invalid or empty geometries")
    if not frame.geometry.geom_type.isin(["Polygon", "MultiPolygon"]).all():
        raise ValueError("AOI must contain polygonal geometry")
    projected = frame.to_crs("EPSG:26910")
    merged = projected.geometry.union_all()
    area_km2 = float(merged.area / 1_000_000)
    if area_km2 <= 0:
        raise ValueError("AOI has non-positive area")
    bounds = [round(float(v), 7) for v in frame.to_crs("EPSG:4326").total_bounds]
    result = {
        "aoi_file": str(aoi_file),
        "aoi_sha256": hashlib.sha256(aoi_file.read_bytes()).hexdigest(),
        "aoi_input_crs": str(frame.crs),
        "aoi_features": len(frame),
        "aoi_bounds_wgs84_west_south_east_north": bounds,
        "aoi_unique_area_km2_metric_epsg26910": round(area_km2, 6),
        "step34_area_comparison": None,
    }
    if expected_quality_file is not None and expected_quality_file.exists():
        qa = json.loads(expected_quality_file.read_text(encoding="utf-8"))
        expected = qa.get("unique_terrain_data_search_area_km2")
        if expected is None:
            raise ValueError("Missing terrain area in Step 3.4 QA report")
        deviation = abs(area_km2 - float(expected))
        result["step34_area_comparison"] = {
            "reported_km2": expected,
            "difference_km2": round(deviation, 6),
            "matches_0_01_km2_tolerance": deviation <= 0.01,
        }
        if deviation > 0.01:
            raise ValueError("AOI area is inconsistent with Step 3.4 QA")
    return result


def produce_audit(catalog: Path, aoi: Path, qa: Path | None) -> dict:
    items = read_and_validate_catalog(catalog)
    geometry = audit_aoi(aoi, qa)
    return {
        "step": "4.1",
        "status": "SOURCE_REGISTRY_AND_AOI_QA_PASSED_COVERAGE_NOT_VERIFIED",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "catalog_file": str(catalog),
        "catalog_sha256": hashlib.sha256(catalog.read_bytes()).hexdigest(),
        "source_count": len(items),
        "sources_by_priority": dict(sorted(Counter(item["priority"] for item in items).items())),
        "sources_with_local_corridor_coverage_verified": sum(
            item["highway99_coverage_verified"] == "YES" for item in items
        ),
        "aoi": geometry,
        "required_next_steps": [
            "Inspect NRCan HRDEM footprints and MRDEM data access within the exact AOI",
            "Obtain coverage and local event counts for the landslide inventory",
            "Check BC geological/TSM/DRA extent and legal reuse terms",
            "Download a clipped MRDEM DTM in Step 4.2 only after source and coverage QA",
        ],
        "scientific_boundary": (
            "This is an OFFLINE catalog/AOI consistency check. No third-party dataset "
            "was downloaded, inspected, or certified for Highway 99 coverage. "
            "The 5 km terrain search AOI is not a landslide hazard or runout zone."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=ROOT / "data/catalog/step4_1_source_registry.csv")
    parser.add_argument("--aoi", type=Path, default=ROOT / "data/processed/study_area/highway99_terrain_search_5000m.geojson")
    parser.add_argument("--step34-qa", type=Path, default=ROOT / "data/processed/study_area/highway99_spatial_preparation_quality.json")
    parser.add_argument("--output", type=Path, default=ROOT / "results/diagnostics/step4_1_source_discovery_audit.json")
    args = parser.parse_args()
    audit = produce_audit(args.catalog, args.aoi, args.step34_qa)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("STEP 4.1 — OFFLINE SOURCE REGISTRY AND AOI AUDIT: PASSED")
    print(f"Sources documented: {audit['source_count']}")
    print(f"AOI area (km2): {audit['aoi']['aoi_unique_area_km2_metric_epsg26910']}")
    print(f"Extent WGS84: {audit['aoi']['aoi_bounds_wgs84_west_south_east_north']}")
    print("NOTE: Actual source coverage has NOT been spatially verified or downloaded.")
    print("OUTPUT:", args.output)


if __name__ == "__main__":
    main()
