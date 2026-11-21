"""Small read-only package availability check; does not simulate hazards or losses."""

from importlib import metadata, util
import sys

PACKAGES = {
    "climada": "climada",
    "climada-petals": "climada_petals",
    "geopandas": "geopandas",
    "rasterio": "rasterio",
    "networkx": "networkx",
    "osmnx": "osmnx",
}


def main() -> int:
    print(f"Python: {sys.version.split()[0]}")
    missing = []
    for package, module in PACKAGES.items():
        found = util.find_spec(module) is not None
        try:
            version = metadata.version(package)
        except metadata.PackageNotFoundError:
            version = "unknown/not installed"
        print(f"{package}: {'available' if found else 'MISSING'} (version: {version})")
        if not found:
            missing.append(package)

    if missing:
        print("Missing packages:", ", ".join(missing))
        return 1
    print("Package presence check passed. Run the official CLIMADA unit test next.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
