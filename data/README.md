# Dataset organization and responsible storage

Only `data/README.md` and `data/metadata/*` are intended to be tracked in Git by default.

**Local working folders** (ignored by Git):
- `data/raw/`: original downloaded public source datasets (unchanged)
- `data/external/`: third-party datasets with restrictions, if permissions allow local storage
- `data/interim/`: reprojected/clipped/pre-cleaned sources
- `data/processed/`: analysis-ready road segments, rasters, event inventory and tabular data

Fill out `data/metadata/sources.csv` for every incoming dataset, including official URL, license, access date, CRS and limitations. Never upload restricted or personal datasets merely because a repository is private.

For large data, keep an authorized storage location and provide a reproducible download/process script whenever licensing permits.
