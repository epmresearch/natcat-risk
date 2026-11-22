# Development Environment Setup

## 1. Purpose

This document describes the validated local Python environment for the Sea-to-Sky Highway landslide-risk thesis project. The environment is isolated from the system Python installation using Miniforge/Mamba.

## 2. Environment details

- **Operating system used for validation:** Windows
- **Environment manager:** Miniforge / Conda / Mamba
- **Conda environment name:** `sea_to_sky_risk`
- **Python version:** `3.11.18`
- **Validation date:** 2026-11-22

### Verified package versions

| Package | Verified version |
| --- | --- |
| `climada` | `6.1.0` |
| `climada-petals` | `6.2.0` |
| `geopandas` | `1.1.4` |
| `rasterio` | `1.4.3` |
| `networkx` | `3.6.1` |
| `osmnx` | `2.1.2` |

These versions were observed in the local verification output. The project's `environment.yml` remains the primary environment specification.

## 3. Create and activate the environment

From the root of the cloned repository, using Miniforge Prompt:

```bat
mamba env create -f environment.yml
mamba activate sea_to_sky_risk
python --version
```

Expected Python version for the validated setup: `3.11.18`.

Do not install project dependencies into the system-wide Python environment.

## 4. Run verification checks

From the repository root, with `sea_to_sky_risk` activated:

```bat
python src/check_environment.py
python -m unittest climada.engine.test.test_impact
python -c "from climada_petals.hazard.landslide import Landslide; print('Landslide import OK')"
```

### Observed results (2026-11-22)

- `src/check_environment.py`: **passed** (all listed packages detected).
- CLIMADA impact unit tests: **44 tests run, OK**.
- CLIMADA Petals Landslide import: **passed** (`Landslide import OK`).
- Some non-fatal deprecation and user warnings appeared during the CLIMADA tests.

These tests verify installation and selected package functionality; they **do not** validate the thesis landslide hazard model or risk calculations.

## 5. VS Code and Jupyter

1. Open the cloned repository folder in VS Code.
2. Install the Microsoft **Python** and **Jupyter** extensions.
3. Select the Python interpreter from the `sea_to_sky_risk` Conda environment.
4. For notebooks, select the kernel corresponding to the same environment.

A Jupyter cell successfully imported `climada`, `climada_petals`, and `Landslide`, and confirmed that the active Python executable belonged to the `sea_to_sky_risk` Conda environment.

## 6. Repository and data hygiene

- Commit configuration (`environment.yml`), source code, tests, and documentation.
- Do **not** commit the Conda environment directory, caches, credentials, or private datasets.
- Store large raw geospatial datasets outside ordinary Git history, and document their sources and preprocessing steps.
- When dependencies change, update `environment.yml`, rerun the checks, and record any relevant change here.

## 7. Next step

Proceed to defining and preparing the Sea-to-Sky Highway 99 study-area geospatial data.
