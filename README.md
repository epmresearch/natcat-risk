# Sea-to-Sky Highway Landslide Risk and Resilience Study

**Status:** Research starter / Step 1 draft (not yet approved by supervisor)  
**Case study (provisional):** Highway 99, Vancouver–Whistler corridor, British Columbia, Canada  
**Software:** CLIMADA Core + CLIMADA Petals, QGIS, Python  
**Research focus:** Probabilistic landslide risk, direct infrastructure damage, road-closure impacts, indirect economic losses, and mitigation cost–benefit assessment.

## Research question
How can probabilistic landslide scenarios, infrastructure vulnerability, and road-network disruption be combined to estimate the **direct and indirect economic consequences** of closures along the Sea-to-Sky Highway, and to compare risk-reduction measures?

This repository contains the **study's own code and reproducible workflow**. It does not contain a copy of the CLIMADA source code. Do not upload raw geospatial datasets, credentials, personal records, or third-party data without checking permission and licensing.

## Start here
1. Read [`docs/project_scope.md`](docs/project_scope.md) (draft scope, questions, and open decisions).
2. Read [`docs/roadmap.md`](docs/roadmap.md) (18 steps and deliverables).
3. Read [`docs/git_workflow.md`](docs/git_workflow.md) or [`docs/git_upload_guide_fa.md`](docs/git_upload_guide_fa.md) to use your **existing university repository**.
4. Install [Miniforge / Mamba](https://github.com/conda-forge/miniforge) and consult the [official CLIMADA install guide](https://climada-python.readthedocs.io/en/stable/getting-started/install.html).
5. When ready for Step 2, try `mamba env create -f environment.yml` then `mamba activate sea_to_sky_risk`. **The environment specification has not yet been tested on your computer.**
6. Run `python src/check_environment.py` and `python -m unittest climada.engine.test.test_impact` to check the installation.

## Repository layout
| Path | Contents |
| --- | --- |
| `docs/` | Research questions, methodology, decisions, workflow, bibliography |
| `config/` | Human-readable configurations; geographic study bounds pending |
| `data/metadata/` | Dataset manifest, provenance and licenses (**tracked**) |
| `data/raw/`, `data/processed/`, `data/external/` | Working datasets (**ignored by Git**) |
| `notebooks/` | Exploratory and explanatory notebooks |
| `src/` | Reusable scripts for hazard, exposure, loss and network models |
| `tests/` | Sanity checks, unit/integration tests and model diagnostics |
| `results/figures/`, `results/tables/` | Selected publication-ready outputs |
| `thesis/` | Thesis outlines, drafts, figure index and references |

## Scientific caution
- A landslide **susceptibility** map is not automatically an annual **occurrence-probability** map.
- A reported historical landslide point is not necessarily a measured event footprint or frequency.
- A hazard reaching the road, physical road damage, and traffic closure are separate processes.
- Indirect costs are scenario-dependent and must avoid double counting lost tourism spending and traveler costs.
- Calculate AAL only when defensible annual frequencies and event losses are available.
- Monetary units, price base year and assumptions must be stated explicitly before publishing any economic result.

## External software
- [CLIMADA Core](https://github.com/CLIMADA-project/climada_python)
- [CLIMADA Petals](https://github.com/CLIMADA-project/climada_petals)
- [Landslide tutorial](https://climada-petals.readthedocs.io/en/latest/tutorial/climada_hazard_Landslide.html)

## Acknowledgment / governance
Research concept developed for an MSc research project with University of Calgary supervision. This starter is a working draft and does not imply academic approval of scope, datasets or results. Follow the university repository's existing policies; do not replace existing project files without reviewing them. Repository licensing and data-sharing permissions are to be confirmed with the supervisor.
