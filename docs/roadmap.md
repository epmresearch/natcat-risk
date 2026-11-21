# End-to-end research roadmap (18 steps)

Track each step with a small commit / pull request (if team workflow permits). Some writing and validation proceed in parallel, but core dependency ordering should be preserved.

| Step | Work | Version-controlled deliverables |
| ---: | --- | --- |
| 01 | Research definition / questions | `docs/project_scope.md`, `docs/decisions.md` |
| 02 | Environment / Git workflow | `environment.yml`, `src/check_environment.py`, `docs/git_workflow.md` |
| 03 | Study area / road segmentation | `config/study_area.yml`, reproducible GIS scripts, figure |
| 04 | Data inventory and access | `data/metadata/sources.csv`, licensing notes |
| 05 | Cleaning, CRS and preprocessing | Scripts / notebooks, QA checks |
| 06 | Landslide hazard modelling | Model methodology, reproducible event generation |
| 07 | CLIMADA hazard encoding | Hazard pipeline / tests |
| 08 | Exposure model | Highway asset model / valuation assumptions |
| 09 | Impact / vulnerability functions | Functions, references, validations |
| 10 | Direct losses | Event losses / segment maps |
| 11 | Closure/recovery response | Closure duration / probability model |
| 12 | Road network disruption | Routing/detour code and travel metrics |
| 13 | Indirect economic losses | Cost model and non-double-counting checks |
| 14 | Probabilistic aggregation | AAL if validated; uncertainty / event curves |
| 15 | Resilience scenarios | Mitigation assumptions and simulations |
| 16 | CBA, validation and sensitivity | Benefit-cost / uncertainty results |
| 17 | Reproducible final outputs | Figures, tables, scripts, run instructions |
| 18 | Thesis writing / defense | Final draft, references, slide materials as permitted |

**Each step is 'done' only when:** code/docs reviewed, input sources documented, run steps reproducible, outputs examined, and relevant changes committed to the university remote repository.

**Important:** The current files are setup drafts, not completed analyses. The supervisor must approve scientific choices and the remote repository's workflow.
