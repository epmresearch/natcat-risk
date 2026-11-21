# Reusable research code

The package will grow in phases:
- `src/hazard/` — inventory, predictors, hazard event sets
- `src/exposure/` — road segmentation and asset valuation
- `src/impact/` — vulnerability and direct damage
- `src/transport/` — closures, alternative routes, time/trip disruption
- `src/economics/` — indirect economic costs, anti-double-counting
- `src/mitigation/` — strategy scenarios and cost-benefit analysis
- `src/common/` — reproducible paths, CRS and metadata helpers

Modules and implementations will be added when their datasets, scientific assumptions and interfaces have been chosen. The current only executable script is `check_environment.py`.
