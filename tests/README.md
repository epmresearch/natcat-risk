# Testing approach

Initial check (after Step 2 installation):

```bash
python src/check_environment.py
python -m unittest climada.engine.test.test_impact
```

Later test categories:
- CRS consistency and geometry/topology checks
- Hazard frequency and raster/event alignment checks
- Loss non-negativity and unit consistency checks
- Network closure connectivity and route invariants
- Direct/indirect component reconciliation and double-counting checks
- Economic model scenario regression tests
