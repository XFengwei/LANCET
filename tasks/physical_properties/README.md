# Physical properties

This task calculates background surface densities, core optical depth,
opacity-corrected surface density and mass, and mass/number volume densities
for every source in the generated source catalogue.

All physical constants and uncertainty assumptions are declared in
`task.toml`. Monte Carlo propagation uses 10,000 samples and a fixed seed so
the table is deterministic.

```bash
conda run -n py312 python scripts/physical_properties.py --stage inspect
conda run -n py312 python scripts/physical_properties.py --stage build
```

Machine-readable and publication tables are written under
`results/physical_properties/tables/`.
