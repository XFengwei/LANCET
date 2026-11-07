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

The catalogue preserves `Resolvability` for every source. N marks upper size
limits and lower volume-density limits at the adopted mass. The publication
table marks these limits explicitly. `core_property_summary.ecsv` applies
the per-quantity policy in `[statistics]`: size and volume-density summaries
use the 209 resolved cores. Peak-derived, beam-averaged surface density and
optical depth are independent of the size limit and retain all 231 sources,
as do the integrated-mass estimates. Environmental quantities retain all
valid measurements. Numerical limit values remain in the catalogue.

## Temperature-prior diagnostic

The optional `[temperature_sampling]` configuration defines the Paper II
temperature stress test without replacing the adopted physical catalogue.
Run `PYTHONPATH=src conda run -n py312 python scripts/core_temperature_sampling.py`.
It compares uniform temperatures with the author-approved peak-brightness
beta prior, with independent core draws and a shared-regional-percentile
dependence check. Each draw recomputes optical depth with the environment
kept at Tenv. The isolated output directory includes per-core priors,
regional mass draws and conditional temperature-only intervals; these are
not a data-constrained temperature posterior or total mass uncertainties.
Use `scripts/plot_core_temperature_sampling.py --qa-scripts /path/to/nature-figure/scripts`
to replot the mass hierarchy for both temperature models.
Use `--trials` for the two Figure 9 integration layouts. The current regional
upper temperatures are 30/40/50 K (Young/Intermediate/Evolved); the endpoint
runner reads the same configuration. After rerunning both numerical runners,
`scripts/report_core_temperature_tests.py` refreshes the numerical reports
and trial captions from the saved products.
