# Physical properties

This task calculates background surface densities, core optical depth,
opacity-corrected surface density and mass, and mass/number volume densities
for every source in the generated source catalogue.

All physical constants and uncertainty assumptions are declared in
`task.toml`. Monte Carlo propagation uses 10,000 samples and a fixed seed so
the table is deterministic.

Input: `results/source_catalogue/tables/source_catalogue.ecsv`. It uses
`Tenv`, `Ibkg1` (intercore background → `Surface_intercorebkg`), `Ibkg2`
(12m+7m+TP background → `Surface_corebkg` and the background in the core
optical depth), the 12m+7m peak `FXP_BST01` and total flux `FXT_BST01` with
their errors (gain-corrected in `source_catalogue` where gain < 1), and the
deconvolved sizes.

```bash
PYTHONPATH=src conda run -n py312 python tasks/core_physical_properties/run.py --stage inspect
PYTHONPATH=src conda run -n py312 python tasks/core_physical_properties/run.py --stage build --overwrite
```

Machine-readable and publication tables are written under
`results/core_physical_properties/tables/`.

The catalogue preserves `Resolvability` for every source. N marks upper size
limits and lower volume-density limits at the adopted mass. The publication
table marks these limits explicitly. `core_property_summary.ecsv` applies
the per-quantity policy in `[statistics]`: size and volume-density summaries
use the 209 resolved cores. Peak-derived, beam-averaged surface density and
optical depth are independent of the size limit and retain all 231 sources,
as do the integrated-mass estimates. Environmental quantities retain all
valid measurements. Numerical limit values remain in the catalogue.

## Change from Paper II (2026-11-20)

The code reproduces the Paper II tables exactly from the Paper II catalogue
(all columns equal; `.tex` and summary byte-identical). With the reviewed
`source_catalogue` (new 12m+7m+TP background, gain correction):

- `Surface_corebkg`: median 2.73 → 1.39 g cm⁻² (new/old 0.22–1.20).
- Gain = 1 sources (212): mass and volume density unchanged to 1e-3; optical
  depth and core surface density × 0.985–1.000 (smaller background).
- Gain < 1 sources (19): mass, Σ_core, τ and volume density × 1.010–2.17;
  `Surface_intercorebkg` × 1.016–2.16.
- Total core mass 545.7 → 552.6 M☉ (Young 41.9 → 43.6, Intermediate
  255.9 → 261.0, Evolved 247.9 → 248.0); median mass 1.40 M☉, unchanged.

## Temperature-prior diagnostic

The temperature-prior sampling of the core masses (Paper II Fig. 9) moved to
`core_formation_efficiency` (stage `temperature_sampling`, settings in its
`[temperature_sampling]`). The endpoint stress test and its plots
(`scripts/core_temperature_stress_test.py`,
`scripts/plot_core_temperature_*.py`, `scripts/report_core_temperature_tests.py`)
read the upper temperatures from there.
