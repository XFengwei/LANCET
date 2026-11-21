# Tasks

Each folder is one self-contained analysis: `task.toml` (inputs, outputs,
settings), `run.py` (the command) and `README.md` (method, outputs, status).
A task writes only to `results/<task>/`. Run a stage from the repository root:

```bash
PYTHONPATH=src conda run -n py312 python tasks/<task>/run.py --stage <name> --overwrite
```

## Run order

"Needs" lists the tasks whose `results/` a task reads (from its `task.toml`).
Tasks with the same level do not depend on each other and can run in any
order; a task can run once every task it needs has run. A fresh clone runs
them by increasing level.

| Level | Task | Needs | Purpose |
|---|---|---|---|
| 0 | `background_effects` | — | Effect of foreground/background dust on a core's optical depth (model only, no file inputs) |
| 0 | `c18o_moments` | — | C18O cube derivatives and moment maps from the SEDIGISM cube |
| 0 | `construct_freefree_model` | — | 1.3 mm free-free model from MeerKAT and CORNISH-South |
| 0 | `gas_to_dust_ratio` | — | Gas-to-dust ratio at the filament's Galactocentric radius |
| 0 | `planck_polarization` | — | Planck 353 GHz polarization products and field-orientation figure |
| 0 | `vvv_mosaic` | — | VVV J/H/Ks mosaic and colour composite |
| 1 | `synthesize_continuum` | `construct_freefree_model` | 1.3 mm total-power model (dust + free-free) and its SIMBA check |
| 2 | `column_density` | `construct_freefree_model`, `synthesize_continuum` | Column-density and surface-density maps, 18″ dust temperature (Fig. 4) |
| 2 | `multiwavelength_strip` | `synthesize_continuum` | Five-panel multiwavelength strip (Fig. 1) |
| 2 | `plot_allscale_continuum` | `synthesize_continuum` | Views of the 12m+7m+TP continuum (Fig. 3) |
| 3 | `source_catalogue` | `column_density` | Source catalogue from the getsf extraction, gain-corrected, with backgrounds |
| 3 | `volume_density_mapper` | `column_density` | Volume-density map from the all-scale surface density |
| 4 | `combine_continuum` | `source_catalogue`, `synthesize_continuum` | uv-plane combination of 12m+7m and TP (J-comb, feather) and its validation (Fig. 2) |
| 4 | `core_physical_properties` | `source_catalogue` | Core optical depth, surface density, mass and volume density |
| 5 | `core_mass_function` | `core_physical_properties` | Global and subregion core mass functions |
| 5 | `environ_measure` | `c18o_moments`, `column_density`, `core_physical_properties` | Environment of each core in a 15″ aperture |
| 6 | `bonnor_ebert` | `environ_measure` | Bonnor–Ebert stability of the cores |
| 6 | `core_formation_efficiency` | `core_physical_properties`, `environ_measure`, `source_catalogue` | Core and fragment formation efficiency; temperature sampling |
| 6 | `environ_core_correlation` | `environ_measure` | Correlations between core and environment properties |
| 7 | `freefree_core_contamination` | `core_formation_efficiency`, `core_physical_properties`, `source_catalogue` | Free-free contamination of core and fragment masses |

`tests/test_task_layout.py` checks this table against the `task.toml` files:
every task is listed once, "Needs" matches its inputs, and the level is one
above its deepest dependency.

## Delivered products made by a task

`combine_continuum` makes the combined 12m+7m+TP continuum image. Its adopted
product was promoted to `data/delivered/alma/continuum/`
(`IRDC316.8.TP7MTM2TM1.robust0.5.cont*.fits`), and `column_density`,
`source_catalogue`, `multiwavelength_strip` and `plot_allscale_continuum` read
that delivered copy, not `results/combine_continuum/`. Rerunning
`combine_continuum` changes those tasks only after a new product is promoted.

`multiwavelength_strip` also reads the mosaic settings in
`tasks/vvv_mosaic/task.toml`; it builds its own near-infrared panel and does not
need the `vvv_mosaic` results.

## Local-only coordination

`paper_ii_referee_response/` (gitignored) tracks the Paper II referee response.
It runs no analysis and is not part of the run order.
