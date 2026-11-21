# Construct free-free model

This task constructs a 1.3 mm free-free emission model from MeerKAT and
CORNISH-South radio continuum data. The model is evaluated at the ALMA
continuum reference frequency, `output_frequency_ghz = 223.8` in `[fit]`
(file names `ff224_map*.fits`, FITS card `FREQ`).

## Boundaries

The task reads the three declared FITS inputs under
`data/delivered/archive/radio/` and never modifies them. Recreated CORNISH
products, fitted maps, model cubes, figures, and run metadata belong under
`results/construct_freefree_model/`.

Reusable physical calculations and pixel fitting live in
`src/lancet/freefree.py`. The task is run through
`tasks/construct_freefree_model/run.py`; no Python `Task` object is required.

## Stages

1. `inspect` validates the inputs and reports exact MeerKAT frequencies from
   FITS cards `FREQ0001` through `FREQ0014`.
2. `prepare` smooths the 2.5-arcsec CORNISH cutout to 8 arcsec and reprojects
   it onto the MeerKAT grid.
3. `fit` constructs the 15-band SED cube and fits the free-free model at each
   configured electron temperature.
4. `postprocess` selects the configured reference temperature, smooths its
   224 GHz map from 8 to 18 arcsec, and reprojects it to the delivered dust
   grid. This creates the free-free input required by continuum synthesis.
5. `plot` compares the fitted 6000 K and 10000 K maps with the 8000 K
   reference using the absolute fractional 224 GHz flux difference.

Temperature comparison and additional pixel diagnostics can consume the
fitted task outputs independently.

The temperature-comparison plot is written as PDF and PNG under
`results/construct_freefree_model/figures/` with the stable basename
`FreeFreeTemplate_Te_varying_flux_difference`. Its reference temperature,
comparison temperatures, bin count, colors, resolution, and output formats
are declared in `task.toml`.

## Output frequency

Until 2026-11-20 the map was evaluated at a hard-coded 230 GHz
(`ff230_map*.fits`, as used in Paper II). It is now a setting, adopted as
223.8 GHz. The fitted emission measures are unchanged (byte-identical rerun);
the free-free flux density is higher by exactly (223.8/230)^-0.1 = 1.002736
in every pixel of the native, smoothed and dust-grid maps. The old
`ff230_map*.fits` files were removed on 2026-11-20 after `synthesize_continuum`
and `column_density` were repointed to `ff224_map*`.

## Paper II fitting profile

The maintained fitting functions in `src/lancet/freefree.py` preserve the
verified standalone implementation, but the old prototype's example settings
were not the configuration that generated the stored Paper II products. The
legacy `fitmask_map`, `nband_map`, and sampled `EM_map` values establish that
Paper II used:

- a fixed RMS of `3e-4 Jy/beam` for each nominal MeerKAT subband;
- a fixed effective RMS of `1e-2 Jy/beam` for CORNISH;
- a per-band threshold of `S/N > 5`;
- `mfs_map > 5e-3 Jy/beam`;
- `-0.5 <= alpha <= 2.0`; and
- the rounded frequency list used for the Paper II products.

These values are now explicit in `task.toml`. The two RFI-flagged MeerKAT
planes contain only NaNs, so a fitted pixel uses either 12 MeerKAT bands or
those 12 bands plus CORNISH.

The RMS values are fixed effective uncertainties, not values calculated during
the task. In particular, the CORNISH value is deliberately much larger than
its thermal image RMS. Together with the 5-sigma cut, it restricts CORNISH to
pixels brighter than `0.05 Jy/beam`, consistent with the manuscript statement
that ATCA is used only in the central bright region because of missing flux.

Run the safe stages with:

```bash
conda run -n py312 python tasks/construct_freefree_model/run.py --stage inspect
conda run -n py312 python tasks/construct_freefree_model/run.py --stage prepare
conda run -n py312 python tasks/construct_freefree_model/run.py --stage fit
conda run -n py312 python tasks/construct_freefree_model/run.py --stage postprocess
conda run -n py312 python tasks/construct_freefree_model/run.py --stage plot
```

## Core contamination

The estimate of smooth radio emission within the core footprints (Paper II
referee item R7) moved on 2026-11-20 to its own task,
`tasks/freefree_core_contamination/`, because it depends on later tasks
(`source_catalogue`, `core_physical_properties`, `core_formation_efficiency`).
