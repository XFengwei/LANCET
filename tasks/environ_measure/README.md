# Environment measurements

This task measures the local environment of every core inside a 15″
circular aperture. For each map it records the mean, median, standard
deviation, median absolute deviation and valid-pixel fraction. Non-finite
and non-positive pixels are excluded (`minimum_valid_value = 0`).

## Inputs

| Key | File | Column prefix |
|---|---|---|
| `physical_catalogue` | `results/core_physical_properties/tables/core_physical_properties.ecsv` | (copied into the output) |
| `temperature` | `data/delivered/sedfitter/combine/Tdust_18arcsec+36arcsec.fits` | `Temp` |
| `velocity_dispersion` | `results/c18o_moments/moments/c18o.mom2.fits.gz` | `Vsig` |
| `intercore_surface_density` | `results/column_density/maps/surface_density_intercore.fits` | `Interco` |
| `surface_density` | `results/column_density/prepared/surface_density_lowres.fits` | `Surface` |

The C18O moment-2 map is a product of the `c18o_moments` task.

```bash
PYTHONPATH=src conda run -n py312 python tasks/environ_measure/run.py --stage inspect
PYTHONPATH=src conda run -n py312 python tasks/environ_measure/run.py --stage build --overwrite
```

Output: `results/environ_measure/tables/core_environment_catalogue.ecsv`,
the physical-property catalogue plus the `*_env_*` columns. It is the input
of `environ_core_correlation` and `bonnor_ebert`.

## Change from Paper II (2026-11-20)

- Before the gain correction below, all `*_env_*` measurements were
  unchanged: `Temp`, `Vsig` and `Interco` byte-identical, `Surface` within 1e-6 (float32 rounding of the regenerated
  low-resolution map). The intercore map now has NaN where the intensity is
  negative; those pixels were already excluded as non-positive.
- The copied physical-property columns follow the reviewed
  `core_physical_properties` (new 12m+7m+TP background, gain-corrected cores).

## Gain correction of the intercore map (2026-11-20)

The intercore surface-density map is built from the getsf background of the
12m+7m clean image (sky × gain). Since 2026-11-20 `column_density` divides it
by the gain where gain < 1, matching the per-core `Ibkg1` of
`source_catalogue` (author). The 15″ apertures of 108 cores include
gain < 1 pixels (more than 25 % for 43). `Interco_env_mean` rises for 107
cores (median of those × 1.018, at most × 1.448) and `Interco_env_median` for
101 (median × 1.036, at most × 1.342); valid fractions and all other columns
are unchanged.
