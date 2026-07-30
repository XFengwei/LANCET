# Environment measurements

This task measures the local environment of every generated core inside a
15-arcsec circular aperture. For each map it records the mean, median,
standard deviation, median absolute deviation, and valid-pixel fraction.
Non-finite and non-positive pixels are excluded.

The inputs are the generated physical-properties catalogue, the delivered
dust-temperature result, the generated C18O moment-2 map, and the generated
surface-density products. Using the current intercore surface-density map
keeps the environment measurements consistent with the opacity adopted by
the physical-properties and column-density tasks. Its measurements are 1.8
times the retained older values because the mass scale changes by 0.9/0.5.

```bash
conda run -n py312 python scripts/environment_measurements.py --stage inspect
conda run -n py312 python scripts/environment_measurements.py --stage build
```

The augmented machine-readable catalogue is written to
`results/environment_measurements/tables/core_environment_catalogue.ecsv`.
It is the input catalogue for subsequent correlation analyses.
