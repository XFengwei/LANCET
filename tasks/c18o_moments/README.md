# C18O moments

This task reconstructs the C18O cube derivatives and moment products from the
single delivered SEDIGISM cube, its only input.

## Stages

- `inspect` reports the delivered cube and the mask settings.
- `prepare` creates a Jy/beam cube with a km/s spectral axis and a
  brightness-temperature cube in K. Both use float32 storage.
- `moments` applies the declared 9-sigma to 2-sigma dilated mask and writes
  moments, uncertainties, noise, mask, and flux products.
- `plot` creates the stacked moment-map PDF and PNG.
- `all` runs every stage in dependency order.

All outputs are written under `results/c18o_moments/`. Nothing is written to
`data/`.

```bash
PYTHONPATH=src conda run -n py312 python tasks/c18o_moments/run.py --stage inspect
PYTHONPATH=src conda run -n py312 python tasks/c18o_moments/run.py --stage all
```

Pass `--overwrite` only when intentionally replacing task products.

## Readers

`environment_measurements` (`moments/c18o.mom2.fits.gz` → `Vsig`).

## Change (2026-11-21)

Outputs moved from `results/moment_maps/` to `results/c18o_moments/` (D5).
Rerun: every cube, moment, mask and flux table is data-identical to the
previous products. The `reproject` stage (moment 2 on the 0.1″ ALMA grid)
was removed with its `alma_grid` input and `lancet.moments`, because nothing
read it (author); `environment_measurements` uses the native moment-2 map.
