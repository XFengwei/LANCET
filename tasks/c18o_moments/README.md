# C18O moments

This task reconstructs the C18O cube derivatives and moment products from the
single delivered SEDIGISM cube. The delivered ALMA continuum is used only as a
target WCS for the velocity-dispersion reprojection.

## Stages

- `inspect` reports the delivered cube, ALMA grid, and mask settings.
- `prepare` creates a Jy/beam cube with a km/s spectral axis and a
  brightness-temperature cube in K. Both use float32 storage.
- `moments` applies the declared 9-sigma to 2-sigma dilated mask and writes
  moments, uncertainties, noise, mask, and flux products.
- `reproject` maps moment-2 and its footprint to the ALMA grid as float32.
- `plot` creates the stacked moment-map PDF and PNG.
- `all` runs every stage in dependency order.

All outputs are written under `results/moment_maps/`. Nothing is written to
`data/`.

```bash
conda run -n py312 python scripts/c18o_moments.py --stage inspect
conda run -n py312 python scripts/c18o_moments.py --stage all
```

Pass `--overwrite` only when intentionally replacing task products.
