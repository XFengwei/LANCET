# Column density

This task converts the all-scale 1.3 mm dust emission and the extracted
intercore background into H2 column-density and gas surface-density maps.

The configured finite-opacity inversion uses 223.8 GHz, a total-gas opacity
of 0.005 cm2/g, and mean molecular weight 2.8. These values are explicit task
settings so later analyses use the same physical assumptions.

## Stages

- `inspect` reports input readiness and model parameters.
- `prepare` regrids the lower-resolution dust temperature and converts the
  retained lower-resolution column-density maps to surface density.
- `maps` creates all-scale and intercore column- and surface-density maps.
- `plot` creates the full-field, free-free-overlay, region-box, and three
  subregion surface-density figures. Full-field views include the adopted
  12m+7m primary-beam response contour at 0.2.
- `dense` measures mass above 1 g/cm2 and its fraction of the adopted cloud
  mass in each evolutionary subregion.
- `all` runs the stages in dependency order.

All generated FITS maps use float32 storage. All outputs remain under
`results/column_density/`; the input data and external extraction products are
read-only. The adopted temperature and lower-resolution column-density maps
are declared under `data/delivered/sedfitter/combine/`.

```bash
conda run -n py312 python scripts/column_density.py --stage inspect
conda run -n py312 python scripts/column_density.py --stage all
```

Pass `--overwrite` only when intentionally replacing task products.
