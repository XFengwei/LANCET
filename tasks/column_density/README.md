# Column density

This task converts the all-scale 1.3 mm dust emission and the getsf intercore
background into H2 column-density and gas surface-density maps (Paper II
Fig. 4, Eq. for Σ). It also prepares the 18″ dust temperature on the ALMA
grid and the Paper I column densities as surface densities.

## Inputs

| Key | File | Notes |
|---|---|---|
| `grid_reference` | `data/delivered/alma/continuum/IRDC316.8.7MTM2TM1.robust0.5.cont.clean.fits` | 12m+7m MIRIAD image; only its 2-D grid is used |
| `dust_continuum_allscale` | `results/combine_continuum/products/jcomb_dust_224GHz_extended.fits` | J-comb 12m+7m + dust-only TP (TP factor 1.0), Jy/beam, converted to MJy/sr here. The extended product continues into the TP-only area beyond the mosaic. Not yet promoted to `data/delivered/`; only a delivered file or a `combine_continuum` product is accepted |
| `intercore_background` | `extraction/getsf/triple/results/+components/G316.sw.1302.obs.sbackground.fits` | getsf background of the 12m+7m image, MJy/sr; getsf is not rerun (D7) |
| `temperature_lowres`, `column_density_lowres`, `column_density_multiscale` | `data/delivered/sedfitter/combine/` | Paper I SED products (unscaled originals) |
| `freefree_smoothed` | `results/construct_freefree_model/Te_8000K/ff224_map_smoothed_18arcsec.fits` | contours on the overlay figure |
| `primary_beam` | `data/delivered/alma/continuum/IRDC316.8.7MTM2TM1.robust0.5.cont.pb.fits` | MIRIAD mosaic pb; its 0.2 contour is drawn |
| region files | `region/` | subregion boxes and the dense-gas regions |

## Method

Finite-opacity inversion at 223.8 GHz with a total-gas opacity of
0.005 cm² g⁻¹ (dust κ = 0.5 cm² g⁻¹, R_gd ≈ 100) and μ_H2 = 2.8. The
temperature is the Paper I 18″+36″ dust temperature, regridded to the ALMA
grid. Pixels with intensity ≤ 0 get column density 0
(`zero_nonpositive_intensity`); NaN stays NaN. All settings are in `task.toml`.

## Stages

- `inspect`: input readiness and model parameters; writes nothing.
- `prepare`: regrids the temperature and converts the lower-resolution column
  densities to surface density (`prepared/`).
- `maps`: all-scale and intercore column- and surface-density maps (`maps/`).
- `plot`: full-field, free-free overlay, subregion boxes and three subregion
  figures, PDF and PNG (`figures/`).
- `dense`: mass above 1 g cm⁻² and its fraction of the configured cloud masses
  in the three evolutionary regions (`measurements/dense_gas.json`). No other
  task reads it.
- `all`: all stages in order.

```bash
PYTHONPATH=src conda run -n py312 python tasks/column_density/run.py --stage inspect
PYTHONPATH=src conda run -n py312 python tasks/column_density/run.py --stage all --overwrite
```

Existing products are kept unless `--overwrite` is given. All FITS products
are float32.

## Readers

`volume_density_mapper` (all-scale and low-resolution surface density),
`environment_measurements` (intercore and low-resolution surface density),
`source_catalogue` (temperature on the ALMA grid).

## Change from Paper II (2026-11-20)

Paper II used the delivered `allscale_combined_v3_dust` map. Its TP
contribution is about 3.1 times the delivered TP model (v3 / new = 3.15 in the
TP-only area and 3.09 where pb ≥ 0.5). The new input uses the TP as delivered
(factor 1.0), so the all-scale surface density is about 3 times lower in
diffuse gas and close to unchanged at the brightest cores:

- median Σ where pb ≥ 0.3: 1.33 → 0.43 g cm⁻²; fraction above 1 g cm⁻²:
  77 % → 7 %;
- dense mass (young, intermediate, evolved): 15 732, 24 246, 13 215 →
  536, 1714, 1849 M☉; fraction of the cloud mass 1.43, 1.43, 1.65 →
  0.05, 0.10, 0.23.

The figure colour scale follows the lower values: grey 0–1 g cm⁻² (linear),
colour 1–3 g cm⁻² (square root), against 0–2 and 2–10 in Paper II
(`[plot_style]`, author 2026-11-20).

The temperature and the intercore maps are byte-identical to the Paper II
versions. `diagnostics/2026-11-13/` records the ×2.8 normalisation variants of
the Paper I column-density maps; they are not task inputs.

## Known limitations

- The combined map's absolute scale is marked "pending" in `combine_continuum`;
  rerun this task once the product is promoted.
- The intercore background comes from the old (CASA) 12m+7m image used by getsf.
