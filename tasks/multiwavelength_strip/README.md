# Multiwavelength strip

Paper II Fig. 1: five panels on one rotated TAN grid that follows the filament
spine (14′ × 4′, 0.3″ pixels). From top to bottom: VVV J/H/Ks, Spitzer 3.6 and
5.8 µm with Herschel 70 µm, Herschel 250 µm, the ALMA 12m+7m+TP 1.3 mm
continuum, and MeerKAT 1.38 GHz. Geometry, inputs, stretches, annotations and
figure settings are in `task.toml`.

## Inputs

| Panel | Input |
|---|---|
| `nir` | VVV tiles, mosaicked with the `vvv_mosaic` settings (`tasks/vvv_mosaic/task.toml`) |
| `mir` | `data/delivered/archive/infrared/image_HHLI-PACS070.fits`, `image_GLIMPSE-IRAC3.fits`, `image_GLIMPSE-IRAC1.fits` |
| `fir` | `data/delivered/archive/infrared/image_HHLI-SPIRE250.fits` (saturated holes filled) |
| `mm` | `data/delivered/alma/continuum/IRDC316.8.TP7MTM2TM1.robust0.5.cont.image.pbcor.extended36.fits`, Jy/beam converted to MJy/sr (`to_mjysr`); beyond the 12m+7m grid filled from `results/synthesize_continuum/model/continuum_224GHz_fill36.fits` (`gap_fill_input`) |
| `radio` | `data/delivered/archive/radio/MeerKAT_subimage_1.38GHz_reprojected.fits` |

The `mm` product is the J-comb 12m+7m + dust + free-free TP of
`combine_continuum` (TP factor 1.0). Beyond the mosaic it is the TP model
alone, so the gap fill uses the same model without a fit. It is shown in
inferno with a square-root stretch from 0 to 1200 MJy/sr (`style = "sqrt"`),
the inferno view of Fig. 3 (`plot_allscale_continuum`, `[styles.sqrt]`).

## Stages

- `reproject`: every input on the strip grid, `results/multiwavelength_strip/grid/`.
- `plot`: `results/multiwavelength_strip/lancet_multiwavelength_strip.{png,pdf}`.
- `all` (default): both.

```bash
PYTHONPATH=src conda run -n py312 python tasks/multiwavelength_strip/run.py --stage all
```

The stages rewrite their products on every run (there is no `--overwrite`
guard).

## Change from Paper II (2026-11-21)

- Code check: with the removed v3 map rebuilt from the archived copy, the
  reprojection reproduces all eight infrared and radio grids exactly and the
  `mm` grid within 9e-4 (only in the filled area; Paper II filled with the
  earlier 230 GHz model).
- The `mm` panel now shows the delivered J-comb product. Paper II showed the
  v3 map, whose TP part was about 3 × the TP model: the old fill matched the
  model to v3 with a factor 2.92. The fitted model fill (`fill_input` and its
  settings) is replaced by the direct `gap_fill_input`; the colour map is
  inferno with a square-root stretch (author), replacing the grey/mako
  dual-colour scale.
- Display cross-fade (`blend_pb_input`, `blend_pb_range = [0.2, 0.5]`): the
  panel shows w·image + (1 − w)·TP model with a cosine weight w rising from 0
  at mosaic pb 0.2 to 1 at pb 0.5; outside the mosaic it is the model alone.
  This removes the abrupt rim and the seam at the grid edge. It affects the
  figure only. Negative noise pixels in the faintest parts of the mosaic
  (12m+7m rms 7.9 MJy/sr) show as black specks at the bottom of the scale.
- `reproduction_20261106/` in the results folder is the earlier reproduction
  record of the Paper II figure and is kept as it is.
