# VVV mosaic

Seamless VVV J, H and Ks mosaic of the field and its colour composite. The
mosaic settings are also used by the near-infrared panel of
`multiwavelength_strip` (Paper II Fig. 1, `vvv_config`).

## Inputs

`data/delivered/archive/infrared/VVV/{d053,d054,d091,d092}_{j,h,ks}.fits`
(four tiles per band). The tile directory must be under `data/delivered/`.

## Method

- Each tile is rescaled to a common zero point (25.0) and its padded,
  under-exposed border (100 px) is blanked.
- Tile backgrounds are compared on 16-pixel median cells; smooth offset
  surfaces make overlapping tiles agree, fading into the tile interiors over
  240″, and each tile's weight ramps up over 40″ from its edge.
- The bands are coadded on a north-up ICRS TAN grid (0.34″ pixels,
  2725 × 3523).
- The colour composite maps Ks, H, J to red, green, blue with a clipped-noise
  black point, a 99.5 % white point, an asinh stretch and a saturation boost.

All settings are in `task.toml`.

## Stages

- `mosaic`: `results/vvv_mosaic/vvv_{j,h,ks}_mosaic.fits`.
- `rgb`: `g316.8_vvv_jhks_full.png` and a preview.
- `all` (default): both.

```bash
PYTHONPATH=src conda run -n py312 python tasks/vvv_mosaic/run.py --stage all
```

Products are rewritten on every run (no `--overwrite` guard).

## Review (2026-11-21)

Rerun: the three band mosaics and the colour PNG are byte-identical to the
previous products. Script moved from `scripts/`; README and tests added.
