# Synthesize continuum

This task builds the 1.3 mm total-power (TP) continuum model: the Paper I dust
model plus the free-free model, both at 223.8 GHz. It also makes a 250 GHz
version, which is compared with the archival SIMBA 1.2 mm map. The 223.8 GHz
products are the TP inputs of `linear_combine_continuum`.

## Inputs

| Key | File | Notes |
|---|---|---|
| `dust` | `data/delivered/sedfitter/Dustmodel_223.8GHz.fits` | Paper I dust model, 18.1″ beam, 6″ grid; byte-identical to the published `TP_SEDsynthetic230GHz_cont.fits` (SHA-256 `fb3c3de0…`), whose "230GHz" is only a label |
| `freefree` | `results/construct_freefree_model/Te_8000K/ff224_map_smoothed_18arcsec_regrid2dust.fits` | 8000 K free-free model at 223.8 GHz, smoothed to 18″ and put on the dust grid |
| `simba_250ghz` | `data/delivered/archive/infrared/SIMBA_1.2mm.fits` | archival SIMBA map, for validation only |

## Model

Both inputs are at the model frequency `reference_frequency_ghz = 223.8` (the
ALMA continuum reference frequency). For a target frequency `nu`:

```text
dust(nu)     = dust(223.8 GHz)     * (nu / 223.8 GHz)^3.8
freefree(nu) = freefree(223.8 GHz) * (nu / 223.8 GHz)^-0.1
continuum    = dust(nu) + freefree(nu)
```

Missing free-free pixels are set to zero before addition; missing dust pixels
stay missing. At 223.8 GHz the model is the plain sum of the two inputs.

### Change from Paper II (2026-11-20)

Paper II assumed both inputs were at 230 GHz.

- **1.3 mm model:** unchanged apart from the free-free part (+0.27%,
  from `construct_freefree_model`). The total rises by 0.02%.
- **250 GHz model:** 10.9% brighter (dust ×(230/223.8)^3.8 = 1.109; total
  flux +10.3%).
- **SIMBA bright peaks:** model/SIMBA moved from 0.80, 1.07, 1.23, 1.28 to
  0.86, 1.18, 1.36, 1.39 (`validation/simba_peak_measurements.csv`).

The superseded `SEDsynthetic_Dust_FF_{230GHz,250GHz,250GHz_smooth24arcsec}.fits`
in `data/delivered/sedfitter/` were this task's Paper II outputs (equal to
float32 rounding). They were removed, as were the old `*_230GHz` products.

## Outputs

`results/synthesize_continuum/`:

- `model/`: `dust_`, `freefree_` and `continuum_` maps at `224GHz` and
  `250GHz`. File labels are rounded; the exact frequency is in `FREQGHZ`.
- `validation/`: the 250 GHz model smoothed to 24″, the calibrated SIMBA map,
  and SIMBA on the model grid.
- `figures/`:
  - `continuum_components_224GHz`: dust, free-free and sum (Paper II Fig. C2).
  - `simba_comparison_250GHz`: SIMBA and the beam-matched model (Fig. C1, top).
  - `sed_demonstration_224GHz`: illustrative SED.

Crop, intensity range, contours, colour map, scale bar, resolution and
formats are set in the `[plot]` tables.

## Stages

- `inspect`: checks inputs, grids, beams, units and model readiness; writes nothing.
- `model`: component and combined maps.
- `validate`: SIMBA comparison maps (after `model`).
- `plot`: the three figures (after `validate`).
- `all`: all four in order.

```bash
PYTHONPATH=src conda run -n py312 python tasks/synthesize_continuum/run.py --stage all --overwrite
```

Existing products are kept unless `--overwrite` is given.

## SIMBA comparisons (Paper II R10 and Fig. C1)

Three scripts in this folder read the `validation/` products. Run them after
`validate`:

```bash
PYTHONPATH=src conda run -n py312 python tasks/synthesize_continuum/simba_pixel_comparison.py
PYTHONPATH=src conda run -n py312 python tasks/synthesize_continuum/simba_peak_comparison.py --overwrite
PYTHONPATH=src conda run -n py312 python tasks/synthesize_continuum/simba_fourier_comparison.py --overwrite
```

- **Pixel comparison** (Fig. C1, bottom): SIMBA and the model pixel by pixel,
  above the 0.25 Jy/beam SIMBA threshold. It has no argument parsing and
  always rewrites its outputs.
- **Peak comparison:** SIMBA local maxima ≥ 1 Jy/beam within the display crop.
  Their brightness is compared with the 24″ model at the same positions, both
  directly and after subtracting the median background in a 36–60″ annulus.
- **Fourier comparison:** Fourier amplitudes over the common footprint. See the
  R10 diagnostic report for its window, beam and noise limits.

The 0.20 SIMBA flux calibration factor was confirmed by the author. Annular
contrasts are affected by neighbouring emission and negative bowls. They do
not independently recover missing spatial scales.

These three scripts still contain fixed paths and settings outside
`task.toml`. Folding them into `run.py` as stages is left for later.
