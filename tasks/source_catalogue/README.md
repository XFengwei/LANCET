# Source catalogue and flux diagnostics

This task turns the getsf extraction into the project source catalogue. It
applies the reliability and compactness cuts, corrects the 12m+7m
measurements for the mosaic gain, samples the intercore background, the
12m+7m+TP background and the dust temperature at each source, calculates
effective deconvolved sizes, assigns stable IDs and tags the three
evolutionary subregions. It also draws the two getsf flux diagnostics
(Paper II Fig. 8).

getsf is not rerun (D7). In the getsf catalogue, suffix `01` is the 12m+7m
measurement and `02` the 12m+7m+TP measurement of the original extraction.

## Inputs

| Key | File | Notes |
|---|---|---|
| `getsf_catalogue` | `data/extraction/getsf/triple/results/G316.sw.sources.ok.cat` | getsf catalogue |
| `intercore_background` | `data/extraction/getsf/triple/results/+components/G316.sw.1302.obs.sbackground.fits` | getsf background of the 12m+7m image (MJy/sr) |
| `source_model` | `…/+components/G316.sw.1302.obs.sources.fits` | getsf source model of the 12m+7m image (MJy/sr) |
| `combined_continuum` | `data/delivered/alma/continuum/IRDC316.8.TP7MTM2TM1.robust0.5.cont.image.pbcor.fits` | J-comb 12m+7m + dust + free-free TP (`combine_continuum`), true sky, Jy/beam |
| `gain` | `data/delivered/alma/continuum/IRDC316.8.7MTM2TM1.robust0.5.cont.gain.fits` | MIRIAD mosaic gain |
| `temperature` | `results/column_density/prepared/temperature_allscale.fits` | 18″ dust temperature on the ALMA grid |
| region files | `region/trichotomy/*_shrink.reg` | subregion membership |

The getsf 12m+7m image (`data/extraction/getsf/triple/images/G316_12m7m.fits`)
is the MIRIAD clean image in MJy/sr (ratio 1.000), i.e. sky × gain. Every
12m+7m getsf quantity is therefore attenuated where gain < 1.

## Method

**12m+7m+TP background (stage `background`).** The source model is converted
to true sky and removed from the combined continuum:

```text
background_12m7mTP = combined_continuum [MJy/sr] − source_model / gain
```

Pixels without sources keep the combined value; NaN stays NaN. Written to
`maps/background_12m7mTP.fits` (float32, MJy/sr). It replaces the getsf
background of the old v3 12m+7m+TP image (`sw.1303`).

**Sampling.** `Ibkg1` (intercore), `Ibkg2` (12m+7m+TP background), `Tenv`
and `Gain` are read at the getsf peak pixel (`sample_at_getsf_peaks`).

**Gain correction.** `FXP_BST01`, `FXP_ERR01`, `FXT_BST01`, `FXT_ERR01` and
`Ibkg1` are divided by `Gain` where it is below 1 (19 of 231 sources,
gain 0.462–0.984). The raw getsf values are kept as `*_RAW`. The selection
cuts and the flux diagnostics use the raw getsf measurements.

## Stages

- `inspect`: inputs and cut counts; writes nothing.
- `background`: the 12m+7m+TP background map.
- `build`: `tables/source_catalogue.ecsv` and the publication table
  `source_catalogue.tex` (needs `background`).
- `diagnostics`: `GetsfResult_FluxMeasure` and `GetsfResult_FluxRatio`
  (PDF, PNG) and their tables; peak/total ratio and
  ΔS/S = (S₁₂ₘ₊₇ₘ₊TP − S₁₂ₘ₊₇ₘ)/S₁₂ₘ₊₇ₘ above `peak_snr_min`. Colours are
  set by `reference_color` and `comparison_color` in `[flux_diagnostics]`.
- `all`: the four in order.

```bash
PYTHONPATH=src conda run -n py312 python tasks/source_catalogue/run.py --stage inspect
PYTHONPATH=src conda run -n py312 python tasks/source_catalogue/run.py --stage all --overwrite
```

Existing products are kept unless `--overwrite` is given.

## Readers

`physical_properties` (`Tenv`, `Ibkg1`, `Ibkg2`, `FXT_BST01`, `FXP_BST01`
and errors), `core_formation`, `combine_continuum` (compact-source check),
`freefree_core_contamination`.

## Change from Paper II (2026-11-20)

- `Ibkg2` from the new background: median 309 → 148 MJy/sr (new/old median
  0.53, range 0.22–1.20). The new combined map uses TP factor 1.0; the old
  getsf input was the v3 map. `Ibkg2 > Ibkg1` for 230 of 231 sources.
- 19 sources with gain < 1: fluxes, errors and `Ibkg1` rise by 1.016–2.166.
  All other sources keep their values; IDs and order are unchanged.
- Flux diagnostics: identical.
