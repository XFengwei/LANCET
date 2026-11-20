# Plot all-scale continuum

Publication and clean views of the 12m+7m+TP 1.3 mm continuum (Paper II
Fig. 3, `ImageCombine_Continuum_Dust_FF_TP7MTM2TM1_DualColor`), its dust-only
version, and the TP model.

## Inputs

| Key | File | Notes |
|---|---|---|
| `combined_continuum` | `data/delivered/alma/continuum/IRDC316.8.TP7MTM2TM1.robust0.5.cont.image.pbcor.extended.fits` | J-comb 12m+7m + dust + free-free TP (TP factor 1.0), continued into the TP-only area |
| `dust_continuum` | `…TP7MTM2TM1.robust0.5.cont_dust.image.pbcor.extended.fits` | as above with the dust-only TP |
| `synthetic_continuum` | `results/synthesize_continuum/model/continuum_224GHz.fits` | 223.8 GHz TP model, 18.1″ |
| `primary_beam` | `…7MTM2TM1.robust0.5.cont.pb.fits` | MIRIAD mosaic pb; the annotated views show its 0.2 contour |

Every image is read with `lancet.imaging.read_mjysr` (Jy/beam converted with
its own beam) and plotted in MJy/sr. The combined images must be delivered;
the TP model must be a task result. All inputs are read-only.

## Stages

- `inspect`: input paths, shapes, precision and units; writes nothing.
- `combined`: annotated (PDF, PNG) and clean (PNG) views of the
  dust + free-free map.
- `sqrt`: **Paper II v3 Fig. 3** (author, 2026-11-20): the dust + free-free
  map on one colour bar
  (inferno, square-root stretch), annotated (PDF, PNG) and clean (PNG):
  `allscale_continuum_dust_freefree_sqrt*`.
- `synthetic`: clean view of the TP model.
- `dust`: annotated and clean views of the dust-only map.
- `all`: the four plotting stages.

```bash
PYTHONPATH=src conda run -n py312 python tasks/plot_allscale_continuum/run.py --stage inspect
PYTHONPATH=src conda run -n py312 python tasks/plot_allscale_continuum/run.py --stage all --overwrite
```

Figures go to `results/plot_allscale_continuum/figures/`. Existing figures
are kept unless `--overwrite` is given.

## Colour scales (`[styles]`, MJy/sr)

- `allscale`: −20 / 80 / 1200, confirmed by the author 2026-11-20 (grey linear up to the transition, colour
  above). Paper II used −21 / 233 / 1500 on the v3 map, whose TP contribution
  was about 3.1 × the delivered TP model. The new transition and maximum sit at
  the same percentiles inside pb ≥ 0.2 (233 → 78, 1500 → 1180); −20 is about
  −2σ (σ = 9.2 MJy/sr).
- `sqrt`: 0 to 1200 MJy/sr, inferno, square-root stretch, one colour bar.
  Values below 0 take the lowest colour; blank pixels (outside the TP model)
  are black (`blank_colour`). The pb 0.2 contour is cyan and the labels white.
- `synthetic`: −2.292 / 57.31 / 275.07, the former Jy/beam limits
  (−0.02 / 0.5 / 2.4) converted with the 18.1″ beam, so that view is unchanged.

## Review check (2026-11-20)

With the archived v3 maps (`allscale-archived/`) converted to MJy/sr, the
code reproduces the stored Fig. 3 clean view (45 of 15.7 M pixels differ).
The stored dust-only clean view (29 July) predated later plotting changes.
