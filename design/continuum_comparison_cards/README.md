# Continuum comparison cards

Three 1.3 mm views, in this order:

1. **Total-power**: the SED-based TP model (`synthesize_continuum`, filled
   model `continuum_224GHz_fill36.fits`), shown as it is (factor 1.0).
2. **ALMA Interferometry**: the delivered 12m+7m true-sky image, blank below
   mosaic primary beam 0.2
   (`IRDC316.8.7MTM2TM1.robust0.5.cont.clean.truesky.pb0.2.fits`).
3. **Total-power + Interferometry**: the grid of the ALMA multiwavelength
   card (`design/multiwavelength_cards/grid/mm_i.fits`): the delivered J-comb
   `extended36` product, the TP model beyond the 12m+7m grid, and a display
   cross-fade to the TP model over mosaic pb 0.2–0.5.

All three use the sky frame and pixel shape of `design/multiwavelength_cards`
(16.56 × 10.445 arcmin, north up, east left) and one brightness mapping:
inferno, square-root stretch, 0–1200 MJy/sr, the scale of the ALMA
multiwavelength card and of Paper II Fig. 1 and Fig. 3. Images in Jy/beam are
converted to MJy/sr with their own beams. Blank pixels display as zero. No
colourbar or beam text is shown.

The first card has the N/E compass, the third the lower-left 1 pc scale bar.
All three carry the same lower-right contact details. Each artwork area is
85.6 × 54.0 mm with a 5 mm black border (95.6 × 64.0 mm in total). PNG and LZW
TIFF exports use 1200 dpi; PDF page sizes are exact.

## One value framework

The three grids are in MJy/sr on the same 3312 × 2089 frame and are drawn
with the same mapping, so equal colours mean equal surface brightness. Unit
conversion uses each image's own beam: 1 Jy/beam = 114.6 MJy/sr for the 18.1″
TP model and 8.695 × 10⁴ MJy/sr for the 0.71″ × 0.61″ images.

| Card | Coverage | Median | 99th percentile | Maximum (MJy/sr) |
|---|---|---|---|---|
| Total-power | 100 % | 8.6 | 95 | 265 |
| ALMA Interferometry | 7 % (mosaic) | 13.0 | 229 | 5374 |
| Total-power + Interferometry | 100 % | 8.6 | 95 | 5373 |

Cross-checks (asserted by `verify_cards.py`): the combined map smoothed to
18.1″ over the TP model has a median of 1.009 (16–84 %: 0.973–1.043); the
combined map over the 12m+7m image at bright pixels (> 300 MJy/sr) has a
median of 1.089, the extended emission the TP adds under the compact
structure. Inside the mosaic the medians are 13.0 (12m+7m), 41.9 (combined)
and 41.4 MJy/sr (TP): the interferometer alone misses the extended level.
About 19 % of the 12m+7m pixels are negative and 0.02 % of the combined
pixels exceed 1200 MJy/sr; both clip to the ends of the scale.

## Reproduce

From the repository root:

```sh
PYTHONPATH=src conda run --no-capture-output -n py312 python design/continuum_comparison_cards/make_cards.py --refresh
PYTHONPATH=src conda run --no-capture-output -n py312 python design/continuum_comparison_cards/verify_cards.py
```

`--refresh` reprojects the inputs (`prepare_inputs.py`). `verify_cards.py`
checks the shared WCS, units, image and border dimensions, and packages the
three card exports into `continuum_comparison_cards.zip`.

## History

Redrawn on 2026-11-21 and moved from `results/continuum_comparison_cards/`
to `design/` (author). The earlier version showed the v3 combined map, with
the TP card multiplied by 3.05 to match it and a fitted TP fill outside the
field; with the new combination (TP factor 1.0) no display factor or fit is
needed. `provenance/make_dual.py` is an earlier script, kept as written. Only
the scripts, this README and the JSON settings are tracked in git; the
rendered files and the `grid/` cache stay local.
