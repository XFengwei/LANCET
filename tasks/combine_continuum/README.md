# Combine continuum

This task makes and validates the 1.3 mm continuum image that combines the
MIRIAD 12m+7m mosaic with the SED-based total-power (TP) model, in the uv plane.
It has two phases: the combination (J-comb, adopted; feather, cross-check) with
the checks of the TP factor, and the diagnostics that show how complete the
interferometric and combined images are in the uv plane (Paper II Fig. 2).

It replaces `linear_combine_continuum` and `combination_diagnostic`, which were
retired on 2026-11-20 after the author approved the comparison (see "Replacing
the old tasks"). Nothing is promoted to `data/delivered/` yet.

## Run

```bash
PYTHONPATH=src conda run -n py312 python tasks/combine_continuum/run.py --stage combine --overwrite
PYTHONPATH=src conda run -n py312 python tasks/combine_continuum/run.py --stage diagnostics --overwrite
PYTHONPATH=src conda run -n py312 python tasks/combine_continuum/run.py --stage all --overwrite
```

`combine` runs phase 1, `diagnostics` phase 2, `all` both in order. Any single
stage can be run alone with `--stage <name>`. Existing outputs are kept unless
`--overwrite` is given. A stage whose input from another stage is missing stops
with a message naming that stage.

Resources: J-comb and feather use a padded FFT grid of about 50 M pixels
(complex64, about 2 GB peak memory, a few minutes per product, about 150 MB per
product image). The diagnostic phase takes 10 to 20 minutes and only reads the
measurement sets; run it in the background.

### Stages

Phase 1, combination. Order: `inspect`, `prepare`, `tp_factor`, `jcomb`,
`feather`, `box_spectrum`, `planck_calibration`. `box_spectrum` needs the J-comb
product for its mock correction and stops if it is missing.

| Stage | Does | Writes (under `results/combine_continuum/`) |
|---|---|---|
| `inspect` | inputs, grids, beams, units, frequencies, provenance | nothing |
| `prepare` | 12m+7m true sky in MJy/sr with its mask; both TP models regridded to its grid with the edge taper; checks the TP footprint margin and the regridding power change | `prepared/` |
| `tp_factor` | global overlap factor on the prepared inputs (documented, not applied) | `checks/tp_factor/` |
| `jcomb` | J-comb with dust + free-free and with dust only, standard and extended products, and their measurements | `products/`, `checks/jcomb/`, `checks/zero_spacing/` |
| `feather` | feather, same inputs, products and measurements | `products/`, `checks/feather/`, `checks/zero_spacing/` |
| `box_spectrum` | power spectra in the middle and south boxes, with the mock correction | `checks/box_spectrum/` |
| `planck_calibration` | TP against Planck 217 GHz; nothing is applied | `checks/planck_calibration/` |

Phase 2, diagnostics. Order: `uv`, `overlap`, `full_sampling`,
`amplitude_alignment`, `figure2`, `compare_images`, `dirty_beam`,
`aca_forward_test`, `joint_images`, `pointing_coverage`. It needs the phase 1
products. Outputs go under `results/combine_continuum/diagnostics/`.

| Stage | Does | Writes (under `diagnostics/`) |
|---|---|---|
| `uv` | uv coverage of the measurement sets | `uv_coordinates/`, `figures/` |
| `overlap` | row cache and matched-cell statistics | `overlap_tests/` |
| `full_sampling` | full-sampling profiles | `full_sampling/` |
| `amplitude_alignment` | 7m moved to the 12m primary beam with the 12m+7m image as sky model; TP through the same beam, with and without its own beam; the image curves of Fig. 2 (c) | `amplitude_alignment/` |
| `figure2` | Paper II Fig. 2 | `figures/` |
| `compare_images` | 7m and 12m-only pbcor images against the MIRIAD true-sky image | `figures/` |
| `dirty_beam` | delivered 12m+7m PSF and its restoring-beam fit | `dirty_beam/` |
| `aca_forward_test` | sky images through the native 7m primary beam at the measured ACA uv points | `aca_forward_test/` |
| `joint_images` | CASA (archived) against MIRIAD true-sky: central and south cutouts | `joint_image_comparison/` |
| `pointing_coverage` | pointing centres on the MIRIAD gain and the CASA pb | `pointing_coverage/` |

`manifest.json` at the output root lists every stage run with its settings and
the SHA-256 of its inputs.

## Inputs

All inputs are declared in `task.toml`; the code checks each key's source.

| Key | File | Source | Used by |
|---|---|---|---|
| `interferometer` | `data/delivered/alma/continuum/IRDC316.8.7MTM2TM1.robust0.5.cont.clean.truesky.pb0.2.fits` | delivered | prepare, figure2, compare_images, joint_images, amplitude_alignment |
| `gain` | `…7MTM2TM1.robust0.5.cont.gain.fits` | delivered | prepare, measurement region, pointing_coverage |
| `pb` | `…7MTM2TM1.robust0.5.cont.pb.fits` | delivered | extended products (taper weight) |
| `psf` (+ `psf_imfit_log`) | `…7MTM2TM1.robust0.5.cont.psf.fits` (+ `.psf.imfit.log`) | delivered | dirty_beam |
| `aca_pbcor` | `…IRDC316.8.7M.robust0.5.cont.image.pbcor.fits` | delivered | compare_images, aca_forward_test |
| `twelve_m_pbcor` | `…IRDC316.8.TM2TM1.robust0.5.cont.image.pbcor.fits` | delivered | compare_images |
| `casa_image`, `casa_pb` | `…/casa-archived/IRDC316.8.7MTM2TM1.cont.image.pbcor.fits`, `….cont.pb.fits` | delivered | joint_images, pointing_coverage |
| `tp_dust`, `tp_dust_ff` | `results/synthesize_continuum/model/{dust,continuum}_224GHz.fits` | `synthesize_continuum` product | all TP use |
| `planck_217`, `planck_co`, `planck_rimo`, `planck_fetch_record` | `data/delivered/archive/planck/hfi_cutouts/` | delivered | planck_calibration |
| `atlasgal_870` | `data/delivered/archive/infrared/atlasgal+planck.fits` | delivered | one Planck fit variant |
| `core_catalogue` | `results/source_catalogue/tables/source_catalogue.ecsv` | task product | compact-source check |
| `[measurement_sets]` | the four measurement sets in `data/reduction/alma/7MTM2TM1/cont/ms/` (12m combined, C5, C2, 7m) | declared read-only intermediate | uv, overlap, full_sampling, amplitude_alignment, aca_forward_test, pointing_coverage |

- Phase 1 reads the 12m+7m image and the TP models with `read_plane`, which
  accepts only Jy/beam with a beam; `prepare` converts them to MJy/sr with each
  image's own beam. The gain and pb maps are read with `read_gain` (no unit or
  beam expected; the shape must match the 12m+7m image). Phase 2 reads images with `read_mjysr` or, for the sampled sky
  models, its own reader: Jy/beam is converted with the image's own beam,
  MJy/sr is passed through, anything else is refused.
- The TP frequency is read from the `FREQGHZ` card of the TP models (223.8 GHz;
  the "230 GHz" in the published model name is a label). Both models must carry
  the same value; it sets the product file label. The primary-beam frequency of
  `pointing_coverage` is the same value.
- MIRIAD products: `truesky` is clean / gain, gain = min(1, pb / 0.508), and
  the `pb0.2` image is blank below pb 0.2.

## Method

### Combination

For same-grid, same-unit, padded planes `low` (TP) and `high` (12m+7m):

```text
C(k) = W_low(k) * R(k) * f * FT(low) + [1 - W_low(k)] * FT(high)
combined = FT^-1(C)
```

`k` is the uv distance in kλ.

- **J-comb** (Jiao et al. 2022): `R = G_high / G_low`, with `G_high` the
  elliptical restoring beam of the 12m+7m image (0.712″ × 0.606″, PA 47.4°)
  and `G_low` the circular 18.1″ TP beam. `W_low` is 1 below `k1 = 5 kλ`, 0
  above `k2 = 9 kλ`, and a Butterworth taper of order 5 in between,
  `b(k) = 1 / (1 + (k/k0)^(2n))` with `k0 = (k1 + k2)/2`, rescaled to run from
  1 to 0. **Deviation from the paper:** the printed form diverges at `k0`; the
  rescaled Butterworth is used instead. The TP is used only where
  `G_low >= 0.1`.
- **Feather** (cross-check): `W_low = G_low`, `R = 1`, so
  `C = f * FT(TP) + (1 - G_low) * FT(12m+7m)`.
- **Preparation:** the true-sky image (NaN set to 0 for the FFT; its finite
  pixels define the mosaic mask). The TP is continued past its footprint with
  an 18″ cosine taper and regridded with a cubic spline; this may change its
  power by at most 2 % in 1 to 9 kλ (checked in `prepare`). The TP footprint
  must extend at least 18.1″ beyond the mosaic mask, or `prepare` stops.
  Padding is 120″ (TP mirrored with a cosine roll-off, 12m+7m zeros).

### TP factor

`f = 1.0`: the TP is used as delivered. One factor serves both TP models; it is
determined from the dust + free-free model only, and the dust-only model takes
the same value. The power-spectrum matched value 1.4 of the earlier work is no
longer in the configuration. Three checks support `f = 1.0` and none changes
it:

1. **`tp_factor`**: matches the TP to the 12m+7m power spectrum over 20 to 50″.
2. **`box_spectrum`**: the same comparison in two boxes, with a mock of known
   factor 1.2 to calibrate the estimator.
3. **`planck_calibration`**: TP against Planck 217 GHz (CO removed): linear
   fits and the brightest Planck beam.

Values from the retired `linear_combine_continuum` run (Paper II inputs), recorded
in `results/combine_continuum/comparison/old_vs_new.md` next to this task's values:

| Check | Result |
|---|---|
| `tp_factor`, region as originally determined (includes the mosaic rim, bilinear regridding) | 1.39 |
| `tp_factor`, prepared inputs (gain defined, spline regridding) | 1.07 |
| Mock of true factor 1.2: power spectrum, gain-defined region / gain ≥ 0.5 | 1.00 / 0.87 |
| Mock of true factor 1.2: cross-spectrum at 12m+7m resolution; real data | 1.20; 1.18 |
| `box_spectrum`, power spectra, lowest TP point: middle box; south box | 0.99 (mock 1.03); 1.20 (mock 1.09) |
| `box_spectrum`, amplitude median: middle; south | 1.08 (mock 1.24); 1.21 (mock 1.03) |
| `box_spectrum`, cross-spectrum: middle; south | 0.91 (mock 1.16); 1.08 (mock 1.01) |
| `box_spectrum`, corrected by the mock (× 1.2 / mock): middle; south | 0.94 to 1.15; 1.29 to 1.41 |
| `planck_calibration`, brightest beam, direct ratio TP / Planck (disc mean; peak) | 0.97; 1.02 (implied factor 1.03; 0.98) |
| `planck_calibration`, brightest beam, fitted background removed | 1.17; 1.22 (implied factor 0.85; 0.82) |
| `planck_calibration`, brightest beam, differential against the faintest 20 % | 1.29; 1.36 (implied factor 0.77; 0.74) |
| `planck_calibration`, brightest beam, no CO removal | 0.89; 0.93 (implied factor 1.12; 1.07) |
| `planck_calibration`, linear fit TP = gain × Planck + offset (gain, coverage ≥ 0.9) | 1.15 |

The estimates of the factor range from about 1.0 to 1.4 with method and
region. If the 12m+7m image is incomplete across 5 to 9 kλ, every estimator
returns too low a value. The brightest-beam Planck ratio depends on how much of
Planck's signal is assigned to a background the TP model lacks, and rests on
under one independent beam.

### Extended products

Each stage also writes, per TP model, an extended product that continues the
combination into the TP-only area beyond the mosaic. The 12m+7m plane is
multiplied by a weight in the mosaic primary beam (`pb`):
`w = 0` for pb ≤ pb₀, `w = 1` for pb ≥ pb₁ and
`w = ½[1 − cos(π (pb − pb₀)/(pb₁ − pb₀))]` between, with
`[combine] extended_taper_pb = [0.2, 0.3]`. The small scales fade out across
that rim instead of stopping at pb 0.2, where the true-sky noise is 2.5 times
the centre value. The TP, its factor and all settings are those of the standard
product. The extended product is blank only outside the TP footprint; beyond
the mosaic it is the TP part of the combination alone. Where pb ≥ pb₁ it equals
the standard product apart from the small large-scale change that the rim
weight makes through the FFT. The two combinations run one after the other, so the FFT planes of both are never held at once. The diagnostics keep using
the standard, blanked products.

## Outputs

`results/combine_continuum/products/` (float32, Jy/beam of the 12m+7m
restoring beam, on the 12m+7m grid):

| File | Content |
|---|---|
| `jcomb_continuum_224GHz.fits` | 12m+7m + dust + free-free TP (adopted) |
| `jcomb_dust_224GHz.fits` | 12m+7m + dust-only TP |
| `feather_continuum_224GHz.fits`, `feather_dust_224GHz.fits` | cross-check |
| `{jcomb,feather}_{continuum,dust}_224GHz_extended.fits` | the four products continued into the TP-only area (see "Extended products") |

- Products are blank outside the input mosaic mask (pb ≥ 0.2,
  `fill_outside_mosaic = false`). The previous task kept TP-only pixels there.
- Header: `COMB_*` cards (method, k1, k2, order, minimum TP response, TP factor,
  TP flavour, TP frequency from `FREQGHZ`, padding, absolute scale "pending")
  and the SHA-256 of every input. The extended products add `COMB_EXTENDED = T`,
  `COMB_TAPER_PB0` and `COMB_TAPER_PB1`, and carry `COMB_FILL_OUTSIDE_MOSAIC = T`.
- Extended checks: `checks/zero_spacing/{method}_{kind}_224GHz_extended_at_tp_beam.fits`
  and `zero_spacing_ratio_{method}_{kind}_224GHz_extended.fits`; in
  `checks/<method>/measurements.json`, per TP model, `extended` holds the
  zero-spacing check over `mosaic`, `tp_only` and `all`,
  `tp_only_ratio_to_scaled_tp`, and `ratio_to_standard_pb_ge_0.35` (median
  extended / standard where pb ≥ 0.35 and the standard is above 3σ). All other
  measurements stay on the standard products.

Checks, written as JSON plus figures (PDF and PNG):

- **Mock-sky gate** (`jcomb`): the known-answer test (unit TP factor,
  interferometer high-passed at 6.6 kλ) must recover flux within 2 % and rms
  within 1 %.
- **Zero spacing**: the product smoothed to 18.1″ against the TP; total flux
  inside the mosaic within 5 %.
- **Overlap continuity**: combined / 12m+7m within 2 % from k2 to 20 kλ.
- **Compact sources**: background-subtracted peaks at the catalogue cores,
  product / 12m+7m within 10 %.
- **Negative bowls**: fraction of pixels below −3σ (σ from the noise box).

Measurements use the region with gain ≥ 0.5 of the delivered gain map.

## Paper II Fig. 2

`diagnostics/figures/figure2_combination_diagnostics.{pdf,svg,png}`, with its
manifest and source CSVs. Three panels: (a) pointings, (b) uv coverage with the
short-baseline inset, (c) amplitudes and visibility counts.

Panel (c) shows the 7m, C2 and C5 RMS in matched uv and orientation cells (7m
moved to the 12m primary beam through the 12m+7m sky model), the TP through the
12m beam with and without its own beam, and two image curves, chosen by `image_curves_shown` and
`image_curve_labels` in `[figure2]`: 12m+7m (dashed crimson) and J-comb (solid
green, labelled "combine"). The feather curve is computed and kept in
`image_curves.csv` but not drawn.

Image curves are computed in the `amplitude_alignment` stage with the TP-curve
method:

- Each image is area-binned to the 0.5″ model grid, multiplied by the actual
  12m pointing primary beam at both basebands (216.74 and 230.95 GHz), with the
  same pointing, band and orientation population as the TP curve, and sampled
  along a continuous radius from 1 kλ.
- The image's restoring beam (circular Gaussian of the geometric mean of BMAJ
  and BMIN) and the 0.5″ binning response are divided out. For the 12m+7m beam
  the residual of the elliptical shape stays within 0.5 %.
- The curve stops at the lowest of three limits: beam response 0.1, half the
  model-grid Nyquist limit (about 103 kλ), and the half-pixel bilinear
  interpolation limit (2 % loss, about 26 kλ). The last applies in practice, so
  the curves cover 1 to 26 kλ.
- Below about 7 kλ the curves are smoothed by the 12m primary beam. Each sample
  is the image multiplied by that beam, so it samples the image visibilities
  convolved with the beam's uv footprint. For a Gaussian beam of FWHM θ the
  footprint is a Gaussian of FWHM 4 ln2 / (π θ): with `pb_12m_fwhm_arcsec =
  26.0` at `pb_reference_frequency_ghz = 225.0`, θ = 26.99″ and 25.33″ at
  216.74 and 230.95 GHz, and the footprint FWHM is 6.74 and 7.19 kλ. Below
  about 5 kλ each sample therefore mixes in power from 6 to 10 kλ. The curves
  are kept from 1 kλ, but the 12m+7m curve does not fall away there: in the
  first full run it flattens at roughly half the beam-corrected TP curve
  (ratio 0.49 to 0.56 over 1 to 4 kλ, 0.63 at 5 kλ), while J-comb and feather
  stay at 0.83 to 0.94 of it over 1 to 5 kλ
  (`diagnostics/amplitude_alignment/image_curves.csv` against
  `tp_corrected_common_mjy` in `tp_continuous.csv`). The image-curve
  description in the `amplitude_alignment` and Fig. 2 manifests states the
  same smoothing.

### ACA forward test

Sky models: 12m+7m (no TP), J-comb, feather and the ACA-only image. The TP
model itself is always predicted. For each sky model and uv bin the stage
reports the ratio of predicted to measured ACA amplitude (`<model>_to_aca`, with
pointing-pair bootstrap percentiles). No TP-scale estimate is derived from it:
the 12m+7m image is not the f = 0 member of J-comb or feather, so a scale
interpolated between them has no meaning. A TP-only sky model is not included,
because the forward-test code cannot take one.

## Limitations

- The TP is an SED model, not a measurement. Its scale rests on the three
  checks above, which disagree between about 1.0 and 1.4.
- The dust-only images are not free of free-free emission at small scales, and
  show negative bowls around bright sources.
- The 12m+7m image is assumed complete across 5 to 9 kλ; its shortest spacing is
  6.6 kλ.
- CASA 6.7.2 `feather` was tested and not used: on these images it returns the
  sum of its inputs. The feather here is the task's own implementation.
- The feather `high_k_join` check is flagged by construction: feather weights
  the 12m+7m image by 1 − G_low, which is 0.89 at k2 = 9 kλ and 0.98 at 12 kλ,
  so combined / 12m+7m from k2 to 20 kλ falls short of 1 by more than the 2 %
  tolerance (0.975 for both TP models in the first full run,
  `checks/feather/measurements.json`). It is a property of the method, not a
  failure of the product.
- `[absolute_scale]` is pending: the products are on the scale set by the 12m+7m
  calibration and `f`. The `planck_calibration` stage gives the comparison the
  decision rests on.
- The Paper II numbers quoted above come from the previous task and its inputs;
  they will change with the true-sky input and the 223.8 GHz TP.

## Promotion

Done only on the author's go-ahead, after the TP factor and the absolute scale
are decided:

- copy `jcomb_continuum_224GHz.fits` and `jcomb_dust_224GHz.fits` to
  `data/delivered/alma/continuum/` as
  `IRDC316.8.TP7MTM2TM1.robust0.5.cont.jcomb.fits` and
  `IRDC316.8.TP7MTM2TM1.robust0.5.cont.jcomb.dust.fits`;
- add README rows there (source run, settings, checksum);
- then repoint the downstream tasks that read a combined image
  (`core_formation`, `column_density` and others).

The feather products are not promoted.

## Replacing the old tasks

After the first full run, J-comb inside gain = 1 was compared with the results
of the retired `linear_combine_continuum` run (expected agreement about 1 %;
differences come from the true-sky rim and the +0.27 % free-free change), the
factor numbers were tabulated old against new, and Fig. 2 and the kept
diagnostic figures were compared with those of the retired
`combination_diagnostic` run. The comparison is
`results/combine_continuum/comparison/old_vs_new.md`. After the author approved
it, both old tasks, their scripts and their results were retired on
2026-11-20.
