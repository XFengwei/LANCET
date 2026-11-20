# Design: continuum combination task (`combine_continuum`)

Date: 2026-11-20 (replaces the 2026-11-15 version). Status: design approved
section by section in conversation; this written spec awaits the author's
review before an implementation plan is written.

## 1. Purpose

One task produces and validates the 1.3 mm continuum image that combines the
MIRIAD 12m+7m mosaic with the SED-based total-power (TP) model. It replaces two
tasks that grew side by side:

- `linear_combine_continuum` (J-comb and feather products, TP-factor checks,
  Planck check), and
- `combination_diagnostic` (visibility and image diagnostics, Paper II Fig. 2).

The task has two phases:

1. **Combination**: prepare the inputs, check the TP factor three ways, and
   write the J-comb (adopted) and feather (cross-check) products.
2. **Diagnostics**: show how complete the interferometric and combined images
   are in the uv plane (Paper II Fig. 2) and keep the image diagnostics the
   paper and the referee response rely on.

Decisions taken with the author in the stage-by-stage review (2026-11-20):

- Interferometer input: the MIRIAD true-sky image cut at pb 0.2
  (`truesky.pb0.2`). The old `prepare` stage added the TP (true sky) to
  sky × gain in the mosaic rim; using true sky removes that inconsistency.
  The attenuated (`pbant`) input and the gain derivation are dropped.
- TP models: the 223.8 GHz products of `synthesize_continuum`.
- TP factor f = 1.0 (TP as delivered). Three checks support it and are kept:
  the global overlap factor (`tp_factor`), box power spectra
  (`box_spectrum`) and Planck (`planck_calibration`). The power-spectrum
  matched value 1.4 (`psd_matched`) is removed from the configuration.
- J-comb and feather settings are unchanged.
- Paper II Fig. 2 (three panels) is kept; panel (c) gains 12m+7m, J-comb and
  feather curves.
- Dropped stages (results archived to the Trash after the new task is
  verified): four-panel `figure1`, `visibilities`, `pb-curve`,
  `tp-beam-test`, `pointing-profiles`, `pointing-range`, `total-power`,
  `beam-match`.

Facts this design relies on:

- The Paper I dust model and the ALMA continuum are at 223.8 GHz; the "230GHz"
  in the published model name is a label. `synthesize_continuum` writes
  `FREQGHZ = 223.8`.
- MIRIAD products: `clean` is true sky only where gain = 1 (sky × gain
  beyond); `truesky` = clean / gain; gain = min(1, pb / 0.508); the `pb0.2`
  versions are blank below pb 0.2.
- The TP beam (18.1″) transmits 50 % at 5.0 kλ and 10 % at 9.2 kλ; the 12m+7m
  image's shortest spacing is 6.6 kλ.
- CASA 6.7.2 `feather` returns the sum of the two inputs on these images; the
  feather here is the task's own implementation.

## 2. Layout

- `tasks/combine_continuum/run.py`: loads `task.toml`, resolves paths through
  `ProjectPaths`, exposes `--config`, `--stage <name>` and `--overwrite`.
  Existing outputs are never replaced without `--overwrite`.
- `tasks/combine_continuum/task.toml`: all paths and settings.
- `tasks/combine_continuum/README.md`: run commands, inputs, method, outputs,
  limitations.
- Calculations stay in the existing `src/lancet/` modules
  (`combination`, `combination_prep`, `combination_validation`,
  `combination_figures`, `planck_calibration`, `combination_tests`,
  `full_sampling`, `amplitude_alignment`, `combination_figure`,
  `aca_forward_test`, `joint_image_comparison`, `dirty_beam`,
  `interferometry`, `visibility_diagnostics`, `pointing_coverage`).
  Code used only by dropped stages is removed; `smooth_to_circular_beam`,
  `_predict`, `sky_model` and `beam_match_and_reproject` stay because kept
  stages use them, and `C_METRES_PER_SECOND` moves to the module that still
  needs it.
- Outputs only under `results/combine_continuum/`.

### Stages

`--stage combine` runs phase 1, `--stage diagnostics` phase 2, `--stage all`
both in order. Each stage can also be run alone.

Phase 1, combination:

| Stage | Does | Writes |
|---|---|---|
| `inspect` | inputs, grids, beams, units, frequencies, provenance; footprint margin | nothing |
| `prepare` | 12m+7m true-sky in MJy/sr, mask from its finite pixels; both TP models regridded to its grid with the edge taper | `prepared/` |
| `tp_factor` | global overlap factor on the prepared inputs (documented, not applied) | `checks/tp_factor/` |
| `jcomb` | J-comb, dust and dust + free-free, plus its measurements | `products/`, `checks/jcomb/`, `checks/zero_spacing/` |
| `feather` | feather, same inputs and measurements | `products/`, `checks/feather/`, `checks/zero_spacing/` |
| `box_spectrum` | power spectra in the middle and south boxes, with the mock correction | `checks/box_spectrum/` |
| `planck_calibration` | TP against Planck 217 GHz: linear fits and the brightest beam; nothing applied | `checks/planck_calibration/` |

Phase 2, diagnostics (needs the phase 1 products). Its outputs go under
`results/combine_continuum/diagnostics/`, in the subfolders the modules already
use (`uv_coordinates/`, `overlap_tests/`, `full_sampling/`,
`amplitude_alignment/`, `aca_forward_test/`, `dirty_beam/`,
`joint_image_comparison/`, `pointing_coverage/`, and `figures/` for the figure
files); the "Writes" column below is relative to that folder.

| Stage | Does | Writes |
|---|---|---|
| `uv` | uv coverage of the measurement sets | `uv_coordinates/`, `figures/` |
| `overlap` | row cache and matched-cell statistics (old `overlap-tests`) | `overlap_tests/` |
| `full_sampling` | full-sampling profiles | `full_sampling/` |
| `amplitude_alignment` | 7m moved to the 12m PB with the 12m+7m image as sky model; TP through the same PB, with and without its beam; also redraws its own figure; also the image curves of Fig. 2 (c) | `amplitude_alignment/` |
| `figure2` | Paper II Fig. 2, three panels (section 6) | `figures/` |
| `compare_images` | 7M and TM2TM1 pbcor against the MIRIAD true-sky image | `figures/` |
| `dirty_beam` | the delivered 12m+7m PSF and its restoring-beam fit | `dirty_beam/` |
| `aca_forward_test` | sky images through the native 7m PB at measured ACA uv points | `aca_forward_test/` |
| `joint_images` | CASA (archived) against MIRIAD true-sky: central and south cutouts (response-letter figure) | `joint_image_comparison/` |
| `pointing_coverage` | pointing centres on the MIRIAD gain and the CASA pb | `pointing_coverage/` |

A stage that needs another stage's output and does not find it stops with a
message naming that stage.

## 3. Inputs

All declared in `task.toml`; the code checks each key's source.

| Key | Path | Source rule | Used by |
|---|---|---|---|
| `interferometer` | `data/delivered/alma/continuum/IRDC316.8.7MTM2TM1.robust0.5.cont.clean.truesky.pb0.2.fits` | delivered | prepare, figure2, compare_images, joint_images, amplitude_alignment sky model |
| `gain` | `data/delivered/alma/continuum/IRDC316.8.7MTM2TM1.robust0.5.cont.gain.fits` | delivered | prepare, measurement region, pointing_coverage |
| `pb` | `data/delivered/alma/continuum/IRDC316.8.7MTM2TM1.robust0.5.cont.pb.fits` | delivered | extended products |
| `psf` | `data/delivered/alma/continuum/IRDC316.8.7MTM2TM1.robust0.5.cont.psf.fits` (+ `.psf.imfit.log`) | delivered | dirty_beam |
| `aca_pbcor` | `data/delivered/alma/continuum/IRDC316.8.7M.robust0.5.cont.image.pbcor.fits` | delivered | compare_images, aca_forward_test |
| `twelve_m_pbcor` | `data/delivered/alma/continuum/IRDC316.8.TM2TM1.robust0.5.cont.image.pbcor.fits` | delivered | compare_images |
| `casa_image`, `casa_pb` | `data/delivered/alma/continuum/casa-archived/IRDC316.8.7MTM2TM1.cont.image.pbcor.fits`, `….cont.pb.fits` | delivered | joint_images, pointing_coverage |
| `tp_dust`, `tp_dust_ff` | `results/synthesize_continuum/model/{dust,continuum}_224GHz.fits` | `synthesize_continuum` product | all TP use |
| `planck_217`, `planck_co`, `planck_rimo`, `planck_fetch_record` | `data/delivered/archive/planck/hfi_cutouts/` | delivered | planck_calibration |
| `atlasgal_870` | `data/delivered/archive/infrared/atlasgal+planck.fits` | delivered | one Planck fit variant |
| `core_catalogue` | `results/source_catalogue/tables/source_catalogue.ecsv` | task product | compact-source check |
| `[reduction_inputs]` | the four measurement sets in `data/reduction/alma/7MTM2TM1/cont/ms/` (12m combined, C5, C2, 7m) | declared read-only intermediate | uv, overlap, full_sampling, amplitude_alignment, aca_forward_test, pointing_coverage |

- Phase 1 reads Jy/beam images with `read_plane` (a beam is required); the
  gain and pb maps are read with `read_gain` (no unit or beam). Phase 2 reads
  images with `read_mjysr` (Jy/beam converted with its own beam; MJy/sr passed
  through; anything else refused).
- The TP frequency is read from the TP header (`FREQGHZ`), not set in the
  configuration. The old `tp_frequency_ghz = 224.0` and the
  `pointing_coverage` `frequency_ghz = 224.4` are removed; the primary-beam
  frequency for pointing coverage is the same `FREQGHZ` (223.8 GHz, the ALMA
  continuum reference frequency). The MIRIAD header carries only
  `RESTFREQ = 224 GHz` on a velocity axis and is not used for this.
- No input is read from `results/prepare_continuum_images/` (removed, D8) or
  from another task's diagnostics folder.

## 4. Method

### 4.1 Core

For same-grid, same-unit, padded planes `low` (TP) and `high` (12m+7m):

    C(k) = W_low(k) · R(k) · f · FT(low) + [1 − W_low(k)] · FT(high)
    combined = FT⁻¹(C)

k is the uv distance in kλ from `fftfreq` and the pixel size in radians.

### 4.2 Weights

- **J-comb** (Jiao et al. 2022):
  - R(k) = G_high(k) / G_low(k): G_high the elliptical 12m+7m restoring beam,
    G_low the circular 18.1″ TP beam, both 1 at k = 0.
  - W_low = 1 below k₁, 0 above k₂, Butterworth in between
    (b(k) = 1 / (1 + (k/k₀)^(2n)), k₀ = (k₁ + k₂)/2, rescaled to run from 1
    to 0). The paper's printed form diverges at k₀; the README states the
    deviation.
  - k₁ = 5 kλ, k₂ = 9 kλ, n = 5; TP used only where G_low ≥ 0.1.
- **Feather**: W_low = G_low, R = 1, so
  C = f · FT(TP) + (1 − G_low) · FT(12m+7m).
- **TP factor** f = 1.0 for both TP models; it is determined from the
  dust + free-free model only, and the dust-only model takes the same value.

### 4.3 Preparation

- Interferometer: `truesky.pb0.2` in MJy/sr; NaN → 0 for the FFT; its finite
  pixels define the mosaic mask.
- TP: continued past its footprint with an 18″ cosine taper, regridded with a
  cubic spline (order 3) to the interferometer grid; the regridding may change
  its power by at most 2 % in 1–9 kλ (checked in `prepare`).
- Margin: the TP footprint must extend at least 18.1″ beyond the mosaic mask,
  or `prepare` stops.
- Padding: 120″; TP mirrored with a cosine roll-off, 12m+7m zeros.

### 4.4 Settings carried over

From `linear_combine_continuum/task.toml`: `[prepare]` (without
`tp_frequency_ghz` and `gain_maximum`), `[tp_factor]` (without
`psd_matched`), `[jcomb]`, `[measure]`, `[box_spectrum]`,
`[planck_calibration]`, `[absolute_scale] status = "pending"`.
`[combine] fill_outside_mosaic` becomes `false` (section 5).

From `combination_diagnostic/task.toml`: `[uv_plot]`, `[overlap_tests]`,
`[image_comparison_plot]`, `[dirty_beam]` (repointed to the delivered PSF),
`[joint_image_comparison]` (repointed), `[aca_forward_test]` (new sky models,
section 6), `[amplitude_alignment]` (sky model = `interferometer`),
`[pointing_coverage]` (repointed, no frequency), and `[figure1]` renamed
`[figure2]` without `tp_display_scale` (the modules that read `figure1` are
updated). `[beam_match]` and `[pointing_range]` are dropped with their
stages; `[visibility_profile]` keeps only `reference_frequency_ghz` and
`uv_max_klambda`, which kept code reads.

## 5. Products and promotion

`results/combine_continuum/products/`:

| File | Content |
|---|---|
| `jcomb_continuum_224GHz.fits` | 12m+7m + dust + free-free TP (adopted) |
| `jcomb_dust_224GHz.fits` | 12m+7m + dust-only TP |
| `feather_continuum_224GHz.fits`, `feather_dust_224GHz.fits` | cross-check |

- float32, Jy/beam of the 12m+7m restoring beam (0.712″ × 0.606″, PA 47.4°),
  on the 12m+7m grid.
- Blank outside the input mosaic mask (pb ≥ 0.2). This changes the current
  `fill_outside_mosaic = true`, which keeps TP-only pixels outside the mosaic.
- Header: `COMB_*` cards (method, k₁, k₂, order, minimum TP response, TP
  factor, TP flavour, TP frequency from `FREQGHZ`, padding, absolute scale
  "pending") and the SHA-256 of every input.
- `manifest.json` at the output root lists every stage run with its settings
  and input hashes.

Promotion is a separate step, done only on the author's go-ahead after the TP
factor and the absolute scale are decided:

- copy `jcomb_continuum_224GHz.fits` and `jcomb_dust_224GHz.fits` to
  `data/delivered/alma/continuum/` as
  `IRDC316.8.TP7MTM2TM1.robust0.5.cont.jcomb.fits` and
  `IRDC316.8.TP7MTM2TM1.robust0.5.cont.jcomb.dust.fits`;
- add README rows (source run, settings, checksum);
- then repoint the downstream tasks (`core_formation`, `column_density` and
  the others that read a combined image).

The feather products are not promoted.

### 5.1 Extended products (added 2026-11-20 at the author's request)

Next to each product above, the `jcomb` and `feather` stages write an
extended version that continues smoothly into the TP-only area beyond the
mosaic: `products/{jcomb,feather}_{continuum,dust}_224GHz_extended.fits`.

- Interferometer input: the prepared true-sky plane multiplied by a weight in
  the mosaic primary beam, read from the delivered
  `IRDC316.8.7MTM2TM1.robust0.5.cont.pb.fits` (new input `pb`):
  w = 0 for pb ≤ pb₀, w = 1 for pb ≥ pb₁, and
  w = ½[1 − cos(π (pb − pb₀)/(pb₁ − pb₀))] in between, with
  `[combine] extended_taper_pb = [0.2, 0.3]`. The small scales of the
  12m+7m image fade out across that rim instead of stopping at pb 0.2, where
  the true-sky noise is 2.5 times the centre value.
- TP input, weights, TP factor and settings are those of the standard product.
- Blank only outside the TP footprint. Beyond the mosaic the image is the TP
  part of the combination alone (J-comb: W_low · R · f · TP; feather: f · TP).
- Where pb ≥ pb₁ the extended product equals the standard one apart from the
  small large-scale change that the rim weight makes through the FFT; the
  stage records the median extended/standard ratio for pb ≥ 0.35.
- Header: the standard `COMB_*` cards plus `COMB_EXTENDED = T` and
  `COMB_TAPER_PB0`, `COMB_TAPER_PB1`; `COMB_FILL_OUTSIDE_MOSAIC = T`.
- Checks: the zero-spacing check and the TP-only ratio are reported for the
  extended products over the mosaic, the TP-only area and both together.
  The other measurements stay on the standard products.
- The diagnostics (phase 2) keep using the standard, blanked products.
- Promotion of an extended version is decided with the rest of section 5.

## 6. Paper II Fig. 2

The three-panel figure is the only layout; the `three_panel` switch and the
four-panel code are removed. Output:
`diagnostics/figures/figure2_combination_diagnostics.{pdf,svg,png}` with its
manifest and source CSVs.

- Panel (a): pointings. Unchanged.
- Panel (b): uv coverage with the short-baseline inset. Unchanged.
- Panel (c): amplitudes and visibility counts.
  - Kept: 7m, C2 and C5 RMS in matched uv and orientation cells (7m moved to
    the 12m PB through the 12m+7m sky model); TP through the 12m PB, with its
    beam and with its beam divided out.
  - Added: 12m+7m (true-sky), J-comb and feather (dust + free-free).

Image curves use the TP-curve method (`continuous_rms` in the
`amplitude_alignment` stage, which already computes the TP curves):

- each image is area-binned to the 0.5″ model grid, multiplied by the actual
  12m pointing PB at both basebands (216.74 and 230.95 GHz), with the same
  pointing, band and orientation population as the TP curve, and sampled
  along a continuous radius from 1 kλ;
- the image's restoring beam and the 0.5″ pixel-binning response are divided
  out. The beam is corrected as the circular Gaussian of the geometric mean
  of BMAJ and BMIN; for the 12m+7m beam (0.712″ × 0.606″) the residual
  depends on orientation and stays within 0.5 % at 26 kλ;
- the curve stops at the lowest of three limits: the beam response falling
  to 0.1 (the rule the TP curve uses); half the model-grid Nyquist limit
  (about 103 kλ); and the interpolation limit. The sampling patch is
  interpolated linearly from the binned grid, which is not corrected: a
  pointing half a pixel off the grid multiplies each axis by |cos(π u p)|.
  The curve stops where that falls to 1 − `image_curve_max_interpolation_loss`
  (0.02), about 26 kλ at 0.5″. In practice this limit applies, so the image
  curves cover 1–26 kλ.

Expected: each image sample is the image multiplied by the 12m PB, whose uv
footprint is a Gaussian of FWHM 4 ln2 / (π θ) for a PB of FWHM θ: with
`pb_12m_fwhm_arcsec = 26.0` at 225 GHz, θ = 26.99″ and 25.33″ at 216.74 and
230.95 GHz, and the footprint FWHM is 6.74 and 7.19 kλ. Below about 5 kλ each
sample therefore mixes in power from 6–10 kλ, so the 12m+7m curve does not fall
away below about 7 kλ: it flattens at roughly half the TP level below about
5 kλ. J-comb and feather follow the TP curve there and meet the 7m points
across 5–9 kλ. The curves are kept from 1 kλ, and the manifests and the README
state this smoothing.

`tp_display_scale` is removed (it is 1.0). The `aca_forward_test` sky models
become 12m+7m true-sky, J-comb, feather and TP; the v3 and CASA-feather
models go.

## 7. Checks

Written by the combination phase as JSON plus figures (PDF + PNG; PNGs
inspected before results are reported):

- **Mock-sky gate** (`jcomb`): the existing known-answer test (unit TP factor,
  interferometer high-passed at 6.6 kλ) must recover flux within 2 % and rms
  within 1 %.
- **Zero spacing**: product smoothed to 18.1″ against the TP; total flux
  inside the mosaic within 5 %; ratio map and median in
  `checks/zero_spacing/`.
- **Overlap continuity**: product amplitude joins TP below k₁ and 12m+7m above
  k₂; combined / 12m+7m within 2 % from k₂ to 20 kλ.
- **Compact sources**: background-subtracted peaks at the catalogue cores,
  product / 12m+7m within 10 %.
- **Negative bowls**: fraction of pixels below −3σ (σ from the noise box);
  the README notes the dust-only negative bowls.
- **TP factor**: `tp_factor`, `box_spectrum` and `planck_calibration` as in
  section 2; none changes f.

Measurements use the region with gain ≥ 0.5 from the delivered gain map.

## 8. Tests

`tests/test_combine_continuum.py` for the task script; the existing module
tests stay.

- Inputs: each key is refused when it comes from the wrong source
  (delivered, `synthesize_continuum`, declared reduction input); the TP
  frequency comes from `FREQGHZ`.
- Products: float32, `COMB_*` cards and input hashes present, blank outside
  the input mask.
- Image amplitude curve: on a Gaussian test image, dividing out the beam
  gives the analytic Fourier amplitude; the curve stops at beam response 0.1.
- Stages: the stage list and order; `--overwrite` protection; `figure2`
  refuses to run without the J-comb and feather products.
- Tests of `linear_combine_continuum` and `combination_diagnostic` that cover
  kept code move here; tests of dropped stages go with them.

`PYTHONPATH=src conda run -n py312 python -m pytest` runs before work is
reported as done.

## 9. Retiring the old tasks

After the first full run of `combine_continuum`:

1. Compare with `results/linear_combine_continuum/`: J-comb inside gain = 1
   should agree within about 1 % (differences come from the true-sky rim and
   the +0.27 % free-free change). The factor numbers (overlap, boxes, Planck)
   are tabulated old against new. If the agreement is worse, stop and report.
2. Compare Fig. 2 and the kept diagnostic figures with
   `results/combination_diagnostic/`.
3. Only after the author has seen the comparison: remove
   `tasks/linear_combine_continuum/`, `tasks/combination_diagnostic/`,
   `scripts/linear_combine_continuum.py`, `scripts/combination_diagnostic.py`
   and their tests, and move both old results folders to the Trash.
   Update the root `README.md`, `docs/project_reorganization_plan.md` and
   `tasks/paper_ii_referee_response/task.toml` references.

## 10. Resources

- J-comb and feather: padded FFT grid of about 50 M pixels, complex64, about
  2 GB peak memory, a few minutes per product; about 150 MB per product image.
- Diagnostic phase: 10–20 minutes, run in the background; the measurement sets
  are only read.

## 11. Out of scope

- The absolute flux scale against Planck (`[absolute_scale]` stays pending).
- Paper II text and figure updates for the 223.8 GHz and 250 GHz changes.
- The CASA-pb inputs of `core_formation` and `volume_density_mapper` (fixed
  at their own reviews).
- Line cubes.

## 12. Dependencies

numpy, scipy, astropy, reproject, matplotlib and `casatools` (`table`, to
read the measurement sets), all already used. No new dependencies.
