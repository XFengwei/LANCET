# Design: reproducible continuum combination (`combine_continuum`)

Date: 2026-11-15. Status: design approved in conversation; awaiting review of
this written spec before an implementation plan is written.

## 1. Purpose

Paper II uses a 12m+7m+TP 1.3 mm continuum image (`allscale_combined_v3_*`)
produced with the J-comb algorithm by Sihan Jiao; the scripts were not shared.
This task replaces those products with a reproducible combination built inside
LANCET from declared inputs, with documented settings and validation.

Decisions taken with the author:

- Method: J-comb (Jiao et al. 2022, Sci. China PMA 65, 299511; arXiv:2208.00588)
  is the adopted product, so the Paper II Methods text keeps naming J-comb.
  Feather is computed by the same code as a cross-check, and CASA `feather`
  anchors that code.
- Scope: continuum only (dust and dust + free-free TP models). The core works
  on any 2-D plane so that a later cube task can reuse it.
- Parameters are fixed before results are seen: transition thresholds from the
  measured uv overlap, TP flux factor 1.0. A small sensitivity grid is reported
  as a systematic, never used for tuning.
- Architecture: one Python uv-plane core in which J-comb and feather differ only
  in their weights (approach A).

Facts this design relies on (established 2026-11-14/15, R4c diagnostics):

- The "230 GHz" TP models are evaluated at ~224.4 GHz (dust; free-free at
  230 GHz, 0.3 % effect), the ALMA continuum frequency. No frequency scaling
  is applied.
- The TP flux scale agrees with Planck 217 GHz: TP = 1.15 × Planck − 4.7 MJy/sr
  at the 5′ beam; the intercept is line-of-sight fore/background and is not
  added to the combination
  (`results/paper_ii_referee_response/diagnostics/R4c/2026-11-15/report.md`).
- The ×2.8 TP factor in v3 came from mistaking the TP beam's uv attenuation
  for a calibration error. The TP beam (18.1″) transmits 50 % at 5.0 kλ and
  10 % at 9.2 kλ; the ACA–TP overlap measured in `combination_diagnostic` is
  4.8–9.2 kλ.

## 2. Task layout

- `tasks/combine_continuum/task.toml`, `tasks/combine_continuum/README.md`
- `scripts/combine_continuum.py` with `--config`, `--stage <name>`,
  `--overwrite`; existing outputs are never replaced without `--overwrite`.
- `src/lancet/combination.py` (calculations), `tests/test_combination.py`.
- Outputs under `results/combine_continuum/` only.

### Inputs (read-only, declared in `[inputs]`)

| Key | Path | Notes |
|---|---|---|
| `interferometer` | `data/delivered/alma/7MTM2TM1/cont/miriad/IRDC316.8.7MTM2TM1.robust0.5.clean.fits` | MIRIAD 12m+7m joint mosaic, PB-corrected, 0.712″×0.606″ PA 47.4°, 0.1″ px, 7096×5127. Byte-identical to Sihan's `interferometer_dust_ff.fits`. |
| `interferometer_attenuated` | `data/delivered/alma/7MTM2TM1/cont/miriad/IRDC316.8.7MTM2TM1.robust0.5.clean.pbant.fits` | Same image multiplied by the mosaic gain; gain = pbant / clean. |
| `tp_dust` | `data/delivered/sedfitter/Dustmodel_230GHz.fits` | Dust model, 18.1″, 6″ px, Jy/beam, ~224 GHz; byte-identical to the Paper I model `data/delivered/paper1/TP_SEDsynthetic230GHz_cont.fits`. |
| `tp_dust_ff` | `data/delivered/sedfitter/SEDsynthetic_Dust_FF_230GHz.fits` | Dust + free-free model, same grid. |
| `v3_dust`, `v3_dust_ff` | `data/delivered/alma/TP7MTM2TM1/cont/allscale_combined_v3_dust{,_ff}.fits` | Comparison only. |
| `core_catalogue` | `results/source_catalogue/tables/source_catalogue.ecsv` | Core positions for the compact-source check (an existing task product, read-only). |

### Stages

1. `inspect` — writes nothing. Reports grids, units, beams, BUNIT, frequency
   keywords and input hashes. Derives the gain map and checks it (0 < gain ≤ 1.01 where
   defined, smooth, and zero exactly where the clean image is blank). Measures the margin between the ALMA
   mosaic (gain > 0) and the TP footprint edge; if it is below one TP FWHM
   (18″), the stage fails and the author decides, rather than the code
   inventing a fill. Checks that the spline-regridded TP keeps its azimuthal
   power within 2 % of the native TP out to k₂.
2. `prepare` — puts each TP on the 0.1″ interferometer grid, converts both
   images to MJy/sr, applies the TP edge taper, writes the prepared planes and
   the gain map (float32) under `prepared/`.
3. `combine` — J-comb (adopted) and Python feather for both TP flavours.
4. `reproduce_v3` — attempts to reproduce the delivered v3 images (see 4.C).
5. `casa_feather` — *dropped on 2026-11-16 (CASA 6.7.2 did not reproduce a textbook feather on this grid; see the task README)*. CASA `feather` on the prepared inputs
   (`casa --nologger --nogui -c`), compared with the Python feather.
6. `sensitivity` — the parameter grid; writes measurements only.
7. `validate` — measurements and figures.
8. `all` — stages 1–7 in order.

### Outputs

- `products/jcomb_dust.fits`, `products/jcomb_dust_ff.fits` (adopted),
  `products/feather_dust.fits`, `products/feather_dust_ff.fits`:
  float32, Jy/beam of the 12m+7m restoring beam, on the 12m+7m grid, NaN
  outside the mosaic (gain = 0), as v3. Header records method, k₁, k₂, n,
  TP factor, TP frequency, working units, and SHA-256 of every input.
- `prepared/`, `casa/`, `validation/` (CSV + JSON with methods text),
  `figures/` (PDF + PNG), `manifest.json`.

Promotion to `data/delivered/alma/TP7MTM2TM1/cont/` is a separate,
author-approved step, under new names; v3 files are never modified.

### Resources

Padded FFT grid ≈ 50 M pixels; complex64 working arrays, ≈ 2 GB peak memory;
≈ 150 MB per product image; ≈ 1.5 GB for all products, prepared planes and
CASA images. Sensitivity runs keep measurements only.

## 3. Algorithm

### 3.1 Core

`combine_plane(low, high, pixel_arcsec, weight_low, beam_ratio, tp_factor)`:

    C(k) = W_low(k) · R(k) · f · FT(low) + [1 − W_low(k)] · FT(high)
    combined = FT⁻¹(C)

`low` (TP) and `high` (12m+7m) are same-grid, same-unit, padded planes; k is
the uv distance in kλ computed from `fftfreq` with the pixel size in radians,
in the unshifted FFT layout; real FFTs in complex64.

### 3.2 Weights

- **J-comb** (Jiao et al. 2022, Eqs. 17–19):
  - R(k) = G_high(k) / G_low(k), with G_high the elliptical 12m+7m restoring
    beam (bmaj, bmin, bpa) and G_low the circular 18.1″ TP beam, both
    normalised to 1 at k = 0. R is evaluated only where W_low > 0.
  - W_low = 1 for k ≤ k₁, 0 for k ≥ k₂, and in between
    W_low = [b(k) − b(k₂)] / [b(k₁) − b(k₂)], b(k) = 1 / (1 + (k/k₀)^(2n)),
    k₀ = (k₁ + k₂)/2.
  - Deviation: the paper prints f(k) = 1/(1 − (k₀/k)^(2n)), which diverges at
    k₀; the authors' code uses the well-behaved Butterworth above. We use the
    latter and state this in the README.
  - Settings: k₁ = 5 kλ, k₂ = 9 kλ, n = 5. Guard: the code refuses
    k₂ where G_low(k₂) < 0.1 (maximum TP amplification 10×).
  - In the noise-free limit the combined beam is exactly G_high.
- **Feather**: W_low = G_low, R = 1 (the form of CASA feather with
  sdfactor = f). The combined beam is not Gaussian; validation measures it.
- **TP factor** f multiplies the TP; default 1.0. The interferometer is taken
  as correctly calibrated.

### 3.3 Preparation

- Units: Jy/beam → MJy/sr with Ω = π·bmaj·bmin/(4 ln 2) for each image; the
  product is converted back to Jy per 12m+7m beam.
- TP regrid: cubic spline (order 3) from 6″ to 0.1″ pixels in the
  interferometer's WCS; the power check in `inspect` guards against
  interpolation smoothing.
- TP edges: outside its footprint the TP is set to zero behind a cosine taper
  of width 18″ starting at the footprint edge (the `inspect` margin check
  keeps the taper outside the mosaic).
- Interferometer: NaN outside the mosaic → 0 for the FFT; the product is
  blanked there afterwards.
- Padding: ≥ 1′ of zeros on each side, rounded up to a fast FFT size
  (`scipy.fft.next_fast_len`).
- Primary-beam frame: J-comb and the Python feather combine the PB-corrected
  images, as the Paper II Methods state. The CASA check uses CASA's
  convention: TP × gain feathered with the attenuated (`pbant`) image, then
  divided by the gain. A difference between the two frames inside gain ≥ 0.5
  is reported.

## 4. Validation

### A. Gates (must pass before the product is used)

1. Synthetic sky on the real grid. Mock sky = upsampled TP model plus
   compact Gaussian sources. Mock TP = sky ⊛ 18.1″; mock interferometer = sky
   ⊛ 12m+7m beam, high-pass filtered by a Gaussian with half power at the
   ACA's shortest spacing (≈ 6.6 kλ). J-comb must return sky ⊛ 12m+7m beam
   with total flux within 2 % and residual rms < 1 % of the peak, inside
   gain ≥ 0.5. This exercises regridding, taper and padding together.
2. Python feather vs CASA feather: median ratio within 1 %, rms difference
   < 2 %, inside gain ≥ 0.5.

### B. Measurements on the adopted J-comb products (reported, flagged if exceeded)

- Large scales: product smoothed to 18.1″ / TP: median within 1 ± 0.05 inside
  gain ≥ 0.5. Through the TP this also places the product on the Planck scale.
- uv amplitude: azimuthal |FT| of the product joins the TP (×R) below k₁ and
  the interferometer above k₂, with no step > 10 % across k₁–k₂.
- Compact sources: at catalogued core positions, the product peak equals the
  12m+7m peak plus the smooth TP level there, within 10 %.
- Negative bowls: fraction of pixels below −3σ (σ from an emission-free box)
  in the product and in the 12m+7m image.
- Effective beam: J-comb is G_high by construction; for feather, the fraction
  of beam volume outside the 12m+7m Gaussian.

### C. Context

- `reproduce_v3`: run our core with f = 2.8 and each kernel available in the
  J-comb code (hard step, step with G_high/G_low, feather), thresholds from
  the Paper I overlap (18.2″–28″) and from the code's defaults. Report which,
  if any, reproduces v3 (median within 2 %, rms within 5 %, gain ≥ 0.5).
  If none does, say so; this is not a gate for the new product.
- v3 vs new J-comb: ratio maps at 2″ and 18″ (expected ≈ 2.8 at large scales,
  ≈ 1 in compact structure).
- Sensitivity grid: (k₁, k₂) ∈ {(4, 8), (5, 9), (6, 9)} kλ × f ∈ {1.0, 1.15};
  B measurements per grid point as the combination systematic.

### Figures

- Summary: maps (12m+7m, TP, J-comb), uv amplitude profile with k₁/k₂,
  large-scale ratio map, sensitivity spread.
- v3 comparison.

Both as PDF and PNG; the PNGs are inspected before results are reported.

## 5. Unit tests (`tests/test_combination.py`)

Small synthetic arrays, no data files:

- Taper: W_low = 1 below k₁, 0 above k₂, monotonic and continuous; invalid
  k₁ ≥ k₂ or n < 1 rejected.
- Beam transfer: elliptical Gaussian transfer matches the analytic form,
  including the position angle; the k grid matches `numpy.fft.fftfreq`.
- Fully sampled sky: low = sky ⊛ B_low and high = sky ⊛ B_high give
  sky ⊛ B_high to numerical precision for any taper; feather keeps the TP's
  total flux.
- Guard: G_low(k₂) < 0.1 is refused.
- TP factor scales only the low-k part.
- Jy/beam → MJy/sr → Jy/beam round-trip.
- Blank regions stay blank; an edge source does not wrap around the padding.
- Products are float32 and carry the settings in the header.

`conda run -n py312 python -m pytest` runs before work is reported as done.

## 6. Documentation and paper handoff

- `tasks/combine_continuum/README.md`: run command, inputs, method (with the
  Butterworth deviation), outputs, limitations: the CLEAN image mixes restored
  components and residuals with different beams; the TP is an SED model, not
  an observation; the Planck background is deliberately excluded; the
  thresholds are fixed a priori.
- Root `README.md` lists the task.
- Handoff, each step only on the author's go-ahead:
  1. Promote the J-comb products to `data/delivered/alma/TP7MTM2TM1/cont/`
     under new names, with provenance; v3 is kept.
  2. Repoint `prepare_continuum_images`, `column_density`, `c18o_moments`,
     `multiwavelength_strip` and `plot_allscale_continuum`, then rerun them.
  3. Paper II: one Methods sentence on the Planck check of the TP flux scale,
     plus the thresholds and TP factor; the v3 comparison goes to the referee
     response if needed.

## 7. Out of scope

- Line cubes (the core is reusable).
- `synthesize_continuum` declaring 230 GHz (its 250 GHz SIMBA model is 9 %
  faint).
- The 1.6× scale of the ATLASGAL+Planck 870 µm map.
- Renaming the `Figure1_*` files.
- Adding the Planck background to any product.

## 8. Dependencies

numpy, scipy (`scipy.fft`, `scipy.ndimage`), astropy, reproject and
matplotlib, all already used in LANCET; CASA 6.7.2 for the `casa_feather`
stage only. No new dependencies.
