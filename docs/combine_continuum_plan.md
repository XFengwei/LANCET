# Continuum combination (`combine_continuum`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A reproducible LANCET task that combines the MIRIAD 12m+7m continuum image with the SED-based TP models using J-comb (adopted) and feather (cross-check), validates the result, and replaces the unshared `allscale_combined_v3_*` products.

**Architecture:** One uv-plane core (`combine_plane`) in which J-comb and feather differ only in their weight functions. Preparation (units, regridding, edge taper, padding) and validation measurements live in separate focused modules. A task script runs the stages `inspect`, `prepare`, `combine`, `reproduce_v3`, `casa_feather`, `sensitivity`, `validate`, `all`, with every setting in `tasks/combine_continuum/task.toml`.

**Tech Stack:** Python 3.12 (conda env `py312`), numpy, scipy (`scipy.fft`, `scipy.ndimage`), astropy, matplotlib, pytest; CASA 6.7.2 for one stage.

**Spec:** `docs/combine_continuum_design.md` (read it first; this plan argues from it).

## Global Constraints

- Run everything from the repository root: `PYTHONPATH=src conda run -n py312 python ...`; tests with `PYTHONPATH=src conda run -n py312 python -m pytest`.
- Known pre-existing failures to ignore (exactly these three): `test_imaging.py::test_generated_floating_point_fits_use_float32_storage`, `test_interferometry.py::test_array_comparison_inspection_is_read_only`, `test_task_independence.py::test_active_tasks_do_not_name_exploratory_source_documents`. Any other failure is yours.
- No new dependencies. No paths or scientific parameters hard-coded in code; all come from `task.toml`.
- Analyses never write inside `data/`; all outputs go under `results/combine_continuum/` via `ProjectPaths.validate_output()`. Delivered files are never modified.
- Existing outputs are not replaced without `--overwrite`.
- FITS products are float32 (`BITPIX = -32`).
- Figures are written as PDF and PNG; inspect the PNG before reporting.
- Task files must not name exploratory notebooks or diagnostic documents.
- Commits: none unless the author asks. Each task ends with a checkpoint (tests pass, files listed) instead of a commit.
- Fixed settings (spec §3.2, §4): `k1_klambda = 5.0`, `k2_klambda = 9.0`, `butterworth_order = 5`, `tp_factor = 1.0`, `minimum_tp_response = 0.1`, `tp_fwhm_arcsec = 18.1`, `edge_taper_arcsec = 18.0`, `minimum_margin_arcsec = 18.1`, `padding_arcsec = 60.0`, `regrid_power_tolerance = 0.02`, `gain_maximum = 1.01`.
- The TP models are treated as 224 GHz; no frequency scaling anywhere.

## Review Focus

1. MIRIAD FITS files are 3-D/4-D (`NAXIS3 = 1`, `CTYPE3 = 'VELO-LSR'`) with `BUNIT = 'JY/BEAM'` in upper case. Reading must squeeze degenerate axes, use the celestial WCS and accept any case of `Jy/beam`. Test owned by Task 2.
2. The images have odd sizes (7096 × 5127 → padded) and `irfft2` silently returns the wrong shape unless `s=` is given. `combine_plane` on odd shapes must return the input shape and keep the fully-sampled identity. Test owned by Task 1.
3. Beam position-angle convention. `BPA` is measured east of north, and `CDELT1 < 0` (east is −x in pixels). A wrong sign would rotate the beam ratio and leave a quadrupole ripple. Test owned by Task 1.
4. A mosaic pixel outside the TP footprint (TP NaN inside the ALMA region) must stop `inspect` with a message naming the margin, not be filled silently. Test owned by Task 2.
5. Re-running a stage over existing outputs must refuse without `--overwrite`, and a configured output root inside `data/` must be rejected. Test owned by Task 3.

---

### Task 1: uv-plane core

**Files:**
- Create: `src/lancet/combination.py`
- Test: `tests/test_combination.py`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `uv_distance_klambda(shape: tuple[int, int], pixel_arcsec: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]`. Returns `(u, v, k)` in kλ for the **rfft2 layout** (shape `(ny, nx//2+1)`), from `np.fft.fftfreq(ny, d)` and `np.fft.rfftfreq(nx, d)`, where `d` is the pixel size in radians; `u` is along x, `v` along y.
  - `gaussian_transfer(u, v, bmaj_arcsec, bmin_arcsec, bpa_deg, cdelt1_sign=-1) -> np.ndarray`, normalised to 1 at k = 0. Convention: BPA east of north; with `cdelt1_sign = -1`, east is −x.
  - `butterworth_taper(k, k1, k2, order) -> np.ndarray`. Implements spec §3.2: 1 for k ≤ k1, 0 for k ≥ k2, rescaled `1/(1+(k/k0)**(2n))` in between with k0 = (k1+k2)/2. Raises `ValueError` if `k1 >= k2` or `order < 1`.
  - `step_taper(k, k_cut) -> np.ndarray`: 1 for k ≤ k_cut, else 0 (used only by `reproduce_v3`).
  - `jcomb_weights(u, v, low_beam: tuple[float, float, float], high_beam: tuple[float, float, float], k1, k2, order, minimum_tp_response=0.1, taper="butterworth", allow_low_response=False) -> tuple[np.ndarray, np.ndarray]`. Returns `(w_low, ratio)`, where `ratio = G_high/G_low` is set where `w_low > 0` and 0 elsewhere. Raises `ValueError` if `G_low < minimum_tp_response` anywhere `w_low > 0` and `allow_low_response` is False. With `taper="step"`, `k1` is the cut and `k2` is ignored.
  - `feather_weights(u, v, low_beam) -> tuple[np.ndarray, np.ndarray]`. Returns `(G_low, ones)`.
  - `combine_plane(low: np.ndarray, high: np.ndarray, w_low: np.ndarray, ratio: np.ndarray, tp_factor: float = 1.0) -> np.ndarray`. Inputs are finite float arrays of equal shape. Computes C = w_low·ratio·tp_factor·FT(low) + (1−w_low)·FT(high) with `scipy.fft.rfft2`/`irfft2(..., s=shape)` in complex64, and returns float32.

- [ ] **Step 1: Write the failing tests** in `tests/test_combination.py`:
  - `test_uv_grid_matches_fftfreq`: `shape=(64, 81)`, `pixel_arcsec=0.5`. `u[0, :]` equals `np.fft.rfftfreq(81, d)/1e3` and `v[:, 0]` equals `np.fft.fftfreq(64, d)/1e3`, with `d = 0.5/206264.806` (rtol 1e-12).
  - `test_gaussian_transfer_matches_numerical_fft` (Review Focus 3): build an image-plane elliptical Gaussian (bmaj 6″, bmin 3″, BPA 30°, pixel 0.5″, 128×129 grid) with its major axis rotated from +y toward −x by BPA. Its normalised `|rfft2|` matches `gaussian_transfer(..., 6, 3, 30)` to within 1e-3 absolute. The same comparison with BPA −30° fails by more than 0.05, which proves the sign matters.
  - `test_butterworth_taper_limits`: for `k = np.linspace(0, 12, 1201)`, `k1=5`, `k2=9`, `order=5`, the taper is exactly 1 on `k <= 5` and exactly 0 on `k >= 9`, non-increasing, with maximum adjacent jump < 0.02. `butterworth_taper(k, 9, 5, 5)` and `butterworth_taper(k, 5, 9, 0)` raise `ValueError`.
  - `test_guard_refuses_weak_tp_response`: low beam (18.1, 18.1, 0), high beam (0.71, 0.61, 47.4), k1 = 5, k2 = 12 → `ValueError`; with `allow_low_response=True` it returns arrays.
  - `test_fully_sampled_identity_jcomb` (Review Focus 2): random positive sky on a `(255, 301)` grid, pixel 0.5″. Low = sky ⊛ (18.1″) and high = sky ⊛ (2″×1.5″, PA 20°), both convolved in Fourier space with the same `gaussian_transfer`. `combine_plane` with J-comb weights (k1 = 5, k2 = 9) returns high to within 1e-4 relative to `high.max()`. The output shape is `(255, 301)` and the dtype float32.
  - `test_feather_keeps_tp_total_flux`: same inputs with high replaced by a high-passed version (its k = 0 term set to 0); the feather output sum equals the low sum within 1e-4 relative.
  - `test_tp_factor_scales_low_k_only`: J-comb with `tp_factor=2` against 1. The difference image's FT is zero (< 1e-5 relative) where `w_low == 0`, and equals `FT(low)·ratio·w_low` where `w_low > 0`.

- [ ] **Step 2: Run to verify they fail**
  Run: `PYTHONPATH=src conda run -n py312 python -m pytest tests/test_combination.py -v`
  Expected: errors with `ModuleNotFoundError: lancet.combination.uvplane`.

- [ ] **Step 3: Implement the functions above in `src/lancet/combination.py`.**
  Rotated-Gaussian transfer: with σ = FWHM/√(8 ln 2) in radians and θ = BPA, the transfer is
  `exp(-2π²[σ_maj² (u·sinθ·s + v·cosθ)² + σ_min² (u·cosθ − v·sinθ·s)²])` with `s = -cdelt1_sign` (so east is −x when `cdelt1_sign = -1`), and u, v in λ. Adjust only if the Step 1 test says otherwise, and keep the test as the authority.

- [ ] **Step 4: Run to verify they pass** (same command; all 7 pass).

- [ ] **Step 5: Checkpoint.** No commit unless the author asks.

### Task 2: preparation helpers

**Files:**
- Create: `src/lancet/combination_prep.py`
- Test: `tests/test_combination_prep.py`

**Interfaces:**
- Consumes: nothing from Task 1.
- Produces:
  - `read_plane(path) -> tuple[np.ndarray, fits.Header, WCS]`. Squeezes degenerate axes, returns a float32 2-D array, the original header and the celestial WCS. Raises `ValueError` unless `BUNIT.lower() == "jy/beam"` and `BMAJ`/`BMIN` exist.
  - `beam_parameters(header) -> tuple[float, float, float]`: (bmaj″, bmin″, bpa°), with `BPA` defaulting to 0.
  - `jybeam_to_mjysr(data, bmaj_arcsec, bmin_arcsec)` and `mjysr_to_jybeam(data, bmaj_arcsec, bmin_arcsec)`, using Ω = π·bmaj·bmin/(4 ln 2).
  - `derive_gain(clean, attenuated, maximum=1.01) -> np.ndarray`: `attenuated/clean` where clean is finite and nonzero, 0 where clean is blank. Raises `ValueError` if any value lies outside `[0, maximum]` (ignoring pixels where `|clean|` is below 1e-3 of its maximum, which are set by interpolating the median gain of their 3×3 neighbours).
  - `regrid_to(data, source_wcs, target_wcs, target_shape, order=3) -> np.ndarray`. Uses `scipy.ndimage.map_coordinates` (spline order 3, `cval=np.nan`, NaNs set to 0 for interpolation and re-masked from a nearest-neighbour footprint), in blocks of rows to bound memory.
  - `footprint_margin_arcsec(tp_footprint: np.ndarray, mosaic: np.ndarray, pixel_arcsec) -> float`. Returns the minimum distance from any mosaic pixel to the TP footprint edge (`scipy.ndimage.distance_transform_edt` on the footprint), negative if a mosaic pixel lies outside the footprint.
  - `taper_edge(data, footprint, width_arcsec, pixel_arcsec) -> np.ndarray`: inside the footprint unchanged; outside, a cosine falloff over `width_arcsec` from the edge value, then 0.
  - `pad(data, pad_pixels) -> tuple[np.ndarray, tuple[slice, slice]]`: zero-pads to `scipy.fft.next_fast_len` sizes ≥ shape + 2·pad. Returns the padded array and the slices that undo it.
  - `azimuthal_power(data, pixel_arcsec, edges_klambda) -> np.ndarray`: mean |FT|² in annuli.

- [ ] **Step 1: Write the failing tests** in `tests/test_combination_prep.py`:
  - `test_read_plane_miriad_layout` (Review Focus 1): write a tmp FITS of shape `(1, 1, 20, 30)` with `CTYPE3='VELO-LSR'`, `CTYPE4='STOKES'` and `BUNIT='JY/BEAM'`. `read_plane` returns shape `(20, 30)`, dtype float32 and `wcs.naxis == 2`. With `BUNIT='K'` it raises `ValueError`.
  - `test_unit_round_trip`: `mjysr_to_jybeam(jybeam_to_mjysr(x, 0.712, 0.606), 0.712, 0.606)` equals x (rtol 1e-6), and 1 Jy/beam for an 18.1″ beam is 2.6455 MJy/sr (rtol 1e-4; Ω = 3.7799e-8 sr).
  - `test_derive_gain`: clean = 2 everywhere except a NaN block, attenuated = 2·g with g ∈ [0.1, 1]. The result recovers g within 1e-6 and is 0 on the NaN block. An attenuated value of 2.2·clean raises `ValueError`.
  - `test_regrid_preserves_power`: an 18.1″ Gaussian blob sampled at 6″ and regridded to 1″. Its `azimuthal_power` agrees with the analytic Gaussian's within 2% for every annulus up to 9 kλ.
  - `test_margin_negative_when_mosaic_outside` (Review Focus 4): a footprint disc of radius 50 px; a mosaic square fully inside gives a margin > 0; one extending 3 px beyond gives a margin ≈ −3·pixel (± 1 px).
  - `test_taper_edge_shape`: inside the footprint values are unchanged; at `width` beyond the edge, 0; monotonic in between.
  - `test_pad_no_wraparound`: a delta at pixel (0, 0) convolved through `pad` → Gaussian FFT → unpad has no flux on the opposite edge (< 1e-6 of peak) when `pad_pixels` ≥ 3σ of the kernel.

- [ ] **Step 2: Run to verify they fail:** `PYTHONPATH=src conda run -n py312 python -m pytest tests/test_combination_prep.py -v` → `ModuleNotFoundError`.

- [ ] **Step 3: Implement the functions above.**

- [ ] **Step 4: Run to verify they pass** (all 7).

- [ ] **Step 5: Checkpoint.**

### Task 3: task configuration, script skeleton and `inspect`

**Files:**
- Create: `tasks/combine_continuum/task.toml`, `tasks/combine_continuum/README.md` (skeleton: run command, inputs, stages; completed in Task 8), `scripts/combine_continuum.py`
- Test: `tests/test_combine_continuum_script.py`

**Interfaces:**
- Consumes: Task 2 `read_plane`, `beam_parameters`, `derive_gain`, `regrid_to`, `footprint_margin_arcsec`, `azimuthal_power`.
- Produces (in `scripts/combine_continuum.py`, importable by tests):
  - `load_task_config(path) -> dict`
  - `resolve_task_paths(project, config) -> dict[str, Path]`. Requires `interferometer`, `interferometer_attenuated`, `tp_dust`, `tp_dust_ff`, `v3_dust` and `v3_dust_ff` to be under `data/delivered`. `core_catalogue` may be under `results/`. The output root goes through `validate_output`.
  - `guard_outputs(paths: list[Path], overwrite: bool) -> None`: raises `FileExistsError` listing the existing files when `overwrite` is False.
  - `run_inspect(paths, config) -> dict`: returns and prints a report covering shapes, beams, units, input SHA-256, gain range, footprint margin (″) and regrid power deviation. Raises `RuntimeError` if the margin is < `minimum_margin_arcsec` or the power deviation is > `regrid_power_tolerance`.
  - `main(argv=None) -> int` with `--config` (default `tasks/combine_continuum/task.toml`), `--stage {inspect,prepare,combine,reproduce_v3,casa_feather,sensitivity,validate,all}` and `--overwrite`. Unimplemented stages raise `NotImplementedError` until their task lands.

`task.toml` tables:
- `[task] name = "combine_continuum"`.
- `[inputs]`: the six delivered paths from spec §2 and `core_catalogue = "results/source_catalogue/tables/source_catalogue.ecsv"`.
- `[outputs] root = "results/combine_continuum"`.
- `[prepare]`: `tp_fwhm_arcsec = 18.1`, `edge_taper_arcsec = 18.0`, `minimum_margin_arcsec = 18.1`, `padding_arcsec = 60.0`, `regrid_order = 3`, `regrid_power_tolerance = 0.02`, `regrid_check_edges_klambda = [1, 2, 3, 4, 5, 6, 7, 8, 9]`, `gain_maximum = 1.01`.
- `[jcomb]`: `k1_klambda = 5.0`, `k2_klambda = 9.0`, `butterworth_order = 5`, `tp_factor = 1.0`, `minimum_tp_response = 0.1`.
- `[feather]`: `tp_factor = 1.0`.
- `[reproduce_v3]`, `[casa_feather]`, `[sensitivity]` and `[validate]` are filled by later tasks.

- [ ] **Step 1: Write the failing tests:**
  - `test_output_root_inside_data_rejected` (Review Focus 5): a config with `root = "data/delivered/x"` → `resolve_task_paths` raises `DataOutputError`.
  - `test_guard_outputs_refuses_existing` (Review Focus 5): an existing tmp file → `FileExistsError`; `overwrite=True` passes.
  - `test_tp_inputs_must_be_delivered`: `tp_dust` pointing under `results/` → `ValueError`.
  - `test_inspect_fails_on_small_margin`: tiny synthetic FITS set (interferometer 60×60 at 1″ with mosaic NaN border; TP 20×20 at 6″ with a footprint ending 5″ from the mosaic) → `run_inspect` raises `RuntimeError` mentioning `margin`.

- [ ] **Step 2: Run to verify they fail:** `PYTHONPATH=src conda run -n py312 python -m pytest tests/test_combine_continuum_script.py -v`.

- [ ] **Step 3: Implement**, following `scripts/synthesize_continuum.py` for config loading and path resolution.

- [ ] **Step 4: Run to verify they pass.**

- [ ] **Step 5: Run `inspect` on the real data:**
  `PYTHONPATH=src conda run -n py312 python scripts/combine_continuum.py --stage inspect`
  Expected: a report with 12m+7m beam 0.712″ × 0.606″ PA 47.4°, grid 7096 × 5127, TP beam 18.1″, gain max ≤ 1.01, margin ≥ 18.1″, power deviation ≤ 2%. **If the margin or the power check fails, stop and report to the author** (spec §2, stage 1); do not change the thresholds.

- [ ] **Step 6: Checkpoint.**

### Task 4: `prepare` and `combine` stages, product writer

**Files:**
- Modify: `scripts/combine_continuum.py`
- Modify: `src/lancet/combination_prep.py` (add `write_product`, `combine_prepared`)
- Test: `tests/test_combine_continuum_script.py`, `tests/test_combination_prep.py`

**Interfaces:**
- Consumes: Task 1 (`uv_distance_klambda`, `jcomb_weights`, `feather_weights`, `combine_plane`) and Task 2 helpers.
- Produces:
  - `write_product(path, data_jybeam, template_header, settings: dict[str, str | float], input_hashes: dict[str, str]) -> None` in `combination_prep.py`. Writes float32 with the 2-D celestial header from the 12m+7m image (BMAJ/BMIN/BPA and BUNIT `Jy/beam` kept), plus `HIERARCH COMB_*` keywords for method, k1, k2, n, tp_factor, `TPFREQ = 224 GHz` and working unit, and one `HISTORY` line per input hash.
  - `run_prepare(paths, config) -> dict[str, Path]`. Writes `prepared/interferometer_mjysr.fits`, `prepared/tp_dust_mjysr.fits`, `prepared/tp_dust_ff_mjysr.fits` (regridded, tapered, unpadded, float32) and `prepared/gain.fits`.
  - `combine_prepared(low, high, pixel_arcsec, high_beam, low_beam, method: str, settings: dict) -> np.ndarray` in `src/lancet/combination_prep.py`, so that `combination_validation` can import it. It pads, builds weights, calls `combine_plane`, unpads, and converts MJy/sr → Jy/beam of the high beam. Methods:
    - `jcomb` → `jcomb_weights(taper="butterworth")`;
    - `feather` → `feather_weights`;
    - `step` → (`step_taper(k, settings["k1"])`, ratio = 1);
    - `step_ratio` → `jcomb_weights(taper="step", allow_low_response=settings.get("allow_low_response", False))`.

    `settings` keys: `k1`, `k2`, `order`, `tp_factor`, `minimum_tp_response`, `padding_arcsec`. Used by `combine`, `reproduce_v3`, `sensitivity` and `synthetic_gate`.
  - `run_combine(paths, config) -> dict[str, Path]`. Writes `products/{jcomb,feather}_{dust,dust_ff}.fits`, NaN where gain = 0.

- [ ] **Step 1: Write the failing tests:**
  - `test_write_product_header`: the written file has `BITPIX == -32`, `BUNIT == 'Jy/beam'`, `COMB_METHOD == 'jcomb'`, `COMB_K1 == 5.0` and a HISTORY line containing a given hash.
  - `test_combine_prepared_end_to_end_small`: synthetic prepared planes (as Task 1's identity test, on a 200×240 grid with a NaN mosaic border converted to 0). The J-comb output equals the high-resolution sky in Jy/beam within 1e-3 relative inside the mosaic, and is NaN wherever the gain is 0.

- [ ] **Step 2: Run to verify they fail.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run to verify they pass.**
- [ ] **Step 5: Run on the real data.** Expected runtime is under 10 min and disk about 1.2 GB.
  `PYTHONPATH=src conda run -n py312 python scripts/combine_continuum.py --stage prepare`, then `--stage combine`.
  Expected: seven prepared files and four products, all float32. Peak of `jcomb_dust_ff` within a few per cent of the 12m+7m peak (0.0806 Jy/beam) plus the local TP level.
- [ ] **Step 6: Checkpoint.**

### Task 5: validation measurements and the synthetic-sky gate

**Files:**
- Create: `src/lancet/combination_validation.py`
- Test: `tests/test_combination_validation.py`

**Interfaces:**
- Consumes: Task 1 `gaussian_transfer`, `uv_distance_klambda`; Task 2 `azimuthal_power`.
- Produces:
  - `compare_planes(a, b, mask) -> dict`: `median_ratio`, `rms_difference_fraction` (rms(a−b)/rms(b)) and `n`.
  - `smooth_to(data, pixel_arcsec, from_beam, to_fwhm_arcsec) -> np.ndarray` (Fourier-domain Gaussian).
  - `large_scale_ratio(product, tp, pixel_arcsec, product_beam, tp_fwhm, mask) -> dict`: median and p16/p84 of smooth_to(product, tp_fwhm)/tp.
  - `amplitude_profile(data, pixel_arcsec, edges_klambda) -> np.ndarray`: azimuthal mean |FT|.
  - `band_step(profile_product, profile_reference, edges, k1, k2) -> float`: maximum relative jump of product/reference between adjacent bins inside [k1, k2].
  - `negative_fraction(data, sigma, mask, threshold=-3.0) -> float`.
  - `core_peak_check(product, interferometer, tp_smooth_level, positions_pix, box=5) -> dict`: median and max of |product_peak − (interf_peak + tp_level)| / product_peak.
  - `beam_volume_outside(method_transfer, high_transfer) -> float`: the fraction of the effective PSF volume not explained by the high-res Gaussian (`1 − Σ(min(eff, high))/Σ eff` in the image plane).
  - `synthetic_gate(gain, pixel_arcsec, high_beam, low_beam, settings, rng_seed=3168) -> dict`. Builds the spec §4.A.1 mock (upsampled TP-like smooth field plus 30 compact Gaussians; mock interferometer high-passed by a Gaussian with half power at `aca_shortest_klambda`), runs `combine_prepared(method="jcomb")`, and returns `flux_error`, `rms_fraction_of_peak` and `passed`.

- [ ] **Step 1: Write the failing tests:**
  - `test_compare_planes_identity_and_scale`: a vs a gives ratio 1, rms 0; 1.1·a vs a gives ratio 1.1.
  - `test_large_scale_ratio_of_scaled_tp`: product = 1.03·smooth sky at high res and tp = sky smoothed to 18.1″ gives a median of 1.03 ± 0.005.
  - `test_band_step_detects_jump`: a reference profile of ones and a product with a 15% jump at 7 kλ gives a step ≥ 0.15; a smooth product gives < 0.01.
  - `test_beam_volume_outside_jcomb_zero`: an effective transfer equal to the high transfer gives 0 (< 1e-6).
  - `test_synthetic_gate_passes_on_small_grid`: a 400×400 grid, gain = 1 inside a disc. `passed` is True, `flux_error < 0.02` and `rms_fraction_of_peak < 0.01`.

- [ ] **Step 2: Run to verify they fail.**
- [ ] **Step 3: Implement.** Add `[validate]` to `task.toml`: `aca_shortest_klambda = 6.6`, `gate_flux_tolerance = 0.02`, `gate_rms_tolerance = 0.01`, `gate_min_gain = 0.5`, `large_scale_tolerance = 0.05`, `band_step_tolerance = 0.10`, `core_peak_tolerance = 0.10`, `negative_sigma = -3.0`, `noise_box_pix = [x0, y0, x1, y1]` (from an emission-free region identified during `inspect`, recorded with its reason), `profile_edges_klambda = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30, 50, 100, 200, 400]`.
- [ ] **Step 4: Run to verify they pass.**
- [ ] **Step 5: Checkpoint.**

### Task 6: `casa_feather` stage

**Files:**
- Modify: `scripts/combine_continuum.py`
- Create: `src/lancet/casa_feather.py`
- Test: `tests/test_casa_feather.py`

**Interfaces:**
- Consumes: Task 4 prepared files; Task 5 `compare_planes`.
- Produces:
  - `casa_feather_script(high_attenuated_fits, low_fits, gain_fits, out_fits, sdfactor: float, workdir) -> str`. Returns CASA Python text that:
    1. `importfits` all three inputs;
    2. multiplies the low image by the gain with `immath`;
    3. runs `feather(imagename, highres, lowres, sdfactor=sdfactor)`;
    4. divides by the gain where gain ≥ `minimum_gain`;
    5. `exportfits` with `bitpix=-32`.
  - `run_casa_feather(paths, config) -> dict`. Writes the script to `casa/run_feather.py` and runs `[casa_executable, "--nologger", "--nogui", "-c", script]` via `subprocess.run(check=True)`, with the log to `casa/casa.log`. Compares `casa/feather_{dust,dust_ff}.fits` with the Python feather products using `compare_planes` where gain ≥ 0.5, and writes `validation/casa_vs_python_feather.json`.
  - `[casa_feather]` in `task.toml`: `casa_executable = "/Applications/CASA-6.7.2.app/Contents/MacOS/casa"`, `minimum_gain = 0.2`, `median_tolerance = 0.01`, `rms_tolerance = 0.02`.

- [ ] **Step 1: Write the failing test** `test_casa_script_contents`: the returned text contains `feather(`, `sdfactor=1.0`, `bitpix=-32`, an `immath` multiply by the gain image, and no path under `data/`.
- [ ] **Step 2: Run to verify it fails.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run to verify it passes.**
- [ ] **Step 5: Run on the real data** (expected a few minutes; CASA images about 0.6 GB in `casa/`):
  `PYTHONPATH=src conda run -n py312 python scripts/combine_continuum.py --stage casa_feather`
  Expected: the JSON shows median ratio within 1% and rms within 2%. **If not, stop and report to the author with both numbers.** A mismatch means the Python core and CASA disagree, and that must be understood before J-comb results are used.
- [ ] **Step 6: Checkpoint.**

### Task 7: `reproduce_v3` and `sensitivity` stages

**Files:**
- Modify: `scripts/combine_continuum.py`, `tasks/combine_continuum/task.toml`
- Test: `tests/test_combine_continuum_script.py`

**Interfaces:**
- Consumes: Task 4 `combine_prepared`; Task 5 `compare_planes`, `large_scale_ratio`, `amplitude_profile`, `band_step`, `negative_fraction`, `core_peak_check`.
- Produces:
  - `expand_variants(section: dict) -> list[dict]`: the explicit list of `{name, method, k1, k2, tp_factor}` from the TOML tables below.
  - `run_reproduce_v3(paths, config) -> Path`. Runs each variant for `dust_ff` and compares it with `v3_dust_ff` where gain ≥ 0.5. Writes `validation/reproduce_v3.csv` (one row per variant: `median_ratio`, `rms_difference_fraction`, `matched`) and states plainly in `validation/reproduce_v3.json` whether any variant matched (median within 0.02 and rms within 0.05). No full images are kept.
  - `run_sensitivity(paths, config) -> Path`. Writes `validation/sensitivity.csv` with the spec §4.B measurements for each grid point; no images.

TOML:
- `[reproduce_v3]`:
  - `tp_factor = 2.8`;
  - `variants`:
    - `{method = "step", k_cut = [7.37, 8.06, 11.33]}` (`k_cut` becomes `settings["k1"]`);
    - `{method = "step_ratio", k_cut = [7.37, 8.06, 11.33], allow_low_response = true}`;
    - `{method = "feather"}`.
  - Where the cuts come from: 7.37 and 11.33 kλ are the Paper I overlap scales 28″ and 18.2″; 8.06 kλ is the J-comb code's `step` cut at √2 × 18.1″.
  - Also `median_tolerance = 0.02` and `rms_tolerance = 0.05`.
- `[sensitivity]`: `thresholds_klambda = [[4.0, 8.0], [5.0, 9.0], [6.0, 9.0]]`, `tp_factors = [1.0, 1.15]`.

- [ ] **Step 1: Write the failing test** `test_expand_variants_counts`: the `[reproduce_v3]` table above expands to 7 variants, and `[sensitivity]` to 6 grid points, each carrying k1/k2/tp_factor.
- [ ] **Step 2: Run to verify it fails.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run to verify it passes.**
- [ ] **Step 5: Run on the real data:** `--stage reproduce_v3`, then `--stage sensitivity`. Expected under 30 min and only a few MB of CSV/JSON. Record whether v3 was reproduced; either outcome is acceptable, but report it plainly.
- [ ] **Step 6: Checkpoint.**

### Task 8: `validate` stage, figures, documentation, full run

**Files:**
- Modify: `scripts/combine_continuum.py`, `tasks/combine_continuum/README.md`, `README.md` (task list)
- Create: `src/lancet/combination_figures.py`
- Test: `tests/test_combination_figures.py`

**Interfaces:**
- Consumes: everything above.
- Produces:
  - `run_validate(paths, config) -> Path`. Runs `synthetic_gate` on the real grid (spec §4.A.1) and the §4.B measurements for `jcomb_{dust,dust_ff}`. It also produces the v3 ratio maps at 2″ and 18″ (written as float32 FITS in `validation/`) and writes `validation/summary.json`, with each measurement's value, tolerance and `flag` (true if exceeded), plus a `methods` string.
  - `plot_summary(...)` and `plot_v3_comparison(...)` in `combination_figures.py`, writing `figures/combination_summary.{pdf,png}` and `figures/v3_comparison.{pdf,png}`. Panels follow spec §4 "Figures". Mark k1/k2 on the uv profile, and show the sensitivity range as a band.
  - `run_all` calls the stages in order.
  - Every stage appends to `manifest.json` in the output root: stage name, UTC time, settings used, input SHA-256 and output paths.

- [ ] **Step 1: Write the failing test** `test_plot_summary_writes_both_formats`: small synthetic inputs → both files exist and are non-empty.
- [ ] **Step 2: Run to verify it fails.**
- [ ] **Step 3: Implement.** Then complete `tasks/combine_continuum/README.md`:
  - run command, inputs, method and stages;
  - the Butterworth deviation (spec §3.2);
  - outputs;
  - limitations from spec §6;
  - the promotion step, described as author-approved.

  Add the task to the root `README.md` task list.
- [ ] **Step 4: Run to verify it passes.**
- [ ] **Step 5: Full run and checks:**
  `PYTHONPATH=src conda run -n py312 python scripts/combine_continuum.py --stage validate`
  Expected:
  - `summary.json` shows the synthetic gate `passed: true`;
  - large-scale ratio 1 ± 0.05;
  - band step < 0.10;
  - flags listed honestly.

  Inspect both PNGs for label collisions, clipping and plausibility. Then run the full suite: `PYTHONPATH=src conda run -n py312 python -m pytest`. Only the three known failures are allowed.
- [ ] **Step 6: Report to the author:**
  - gate results, flags, the v3 reproduction outcome and the sensitivity spread;
  - figure paths;
  - disk used under `results/combine_continuum/`.

  Promotion, downstream repointing and paper text wait for the author's go-ahead (spec §6).
