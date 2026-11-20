# combine_continuum Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** One task, `combine_continuum`, that prepares the 12m+7m true-sky
image and the 223.8 GHz TP models, writes the J-comb and feather products with
their checks, and draws the diagnostics including Paper II Fig. 2. It replaces
`linear_combine_continuum` and `combination_diagnostic`.

**Architecture:** `tasks/combine_continuum/run.py` holds the CLI, the
configuration and path resolution, and stage dispatch. The phase 1 stage
runners move from `scripts/linear_combine_continuum.py` into
`src/lancet/combination_stages.py`, and the phase 2 runners from
`scripts/combination_diagnostic.py` into
`src/lancet/combination_diagnostics.py`. The existing calculation modules stay
and are adjusted. The old scripts keep working, untouched, until Task 9
removes them.

**Tech Stack:** Python 3.12 (conda env `py312`), numpy, scipy, astropy,
reproject, matplotlib, `casatools.table`; pytest.

**Spec:** `docs/combine_continuum_design.md`

## Global Constraints

- Run tests with `PYTHONPATH=src conda run -n py312 python -m pytest`. Subprocess tests need `PYTHONPATH=src`.
- FITS products are float32 (`BITPIX = -32`).
- Figures are PDF + PNG (+ SVG where the settings list it). Inspect the PNG before reporting a figure.
- Analyses never write inside `data/`. All outputs go under `results/combine_continuum/`, through `ProjectPaths.validate_output()`.
- Existing outputs are never replaced without `--overwrite`.
- No paths or parameters are hard-coded in code; everything comes from `tasks/combine_continuum/task.toml`.
- No new dependencies.
- Do not edit, move or delete anything under `data/`, `results/linear_combine_continuum/` or `results/combination_diagnostic/` before Task 9, and do Task 9 only after the author approves the Task 8 comparison.
- Commit only when the author asks. The commit steps below are proposals: stage only the files the task names, because other sessions have uncommitted edits in this checkout (`tasks/linear_combine_continuum/task.toml`, `tasks/combination_diagnostic/*`, `scripts/combination_diagnostic.py`, `src/lancet/pointing_coverage.py`, `tests/test_pointing_coverage.py`).
- The phase 2 output root is `results/combine_continuum/diagnostics/`. Modules keep their own subfolder names under it.
- TP frequency: read from the TP header card `FREQGHZ` (223.8). File labels use the rounded value: `224GHz`.

## Review Focus

1. **A TP file without `FREQGHZ`** (for example an old `Dustmodel_*` path put back in the config): `inspect` and `prepare` must stop with a message naming the file and the missing card, not fall back to 224.0. Test in Task 1.
2. **The true-sky input is 4-D (`NAXIS = 4`, 7096 × 5127 × 1 × 1) with `BUNIT = 'JY/BEAM'` in capitals**: every reader must squeeze it and accept the unit in any case. Test in Task 2 with a 4-D synthetic file.
3. **A diagnostic stage run before `jcomb`/`feather`**: `figure2`, `overlap`, `full_sampling`, `amplitude_alignment` and `aca_forward_test` must stop with a message that names `jcomb` or `feather`, before any long measurement-set read. Test in Task 4.
4. **Rerun without `--overwrite`**: every stage refuses before it computes anything. That includes the stages from `combination_diagnostic`, whose modules check outputs at different points. Test in Task 4 (`figure2`) and Task 3 (`jcomb`).
5. **Image curves beyond their valid range**: the 12m+7m, J-comb and feather curves in Fig. 2 (c) must stop at min(beam response 0.1, half the 0.5″ model-grid Nyquist limit). A curve plotted past that would show binning loss as a drop in amplitude. Test in Task 5.

---

### Task 1: Task folder, configuration and path resolution

**Files:**
- Create: `tasks/combine_continuum/task.toml`
- Create: `tasks/combine_continuum/run.py`
- Modify: `src/lancet/imaging.py` (add `read_frequency_ghz`, `frequency_label`)
- Test: `tests/test_combine_continuum.py`, `tests/test_imaging.py`

**Interfaces:**
- Produces, in `src/lancet/imaging.py`:
  - `read_frequency_ghz(path) -> float` returns the `FREQGHZ` card of a FITS file. It raises `ValueError(f"{path}: no FREQGHZ card")` when the card is absent.
  - `frequency_label(frequency_ghz: float) -> str` returns `f"{round(frequency_ghz):d}GHz"`.
- Produces, in `tasks/combine_continuum/run.py`:
  - `DEFAULT_CONFIG = Path("tasks/combine_continuum/task.toml")`.
  - `DELIVERED_INPUTS = ("interferometer", "gain", "psf", "psf_imfit_log", "aca_pbcor", "twelve_m_pbcor", "casa_image", "casa_pb", "planck_217", "planck_co", "planck_rimo", "planck_fetch_record", "atlasgal_870")`.
  - `SYNTHESIS_INPUTS = ("tp_dust", "tp_dust_ff")`.
  - `COMBINE_STAGES = ("inspect", "prepare", "tp_factor", "box_spectrum", "planck_calibration", "jcomb", "feather")`.
  - `DIAGNOSTIC_STAGES = ("uv", "overlap", "full_sampling", "amplitude_alignment", "figure2", "compare_images", "dirty_beam", "aca_forward_test", "joint_images", "pointing_coverage")`.
  - `load_task_config(path: Path) -> dict`.
  - `resolve_task_paths(project: ProjectPaths, config: Mapping) -> dict` returns:
    - every `[inputs]` key at the top level as a `Path`, and the same mapping under `"inputs"`;
    - `"measurement_sets"`: name → `Path`;
    - `"regions"`: name → `Path`;
    - `"products"`: `product_paths(output_root)`;
    - `"output_root"`: `results/combine_continuum`;
    - `"diagnostics_root"`: `output_root / "diagnostics"`.
  - `product_paths(output_root: Path, label: str = "224GHz") -> dict[str, Path]` has keys `jcomb_continuum`, `jcomb_dust`, `feather_continuum`, `feather_dust` → `output_root/"products"/f"{method}_{kind}_{label}.fits"`.
  - `stages_for(selection: str) -> list[str]`: `"combine"`, `"diagnostics"`, `"all"` or one stage name. An unknown name raises `ValueError`.
  - `require_outputs(paths: Sequence[Path], stage: str) -> None` raises `FileNotFoundError` naming the stage and the first missing path.
  - `STAGES: dict[str, Callable[[dict, dict, bool], object]]`. In this task only `"inspect"` is registered; later tasks add the others.
  - `main(argv=None) -> int` with `--config`, `--stage` (choices: the stage names plus `combine`, `diagnostics`, `all`) and `--overwrite`. It runs the stages in order and appends one manifest entry per stage (copy `append_manifest`, `input_hashes`, `sha256`, `write_json`, `output_paths` from `scripts/linear_combine_continuum.py`).

`task.toml` contents:
- `[task] name = "combine_continuum"`; `[outputs] root = "results/combine_continuum"`.
- `[inputs]`: the paths of spec section 3, with these exact keys:
  - `interferometer`, `gain`, `psf`, `psf_imfit_log`, `aca_pbcor`, `twelve_m_pbcor`;
  - `casa_image` = `data/delivered/alma/continuum/casa-archived/IRDC316.8.7MTM2TM1.cont.image.pbcor.fits`;
  - `casa_pb` = `…/casa-archived/IRDC316.8.7MTM2TM1.cont.pb.fits`;
  - `tp_dust`, `tp_dust_ff`;
  - `planck_217`, `planck_co`, `planck_rimo`, `planck_fetch_record` (from `data/delivered/archive/planck/hfi_cutouts/`);
  - `atlasgal_870`, `core_catalogue`.
- `[measurement_sets]`: the four entries of `tasks/combination_diagnostic/task.toml` (`combined_12m`, `configuration_c5`, `configuration_c2`, `aca_7m`).
- `[regions] successful_combinations = "region/success.reg"`.
- Copied verbatim from `tasks/linear_combine_continuum/task.toml`: `[prepare]`, `[tp_factor]`, `[combine]`, `[jcomb]`, `[measure]`, `[box_spectrum]`, `[planck_calibration]`, `[absolute_scale]`, with three deletions and one change:
  - delete `tp_frequency_ghz` and `gain_maximum` from `[prepare]`;
  - delete `psd_matched` from `[tp_factor]`;
  - set `[combine] fill_outside_mosaic = false`.
- Copied from `tasks/combination_diagnostic/task.toml`:
  - unchanged: `[uv_plot]`, `[overlap_tests]`, `[image_comparison_plot]`, `[dirty_beam]`, `[joint_image_comparison]`, `[aca_forward_test]`, `[amplitude_alignment]`, `[pointing_coverage]`;
  - `[visibility_profile]` with only `reference_frequency_ghz` and `uv_max_klambda`;
  - `[figure1]` renamed `[figure2]`, without `tp_display_scale`.
- Tasks 4–6 adjust the copied path keys inside those tables.

- [ ] **Step 1: Write the failing tests**

In `tests/test_imaging.py`:

```python
def test_read_frequency_ghz(tmp_path):
    path = tmp_path / "tp.fits"
    fits.PrimaryHDU(np.zeros((2, 2), np.float32)).writeto(path)
    with pytest.raises(ValueError, match="FREQGHZ"):
        read_frequency_ghz(path)
    fits.setval(path, "FREQGHZ", value=223.8)
    assert read_frequency_ghz(path) == 223.8
    assert frequency_label(223.8) == "224GHz"
```

In `tests/test_combine_continuum.py`:
- Load `run.py` with `spec_from_file_location` as in `tests/test_continuum.py`.
- Copy the `_config(changes)` helper from `tests/test_linear_combine_continuum_script.py`.

Tests:
- `test_delivered_inputs_must_be_delivered`, parametrized over `DELIVERED_INPUTS`: set the key to `"results/synthesize_continuum/model/dust_224GHz.fits"`; expect `ValueError` matching the key.
- `test_tp_inputs_must_come_from_synthesize_continuum`: `tp_dust` = `data/delivered/sedfitter/Dustmodel_223.8GHz.fits` raises `ValueError` matching `"tp_dust"`.
- `test_measurement_sets_must_be_reduction_inputs`: `aca_7m` = `"data/delivered/alma/continuum"` raises `ValueError`.
- `test_config_points_at_decided_inputs`:
  - `inputs.interferometer` ends with `"cont.clean.truesky.pb0.2.fits"`;
  - `tp_dust_ff == "results/synthesize_continuum/model/continuum_224GHz.fits"`;
  - `"tp_frequency_ghz" not in config["prepare"]`;
  - `"psd_matched" not in config["tp_factor"]`;
  - `config["combine"]["fill_outside_mosaic"] is False`;
  - `"tp_display_scale" not in config["figure2"]`.
- `test_resolved_groups`: the real config resolves; `paths["inputs"]["gain"] == paths["gain"]`; `paths["products"]["jcomb_continuum"].name == "jcomb_continuum_224GHz.fits"`; `paths["diagnostics_root"] == paths["output_root"] / "diagnostics"`.
- `test_stage_selection`: `stages_for("combine") == list(COMBINE_STAGES)`, `stages_for("diagnostics") == list(DIAGNOSTIC_STAGES)`, `stages_for("all") == [*COMBINE_STAGES, *DIAGNOSTIC_STAGES]`, `stages_for("jcomb") == ["jcomb"]`, and `stages_for("nope")` raises `ValueError`.
- `test_require_outputs_names_the_stage`: a missing `tmp_path / "x.fits"` raises `FileNotFoundError` matching `"jcomb"`.
- `test_inspect_writes_nothing`: run `[sys.executable, "tasks/combine_continuum/run.py", "--stage", "inspect"]` with `env={**os.environ, "PYTHONPATH": "src"}`. Assert:
  - return code 0;
  - stdout contains `"TP frequency: 223.8 GHz"` and `"truesky.pb0.2"`;
  - the set of files under `results/combine_continuum` is the same before and after (an empty set when the folder is absent).

- [ ] **Step 2: Run them to verify they fail**

Run: `PYTHONPATH=src conda run -n py312 python -m pytest tests/test_combine_continuum.py tests/test_imaging.py -q`
Expected: collection error or failures (`run.py` missing, `read_frequency_ghz` not importable).

- [ ] **Step 3: Implement.**
- Write `task.toml` as specified above, and `read_frequency_ghz` / `frequency_label` in `imaging.py`.
- In `run.py`:
  - `resolve_task_paths` applies the source rules: delivered keys under `project.delivered`; synthesis keys under `results/synthesize_continuum`; measurement sets under `project.reduction` and existing directories; every file input must exist.
  - `inspect` prints, for each image input, its shape, `BUNIT` and beam. For each TP it prints `TP frequency: {read_frequency_ghz(...)} GHz`. It also lists the measurement sets and the output root, and writes nothing.
  - `inspect` does not compute the footprint margin: `prepare` checks it (Task 2), because `inspect` must stay fast.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `PYTHONPATH=src conda run -n py312 python -m pytest tests/test_combine_continuum.py tests/test_imaging.py -q`
Expected: all pass.

- [ ] **Step 5: Commit (on the author's go-ahead)**

```bash
git add tasks/combine_continuum/task.toml tasks/combine_continuum/run.py src/lancet/imaging.py tests/test_combine_continuum.py tests/test_imaging.py
git commit -m "Add combine_continuum task skeleton: config, input provenance, stage selection"
```

---

### Task 2: Phase 1 module and the `prepare` stage on the true-sky input

**Files:**
- Create: `src/lancet/combination_stages.py`
- Modify: `tasks/combine_continuum/run.py` (register `prepare`)
- Test: `tests/test_combination_stages.py`

**Interfaces:**
- Consumes: `read_frequency_ghz`, and `resolve_task_paths`'s flat keys (`interferometer`, `gain`, `tp_dust`, `tp_dust_ff`, `output_root`).
- Produces, in `src/lancet/combination_stages.py`, everything from `scripts/linear_combine_continuum.py` between `load_task_config` and `run_planck_calibration` inclusive, except:
  - `load_task_config`, `_project_path`, `resolve_task_paths`, `main` and the CLI constants are not copied;
  - `interferometer_within_gain` is removed (the input is already blank outside pb 0.2).
- New: `mosaic_gain(gain: np.ndarray, interferometer: np.ndarray) -> np.ndarray` returns the delivered gain where `interferometer` is finite and 0 elsewhere, as float32.
- `run_prepare(paths, config, overwrite=False) -> dict[str, Path]` keeps its outputs (`prepared_paths`). It changes in four ways:
  - it reads `paths["interferometer"]` (true sky) and `paths["gain"]` with `read_plane`;
  - its prepared `gain.fits` is `mosaic_gain(gain, truesky)`;
  - its prepared interferometer is `jybeam_to_mjysr(np.nan_to_num(truesky), *beam[:2])`;
  - the TP metadata `tp_frequency_ghz` is `read_frequency_ghz(paths[f"tp_{flavour}"])`.
- Input hashes use `["interferometer", "gain"]`.
- `measurement_context` takes `clean` as `np.nan_to_num(read_plane(paths["interferometer"])[0])`, without `interferometer_within_gain`.

- [ ] **Step 1: Write the failing tests**

In `tests/test_combination_stages.py`, copy the `_header(n, pixel_arcsec, beam_arcsec)` helper from `tests/test_linear_combine_continuum_script.py`. Tests:
- `test_mosaic_gain_is_zero_where_the_interferometer_is_blank`: with `gain = np.full((4, 4), 0.8)` and an interferometer that is NaN in its first row, the result is 0 in row 0, 0.8 elsewhere, and dtype float32.
- `test_prepare_reads_4d_truesky_and_tp_frequency(tmp_path)`. Synthetic inputs:
  - interferometer: 400 × 400 px at 0.5″, written 4-D with shape (1, 1, 400, 400) and `BUNIT = "JY/BEAM"`. It is finite only in the central 200 × 200 px; the rest is NaN.
  - gain: same grid, all 0.9.
  - two TP models: 61 × 61 px at 6″, 18.1″ beam, `FREQGHZ = 223.8`, smooth Gaussian content.

  Use the config from `tests/test_combine_continuum.py`'s `_config` with `margin_decimation = 1`. Assert:
  - the prepared `gain.fits` is 0 outside the central block and 0.9 inside;
  - the prepared TP header records `tp_frequency_ghz` 223.8, wherever `write_product` stores metadata (read it back the same way `tests/test_linear_combine_continuum_script.py::test_product_metadata_records_padding` does);
  - the prepared interferometer is float32 and 2-D.
- `test_prepare_stops_on_tp_without_freqghz(tmp_path)`: the same inputs, but one TP has no `FREQGHZ`; `ValueError` matching `"FREQGHZ"`.
- `test_prepare_fails_on_small_margin(tmp_path)`: port the existing test of the same name onto the new inputs (interferometer finite to the TP edge); expect `RuntimeError` matching `"margin"`.

- [ ] **Step 2: Run them to verify they fail**

Run: `PYTHONPATH=src conda run -n py312 python -m pytest tests/test_combination_stages.py -q`
Expected: `ModuleNotFoundError: lancet.combination_stages`.

- [ ] **Step 3: Create the module by copying and applying the changes listed under Interfaces.** Register `STAGES["prepare"] = combination_stages.run_prepare` in `run.py`.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `PYTHONPATH=src conda run -n py312 python -m pytest tests/test_combination_stages.py tests/test_combine_continuum.py -q`
Expected: all pass.

- [ ] **Step 5: Commit (on the author's go-ahead)**

```bash
git add src/lancet/combination_stages.py tasks/combine_continuum/run.py tests/test_combination_stages.py
git commit -m "combine_continuum: prepare stage on the MIRIAD true-sky image"
```

---

### Task 3: Phase 1 checks and products

**Files:**
- Modify: `src/lancet/combination_stages.py`, `tasks/combine_continuum/run.py`
- Test: `tests/test_combination_stages.py`

**Interfaces:**
- Consumes: `product_paths`, `frequency_label`, `read_frequency_ghz`, `run_prepare` outputs.
- Produces:
  - `FLAVOUR_KIND = {"dust_ff": "continuum", "dust": "dust"}`.
  - `run_tp_factor`, `run_box_spectrum`, `run_planck_calibration` write under `output_root/"checks"/<stage>/`, with file names unchanged.
  - `run_method(method, paths, config, overwrite=False)` writes:
    - the products at `paths["products"][f"{method}_{FLAVOUR_KIND[flavour]}"]`;
    - measurements and figures under `output_root/"checks"/method/`, with the old names;
    - `{method}_{kind}_{label}_at_tp_beam.fits` and `zero_spacing_ratio_{method}_{kind}_{label}.fits` under `output_root/"checks"/"zero_spacing"/`.
  - `product_metadata(config, method, flavour, settings, tp_frequency_ghz: float) -> dict` takes the frequency as an argument; `run_method` passes `read_frequency_ghz(paths[f"tp_{flavour}"])`.
  - `TP_FACTOR_VARIANTS` keeps only `"prepared"`; `run_tp_factor` drops the `as_determined` measurement and any read of `psd_matched`.
  - `run.py` registers `tp_factor`, `box_spectrum`, `planck_calibration`, `jcomb`, `feather`.

- [ ] **Step 1: Write the failing tests**
- `test_product_metadata_takes_tp_frequency`: `product_metadata(config, "jcomb", "dust_ff", {}, 223.8)["tp_frequency_ghz"] == 223.8`.
- `test_product_blank_outside_mosaic_when_not_filling`: with `gain = [[0, 0.5]]` and footprint all True, `product_blank(gain, footprint, False)` is `[[True, False]]`.
- `test_tp_factor_has_only_the_prepared_variant`: `set(TP_FACTOR_VARIANTS) == {"prepared"}`.
- `test_method_output_layout(tmp_path)`: monkeypatch `measurement_context` and `combine_prepared` so `run_method("jcomb", ...)` runs on 32 × 32 arrays. Assert:
  - `products/jcomb_continuum_224GHz.fits` and `products/jcomb_dust_224GHz.fits` exist, have `BITPIX == -32`, and have `COMB_TP_FREQUENCY_GHZ == 223.8`;
  - `checks/jcomb/measurements.json` exists;
  - `checks/zero_spacing/jcomb_continuum_224GHz_at_tp_beam.fits` exists;
  - a second call without `overwrite` raises `FileExistsError` before the patched `combine_prepared` is called again (count its calls).

  If monkeypatching `measure_product`/`zero_spacing` is also needed to keep the test small, patch them to return minimal dicts/arrays.

- [ ] **Step 2: Run to verify they fail**

Run: `PYTHONPATH=src conda run -n py312 python -m pytest tests/test_combination_stages.py -q`
Expected: failures on the new tests (signature, variant set, paths).

- [ ] **Step 3: Implement the changes listed under Interfaces, and register the five stages in `run.py`.**

- [ ] **Step 4: Run the tests to verify they pass**

Run: `PYTHONPATH=src conda run -n py312 python -m pytest tests/test_combination_stages.py tests/test_combine_continuum.py -q`
Expected: all pass.

- [ ] **Step 5: Commit (on the author's go-ahead)**

```bash
git add src/lancet/combination_stages.py tasks/combine_continuum/run.py tests/test_combination_stages.py
git commit -m "combine_continuum: TP-factor checks, J-comb and feather products"
```

---

### Task 4: Phase 2 runners and their inputs

**Files:**
- Create: `src/lancet/combination_diagnostics.py`
- Modify: `src/lancet/combination_tests.py`, `src/lancet/full_sampling.py`, `src/lancet/amplitude_alignment.py`, `src/lancet/aca_forward_test.py`, `src/lancet/combination_figure.py` (keys only; the figure changes are Task 5), `src/lancet/joint_image_comparison.py`, `src/lancet/dirty_beam.py`, `tasks/combine_continuum/run.py`, `tasks/combine_continuum/task.toml`
- Test: `tests/test_combination_diagnostics.py`; existing module tests adjusted where their fixtures use the old keys.

**Interfaces:**
- Consumes: the `resolve_task_paths` groups (`inputs`, `products`, `measurement_sets`, `regions`, `diagnostics_root`).
- Phase 2 runners receive `diagnostic_paths(paths) -> dict`: a copy of `paths` with `output_root` replaced by `paths["diagnostics_root"]`. It is defined in `run.py`.
- Key changes in the modules (the old scripts are not edited; they stop working for these stages, which Task 9 removes):
  - **The TP input:** `paths["reduction_inputs"]["total_power"]` becomes `paths["inputs"]["tp_dust_ff"]` in `combination_figure`, `amplitude_alignment` and `aca_forward_test`.
  - **Images in `combination_tests.run_overlap_tests` and `full_sampling.build_full_sampling_profiles`:** `paths["images"]["allscale_sky"]` becomes `paths["products"]["jcomb_continuum"]`, and `paths["images"]["twelve_m_pbcor"]` becomes `paths["inputs"]["twelve_m_pbcor"]`. Both call `require_inputs([...], "jcomb")` first; see below.
  - **Configuration key:** `config["figure1"]` becomes `config["figure2"]` in every module that reads it (`full_sampling`, `aca_forward_test`, `combination_figure`).
  - **`combination_tests.binned_sky_model`:** reads through `read_mjysr`, so Jy/beam and 4-D inputs are accepted.
  - **`joint_image_comparison.run_joint_image_comparison(paths, config, *, save_figure, overwrite=False)`:** takes `paths`. It reads `paths["inputs"]["casa_image"]`, `paths["inputs"]["interferometer"]` and `paths["inputs"]["casa_pb"]`, uses `read_mjysr` for the two images, and writes to `paths["output_root"]/"joint_image_comparison"`. Remove the `casa`, `miriad` and `casa_pb` keys from `[joint_image_comparison]`.
  - **`dirty_beam.run_dirty_beam(paths, config, *, save_figure, overwrite=False)`:** takes `paths`. It reads `paths["inputs"]["psf"]`, `paths["inputs"]["psf_imfit_log"]` and, as the paper image, `paths["inputs"]["interferometer"]`, and writes to `paths["output_root"]/"dirty_beam"`. Remove `psf_fits`, `imfit_log` and `paper_image` from `[dirty_beam]`.
  - **Amplitude-alignment sky model and ACA forward test:** `[amplitude_alignment] sky_model = "inputs.interferometer"`. The `[aca_forward_test.sky_models]` table becomes:
    - `interferometer = { source = "inputs.interferometer", label = "12m+7m (no TP)", tp_scale = 0.0, tp_beam_transfer = true }`
    - `jcomb = { source = "products.jcomb_continuum", label = "J-comb", tp_scale = 1.0, tp_beam_transfer = true }`
    - `feather = { source = "products.feather_continuum", label = "Feather", tp_scale = 1.0, tp_beam_transfer = true }`
    - `aca_only = { source = "inputs.aca_pbcor", label = "ACA-only image" }`
  - **TP-only sky model:** add `tp = { source = "inputs.tp_dust_ff", label = "TP model", tp_scale = 1.0 }` only if `run_aca_forward_test` accepts a TP-only model; otherwise leave it out and record a ruling.
- New: `require_inputs(paths: Sequence[Path], stage: str) -> None` in `combination_tests.py` (the lowest module of this group, so there is no import cycle), with the same behaviour as `run.require_outputs`. `full_sampling`, `aca_forward_test`, `amplitude_alignment` and `combination_figure` import it from there; no module imports `run.py`.
- `src/lancet/combination_diagnostics.py` is copied from `scripts/combination_diagnostic.py`: `_prepare_output`, `_figure_outputs`, `_save_figure`, `_uv_cache`, `ARRAY_STYLE`, `_plot_style`, `plot_uv_coverage`, `_recenter_wcs`, `_add_regions`, `plot_image_comparison` and `run_pointing_coverage`. Changes:
  - `plot_image_comparison` reads `paths["inputs"]["aca_pbcor"]`, `["twelve_m_pbcor"]` and `["interferometer"]` with `read_mjysr`, in place of `_read_2d`.
  - `run_pointing_coverage(paths, config, *, overwrite=False)` reads `paths["inputs"]["gain"]` and `["casa_pb"]` and `paths["measurement_sets"]`, writes to `paths["output_root"]/"pointing_coverage"`, and computes the primary beam at `read_frequency_ghz(paths["inputs"]["tp_dust_ff"])`. Remove `gain`, `pb` and `frequency_ghz` from `[pointing_coverage]`.
- `run.py` registers the ten phase 2 stages. Each is called with `diagnostic_paths(paths)`. `amplitude_alignment` runs `run_amplitude_alignment` and then `replot_alignment`.

- [ ] **Step 1: Write the failing tests**

In `tests/test_combination_diagnostics.py`:
- `test_diagnostic_paths_use_the_diagnostics_root`: on the real config, `run.diagnostic_paths(paths)["output_root"] == paths["output_root"] / "diagnostics"`.
- `test_overlap_stops_without_jcomb_products(tmp_path)`: resolved paths with `products` pointing into `tmp_path`, which holds no product. `combination_tests.run_overlap_tests(paths, config, save_figure=..., overwrite=False)` raises `FileNotFoundError` matching `"jcomb"`. Monkeypatch `read_test_rows` to raise `AssertionError`, to prove the check comes first.
- `test_binned_sky_model_accepts_4d_jybeam(tmp_path)`: a 4-D Jy/beam file with a constant 1 Jy/beam and a 1″ beam. The binned result is in MJy/sr and equals `jybeam_to_mjysr(1, 1, 1)` within 1e-6 relative.
- `test_joint_and_dirty_beam_take_paths`: `inspect.signature(run_joint_image_comparison).parameters` begins with `paths`, and so does `run_dirty_beam`'s.
- `test_config_has_no_path_keys_in_stage_tables`: `[joint_image_comparison]`, `[dirty_beam]` and `[pointing_coverage]` contain none of `casa`, `miriad`, `casa_pb`, `psf_fits`, `imfit_log`, `paper_image`, `gain`, `pb`, `frequency_ghz`.
- Run the existing `tests/test_combination_tests.py`, `test_full_sampling.py`, `test_aca_forward_test.py`, `test_dirty_beam.py`, `test_pointing_coverage.py` and `test_combination_figures.py`. Update any fixture that builds `paths` with the old keys to use the new groups. Record each changed assertion in the ledger.

- [ ] **Step 2: Run to verify they fail**

Run: `PYTHONPATH=src conda run -n py312 python -m pytest tests/test_combination_diagnostics.py -q`
Expected: failures (module missing, signatures, config keys).

- [ ] **Step 3: Implement the changes listed under Interfaces.**

- [ ] **Step 4: Run the tests to verify they pass**

Run: `PYTHONPATH=src conda run -n py312 python -m pytest tests/test_combination_diagnostics.py tests/test_combination_tests.py tests/test_full_sampling.py tests/test_aca_forward_test.py tests/test_dirty_beam.py tests/test_pointing_coverage.py tests/test_combine_continuum.py -q`
Expected: all pass.

- [ ] **Step 5: Commit (on the author's go-ahead)**

Stage `src/lancet/pointing_coverage.py` and `tests/test_pointing_coverage.py` only if this task changed them. They carry another session's uncommitted work; if they must change, stop and ask the author first.

```bash
git add src/lancet/combination_diagnostics.py src/lancet/combination_tests.py src/lancet/full_sampling.py src/lancet/amplitude_alignment.py src/lancet/aca_forward_test.py src/lancet/combination_figure.py src/lancet/joint_image_comparison.py src/lancet/dirty_beam.py tasks/combine_continuum/run.py tasks/combine_continuum/task.toml tests/test_combination_diagnostics.py
git commit -m "combine_continuum: diagnostic stages on delivered inputs and the new products"
```

---

### Task 5: Fig. 2 with the image curves

**Files:**
- Modify: `src/lancet/amplitude_alignment.py`, `src/lancet/combination_figure.py`, `tasks/combine_continuum/task.toml`, `tasks/combine_continuum/run.py`
- Test: `tests/test_combination_diagnostics.py`

**Interfaces:**
- Consumes: `continuous_rms`, `sky_model`, `response_limit`, `gaussian_transfer` (all in `amplitude_alignment`/`aca_forward_test`), and the `products` group.
- Produces:
  - `image_curve_limit(beam_fwhm_arcsec: float, minimum_response: float, model_pixel_arcsec: float) -> float` in `amplitude_alignment`: `min(response_limit(beam, minimum_response), 0.5 * 1/(2*model_pixel_arcsec*π/648000)/1000)` in kλ.
  - `pixel_binning_transfer(u_klambda, pixel_arcsec) -> np.ndarray`: `np.sinc(u*1000*pixel_rad)` (numpy's normalised sinc), the response of area-binning to `pixel_arcsec`.
  - No new public curve function. `run_amplitude_alignment` computes the curves:
    - for each entry of the new setting `[amplitude_alignment] image_curves` (`interferometer = "inputs.interferometer"`, `jcomb = "products.jcomb_continuum"`, `feather = "products.feather_continuum"`), it bins the image with `sky_model`;
    - it computes `continuous_rms` through the 12m PB (the `"common"` frame) on its own grid, `np.geomspace(continuous_min_klambda, image_curve_limit(...), image_curve_points)` with the new setting `image_curve_points = 60`;
    - it divides by `gaussian_transfer(grid, geometric-mean beam FWHM)` and `pixel_binning_transfer(grid, model_pixel_arcsec)`.
  - The curves are written to `amplitude_alignment/image_curves.csv` (columns `name, uv_klambda, rms_mjy`); `read_products` also returns them, as a fifth value `image_curves: dict[str, tuple[np.ndarray, np.ndarray]]`. Every caller of `read_products` is updated.
  - `make_combination_figure(paths, config, *, save_figure, overwrite=False)`:
    - the `three_panel` parameter is removed, and the branches for `three_panel=False` are deleted;
    - the stem becomes `"figure2_combination_diagnostics"` (+ `output_suffix`);
    - `tp_display_scale` is removed;
    - panel (c) plots the three image curves (labels "12m+7m", "J-comb", "Feather") in addition to the existing curves, and the image-curve CSV is copied to `figures/figure2_combination_diagnostics_image_curves.csv`;
    - the manifest gets an `"image_curves"` entry with the method text from spec section 6;
    - the function first calls `require_inputs([paths["products"]["jcomb_continuum"], paths["products"]["feather_continuum"]], "jcomb/feather")`.
- Colours: three entries in the order of the module's existing palette (no new hues). The 12m+7m curve is dashed, J-comb solid, feather dotted. Use the `dataviz` skill for the colour check only if the palette needs a new entry.

- [ ] **Step 1: Write the failing tests**
- `test_image_curve_limit_is_the_lower_of_beam_and_nyquist`:
  - `image_curve_limit(0.66, 0.1, 0.5)` equals `0.5 * 206.26…` within 1e-3, because the Nyquist half is lower than the beam limit of about 250 kλ;
  - `image_curve_limit(18.1, 0.1, 0.5)` equals `response_limit(18.1, 0.1)`.
- `test_beam_and_binning_corrected_curve_recovers_a_gaussian(tmp_path)`:
  - Build a 2048 × 2048 px, 0.1″ Jy/beam image holding one circular Gaussian source (FWHM 2″, beam 0.7″, so the restored source has FWHM √(2² + 0.7²)″), at the phase centre of one synthetic pointing pair.
  - Run `sky_model` → `continuous_rms` with one pair, one band and one sector, and a PB large enough to be flat over the source. Divide by `gaussian_transfer(grid, 0.7)` and `pixel_binning_transfer(grid, 0.5)`.
  - Compare with the analytic amplitude of the unrestored 2″ Gaussian times the PB-patch normalisation that `continuous_rms` applies. They must agree within 2 % for u ≤ `image_curve_limit(0.7, 0.1, 0.5)`.
  - Build the minimal `settings` dict from `[overlap_tests]`.
- `test_figure2_requires_products(tmp_path)`: `make_combination_figure` with `products` in an empty `tmp_path` raises `FileNotFoundError` matching `"jcomb"`.
- `test_figure_has_no_three_panel_switch`: `"three_panel" not in inspect.signature(make_combination_figure).parameters`.

- [ ] **Step 2: Run to verify they fail**

Run: `PYTHONPATH=src conda run -n py312 python -m pytest tests/test_combination_diagnostics.py -q -k "image_curve or gaussian or figure"`
Expected: failures (functions missing; `three_panel` still present).

- [ ] **Step 3: Implement the changes listed under Interfaces; add `image_curves` and `image_curve_points` to `[amplitude_alignment]`.**

- [ ] **Step 4: Run the tests to verify they pass**

Run: `PYTHONPATH=src conda run -n py312 python -m pytest tests/test_combination_diagnostics.py tests/test_combination_figures.py -q`
Expected: all pass.

- [ ] **Step 5: Commit (on the author's go-ahead)**

```bash
git add src/lancet/amplitude_alignment.py src/lancet/combination_figure.py tasks/combine_continuum/task.toml tasks/combine_continuum/run.py tests/test_combination_diagnostics.py
git commit -m "combine_continuum: Fig. 2 with 12m+7m, J-comb and feather curves"
```

---

### Task 6: Task README and project references

**Files:**
- Create: `tasks/combine_continuum/README.md`
- Modify: `README.md`, `docs/project_reorganization_plan.md`
- Test: `tests/test_task_independence.py` (existing; must not gain new failures)

**Interfaces:** none.

- [ ] **Step 1: Write the README.** Follow `tasks/synthesize_continuum/README.md` in form. It covers:
  - run commands: `--stage combine`, `--stage diagnostics`, `--stage all`, single stages;
  - the inputs table from spec section 3;
  - the method: J-comb with the Butterworth deviation; feather; TP factor 1.0 and its three checks with the Paper II numbers from `tasks/linear_combine_continuum/README.md`;
  - outputs per stage;
  - Fig. 2 and the image-curve method;
  - limitations: the dust-only negative bowls; the TP is an SED model; CASA 6.7.2 `feather` returns the sum of its inputs; `[absolute_scale]` is pending; the shortest 12m+7m spacing is 6.6 kλ;
  - the promotion rule from spec section 5.

  It must not name exploratory notebooks or documents.
- [ ] **Step 2: Add the task to the task list in the root `README.md`**, and a row to the review log in `docs/project_reorganization_plan.md` ("combine_continuum: merged task created; old tasks retired after comparison").
- [ ] **Step 3: Run the task-independence test**

Run: `PYTHONPATH=src conda run -n py312 python -m pytest tests/test_task_independence.py -q`
Expected: the same single pre-existing failure (`test_active_tasks_do_not_name_exploratory_source_documents`), with no `combine_continuum` file in its message.

- [ ] **Step 4: Commit (on the author's go-ahead)**

```bash
git add tasks/combine_continuum/README.md README.md docs/project_reorganization_plan.md
git commit -m "combine_continuum: README and project references"
```

---

### Task 7: First full run

**Files:** none changed; outputs under `results/combine_continuum/`.

- [ ] **Step 1: Check the disk.** Run `df -h .` and require ≥ 10 GB free: the products are about 0.6 GB, the prepared planes about 0.9 GB, and the diagnostic caches a few GB. If less is free, stop and report.
- [ ] **Step 2: Run phase 1 in the background**

Run (background): `PYTHONPATH=src conda run -n py312 python tasks/combine_continuum/run.py --stage combine > <scratchpad>/combine_phase1.log 2>&1`. The log goes to the session scratchpad, not the repository.
Expected: exit 0; four products in `results/combine_continuum/products/`; `checks/` with `tp_factor`, `box_spectrum`, `planck_calibration`, `jcomb`, `feather`, `zero_spacing`. The mock-sky gate passes.

- [ ] **Step 3: Run phase 2 in the background**

Run (background): `PYTHONPATH=src conda run -n py312 python tasks/combine_continuum/run.py --stage diagnostics > <scratchpad>/combine_phase2.log 2>&1`.
Expected: exit 0; `diagnostics/figures/figure2_combination_diagnostics.{pdf,svg,png}` and the other figures.

- [ ] **Step 4: Inspect every PNG** written by both phases (Read tool). For Fig. 2 (c), check:
  - the 12m+7m curve falls below about 7 kλ;
  - J-comb and feather follow the TP curve below 5 kλ and meet the 7m points across 5–9 kλ.

  Note anything else as a finding; do not tune settings.

---

### Task 8: Comparison with the old results (stop for the author)

**Files:**
- Create: `results/combine_continuum/comparison/old_vs_new.json` and `old_vs_new.md`, written by a scratchpad script; nothing under `src/`.

- [ ] **Step 1: J-comb and feather.** Old products are `results/linear_combine_continuum/{jcomb,feather}/{method}_{dust,dust_ff}.fits`; new ones are `results/combine_continuum/products/`. Compute per pair:
  - the median and 5–95 % range of new/old inside delivered gain = 1, where |old| > 3σ (σ from the noise box);
  - the same in the rim, 0 < gain < 1.

  Expected: inside gain = 1 the median is within 1 % of 1.0. Report the rim separately.
- [ ] **Step 2: Factor checks.** Tabulate old against new:
  - `tp_factor` (prepared variant);
  - `box_spectrum` middle and south, measured and mock-corrected;
  - `planck_calibration`: the fit gain and offset for coverage ≥ 0.9, and the brightest-beam ratios.
- [ ] **Step 3: Figures.** Put old and new Fig. 2 side by side, and the old and new `Interferometry_combine_comparison_vertical`, `CASA_MIRIAD_central_south` and dirty-beam figures. List the visible differences.
- [ ] **Step 4: Stop and report** to the author: the tables, the figure list, and whether the 1 % expectation holds. If it fails, say so plainly and do not continue. Task 9 waits for the author's explicit approval.

---

### Task 9: Retire the old tasks (only after the author approves Task 8)

**Files:**
- Delete (`git rm`): `tasks/linear_combine_continuum/`, `tasks/combination_diagnostic/`, `scripts/linear_combine_continuum.py`, `scripts/combination_diagnostic.py`, `tests/test_linear_combine_continuum_script.py`.
- Delete code used only by dropped stages:
  - `src/lancet/pointing_profiles.py` and `tests/test_pointing_profiles.py`;
  - `run_pb_amplitude_curve` and `_plot_pb_amplitude_curve` in `combination_tests.py`;
  - the parts of `tp_beam_test.py` other than `smooth_to_circular_beam` and `_predict`, which move into `aca_forward_test.py`, then `tp_beam_test.py` itself;
  - the `visibility_diagnostics.py` functions that no kept module imports. Keep `C_METRES_PER_SECOND`, `_table_tool` and `average_visibility_rows`, or move `C_METRES_PER_SECOND` to `combination_tests.py` and update its importers.
- Modify: `tasks/paper_ii_referee_response/task.toml` and every other file that names the old tasks or scripts (find them with `grep -rn "linear_combine_continuum\|combination_diagnostic" --include=*.py --include=*.toml --include=*.md . | grep -v '^./results\|^./data\|reorganization_baseline'`).
- Move to the Trash (not `rm`): `results/linear_combine_continuum/` and `results/combination_diagnostic/`, to `~/.Trash/LANCET_old_combination_results_<date>/`.

- [ ] **Step 1: Write the failing test**

`test_no_reference_to_retired_tasks` in `tests/test_combine_continuum.py`: no file under `tasks/`, `scripts/`, `src/` or `tests/` (except this test file) contains `"linear_combine_continuum"` or `"combination_diagnostic"`.

- [ ] **Step 2: Run to verify it fails**

Run: `PYTHONPATH=src conda run -n py312 python -m pytest tests/test_combine_continuum.py -q -k retired`
Expected: FAIL, listing the old files.

- [ ] **Step 3: Delete, move and repoint as listed. Then grep for imports of every deleted function and fix or remove its last users.**

- [ ] **Step 4: Run the full suite**

Run: `PYTHONPATH=src conda run -n py312 python -m pytest -q`
Expected: no failure beyond the pre-existing ones known at the start of this plan:
- `test_continuum_plotting::test_allscale_plot_inspection_is_read_only`
- `test_core_formation::test_core_formation_inspection_reports_optional_reference`
- `test_interferometry::test_array_comparison_inspection_is_read_only`
- `test_interferometry::test_generated_beam_match_agrees_with_retained_product`
- `test_moments::test_moment_task_inspection_is_read_only`
- `test_moments::test_reprojected_moment2_uses_alma_shape_and_float32`
- `test_task_independence::test_active_tasks_do_not_name_exploratory_source_documents`

- [ ] **Step 5: Commit (on the author's go-ahead)**

```bash
git add -A tasks/combine_continuum tasks/paper_ii_referee_response/task.toml src/lancet tests README.md docs/project_reorganization_plan.md
git commit -m "Retire linear_combine_continuum and combination_diagnostic in favour of combine_continuum"
```

Check `git status` before committing: the `git rm` deletions must be staged, and the other sessions' unrelated files must not be.

---

### Task 10: Extended products (author request, 2026-11-20; runs before Task 9)

Spec section 5.1. The author chose the pb taper and "all four products; diagnostics unchanged".

**Files:**
- Modify: `src/lancet/combination_stages.py`, `tasks/combine_continuum/run.py`, `tasks/combine_continuum/task.toml`, `tasks/combine_continuum/README.md`
- Test: `tests/test_combination_stages.py`, `tests/test_combine_continuum.py`

**Interfaces:**
- Consumes: `run_method`, `combine_prepared`, `measurement_context`, `zero_spacing`, `tp_only_ratio`, `product_metadata`, `product_paths`, `DELIVERED_INPUTS`.
- Produces:
  - `[inputs] pb = "data/delivered/alma/continuum/IRDC316.8.7MTM2TM1.robust0.5.cont.pb.fits"`. It is added to `DELIVERED_INPUTS` and read with `read_gain`, which does not check units; the shape must match the interferometer.
  - `[combine] extended_taper_pb = [0.2, 0.3]`. The existing `fill_outside_mosaic = false` stays and applies to the standard products.
  - `pb_taper(pb: np.ndarray, pb0: float, pb1: float) -> np.ndarray` in `combination_stages.py`: float32, 0 at pb ≤ pb0 and for NaN, 1 at pb ≥ pb1, and ½[1 − cos(π (pb − pb0)/(pb1 − pb0))] in between. It raises `ValueError` unless 0 < pb0 < pb1 ≤ 1.
  - `product_paths(output_root, label)` gains four keys, `jcomb_continuum_extended`, `jcomb_dust_extended`, `feather_continuum_extended` and `feather_dust_extended`. They map to `products/{method}_{kind}_{label}_extended.fits`.
  - `run_method` writes the standard and the extended product for each flavour:
    - The extended product combines `high × pb_taper(pb)` with the same TP plane and settings, and is blank only outside the TP footprint (`product_blank(gain, footprint, True)`).
    - The extended product's header carries the `COMB_*` cards with `fill_outside_mosaic = True`, plus `COMB_EXTENDED = T`, `COMB_TAPER_PB0` and `COMB_TAPER_PB1`.
    - `checks/zero_spacing/` also gets `{method}_{kind}_{label}_extended_at_tp_beam.fits` and `zero_spacing_ratio_{method}_{kind}_{label}_extended.fits`.
    - `measurements.json` gets, per flavour, an `"extended"` entry holding:
      - the zero-spacing check over `mosaic`, `tp_only` and `all`;
      - `tp_only_ratio_to_scaled_tp`;
      - `ratio_to_standard_pb_ge_0.35`: the median of extended/standard where pb ≥ 0.35 and the standard is finite and above 3σ.
    - All outputs, extended included, are guarded by `--overwrite` before any computation.
  - `measure_product` and the other standard measurements are unchanged.

- [ ] **Step 1: Write the failing tests**
  - `test_pb_taper`: with pb = [nan, 0.1, 0.2, 0.25, 0.3, 0.9] and (0.2, 0.3), expect [0, 0, 0, 0.5, 1, 1] within 1e-6. Also check that (0.3, 0.2) and (0, 0.3) raise `ValueError`.
  - `test_extended_product_paths`: `product_paths(root, "224GHz")["jcomb_continuum_extended"].name == "jcomb_continuum_224GHz_extended.fits"`, and the same pattern holds for all four keys.
  - `test_pb_is_a_delivered_input`: parametrize the existing delivered-input test with `pb`, and assert `config["combine"]["extended_taper_pb"] == [0.2, 0.3]`.
  - `test_method_writes_extended_products(tmp_path)`. Extend the monkeypatched `test_method_output_layout` pattern with a synthetic pb whose finite mosaic ends at pb 0.2, with TP extending beyond it. Assert:
    - the four extended FITS exist, have BITPIX −32, and carry `COMB_EXTENDED == True` and `COMB_TAPER_PB0 == 0.2`;
    - the extended product is finite in the TP-only region, where the standard one is NaN;
    - where pb ≥ 0.3 the extended and standard products agree within 1e-3 relative on the synthetic data;
    - `measurements.json` has `dust_ff.extended.zero_spacing.tp_only.n > 0`.
  - `test_extended_outputs_guarded`: a second call without `overwrite` raises `FileExistsError` before any combination. Count `combine_prepared` calls.

- [ ] **Step 2: Run them to verify they fail.**
  Run: `PYTHONPATH=src conda run -n py312 python -m pytest tests/test_combination_stages.py tests/test_combine_continuum.py -q`
  Expected: the new tests fail (`pb_taper` missing, keys missing, no `pb` input).

- [ ] **Step 3: Implement the changes listed under Interfaces.** Then add the extended products, their method and their checks to the task README: the products table, the method section and the outputs.

- [ ] **Step 4: Run the tests, then the full suite.**
  Expected: all new tests pass; the full suite shows only the seven known failures.

- [ ] **Step 5: Commit, staging only the listed files.**
