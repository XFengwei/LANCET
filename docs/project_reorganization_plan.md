# LANCET reorganisation plan (draft, 2026-11-18)

Status: Phase 0 started 2026-11-18. Each later phase needs a go-ahead.

## 1. Goal

The project is built from **task → result** pairs:

- `data/` holds everything that goes in.
- `tasks/<task>/` holds everything that defines one analysis: settings, the runnable script, its README.
- `results/<task>/` holds everything that comes out, and can always be deleted and regenerated.
- `src/` holds only code shared between tasks.

Papers, design, docs, presentations, references and the website stand beside this and do not take part in it.

## 2. Where the project is now (survey of 2026-11-18)

| Item | State |
|---|---|
| Tasks | 24 folders in `tasks/`; 23 have a `task.toml`, 1 (`sa_jepa_continuum_pilot`) uses YAML files |
| Scripts | 37 in `scripts/`; 22 belong to a task, 15 have no task folder |
| Library | 36 modules in `src/lancet/` (8,400 lines); 6 have no test |
| Tests | 37 files, 213 tests, 3 failing for a long time; 19 tests call scripts by path |
| Results | 28 folders (8.4 GB); 6 have no task of the same name |
| Data | `delivered/` 48 GB, `reduction/` 135 GB; three path-tracking files that no code reads |
| Other inputs | `extraction/` 16 GB (getsf, astrodendro); `region/` 0.7 MB, used by 3 tasks |
| Git | 30 uncommitted changes, including the whole of the new combination task |

Problems this plan addresses:

1. **Scripts are separated from their tasks**, and 15 scripts belong to no task.
2. **Results and tasks do not match one-to-one** (`moment_maps`, `starcluster`, `multiwavelength_cards`, `continuum_comparison_cards`, `diagnostics`, `compare_interferometer_arrays`, a symlink to `combination_diagnostic`).
3. **Derived products sit among the inputs.** `data/delivered/` contains images that a task can regenerate (unit conversions, regrids, the old v1–v3 combined images).
4. **Tasks read each other's results without a stated order**, and four tasks read from `results/paper_ii_referee_response/`, which is a coordination log, not a product.
5. **Nothing records which data can be regenerated and which cannot.**

## 3. Target layout

```
LANCET/
  data/
    delivered/        stable inputs (read-only)
    reduction/        reduction workspace
    extraction/       source-extraction products (getsf, astrodendro), permanent input
    README.md         rules
    MANIFEST.tsv      one row per input: origin, and whether and how it can be regenerated
  regions/            DS9 regions (renamed from region/)
  tasks/
    README.md         index of tasks, their order and status
    <task>/
      task.toml       inputs, outputs, settings
      run.py          the runnable script (moved from scripts/<task>.py)
      README.md       method, outputs, limitations
      *.py            helpers used only by this task
  src/                shared code, importable as `lancet` (no lancet/ subfolder)
  tests/              tests of src/ and of the task scripts
  results/<task>/     outputs; one folder per task, nothing else
  papers/ design/ docs/ presentation/ references/ website/     stand-alone
  pyproject.toml  README.md  CLAUDE.md
```

Removed: `scripts/`, `region/`, `src/lancet/`, the three path files in `data/`.

Run command for every task:

```bash
conda run -n py312 python tasks/<task>/run.py --stage <name>
```

## 4. Decisions (author, 2026-11-18)

**D1. Import name.** Keep `lancet`, mapped to the flat `src/` in `pyproject.toml`, so it can later
be imported as a package. The project is installed once with `pip install -e .`.

**D2. Tasks are reviewed one by one.** During the sort-up each task is run and checked by the
author, and may be renamed. Scripts that no current task needs are deleted after that check, not
before. The proposal below is the starting point for the review, not a decision:

| Scripts with no task | Proposal |
|---|---|
| `core_freefree_contamination`, `plot_core_freefree_contamination`, `freefree_mass_efficiency_bias` | into their own task `freefree_core_contamination` (done 2026-11-20; author) |
| `core_temperature_sampling`, `core_temperature_stress_test`, their two plot scripts, `report_core_temperature_tests` | new task `core_temperature_tests` |
| `simba_fourier_comparison`, `simba_peak_comparison`, `simba_pixel_comparison` | new task `simba_comparison` |
| `source_deconvolution` | into `source_catalogue` |
| `star_cluster` | new task `star_cluster` |
| `gas_to_dust_ratio` | new task `gas_to_dust_ratio` |

**D3. `data/delivered/` is sorted properly.** First case: remove `data/delivered/paper1/`.
`TP_SEDsynthetic230GHz_cont.fits` is byte-identical to `sedfitter/Dustmodel_230GHz.fits`, so the
two tasks that read it (`construct_freefree_model`, `synthesize_continuum`) are repointed there.
`PaperI_table2.fit` (read by `core_formation`) moves to `data/delivered/archive/paper-1/`.

**D4. `extraction/` moves into `data/`** as a permanent input (`data/extraction/`). The author
prunes useless and intermediate getsf files; before and after, the files needed to reproduce the
Paper II results are checked to exist (list in `docs/reorganization_baseline/`).

**D5. Results with no task.** As proposed: `moment_maps` → `c18o_moments`; `starcluster` →
`star_cluster`; `multiwavelength_cards` and `continuum_comparison_cards` become tasks; `diagnostics`
is sorted into the owning tasks; the `compare_interferometer_arrays` symlink is deleted (done 2026-11-19).

**D6. Reduction scripts stay in `data/reduction/`** (author, 2026-11-18). The reduction tree is a
stand-alone workspace with its own scripts (for example `cont/script/combine_cont.sh`); its finished
products are delivered to `data/delivered/`, and tasks read those. Tasks should not depend on files
inside the reduction tree; `combination_diagnostic`, which reads 13 of them, is settled at its review.

**D7. getsf and the new combination** (author, 2026-11-18). A getsf run is expensive and is not
repeated. The existing extraction stays as it is for the source catalogue and for the getsf figures
(Paper II Fig. 7 and D1). Quantities measured on the map (source background, optical-depth
correction and what follows) use the updated combined map. This needs one new product: the new
12m+7m+TP map minus the getsf source model of the 12m+7m image (a plain subtraction), replacing the
getsf background component of the old combined map.

**D8. No shared image-preparation task** (author, 2026-11-20). `prepare_continuum_images` is
dissolved: each task does its own Jy/beam→MJy/sr conversion or primary-beam cut on the delivered
images it reads, with the `lancet.imaging` functions and settings in its own `task.toml`. Its
readers absorb the step at their reviews: the merged combination task (replacing
`combination_diagnostic` and `linear_combine_continuum`), `core_formation` and
`volume_density_mapper`. Done 2026-11-20: `core_formation` and `volume_density_mapper` now read
delivered MIRIAD images (`lancet.imaging.read_mjysr`), the task folder was removed and
`results/prepare_continuum_images/` (1.8 GB) moved to the Trash.

## 5. Phases

Each phase ends with the test suite passing and one commit. Moves use `git mv`, so file history is kept.

### Phase 0. Checkpoint

- Commit the current state (the new combination task and its tests are not yet in git).
- Record the baseline: test results, and the output of each task's cheapest stage.
- Remove clutter that is not part of any decision: `.ipynb_checkpoints/`, `__pycache__/`, `.DS_Store`; move `G316_combination_tryall.ipynb` to `references/`.

### Phase 1. Tasks own their scripts

- Move `scripts/<task>.py` → `tasks/<task>/run.py` for the 22 paired scripts.
- Go through the tasks one at a time with the author (D2): run it, check the outputs, settle its name, then
  move its script. Unpaired scripts are placed or deleted in the same review.
- Update the 19 tests that call scripts by path, every README run command, and the default config path inside each script.
- Delete `scripts/`.

- 12m+7m continuum inputs (author, 2026-11-19): `data/delivered/alma/7MTM2TM1/cont/` now holds the MIRIAD
  clean, gain, sen, pb, clean.truesky and clean.flatnoise maps directly (see its README); the CASA products
  are in `casa-archived/` and not used for science; drogdeg, pbcut, MJysr and pbant were removed. Each task
  that still reads `cont/miriad/…` or `cont/casa/…` is repointed at its review, with the MIRIAD `pb` replacing
  the CASA primary beam: `prepare_continuum_images`, `construct_freefree_model`, `column_density`,
  `core_formation`, `volume_density_mapper`, `plot_allscale_continuum`, `linear_combine_continuum` (pbant),
  `combination_diagnostic` (CASA pb.tt0), `plot_core_freefree_contamination`, and the beam-matching test in
  `tests/test_interferometry.py` (used the removed convolve2ACA image as reference).

### Phase 2. Flatten `src/`

- Move `src/lancet/*.py` → `src/*.py` and adjust `pyproject.toml` (D1).
- Sort the modules: one used by a single task and with no wider use moves into that task's folder; shared ones stay. First candidates to move: `combination_figure`, `joint_image_comparison`, `tp_beam_test`, `amplitude_alignment`, `aca_forward_test` (all serve `combination_diagnostic` only).
- Add tests for the shared modules that have none (`infrared_mosaic` is the main one).

### Phase 3. Results match tasks

- Rename and sort the six unmatched result folders (D5).
- Every `task.toml` writes only to `results/<its own name>/`; add a test that checks this for all tasks.
- Replace reads from `results/paper_ii_referee_response/` by proper inputs. The Planck extracts used by `linear_combine_continuum` are promoted to `data/delivered/archive/planck/`.
- Write `tasks/README.md`: the list of tasks, which results each one needs, and the order to run them. The order follows from the `task.toml` inputs (for example `construct_freefree_model` → `column_density` → `source_catalogue` → `physical_properties` → `environment_measurements` → `bonnor_ebert`).

### Phase 4. Data

- Rename `region/` → `regions/` and update the 9 references in three `task.toml` files.
- Move `extraction/` → `data/extraction/` and update the 7 references in three `task.toml` files (D4). Done 2026-11-20 at the `source_catalogue` review (also the 4 paths in the `freefree_core_contamination` scripts); pruning still open.
- Remove `data/delivered/paper1/` as in D3.
- Write `data/MANIFEST.tsv`. One row per file or folder of `data/delivered/`, `data/reduction/` (top two levels) and `extraction/`:

  | Column | Meaning |
  |---|---|
  | path, size | what it is |
  | origin | `observed`, `archive` (public, can be downloaded again), `collaborator`, `reduced` (from `reduction/` by CASA or MIRIAD), `derived` (from other stored data) |
  | reproducible | `no`, `download`, `reduction` (hours, external tools), `task` (name of the task that makes it) |
  | used by | tasks that read it |

- Delete `path_aliases.yaml`, `path_ledger.tsv` and `path_ledger.initial-2026-07-27.tsv`. No code reads them; they are mentioned only in `CLAUDE.md`, `data/README.md`, `data/STATUS.md` and two audits in `docs/`, which are updated.
- Label the derived files in `data/delivered/` (D3). Known cases: `alma/TP7MTM2TM1/` (old v1–v3 combined images, now superseded), the `*.MJysr*`, `*.rp*` and `*.pbcut*` variants of the continuum images, and the smoothed SED products.
- Delete `data/delivered/archive/radio/CORNISH_5p5GHz_8arcsec_{native,on_MeerKATgrid}.fits`: same data as `results/construct_freefree_model/prepared/` (headers differ). First repoint `tests/test_freefree.py::test_paper2_profile_reproduces_reference_13_band_pixel`, which reads the delivered copy. Keep `Tile724_CORNISH_5GHz.fits` and `subcube_G316.hdr` as provenance (source tile and cutout header; no reader).
- Let git track `data/README.md` and `data/MANIFEST.tsv` (two exceptions in `.gitignore`), so a fresh clone states what data it needs.

### Phase 5. Documentation and rules

- Rewrite `README.md` and `CLAUDE.md` for the new layout and run command.
- Remove the two obsolete combination design documents in `docs/`.
- Fix the three long-standing test failures: a float64 FITS in `combination_diagnostic`, the `interferometry` inspection test, and the word "notebook" in the referee progress log.

### Phase 6. Proof

- Fresh clone into a temporary folder, install, link `data/`, `extraction/`: the tests pass and every task's inspection stage runs.
- Delete `results/` for two or three complete chains and regenerate them; compare with the baseline of Phase 0.

## 6. What stays as it is

- `data/delivered/` and `data/reduction/` and their rules.
- The content of every task: no science changes in this reorganisation.
- `tests/`. It is needed: it is the only check that a move did not break anything, and it guards the data rules (outputs never inside `data/`, float32 FITS, tasks independent of notebooks).
- `papers/`, `design/`, `docs/`, `presentation/`, `references/`, `website/`.

## 7. Risks

| Risk | Handling |
|---|---|
| A move breaks a task silently | baseline in Phase 0, tests after every phase, regeneration in Phase 6 |
| Paper II work is interrupted | Phases 1–3 touch only code layout; any phase can wait for a quiet moment |
| Results that cannot be regenerated are lost | nothing in `results/` is deleted before its task has regenerated it; hand-made records (referee diagnostics) are kept and marked as such |
| Deleting the path ledgers loses provenance | the manifest replaces what is still relevant; the ledgers remain in git history only if they are committed first (they are in `data/`, which git ignores, so they are copied to `docs/archive/` before deletion if wanted) |

## 8. Size of the work

About 110 file moves, 15 scripts to place, 6 result folders to sort, one manifest of roughly 200 rows, and documentation. Phases 0–3 are mechanical and can be done in one or two sessions; Phase 4 needs your judgement on individual data products.

## 9. Phase 1 review log

One row per task, filled in as the author checks it.

| Task | Checked | Name | Script | Rerun against stored results | Open points |
|---|---|---|---|---|---|
| `prepare_continuum_images` | 2026-11-19 | kept | `tasks/prepare_continuum_images/run.py` | all 15 outputs byte-identical | 9 of 15 outputs have no reader; 13 outputs duplicate files in `data/delivered/` (3.4 GB); 4 operations still use the old v3 combined images |
|  | 2026-11-20 | **removed** (D8) | — | — | readers repointed: `core_formation` overlay → MIRIAD flatnoise + 7m (converted in task); `volume_density_mapper` footprint → MIRIAD pb ≥ 0.2 (99% overlap with the old CASA cut); `combination_diagnostic` inputs handled in the merge |
| `construct_freefree_model` | 2026-11-20 | kept | `tasks/construct_freefree_model/run.py` | all 28 FITS byte-identical at 230 GHz; output frequency then made a setting, 223.8 GHz (author): `ff224_map*` = old × 1.002736 exactly, all other outputs byte-identical | done: dust grid read from `sedfitter/Dustmodel_230GHz.fits` (identical output); y-label fixed; core contamination split into `freefree_core_contamination` (not reviewed, ALMA inputs outdated; results moved to `results/freefree_core_contamination/`). Open: old `ff230_map*` kept until `synthesize_continuum`/`column_density` review; CORNISH 8″ copies in `delivered/` → Phase 4 |
| `synthesize_continuum` | 2026-11-20 | kept | `tasks/synthesize_continuum/run.py` (+ 3 SIMBA scripts in the task folder) | all 9 FITS byte-identical at the old 230 GHz setting; then both inputs set to 223.8 GHz (author): dust unchanged, free-free +0.27%, 250 GHz model +10.9% (SIMBA peak ratios 0.80–1.28 → 0.86–1.39) | done: `Dustmodel_223.8GHz.fits` (renamed by author), 3 delivered `SEDsynthetic_Dust_FF*` and old `*_230GHz`/`ff230` results moved to Trash; `linear_combine_continuum`, `column_density`, `multiwavelength_strip`, `plot_allscale_continuum` repointed (their results stale until rerun). Open: Paper II text and Figs. C1/C2 (250 GHz benchmark changed); SIMBA scripts not yet stages; `archive/infrared/SIMBA_1.2mm_newhdr*.fits` duplicate validation outputs |
| `column_density` | 2026-11-20 | kept | `tasks/column_density/run.py` | old v3 input rebuilt in scratch: all 7 maps within 1 float32 ulp of the stored ones, dense masses to 1e-10; then rerun on the new inputs: temperature and intercore maps byte-identical, all-scale Σ about 3× lower (v3 TP contribution ≈ 3.1 × delivered TP; new J-comb uses 1.0) | done: all-scale input = `combine_continuum` `jcomb_dust_224GHz_extended` (Jy/beam converted in task; not yet promoted), grid from the 12m+7m MIRIAD image, MIRIAD `pb` contour, PDF+PNG, colour scale 0–1 / 1–5 g cm⁻², negative intensity → NaN in the maps plus `*_signed` surface-density maps that keep negatives; figures show |signed Σ| (mirror fill) (author); new `consistency` stage: all-scale Σ smoothed to 18.1″ = TP-model Σ within 0.5 %, 3–8 % below Paper I Σ; input then the delivered `cont_dust` `extended36` product (whole 12m+7m grid; identical within 7.5e-5 g cm⁻² inside pb ≥ 0.3, 5.84 M pixels added; dense masses unchanged), README rewritten. Open: Paper II Fig. 4 and Σ statements (v3 plan); rerun after promotion; `dense` stage has no reader; downstream `volume_density_mapper` stale |
| `plot_allscale_continuum` | 2026-11-20 | kept | `tasks/plot_allscale_continuum/run.py` | archived v3 maps converted to MJy/sr: Fig. 3 clean view reproduced (45 of 15.7 M pixels differ); the stored dust-only clean view (29 July) predated later plotting changes | done: inputs = delivered J-comb `cont` and `cont_dust` `.extended` products (Jy/beam, converted with `read_mjysr`), MIRIAD `pb` contour, annotated views also PNG, colour scale −20 / 80 / 1200 MJy/sr (percentile-matched to Paper II's −21 / 233 / 1500), TP-model limits converted to MJy/sr (view unchanged), README rewritten. colour scale confirmed by the author; added stage `sqrt`: one inferno square-root scale, 0–1200 MJy/sr, blank pixels black; chosen as Paper II v3 Fig. 3 (author); TP-model view now `tp_continuum.png` in the same inferno style (0–275 MJy/sr); inputs then the filled `extended36` products and `continuum_224GHz_fill36` (no blank pixels; synthesize_continuum `fill`, combine_continuum extended TP inputs, author request). Open: Paper II Fig. 3 and caption (v3 plan) |
| `source_catalogue` | 2026-11-20 | kept | `tasks/source_catalogue/run.py` | IDs, order and 34 of 42 columns unchanged; flux diagnostics identical; `Ibkg2` from the new background (median 309 → 148 MJy/sr); 19 gain < 1 sources: 12m+7m fluxes, errors and `Ibkg1` × 1.016–2.166 | done (author): `extraction/` → `data/extraction/`; new stage `background` = J-comb 12m7mTP − getsf 12m7m source model / gain (replaces `sw.1303` background); gain correction of the 12m+7m quantities below gain 1, raw values kept; README rewritten. Open: `source_deconvolution.py` to place; `physical_properties` and later tasks stale |
| `physical_properties` | 2026-11-20 | kept | `tasks/physical_properties/run.py` | from the Paper II catalogue: all columns equal, `.tex` and summary byte-identical; on the new catalogue: `Surface_corebkg` median 2.73 → 1.39 g cm⁻², 19 gain < 1 cores × 1.01–2.17 in mass/τ/Σ, others unchanged to 1e-3; total core mass 545.7 → 552.6 M☉ | done: script moved, README updated, regression values updated, stray `.ipynb_checkpoints/` in results removed. Open: `[temperature_sampling]` belongs to the planned `core_temperature_tests` task; `core_mass_function` reviewed next (tests updated there) |
| `core_mass_function` | 2026-11-20 | kept | `tasks/core_mass_function/run.py` | from the Paper II properties: tables and posterior byte-identical; on the new properties: global Mmin 1.751 → 1.740 M☉, tail 92 → 94, α_MLE 2.415 → 2.407 (posterior 2.41, unchanged at 0.01); Young 2.32 → 2.27, Intermediate 2.49 → 2.46, Evolved unchanged | done: script moved, README with old/new table, regression values updated. Open: Paper II CMF numbers (v3 plan) |
| `environment_measurements` | 2026-11-20 | kept | `tasks/environment_measurements/run.py` | from the Paper II properties: all `*_env_*` columns equal except `Surface_env` (≤ 1e-6, float32 rounding of the regenerated low-resolution map); on the current inputs the environment measurements are unchanged, only the copied physical columns change | done: script moved, README rewritten. Then (author): `column_density` intercore map gain-corrected below gain 1 (new `gain` input), `Interco_env_mean` × up to 1.45 for 107 cores, nothing else changed. Open: `Vsig` input from `results/moment_maps/` (→ `c18o_moments`, D5) |
|  | 2026-11-19 (cont.) | | | | 7MTM2 dropped (no reader; author): its mask and conversion removed from `task.toml`, outputs deleted, the two tests that reproduced the delivered 7MTM2 files replaced by synthetic tests. Still open: 7MTM2TM1 inputs removed from `data/delivered/` (repoint to MIRIAD `pb`) |
| `combine_continuum` | 2026-11-20 | new (merges `linear_combine_continuum` and `combination_diagnostic`) | `tasks/combine_continuum/run.py` | not yet compared with the old results | merged task created; old tasks retired after comparison (J-comb within about 1 % inside gain = 1; factor numbers and Fig. 2 compared; author approves). Open: first full run, promotion of J-comb, `[absolute_scale]` pending |
| `linear_combine_continuum`, `combination_diagnostic` | 2026-11-20 | **removed** (retired, replaced by `combine_continuum`) | — | compared with the first full `combine_continuum` run in `results/combine_continuum/comparison/old_vs_new.md`; author approved | task folders, scripts and their tests deleted; code used only by dropped stages deleted (`pointing_profiles`, `tp_beam_test`, `run_pb_amplitude_curve`, the radial-profile readers in `visibility_diagnostics`); `paper_ii_referee_response/task.toml` repointed to `tasks/combine_continuum/run.py`; `results/linear_combine_continuum/` and `results/combination_diagnostic/` (including hand-written R4a reports) moved to `~/.Trash/LANCET_old_combination_results_2026-11-20/` |
