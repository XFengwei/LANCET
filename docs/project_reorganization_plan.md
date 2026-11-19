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
2. **Results and tasks do not match one-to-one** (`moment_maps`, `starcluster`, `multiwavelength_cards`, `continuum_comparison_cards`, `diagnostics`, an empty `compare_interferometer_arrays`).
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
| `core_freefree_contamination`, `plot_core_freefree_contamination`, `freefree_mass_efficiency_bias` | into `construct_freefree_model` |
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
is sorted into the owning tasks; the empty `compare_interferometer_arrays` is deleted.

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
- Move `extraction/` → `data/extraction/` and update the 7 references in three `task.toml` files (D4).
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
