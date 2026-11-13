# LANCET — instructions for Claude

LANCET is a Python toolkit and task framework for the analysis of a single
target, the G316 filament (IRDC 316.8). G316 is never repeated in package,
configuration or output names. Paper I is published, Paper II is under
revision (A&A, referee response), and further papers (e.g. APEX spectral-line
work) are developed in the same repository. `README.md` describes the layout
for people; this file records the working rules.

## Repository layout

- `src/lancet/` — reusable code, one module per topic. `paths.py` holds the
  global path policy (`ProjectPaths`).
- `scripts/<task>.py` — one runnable command per task.
- `tasks/<task>/` — `task.toml` (all paths and scientific settings) and
  `README.md` (inputs, method, outputs, status).
- `tests/test_<module>.py` — pytest tests for `src/lancet/`.
- `region/` — small DS9 regions that tasks declare.
- `docs/` — audits and design notes.
- Local only (gitignored, never committed): `data/`, `results/`,
  `extraction/`, `papers/`, `presentation/`, `references/`. A fresh clone has
  code and tests but no data; create these folders locally with the same names.

## Data rules

- `data/delivered/` holds stable, immutable inputs. Analyses read from here,
  through `ProjectPaths.delivered_input()`.
- `data/reduction/` is the mutable reduction workspace (CASA, MIRIAD, export
  scripts). Analyses may read it only when the input is declared as a
  read-only intermediate (`reduction_input(..., declared=True)`), for
  reduction diagnostics or method illustrations.
- Analyses never write inside `data/`; `ProjectPaths.validate_output()`
  enforces this. All products go under `results/<task>/`.
- A reduction product becomes an analysis input only after it is promoted to
  `data/delivered/`. Never modify or overwrite delivered files.
- Reruns or tests of a reduction go in a new, clearly named subfolder of
  `data/reduction/...` (for example `<product>_r4b/`), never over existing
  products. Delete their intermediates when the user has what they need;
  MIRIAD and CASA runs can fill the disk.
- `data/path_aliases.yaml` and `data/path_ledger.tsv` map historical paths to
  the current layout.

## Task pattern

Every analysis is a self-contained task:

1. `tasks/<task>/task.toml` declares `[task]`, `[inputs]` (repository-relative
   paths), `[outputs] root = "results/<task>"`, and one table per stage with
   its scientific settings. No paths or parameters are hard-coded in code.
2. `scripts/<task>.py` loads the TOML, resolves paths through `ProjectPaths`,
   and exposes `--config`, `--stage <name>` and `--overwrite`. Existing
   outputs are not replaced without `--overwrite`.
3. Calculations live in `src/lancet/<module>.py`, with tests in `tests/`.
4. `tasks/<task>/README.md` gives the run command, inputs, method, outputs and
   known limitations. Task files must not depend on or name exploratory
   notebooks or documents.

Small one-off checks (for example the SIMBA and core-temperature tests) are
standalone scripts without a task folder; `README.md` lists them under
"Independent calculations". Use a full task for anything a paper relies on.

Run a stage from the repository root:

```bash
PYTHONPATH=src conda run -n py312 python scripts/<task>.py --stage <name> --overwrite
```

## Code conventions

- Python 3.12 in the conda environment `py312`; install with
  `python -m pip install -e ".[dev]"`.
- Run tests with `conda run -n py312 python -m pytest` before reporting
  work as done.
- FITS products are float32 (`BITPIX = -32`). Writing float64 needs a
  documented reason and a regression test.
- Figures are written as PDF plus PNG (and SVG where a task asks for it).
  Inspect the PNG before reporting a figure.
- Match the existing style: short functions, settings passed in from the TOML,
  no new dependencies without asking.

## External tools

Install paths differ between machines; check before use.

- CASA 6.7.2 (`casa --nologger --nogui -c script.py`) for measurement sets.
- MIRIAD runs in Docker on macOS (amd64 under Rosetta):
  - `miriad80:latest` — MIRIAD 6.0.2 with MAXANT=80, needed for the combined
    ALMA 12 m data (71 antenna rows). Built from
    `miriad-alma-maxant80-u2204.tar.gz` and its `Dockerfile.miriad.alma`.
    Tools are on `PATH`; run scripts with `docker run --rm --platform
    linux/amd64 -v "$PWD:/data" -w /data miriad80:latest bash script.sh`.
  - `my-carma-miriad:latest` — older CARMA build (MAXANT=64); needs
    `bash -lc` to load its environment.
  - Keep Docker Desktop's disk in mind: a full host disk corrupts its VM.
- LaTeX: compile with `lualatex` into a scratch output directory when only
  checking a build, so the user's PDFs are not overwritten.

## Papers

- Each paper lives in `papers/<roman numeral>/` (local only). Its
  `AGENTS.md` names the working source and which manuscript variants may be
  edited; read it before touching any `.tex` file.
- Scientific changes for a paper belong in the owning task under `tasks/`,
  not in a paper-specific pipeline. Copy only verified task figures into a
  paper's figure folder and record their source.
- A referee response is tracked by a coordination task (for Paper II,
  `tasks/paper_ii_referee_response/`): `referee_requests.md` lists the items,
  `progress.md` logs each completed work unit in its fixed entry format, and
  diagnostics go under `results/<that task>/diagnostics/`.
- Response and manuscript statements must match completed, verified results.
  Do not finalize a manuscript or generate reserved variants unless asked.

## Working practice

- Confirm a reproduction against the delivered product before drawing
  conclusions from it, and say plainly when it does not match.
- When a long run (CASA, MIRIAD clean) is needed, state the expected runtime
  and disk use first, and run it in the background.
- Report numbers with their source file, so they can be traced later.
