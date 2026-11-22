# LANCET

LANCET is a Python toolkit and task framework for the analysis of a single
target, the G316 filament (IRDC 316.8). Because there is only one target, its
name is not repeated in package, configuration or output names. LANCET Paper I
is published; Paper II is under revision.

## AI assistance

Development of LANCET and its website has been assisted by OpenAI Codex and
Claude.

## Repository layout

| Folder | Content |
|---|---|
| `src/lancet/` | Reusable code, one module per topic. `lancet.combination` holds the code of the continuum combination; `lancet.paths` holds the path policy (`ProjectPaths`). |
| `tasks/<task>/` | One analysis per folder: `task.toml` (inputs, outputs, settings), `run.py` (the command), `README.md` (method, outputs, limitations). |
| `tests/` | pytest tests of `src/lancet/` and of the task layout. |
| `regions/` | Small DS9 regions that tasks declare. |
| `design/` | Stand-alone artwork (print cards); scripts and settings are tracked, rendered files are not. |
| `docs/` | Audits, design notes and `make_data_manifest.py`. |
| `website/` | Project website. |

Not distributed (gitignored, local only): `data/` (except `data/README.md`
and `data/MANIFEST.tsv`), `results/`, `papers/`, `presentation/` and
`references/`. A fresh clone has the code and tests but no data; create these
folders locally with the same names.

## Installation and tests

```bash
conda run -n py312 python -m pip install -e ".[dev]"
```

`volume_density_mapper` also needs the extra `.[volume-density]`. Run the
tests from the repository root:

```bash
conda run -n py312 python -m pytest
```

## Tasks

Every analysis is a task. Run one stage from the repository root:

```bash
PYTHONPATH=src conda run -n py312 python tasks/<task>/run.py --stage <name> --overwrite
```

Each `run.py` takes `--config`, `--stage` and `--overwrite`; existing outputs
are not replaced without `--overwrite`. A task writes only to
`results/<task>/`. [`tasks/README.md`](tasks/README.md) lists all tasks, what
each needs from the others, and the order in which to run them; each task's
README gives its stages, inputs and outputs.

### Standalone calculations

Small calculations without outputs run as modules:

- `python -m lancet.gas_to_dust --radius-kpc 6.5 --radius-error-kpc 0.2`:
  gas-to-dust ratio from Giannetti et al. (2017, Eq. 2). The Paper II value is
  the task `gas_to_dust_ratio`.
- `lancet.deconvolve_gaussian` and `lancet.deconvolution_size`:
  elliptical-Gaussian source deconvolution and its round-trip check.

## Data

`data/README.md` gives the rules, and `data/MANIFEST.tsv` lists every input
file with its origin and the tasks that read it. In short:

- `data/delivered/` holds stable, immutable inputs; analyses read them through
  `ProjectPaths.delivered_input()`.
- `data/reduction/` is the mutable reduction workspace (CASA, MIRIAD). An
  analysis may read it only as a declared read-only intermediate
  (`reduction_input(..., declared=True)`).
- Analyses never write inside `data/`; `ProjectPaths.validate_output()`
  rejects such paths.

`ProjectPaths` finds the repository root from `pyproject.toml`, so commands do
not depend on the folder name or the working directory.
