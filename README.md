# LANCET

LANCET is a reusable Python toolkit and reproducible task framework for the
G316 analysis.

G316 is the single scientific target of LANCET. It is therefore not repeated in
package names, configuration layers, or output directory names.

## AI assistance

Development of LANCET and its website has been assisted by OpenAI Codex and
Claude.

## Repository layout

- Reusable Python code lives under `src/lancet/`.
- Automated checks live under `tests/`.
- Global path policy lives in `src/lancet/paths.py`.
- Independent workflow settings and documentation live under `tasks/`.
- Runnable scientific commands live under `scripts/`.
- Small DS9 region definitions required by configured tasks live under
  `region/`.
- Each integration step is tested with the `py312` conda environment.

Install LANCET in the active Python 3.12 environment with:

```bash
python -m pip install -e ".[dev]"
```

Run the tests from the repository root with:

```bash
conda run -n py312 python -m pytest
```

## Local scientific materials

Observational data, external extraction products, manuscripts, presentations,
references, and generated results are intentionally not distributed in this
repository. The root `.gitignore` excludes:

- `data/` (including `data/extraction/`, the source-extraction products)
- `papers/`
- `presentation/`
- `references/`
- `results/`

Task configuration still declares the expected local input and output paths.
Stable inputs belong under `data/delivered/`, active reduction products under
`data/reduction/`, and generated task products under `results/`. Cloning the
code alone therefore supports installation and unit testing; running a
data-dependent task additionally requires its declared local inputs.

## Independent calculations

### Paper II referee response

The referee-response coordination task (`tasks/paper_ii_referee_response/`)
is kept locally and is not published (gitignored). It maps the referee items
to the analysis tasks and records diagnostics under
`results/paper_ii_referee_response/`.

### Standalone calculations

Calculate the gas-to-dust mass ratio from Giannetti et al. (2017), A&A, 606,
L12, Eq. 2:

```bash
conda run -n py312 python scripts/gas_to_dust_ratio.py \
  --radius-kpc 6.5 \
  --radius-error-kpc 0.2
```

Add `--json` for machine-readable output. The command keeps the statistical,
systematic, and intrinsic-scatter intervals separate.

Test elliptical-Gaussian source deconvolution with the reference example:

```bash
conda run -n py312 python scripts/source_deconvolution.py \
  --output results/diagnostics/source_deconvolution.png
```

The reusable function is available as `lancet.deconvolve_gaussian`.

Inspect the delivered collaborator star-cluster orbit catalogue and produce
the orbit diagnostics with:

```bash
conda run -n py312 python scripts/star_cluster.py
```

The default figures are written under `results/starcluster/`. Use
`--summary-only` to inspect the selected orbit without producing files.

Inspect the delivered radio inputs for the free-free model task:

```bash
conda run -n py312 python tasks/construct_freefree_model/run.py --stage inspect
```

Recreate the smoothed and reprojected CORNISH inputs under
`results/construct_freefree_model/prepared/` with `--stage prepare`. The Paper
II frequency, RMS, S/N, and mask settings recovered from the legacy products
are recorded in
[`tasks/construct_freefree_model/task.toml`](tasks/construct_freefree_model/task.toml).
The FITS provenance audit is in
[`docs/construct_freefree_model_audit.md`](docs/construct_freefree_model_audit.md).

After fitting, create the 18-arcsec free-free map on the dust grid and
synthesize the 223.8 and 250 GHz continuum products with:

```bash
conda run -n py312 python tasks/construct_freefree_model/run.py --stage postprocess
conda run -n py312 python tasks/construct_freefree_model/run.py --stage plot
conda run -n py312 python tasks/synthesize_continuum/run.py --stage inspect
conda run -n py312 python tasks/synthesize_continuum/run.py --stage all
```

The continuum model configuration, component inputs, frequency scaling, and
SIMBA validation settings are declared in
[`tasks/synthesize_continuum/task.toml`](tasks/synthesize_continuum/task.toml).
All generated maps and publication figures are written under
`results/synthesize_continuum/`. To recreate only the figures from existing
task products, run:

```bash
conda run -n py312 python tasks/synthesize_continuum/run.py \
  --stage plot \
  --overwrite
```

Each task converts the delivered continuum images it reads (Jy/beam to MJy/sr,
`lancet.imaging.read_mjysr`) itself; there is no separate image-preparation
task (decision D8 in `docs/project_reorganization_plan.md`).

Combine the 12m+7m true-sky continuum mosaic with the SED-based TP models in
the uv plane (J-comb adopted, feather as cross-check), check the TP factor and
make the uv-completeness diagnostics (Paper II Fig. 2) with:

```bash
PYTHONPATH=src conda run -n py312 python tasks/combine_continuum/run.py --stage combine --overwrite
PYTHONPATH=src conda run -n py312 python tasks/combine_continuum/run.py --stage diagnostics --overwrite
```

Settings are in
[`tasks/combine_continuum/task.toml`](tasks/combine_continuum/task.toml); inputs,
method, products and limitations in
[`tasks/combine_continuum/README.md`](tasks/combine_continuum/README.md).
Outputs go under `results/combine_continuum/`. The diagnostic stages read the
four declared measurement sets under `data/reduction/alma/7MTM2TM1/cont/ms/` as
read-only intermediates. This task replaced `linear_combine_continuum` and
`combination_diagnostic`, which were retired on 2026-11-20.

Create the all-scale continuum display products with:

```bash
conda run -n py312 python scripts/plot_allscale_continuum.py --stage inspect
conda run -n py312 python scripts/plot_allscale_continuum.py --stage all
```

This task reads the delivered dust and dust-plus-free-free all-scale maps and
the current synthesized 223.8 GHz continuum result. It creates annotated
publication PDFs and clean PNGs under
`results/plot_allscale_continuum/figures/` without writing any FITS products.
Its input paths, view, transition values, and output formats are declared in
[`tasks/plot_allscale_continuum/task.toml`](tasks/plot_allscale_continuum/task.toml).

Test `volume-density-mapper` on the exact 1500-pixel central field and then
build the resolution-aware full-map width and mean-density products with:

```bash
python -m pip install -e ".[volume-density]"
conda run -n py312 python tasks/volume_density_mapper/run.py inspect
conda run -n py312 python tasks/volume_density_mapper/run.py pilot
conda run -n py312 python tasks/volume_density_mapper/run.py full
conda run -n py312 python tasks/volume_density_mapper/run.py mask
```

The task uses the finite ALMA footprint to distinguish the 0.65-arcsec inner
resolution from the 18-arcsec outer map, applies a local beam floor, and flags
the resolution transition. The final `mask` stage writes the recommended
PB-response >= 0.20 products with the outer region removed. The native full
3-D cube is deliberately not
written because the upstream pilot fails projected-column closure; details
and product definitions are in
[`tasks/volume_density_mapper/README.md`](tasks/volume_density_mapper/README.md).

Reconstruct the C18O cube derivatives, masked moments, ALMA-grid
velocity-dispersion map, and moment figure from the single delivered line cube
with:

```bash
conda run -n py312 python tasks/c18o_moments/run.py --stage inspect
conda run -n py312 python tasks/c18o_moments/run.py --stage all
```

The dilated-mask thresholds and supporting parameters are declared in
[`tasks/c18o_moments/task.toml`](tasks/c18o_moments/task.toml). All generated
cubes, moment products, and figures are written under
`results/c18o_moments/`; the delivered data tree remains read-only.

Build the source catalogue and its dependent physical-property table with:

```bash
conda run -n py312 python tasks/source_catalogue/run.py --stage inspect
conda run -n py312 python tasks/source_catalogue/run.py --stage all
conda run -n py312 python tasks/core_physical_properties/run.py --stage inspect
conda run -n py312 python tasks/core_physical_properties/run.py --stage build
```

The source task applies the declared getsf cuts, samples the background and
temperature maps, measures effective deconvolved sizes, and assigns the
evolutionary subregions. Its `diagnostics` stage also compares the 12m+7m and
12m+7m+TP peak and integrated flux measurements and writes
`GetsfResult_FluxMeasure` and `GetsfResult_FluxRatio` under
`results/source_catalogue/figures/`. The physical task applies deterministic
Monte Carlo uncertainty propagation. Machine-readable tables and figures are
written under `results/source_catalogue/` and
`results/core_physical_properties/`.

Calculate the supporting foreground/background opacity sensitivity with:

```bash
conda run -n py312 python scripts/background_effects.py --stage inspect
conda run -n py312 python scripts/background_effects.py --stage all
```

The task evaluates the corrected compact-core optical depth for configured
environmental emission and core brightness grids. Its assumptions are
declared in
[`tasks/background_effects/task.toml`](tasks/background_effects/task.toml);
the table, summary, and PDF/PNG figure are written under
`results/background_effects/`.

Measure the regional cloud-to-fragment-to-core mass efficiencies and create
the compact-core/fragment overlay figures with:

```bash
conda run -n py312 python tasks/core_formation_efficiency/run.py --stage inspect
conda run -n py312 python tasks/core_formation_efficiency/run.py --stage all
```

The generated physical-properties catalogue supplies the core masses and
counts. The cloud and Paper I-scale fragment values are explicit adopted
inputs in
[`tasks/core_formation_efficiency/task.toml`](tasks/core_formation_efficiency/task.toml). The
fragment masses are rescaled from the Paper I opacity of 0.9 cm2/g to the
current 0.5 cm2/g opacity; cloud masses remain on their existing SED scale.
The
delivered Paper I fragment catalogue is recorded as an optional provenance
reference rather than a calculation dependency. Tables and figures are
written under `results/core_formation_efficiency/`.

Measure the temperature, C18O velocity-dispersion, total surface-density, and
intercore surface-density environment around every core with:

```bash
conda run -n py312 python tasks/environ_measure/run.py --stage inspect
conda run -n py312 python tasks/environ_measure/run.py --stage build
```

The task uses 15-arcsec apertures and writes the augmented catalogue under
`results/environ_measure/`. Its intercore measurements use the
current 0.5 cm2/g opacity scale and are therefore 1.8 times the retained
values calculated on the earlier 0.9 cm2/g scale.

Calculate the core–environment correlations and the intercore
surface-density regression with:

```bash
conda run -n py312 python tasks/environ_core_correlation/run.py --stage inspect
conda run -n py312 python tasks/environ_core_correlation/run.py --stage all
```

The task writes a machine-readable pairwise-correlation table, a deterministic
MCMC summary and posterior sample, and four PDF/PNG figure pairs under
`results/environ_core_correlation/`, including the clean single-panel paper figure
`correlation_surface_env-surface_core.pdf`. Variables, labels, uncertainty
profiles, sampler settings, and plot styles are declared in
[`tasks/environ_core_correlation/task.toml`](tasks/environ_core_correlation/task.toml).

Calculate thermal Bonnor–Ebert stability and the global/regional core mass
functions with:

```bash
conda run -n py312 python tasks/bonnor_ebert/run.py --stage inspect
conda run -n py312 python tasks/bonnor_ebert/run.py \
  --stage all --profile radius_thermal
conda run -n py312 python tasks/bonnor_ebert/run.py \
  --stage all --profile measured_external
conda run -n py312 python tasks/core_mass_function/run.py --stage inspect
conda run -n py312 python tasks/core_mass_function/run.py --stage all
```

The Bonnor–Ebert task provides isolated radius and measured-external-pressure
profiles. Quantities may come from per-core catalogue columns or configured
scalars; internal support is kept separate from the environmental density and
velocity terms. The measured-external profile converts
`Interco_env_mean` to volume density with an adopted 0.4-pc depth, combines
`Temp_env_mean` and the thermally corrected `Vsig_env_mean` into the external
pressure, and uses a configurable 0–1 fraction of the environmental
nonthermal dispersion for internal support. The retained fraction is
explicitly recorded as a provisional dissipation assumption. The CMF task
fits
`dN/dM proportional to M^(-alpha)`, selects `Mmin` by its KS diagnostic, and
writes fit tables, posterior samples, a local catalogue export, and PDF/PNG
figures. Their settings and assumptions are declared in
[`tasks/bonnor_ebert/task.toml`](tasks/bonnor_ebert/task.toml) and
[`tasks/core_mass_function/task.toml`](tasks/core_mass_function/task.toml);
all products stay under `results/`.

The magnetic-field analysis is maintained as an explicitly ongoing project
idea:

```bash
conda run -n py312 python scripts/magnetic_fields.py --stage inspect
conda run -n py312 python scripts/magnetic_fields.py --stage all
```

It creates provisional Planck polarization masks, debiased polarization
products, orientation maps, and a diagnostic vector figure. It does not
estimate field strength. It declares `status = "ongoing"` and writes
summaries with `paper_ready = false`. Its assumptions, unresolved decisions,
and continuity requirements are recorded in
[`tasks/magnetic_fields/README.md`](tasks/magnetic_fields/README.md).


## Project paths

The package discovers the repository root from `pyproject.toml`, so code does
not depend on the checkout directory being named `G316` or on the command's
current working directory. `ProjectPaths` defines the fixed global layout:
stable inputs under `data/delivered/`, active reduction under
`data/reduction/`, and generated products under `results/`. Workflow-specific
paths and scientific settings are declared in `tasks/*/task.toml`.

Stable analysis inputs should be resolved with `delivered_input()`. Reading an
intermediate reduction product requires an explicit `declared=True` argument,
and output validation rejects paths inside `data/`.

## Task convention

Runnable tasks must be self-contained scientific workflows. New files under
`scripts/` and `tasks/` must describe their inputs, calculations, and outputs
without depending on or naming exploratory source documents.

Floating-point FITS task products use float32 storage (`BITPIX = -32`) by
default. A task may use float64 internally for calculations, but writing a
float64 FITS product requires a documented scientific precision need and a
dedicated regression test.
