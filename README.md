# LANCET

LANCET is a reusable Python toolkit and reproducible task framework for the
G316 analysis.

G316 is the single scientific target of LANCET. It is therefore not repeated in
package names, configuration layers, or output directory names.

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

- `data/`
- `extraction/`
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
conda run -n py312 python scripts/construct_freefree_model.py --stage inspect
```

Recreate the smoothed and reprojected CORNISH inputs under
`results/construct_freefree_model/prepared/` with `--stage prepare`. The Paper
II frequency, RMS, S/N, and mask settings recovered from the legacy products
are recorded in
[`tasks/construct_freefree_model/task.toml`](tasks/construct_freefree_model/task.toml).
The FITS provenance audit is in
[`docs/construct_freefree_model_audit.md`](docs/construct_freefree_model_audit.md).

After fitting, create the 18-arcsec free-free map on the dust grid and
synthesize the 230 and 250 GHz continuum products with:

```bash
conda run -n py312 python scripts/construct_freefree_model.py --stage postprocess
conda run -n py312 python scripts/construct_freefree_model.py --stage plot
conda run -n py312 python scripts/synthesize_continuum.py --stage inspect
conda run -n py312 python scripts/synthesize_continuum.py --stage all
```

The continuum model configuration, component inputs, frequency scaling, and
SIMBA validation settings are declared in
[`tasks/synthesize_continuum/task.toml`](tasks/synthesize_continuum/task.toml).
All generated maps and publication figures are written under
`results/synthesize_continuum/`. To recreate only the figures from existing
task products, run:

```bash
conda run -n py312 python scripts/synthesize_continuum.py \
  --stage plot \
  --overwrite
```

Prepare delivered continuum images for the multi-array and all-scale
comparison tasks with:

```bash
conda run -n py312 python scripts/prepare_continuum_images.py --stage inspect
conda run -n py312 python scripts/prepare_continuum_images.py \
  --stage mask --only 7MTM2
conda run -n py312 python scripts/prepare_continuum_images.py \
  --stage convert --only 7MTM2_pbcut
```

This applies the declared primary-beam cutoff masks before converting the
configured Jy/beam images to MJy/sr. The full configured output is about
1.8 GiB with float32 FITS storage; omit `--only` only when the complete
products are required. The
source audit and remaining task boundaries are recorded in
[`docs/continuum_image_workflow_audit.md`](docs/continuum_image_workflow_audit.md).

After preparing the continuum images, recreate the multi-array comparison
products with:

```bash
conda run -n py312 python scripts/compare_interferometer_arrays.py --stage inspect
conda run -n py312 python scripts/compare_interferometer_arrays.py --stage uv
conda run -n py312 python scripts/compare_interferometer_arrays.py --stage compare-images
conda run -n py312 python scripts/compare_interferometer_arrays.py --stage beam-match
```

The UV stage reads the four declared measurement sets under
`data/reduction/alma/7MTM2TM1/cont/ms/` without modifying them. Later plot
runs reuse the compact coordinate cache unless `--refresh-uv` is requested.

Create the all-scale continuum display products with:

```bash
conda run -n py312 python scripts/plot_allscale_continuum.py --stage inspect
conda run -n py312 python scripts/plot_allscale_continuum.py --stage all
```

This task reads the delivered dust and dust-plus-free-free all-scale maps and
the current synthesized 230 GHz continuum result. It creates annotated
publication PDFs and clean PNGs under
`results/plot_allscale_continuum/figures/` without writing any FITS products.
Its input paths, view, transition values, and output formats are declared in
[`tasks/plot_allscale_continuum/task.toml`](tasks/plot_allscale_continuum/task.toml).

Reconstruct the C18O cube derivatives, masked moments, ALMA-grid
velocity-dispersion map, and moment figure from the single delivered line cube
with:

```bash
conda run -n py312 python scripts/c18o_moments.py --stage inspect
conda run -n py312 python scripts/c18o_moments.py --stage all
```

The dilated-mask thresholds and supporting parameters are declared in
[`tasks/c18o_moments/task.toml`](tasks/c18o_moments/task.toml). All generated
cubes, moment products, reprojections, and figures are written under
`results/moment_maps/`; the delivered data tree remains read-only.

Build the source catalogue and its dependent physical-property table with:

```bash
conda run -n py312 python scripts/source_catalogue.py --stage inspect
conda run -n py312 python scripts/source_catalogue.py --stage all
conda run -n py312 python scripts/physical_properties.py --stage inspect
conda run -n py312 python scripts/physical_properties.py --stage build
```

The source task applies the declared getsf cuts, samples the background and
temperature maps, measures effective deconvolved sizes, and assigns the
evolutionary subregions. Its `diagnostics` stage also compares the 12m+7m and
12m+7m+TP peak and integrated flux measurements and writes
`GetsfResult_FluxMeasure` and `GetsfResult_FluxRatio` under
`results/source_catalogue/figures/`. The physical task applies deterministic
Monte Carlo uncertainty propagation. Machine-readable tables and figures are
written under `results/source_catalogue/` and
`results/physical_properties/`.

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
conda run -n py312 python scripts/core_formation.py --stage inspect
conda run -n py312 python scripts/core_formation.py --stage all
```

The generated physical-properties catalogue supplies the core masses and
counts. The cloud and Paper I-scale fragment values are explicit adopted
inputs in
[`tasks/core_formation/task.toml`](tasks/core_formation/task.toml). The
fragment masses are rescaled from the Paper I opacity of 0.9 cm2/g to the
current 0.5 cm2/g opacity; cloud masses remain on their existing SED scale.
The
delivered Paper I fragment catalogue is recorded as an optional provenance
reference rather than a calculation dependency. Tables and figures are
written under `results/core_formation/`.

Measure the temperature, C18O velocity-dispersion, total surface-density, and
intercore surface-density environment around every core with:

```bash
conda run -n py312 python scripts/environment_measurements.py --stage inspect
conda run -n py312 python scripts/environment_measurements.py --stage build
```

The task uses 15-arcsec apertures and writes the augmented catalogue under
`results/environment_measurements/`. Its intercore measurements use the
current 0.5 cm2/g opacity scale and are therefore 1.8 times the retained
values calculated on the earlier 0.9 cm2/g scale.

Calculate the core–environment correlations and the intercore
surface-density regression with:

```bash
conda run -n py312 python scripts/correlation.py --stage inspect
conda run -n py312 python scripts/correlation.py --stage all
```

The task writes a machine-readable pairwise-correlation table, a deterministic
MCMC summary and posterior sample, and four PDF/PNG figure pairs under
`results/correlation/`, including the clean single-panel paper figure
`Correlation_Surface_env-Surface_core.pdf`. Variables, labels, uncertainty
profiles, sampler settings, and plot styles are declared in
[`tasks/correlation/task.toml`](tasks/correlation/task.toml).

Calculate thermal Bonnor–Ebert stability and the global/regional core mass
functions with:

```bash
conda run -n py312 python scripts/bonnor_ebert.py --stage inspect
conda run -n py312 python scripts/bonnor_ebert.py \
  --stage all --profile radius_thermal
conda run -n py312 python scripts/bonnor_ebert.py \
  --stage all --profile measured_external
conda run -n py312 python scripts/core_mass_function.py --stage inspect
conda run -n py312 python scripts/core_mass_function.py --stage all
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

Two additional analyses are maintained as explicitly ongoing project ideas:

```bash
conda run -n py312 python scripts/magnetic_fields.py --stage inspect
conda run -n py312 python scripts/magnetic_fields.py --stage all
conda run -n py312 python scripts/n_pdf.py --stage inspect
conda run -n py312 python scripts/n_pdf.py --stage all
```

The magnetic-field task creates provisional Planck polarization masks,
debiased polarization products, orientation maps, and a diagnostic vector
figure. It does not estimate field strength. The N-PDF task creates three
regional cutouts and normalized column-density histograms from the currently
adopted retained all-scale map. Both tasks declare `status = "ongoing"` and
write summaries with `paper_ready = false`. Their assumptions, unresolved
decisions, and continuity requirements are recorded in
[`tasks/magnetic_fields/README.md`](tasks/magnetic_fields/README.md) and
[`tasks/n_pdf/README.md`](tasks/n_pdf/README.md).


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
