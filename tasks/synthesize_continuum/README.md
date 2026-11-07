# Synthesize continuum

This task constructs continuum models by combining matched-resolution dust and
free-free emission maps.

## Inputs

- A stable published 230 GHz dust model under `data/delivered/paper1/`.
- The 8000 K, 230 GHz free-free model produced by
  `construct_freefree_model`, smoothed to 18 arcsec and reprojected to the dust
  grid.
- The archival SIMBA 1.2 mm map for optional 250 GHz validation.

All inputs are read-only.

## Model

For a target frequency `nu`, the task evaluates

```text
dust(nu)     = dust(230 GHz)     * (nu / 230 GHz)^3.8
freefree(nu) = freefree(230 GHz) * (nu / 230 GHz)^-0.1
continuum    = dust(nu) + freefree(nu)
```

Missing free-free pixels are set to zero before addition; missing dust pixels
remain missing. The configured indices and frequencies reproduce the Paper II
230 and 250 GHz products.

## Outputs

Component and combined FITS maps are written under
`results/synthesize_continuum/model/` for each configured frequency:

- `dust_<frequency>.fits`;
- `freefree_<frequency>.fits`; and
- `continuum_<frequency>.fits`.

The optional validation products are written under
`results/synthesize_continuum/validation/`:

- the 250 GHz model smoothed from 18 to 24 arcsec;
- the calibrated SIMBA map; and
- SIMBA reprojected to the model grid.

Publication figures are written as PDF and PNG under
`results/synthesize_continuum/figures/`:

- `continuum_components_230GHz` compares dust, free-free, and their sum; and
- `simba_comparison_250GHz` compares SIMBA with the beam-matched model; and
- `sed_demonstration_230GHz` illustrates the dust and free-free spectral
  components and their combined 230 GHz prediction.

The crop, intensity range, contour levels, color map, scale bar, resolution,
and output formats are declared in the `[plot]` configuration tables.

## Stages

- `inspect` checks the declared inputs, grids, beams, units, and model
  readiness without writing files.
- `model` creates the component and combined 230 and 250 GHz maps.
- `validate` creates the SIMBA comparison products after `model`.
- `plot` creates the two figures after `model` and `validate`.
- `all` runs all four stages in order.

Run from the project root:

```bash
conda run -n py312 python scripts/synthesize_continuum.py --stage inspect
conda run -n py312 python scripts/synthesize_continuum.py --stage all
```

Existing products are protected by default. Add `--overwrite` only when the
declared configuration or an upstream input has intentionally changed.

For the Paper II SIMBA peak comparison, first run the `model` and `validate`
stages, then run:

```bash
conda run -n py312 python scripts/simba_peak_comparison.py
```

For the common-footprint Fourier-amplitude comparison, run:

```bash
conda run -n py312 python scripts/simba_fourier_comparison.py
```

The script writes the figure to `results/synthesize_continuum/figures/` and
the binned data and methods to `results/synthesize_continuum/validation/`.
Use `--overwrite` to replace existing outputs. See the R10 Fourier diagnostic
report for the window, beam and noise limitations.

This separately maintained comparison selects SIMBA local maxima at least
1 Jy/beam in the configured display crop, and compares their brightnesses
with the 250 GHz model at the same positions and 24-arcsec resolution. It
also reports peak contrasts after subtracting median backgrounds in 36--60
arcsec annuli. The CSV and JSON measurements are in `validation/`; the
figure is `figures/simba_bright_peak_comparison.pdf`. The 0.20 SIMBA flux
calibration factor was confirmed by the author. Annular contrasts are
affected by neighbouring emission and negative bowls and do not independently
recover missing spatial scales. Use `--overwrite` to replace these outputs.
