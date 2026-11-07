# Plot all-scale continuum

This task creates the publication and clean-display views of the all-scale
1.3 mm continuum emission.

## Inputs

- The delivered dust-plus-free-free all-scale continuum image.
- The delivered dust-only all-scale continuum image.
- The current 230 GHz synthesized-continuum task result.
- The adopted 12m+7m mosaic primary-beam map. Annotated views show its 0.2
  contour, matching the source-extraction selection.

All inputs are read-only. The task does not use or modify reduction data.

## Stages

- `inspect` reports input paths, shapes, FITS precision, and units without
  reading the full image arrays.
- `combined` creates annotated and clean dust-plus-free-free views.
- `synthetic` creates a clean view of the synthesized 230 GHz continuum.
- `dust` creates annotated and clean dust-only views.
- `all` runs all three plotting stages.

All generated figures are written under
`results/plot_allscale_continuum/figures/`.

Run:

```bash
conda run -n py312 python scripts/plot_allscale_continuum.py --stage inspect
conda run -n py312 python scripts/plot_allscale_continuum.py --stage combined
conda run -n py312 python scripts/plot_allscale_continuum.py --stage synthetic
conda run -n py312 python scripts/plot_allscale_continuum.py --stage dust
conda run -n py312 python scripts/plot_allscale_continuum.py --stage all
```

Pass `--overwrite` only when intentionally replacing existing task figures.
