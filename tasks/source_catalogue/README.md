# Source catalogue and flux diagnostics

This task converts the external getsf extraction into the project source
catalogue. It applies the declared reliability and compactness cuts, samples
the intercore background, core background, and dust temperature at each
source, calculates effective deconvolved sizes, assigns stable IDs, and tags
the three evolutionary subregions.

The same task reproduces the two flux-measurement diagnostics. In the getsf
catalogue, suffix `01` is the 12m+7m measurement and suffix `02` is the
12m+7m+TP measurement. The ratio diagnostic plots

\[
R = \frac{S_{\rm peak}}{S_{\rm integrated}},
\]

with independent peak and integrated-flux errors propagated into \(R\). The
comparison inset uses

\[
\Delta S / S =
\frac{S_{\rm 12m+7m+TP} - S_{\rm 12m+7m}}{S_{\rm 12m+7m}}
\]

for sources whose peak measurements have signal-to-noise ratios above the
configured threshold in both images. Plot settings and selection thresholds
are declared in `task.toml`.

The principal data colors are controlled by `reference_color` and
`comparison_color` in the `[flux_diagnostics]` section of `task.toml`.
Marker-edge and error-bar colors are derived automatically from these two
values.

```bash
conda run -n py312 python scripts/source_catalogue.py --stage inspect
conda run -n py312 python scripts/source_catalogue.py --stage build
conda run -n py312 python scripts/source_catalogue.py --stage diagnostics
conda run -n py312 python scripts/source_catalogue.py --stage all
```

The machine-readable ECSV and publication table are written under
`results/source_catalogue/tables/`. The diagnostic figures are written under
`results/source_catalogue/figures/` as
`GetsfResult_FluxMeasure.pdf` and `GetsfResult_FluxRatio.pdf`, with PNG
companions. External extraction products and temperature maps are read-only.
