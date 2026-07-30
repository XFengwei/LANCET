# Core–environment correlations

This task measures pairwise correlations among four core properties and four
environment properties, then fits the intercore surface-density relation with
a power law. Its input is the generated, opacity-consistent environment
catalogue:

`results/environment_measurements/tables/core_environment_catalogue.ecsv`

The task produces a detailed corner plot, a compact correlation matrix, and a
regression-comparison diagnostic from one declared configuration.

The primary regression follows the transferred analysis definition:

```text
Sigma_core^ic = A (Sigma_env^ic)^alpha
```

It is fitted in log10 space with uncertainties in both axes and an intrinsic
scatter term. The primary profile assigns a 33% uncertainty to the
environment measurement. A 20% alternative is retained with neighboring
choices in an explicit sensitivity grid. For a fractional linear uncertainty
`f`, the implemented log uncertainty is `f / ln(10)`.

The current opacity correction scales both intercore axes together. It changes
the fitted normalization but leaves the correlation coefficients and
power-law slope unchanged.

Run:

```bash
conda run -n py312 python scripts/correlation.py --stage inspect
conda run -n py312 python scripts/correlation.py --stage all
```

Use `--stage matrices` for the pairwise table and two matrix figures, or
`--stage regression` for the MCMC regression only. Add `--overwrite` when
replacing existing task products.

Outputs under `results/correlation/` are:

- `tables/pairwise_correlations.ecsv`;
- `tables/intercore_surface_density_powerlaw.json`;
- `posterior/intercore_surface_density_powerlaw.npz`;
- detailed corner-plot PDF and PNG;
- compact Spearman-matrix PDF and PNG; and
- regression-comparison PDF and PNG; and
- the clean paper figure
  `figures/Correlation_Surface_env-Surface_core.pdf` with a matching PNG.
