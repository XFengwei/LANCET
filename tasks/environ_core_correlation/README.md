# Core–environment correlations

`variables.resolved_only_columns` applies per-quantity censoring. Pairs using
core volume density require `Resolvability=Y` (209 sources before
positive/finite filtering). Beam-averaged surface density, mass and environmental pairs use all valid
sources. The same selection is applied to the data used for statistics,
histograms, scatters, contours and any regression whose axes require it.
The intercore–intercore fit retains its original 228 valid positions.

This task measures pairwise correlations among four core properties and four
environment properties, then fits the intercore surface-density relation with
a power law. Its input is the generated, opacity-consistent environment
catalogue:

`results/environ_measure/tables/core_environment_catalogue.ecsv`

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

The opacity scales both intercore axes together. It changes the fitted
normalization but leaves the correlation coefficients and power-law slope
unchanged.

Run:

```bash
PYTHONPATH=src conda run -n py312 python tasks/environ_core_correlation/run.py --stage inspect
PYTHONPATH=src conda run -n py312 python tasks/environ_core_correlation/run.py --stage all
```

Use `--stage matrices` for the pairwise table and two matrix figures, or
`--stage regression` for the MCMC regression only. Add `--overwrite` when
replacing existing task products.

Outputs under `results/environ_core_correlation/` are:

- `tables/pairwise_correlations.ecsv`;
- `tables/intercore_surface_density_powerlaw.json`;
- `posterior/intercore_surface_density_powerlaw.npz`;
- detailed corner-plot PDF and PNG;
- compact Spearman-matrix PDF and PNG; and
- regression-comparison PDF and PNG; and
- the clean paper figure
  `figures/correlation_surface_env-surface_core.pdf` with a matching PNG.

## Change from Paper II (2026-11-20)

The code reproduces the Paper II products from the Paper II environment
catalogue: the pairwise table is identical, OLS, Theil–Sen and ODR agree to
rounding. The MCMC differs at the 1e-3 level (α median 1.2231 → 1.2233)
because the maximum-likelihood start differs in the flat intrinsic-scatter
direction (ln σ −10.0 → −9.14, library versions); reruns are deterministic.

On the reviewed catalogue (gain-corrected cores and intercore map, new
12m+7m+TP core background; `tables/intercore_surface_density_powerlaw.json`):

- Σ^ic_core–Σ^ic_env: N = 228; Pearson 0.720 → 0.707, Spearman 0.722 → 0.705;
  MCMC α = 1.223 (+0.069/−0.065) → 1.242 (+0.072/−0.069), log₁₀A
  0.577 → 0.579.
- Sensitivity to the x uncertainty (0, 0.1, 0.2, 0.33, 0.5):
  α 1.26, 1.26, 1.26, 1.22, 1.01 → 1.29, 1.29, 1.29, 1.24, 1.01.
- Largest Spearman changes: Σ_core–Σ_env 0.19 → 0.16, M_core–Σ_env
  0.14 → 0.11, Σ^ic_env–Σ_env 0.56 → 0.53, Σ^ic_core–Σ_env 0.37 → 0.34.
