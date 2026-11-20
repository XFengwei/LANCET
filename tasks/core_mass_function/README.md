# Core mass function

This task estimates the global and evolutionary-subregion CMFs from the
generated physical-properties catalogue. Rows are masked once so IDs,
subregions, masses, and auxiliary columns remain aligned.

The fitted convention is
`dN/dM proportional to M^(-alpha)`. For each sample, the lower fitting limit
is selected from observed masses by minimizing the KS distance, subject to a
minimum tail count. No external completeness floor is currently imposed; set
`fit.completeness_mass_msun` to a positive value if a completeness experiment
provides one. The exponent is estimated analytically and with a one-parameter
Metropolis posterior. The bootstrap holds the selected lower limit fixed, so
it does not include uncertainty in `Mmin`.

```bash
PYTHONPATH=src conda run -n py312 python tasks/core_mass_function/run.py --stage inspect
PYTHONPATH=src conda run -n py312 python tasks/core_mass_function/run.py --stage all
```

Use `--stage fit` or `--stage plot` for one product group. Add `--overwrite`
when replacing products. The fit tables, posterior samples, safe catalogue
export, and figures are written under `results/core_mass_function/`; no
external project path is written.

## Change from Paper II (2026-11-20)

The code reproduces the Paper II products exactly from the Paper II
`physical_properties` catalogue (all tables and posterior samples
byte-identical). With the reviewed catalogue (new 12m+7m+TP background,
gain-corrected masses for 19 cores; `tables/powerlaw_fits.ecsv`):

| Sample | Mmin (M☉) | N tail | α (MLE) | α posterior median [16, 84 %] |
|---|---|---|---|---|
| All | 1.751 → 1.740 | 92 → 94 | 2.415 → 2.407 | 2.41 [2.28, 2.56] → 2.41 [2.28, 2.55] |
| Young | 1.332 | 12 | 2.323 → 2.269 | 2.36 [2.01, 2.76] → 2.30 [1.97, 2.68] |
| Intermediate | 1.295 → 1.233 | 63 → 68 | 2.494 → 2.461 | 2.52 [2.34, 2.70] → 2.48 [2.32, 2.66] |
| Evolved | 1.751 | 42 | 2.316 | 2.34 [2.14, 2.56], unchanged |
