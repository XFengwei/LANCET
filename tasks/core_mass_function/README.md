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
conda run -n py312 python scripts/core_mass_function.py --stage inspect
conda run -n py312 python scripts/core_mass_function.py --stage all
```

Use `--stage fit` or `--stage plot` for one product group. Add `--overwrite`
when replacing products. The fit tables, posterior samples, safe catalogue
export, and figures are written under `results/core_mass_function/`; no
external project path is written.
