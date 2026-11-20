# Bonnor–Ebert stability

The radius-based profile sets `radius.resolved_only=true` and uses the 209
resolved core sizes. Upper size limits are excluded from its point-estimate
critical masses, instability statistics and plots. External-pressure
profiles use their own configured environmental inputs and validity masks.

This task provides selectable, isolated calculation profiles. Every profile
writes into its own directory under `results/bonnor_ebert/`, so changing a
formula cannot silently overwrite a different result.

`radius_thermal` uses

```text
M_BE = 2.4 R a_support^2 / G
```

with thermal support and `R = 0.5 FWHM_linear`. The table now also records the
edge density and external pressure implied by treating that observed size as
the critical BE truncation radius. These are model-implied quantities, not
environment measurements.

`measured_external` uses

```text
rho_env = Sigma_ic(env) / L
sigma_nt,env^2 = Vsig_env_mean^2 - k_B T_env / (30 m_H)
P_external = rho_env * (c_s,env^2 + sigma_nt,env^2)

sigma_nt,int = f_turb * sigma_nt,env
a_support,int^2 = c_s,env^2 + sigma_nt,int^2

M_BE = 1.18 a_support,int^4 / (G^(3/2) P_external^(1/2))
```

Here `Sigma_ic(env)` is the aperture-mean `Interco_env_mean`, `L = 0.4 pc`,
and `T_env = Temp_env_mean`. `Vsig_env_mean` is treated as the observed
one-dimensional C18O moment-2 dispersion. The C18O thermal width is removed
before the gas thermal sound speed is added, so thermal broadening is not
counted twice in the pressure.

The internal temperature is currently the same `Temp_env_mean`. Internal
nonthermal support is a configurable fraction `f_turb` of the environmental
nonthermal dispersion. The configured fraction is a provisional dissipation
assumption, not a measurement, and is recorded in every result. Set it
anywhere from 0 (complete dissipation) to 1 (no dissipation) under
`profiles.measured_external.internal_nonthermal`.

The output catalogue records the adopted surface density, path length,
volume-averaged mass and H2 number densities, observed C18O dispersion,
derived environmental nonthermal and total dispersions, internal dispersion,
pressure, and turbulence fraction for every core.

The pressure profile can instead accept a directly measured pressure by
changing `external_pressure_mode` to `direct` and configuring
`direct_external_pressure` as either `pressure` or `pressure_over_k`.
Alternative density inputs remain supported as `number_density_h2` or
`mass_density`; alternative velocity inputs can be `nonthermal_sigma`,
`total_sigma`, `observed_sigma`, or `observed_fwhm`.

For every profile, `Alpha_BE = M_core / M_BE`; values at or above one are
reported as supercritical. Its current uncertainty propagates core-mass
uncertainty only. Uncertainties in temperature, surface density, the adopted
0.4-pc depth, velocity dispersion, and `f_turb` are not yet propagated.

```bash
PYTHONPATH=src conda run -n py312 python tasks/bonnor_ebert/run.py --stage inspect
PYTHONPATH=src conda run -n py312 python tasks/bonnor_ebert/run.py \
  --stage all --profile radius_thermal
PYTHONPATH=src conda run -n py312 python tasks/bonnor_ebert/run.py \
  --stage all --profile measured_external
```

The profile named by `calculation.active_profile` is used when `--profile` is
omitted. Use `--stage catalogue` or `--stage plot` for individual product
groups. Add `--overwrite` when replacing products within the same profile.

## Change from Paper II (2026-11-20)

The code reproduces both profiles from the Paper II environment catalogue:
`radius_thermal` exactly, `measured_external` within 3e-9 (astropy
constants). On the reviewed catalogue (gain-corrected masses and intercore
map, new 12m+7m+TP core background):

| | Paper II | New |
|---|---|---|
| `radius_thermal`: supercritical | 196/209 (0.938) | 198/209 (0.947) |
| `radius_thermal`: median α_BE [16, 84 %] | 2.98 [1.46, 10.1] | 3.00 [1.51, 10.1] |
| `measured_external`: supercritical | 224/231 (0.970) | 224/231 (0.970) |
| `measured_external`: median α_BE [16, 84 %] | 6.48 [2.41, 17.4] | 6.48 [2.53, 17.6] |
| `measured_external`: median P/k (K cm⁻³) | 4.94e7 | 4.95e7 |

By region (`radius_thermal` median α_BE): Young 8.55 → 8.55, Intermediate
2.92 → 3.00, Evolved 2.57 → 2.62.
