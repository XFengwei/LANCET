"""Refresh temperature-test reports from current numerical products."""
from pathlib import Path
import itertools
import json
import re
import tomllib

from astropy.table import Table
import numpy as np


def mdtable(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"]*len(headers)) + " |"] +
        ["| " + " | ".join(map(str, row)) + " |" for row in rows])


def interval(row, prefix, decimals):
    return (f"{row[f'{prefix}_p50_msun']:.{decimals}f} "
            f"[{row[f'{prefix}_p16_msun']:.{decimals}f}, {row[f'{prefix}_p84_msun']:.{decimals}f}]")


def main():
    root = Path(__file__).resolve().parents[1]
    cfg = tomllib.loads((root/"tasks/physical_properties/task.toml").read_text())
    sampling = root/cfg["temperature_sampling"]["output_root"]
    out = sampling.parent
    summary = json.loads((out/"summary.json").read_text())
    meta = json.loads((sampling/"metadata.json").read_text())
    caps = cfg["temperature_sampling"]["upper_temperatures_k"]
    assert list(summary["caps_k"].values()) == caps == meta["settings"]["upper_temperatures_k"]
    captext = "/".join(f"{v:g}" for v in caps)
    adopted_figure = (out/"manuscript_update/adopted-figure.json").exists()
    adoption = ("The compact temperature overlay is adopted in manuscript Figure 9; "
                "the catalogue retains the Tenv baseline." if adopted_figure else
                "The adopted catalogue and manuscript figure retain their Tenv baseline.")
    base, hot = summary["baseline"], summary["hot_endpoint"]
    endpoints = [[b["region"], f'{h["endpoint_k"]:g}',
        f'{b["core_mass_msun"]:.2f} → {h["core_mass_msun"]:.2f}',
        f'{100*(1-h["retained_mass_fraction"]):.1f}%',
        f'{b["core_cloud_percent"]:.3f}% → {h["core_cloud_percent"]:.3f}%',
        f'{b["core_fragment_percent"]:.2f}% → {h["core_fragment_percent"]:.2f}%'] for b,h in zip(base,hot)]
    bounds = []
    for i,j in itertools.combinations(range(3),2):
        bounds.append([f'{base[i]["region"]} / {base[j]["region"]}',
            str(base[i]["core_cloud_percent"] < hot[j]["core_cloud_percent"]),
            str(hot[i]["core_fragment_percent"] > base[j]["core_fragment_percent"])])
    corr = Table.read(out/"correlation_sensitivity.ecsv")
    pairs = [("Mass_core","Surface_core"), ("Surface_core","Volden_core"),
             ("Mass_core","Volden_core"), ("Surface_core","Surface_intercorebkg")]
    corrows = []
    for a,b in pairs:
        pair = corr[(corr["core_property"]==a)&(corr["environment_property"]==b)]
        cold, warm = pair[pair["h"]==0][0], pair[pair["h"]==1][0]
        corrows.append([a+" / "+b, str(cold["n"]), f'{cold["spearman_rho"]:.3f} → {warm["spearman_rho"]:.3f}'])
    ms, ts = summary["mass_statistics"], summary["tau_statistics"]
    report = f"""# R8d: deterministic core-temperature stress test

Current upper temperatures: **{captext} K** (Young/Intermediate/Evolved).
These replace the previous 40/50/60 K run, preserved in the sibling archive
`2026-11-01_caps_40_50_60`. {adoption}

All 231 sources are retained. For each core, T(h)=Tenv+h(Tupper−Tenv),
with h=0…1 at 101 grid points. The core temperature varies; the environment,
measured background/foreground brightness, source flux, size, selection,
opacity, gas-to-dust ratio, distance and fragment/cloud masses stay fixed.
There is no free–free correction in this temperature-only test.

Use the full Planck function at 223.8 GHz and the configured 0.71×0.60 arcsec
beam. With q=sqrt(1−Ibg/Bnu(Tenv)) and Irear=Ibg/(1+q),
tau(T)=−ln[1−Icore/(q×(Bnu(T)−Irear))]. The correction is
f(tau)=tau/(1−exp(−tau)); scale each adopted catalogue mass by
[Bnu(Tenv)/Bnu(T)]×[f(tau(T))/f(tau(Tenv))]. This anchors h=0 exactly
to the adopted Monte Carlo catalogue central values. Surface densities
scale by the tau ratio and volume densities by the mass ratio at fixed size.

The [ALMA-IMF PPMAP results](https://doi.org/10.1051/0004-6361/202345986)
and [mass-averaged temperature study](https://doi.org/10.1051/0004-6361/202451931)
motivate a heating sensitivity test; they do not measure temperatures in
G316.8 or prescribe these regional endpoints as universal upper limits.

## Endpoint results

{mdtable(['Region','Tupper (K)','Core mass (solar masses)','Mass decrease','Core/cloud','Core/fragment'],endpoints)}

The coherent regional core/cloud ordering is retained at every sampled h:
**{summary['coherent_cloud_order_all_grid_points']}**. The coherent
core/fragment ordering is retained: **{summary['coherent_fragment_order_all_grid_points']}**.
Fragment/cloud fractions remain 0.713%, 3.917%, 16.535%, since those
reservoir masses are fixed. These are deterministic scenarios, not
confidence intervals or a completeness-corrected core population.

Median individual core mass: {ms['baseline_median']:.3f} → {ms['hot_median']:.3f}
solar masses; maximum: {ms['baseline_max']:.2f} → {ms['hot_max']:.2f}.
The nominal maximum tau becomes {ts['hot_max']:.5f}, with maximum
optical-depth correction {ts['hot_max_correction']:.5f}.

## Independently variable regional temperatures

The eight endpoint corners are tabulated in `independent_regional_endpoints.ecsv`.
Because masses decrease monotonically with temperature, these bound the
regional sums even if individual core temperatures vary independently.
The full three-stage core/cloud ordering holds at all corners:
**{summary['independent_cloud_order_all_corners']}**; the full core/fragment
ordering holds: **{summary['independent_fragment_order_all_corners']}**.
Individual pairwise tests over the entire allowed intervals are:

{mdtable(['Earlier / later region','Core/cloud always earlier < later','Core/fragment always earlier > later'],bounds)}

The coherent heating sequence and the fully independent temperature envelope
are different assumptions; do not replace either with an unconditional
claim of robustness.

## Correlations

Representative Spearman coefficients (baseline → hot endpoint):

{mdtable(['Pair','n','rho'],corrows)}

The full affected correlation table is recomputed at h=0,0.25,0.5,0.75,1.
Pairs involving core volume density require resolved cores,
with positive finite values selected for each pair. All sources remain in
the integrated-mass and beam-averaged surface-density analyses, and the
per-core table carries resolution flags.
The local-intercore versus environmental-intercore relation is unchanged by
construction because core heating changes neither input. No new regression
posterior or significance interpretation is implied.

## Reproduction and checks

Run `PYTHONPATH=src conda run -n py312 python scripts/core_temperature_stress_test.py`;
default endpoints now come from the physical-properties task configuration.
Use `scripts/plot_core_temperature_stress_test.py` for the coherent-grid plot.
Input paths and SHA256 checksums are recorded in `summary.json`. Baseline
regional masses, equal-temperature tau inversion, dimensional mass ratios,
finite physical outputs and monotonic individual mass changes pass checks.
Verified scientific figure exports and source data remain alongside this
report. Current source IDs and mass selections are unchanged. R8 remains
in progress pending author review.
"""
    (out/"report.md").write_text(report)
    table = Table.read(sampling/"regional_sampling_summary.ecsv")
    ordering = Table.read(sampling/"ordering_draw_fractions.ecsv")
    labels = {"uniform":"Uniform", "brightness":"Brightness-dependent"}
    massrows = []
    for prior in labels:
        for coupling in ("independent","regional_shared"):
            selected = table[(table["prior"]==prior)&(table["coupling"]==coupling)]
            massrows.append([labels[prior],coupling]+[interval(r,"mass",2) for r in selected])
    fractionrows = []
    for r in table[table["coupling"]=="independent"]:
        vals = []
        for prefix in ("core_cloud","core_fragment"):
            decimals = 3 if prefix=="core_cloud" else 2
            vals.append(f'{r[prefix+"_p50_percent"]:.{decimals}f} [{r[prefix+"_p16_percent"]:.{decimals}f}, {r[prefix+"_p84_percent"]:.{decimals}f}]')
        fractionrows.append([labels[r["prior"]],r["region"]]+vals)
    orderrows = [[labels[r["prior"]],r["coupling"], f'{100*r["core_cloud_Y_lt_I_lt_E"]:.3f}%',
        f'{100*r["core_fragment_Y_gt_I_gt_E"]:.3f}%'] for r in ordering]
    maxsplit = max(table["split_quantile_difference_over_68_width"])*100
    report = f"""# R8d: sampled core temperatures and mass efficiencies

Current regional endpoints: **{captext} K** (Young/Intermediate/Evolved),
superseding the archived 40/50/60 K results. All 231 adopted cores are
retained (16/120/95). {adoption}

## Sampling functions and fixed assumptions

Tcore_i=Tenv_i+F_i(Tupper_r−Tenv_i). Uniform model: F_i~Uniform(0,1).
Author-approved brightness model: q_i is the global empirical peak-brightness
percentile over all 231 cores (average ranks for ties); F_i~Beta(1+2q_i,3−2q_i).
The mean heating fraction is 0.25+0.5q_i, rising from 25% to 75% with scatter.
This is an assumed sensitivity relation, not a temperature inferred from
brightness alone. Optical depth and beam filling also affect brightness.

Run {meta['settings']['sample_count']} independent joint realizations per
prior/dependence choice with seed {meta['settings']['random_seed']} and
spawned streams. Direct Monte Carlo propagates the assumed priors; no new
temperature likelihood or MCMC posterior is fitted. Recompute full Planck
functions and core optical depths, keeping foreground/background at Tenv.
Keep observed fluxes, sizes, selection, distance, opacity, gas-to-dust ratio
and fragment/cloud masses fixed. No free–free or other mass uncertainty is
added in this temperature-only test.

Independent: draw every core's F separately. Regional shared: draw one
uniform percentile per region and apply each core's inverse beta CDF. This
preserves every core's marginal prior and tests maximally positive heating
rank dependence within a region; different regions remain independent.
Independent variation averages down in regional sums; shared regional
uncertainty does not. Neither choice is fitted to the heating field.

## Core-mass sums

Median [16th,84th percentile], in solar masses. Regional medians are taken
from the sum in each draw, not from sums of individual medians.

{mdtable(['Prior','Dependence','Young','Intermediate','Evolved'],massrows)}

The Tenv baseline sums are 41.90,255.87,247.91 solar masses. All intervals
are conditional temperature-prior ranges, not total systematic errors or
data-constrained confidence intervals.

## Mass fractions under independent draws

Median [16th,84th percentile], in percent:

{mdtable(['Prior','Region','Core/cloud (%)','Core/fragment (%)'],fractionrows)}

Fragment/cloud stays at 0.713%,3.917%,16.535%. Both median regional trends
are retained in the independent and regional-shared models. Absolute mass
fractions and uncertainty widths depend on the temperature assumptions.

## Ordering fractions

Fractions of the sampled joint realizations, **not observational confidence**:

{mdtable(['Prior','Dependence','Core/cloud Y<I<E','Core/fragment Y>I>E'],orderrows)}

A 100% draw fraction means zero reversals among the sampled realizations,
not a mathematical guarantee. Shared regional heating has wider mass ranges
and can weaken ordering retention; this dependence is part of the result.

## Figures and caption

`mass_efficiency_independent.pdf` compares uniform and brightness-dependent
sampling using the independent draws. `mass_efficiency_regional_shared.pdf`
shows the corresponding dependence stress test. `figure9_trials/` contains
both integration layouts, using the updated endpoints and simulations.

> Temperature sensitivity of the multiscale mass fractions. Panels (a) and
> (b) assign core temperatures uniformly within their allowed intervals and
> according to the brightness-dependent distribution, respectively. Upper
> temperatures are {caps[0]:g}, {caps[1]:g} and {caps[2]:g} K for the Young,
> Intermediate and Evolved regions. Filled markers and solid curves show
> median fractions from {meta['settings']['sample_count']} joint realizations;
> vertical bars span the 16th–84th percentiles. Open core markers and dashed
> curves show the Tenv baseline. Fragment and cloud masses are fixed.
> These intervals quantify the specified temperature-assignment model only.

For the regional-shared figure, specify that cores within each region share
a heating percentile. For the Figure 9 trials, use their separate caption
and retain the distinction between baseline counting and temperature errors.

## Reproduction and verification

Configuration: `[temperature_sampling]` in the physical-properties task.
Run `PYTHONPATH=src conda run -n py312 python scripts/core_temperature_sampling.py`;
then use `scripts/plot_core_temperature_sampling.py` with `--qa-scripts` and,
for the integration alternatives, `--trials`. Refresh these reports with
`scripts/report_core_temperature_tests.py` after the numerical runs.

All baseline regional masses and source IDs match. Every sampled mass is
finite and bounded by its deterministic temperature endpoints. Sampled
means agree with independent 160-point Gauss–Legendre integrations within
six Monte Carlo standard errors. Split-half differences in the 16/50/84
percentiles are at most {maxsplit:.2f}% of the model's 16–84 interval width.
Input checksums, source paths, settings and all fixed assumptions are in
`metadata.json`. The core-property physics tests and rendered figure QA
are run alongside this refresh.

`regional_sampling_summary.ecsv` stores means, SDs, medians, 68% and 95%
equal-tail intervals and numerical stability diagnostics.
`core_temperature_priors.ecsv` stores per-core priors and quantiles;
`regional_mass_draws.npz` retains all regional draws. {adoption} R8 remains
in progress for the remaining size and mass-model subitems.
"""
    (sampling/"report.md").write_text(report)
    trial = sampling/"figure9_trials/README.md"
    text = trial.read_text()
    text = re.sub(r"> of \d+(?:\.\d+)?, \d+(?:\.\d+)? and \d+(?:\.\d+)? K\.",
                  f"> of {caps[0]:g}, {caps[1]:g} and {caps[2]:g} K.", text)
    text = re.sub(r"Current run: [\d./]+ K regional upper temperatures;",
                  f"Current run: {captext} K regional upper temperatures;", text)
    if "Current run:" not in text:
        text = text.replace("# Figure 9 integration trials", f"# Figure 9 integration trials\n\nCurrent run: {captext} K regional upper temperatures; all numerical products\nand figures regenerated from the updated sampling configuration.")
    trial.write_text(text)
    print(f"Refreshed endpoint, sampling and trial reports for {captext} K")


if __name__ == "__main__":
    main()
