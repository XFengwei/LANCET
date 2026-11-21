"""Core-temperature sensitivity of the regional mass fractions.

Stages of ``core_formation_efficiency`` (moved from the R8d referee scripts):

- ``run_stress_test``: deterministic heating from Tenv to the regional upper
  temperatures, keeping the foreground/background at Tenv; ratios are
  anchored to the adopted physical catalogue.
- ``plot_stress_test`` and ``plot_sampling``: the heating-grid figure, the
  prior comparison of the temperature sampling and the two Fig. 9 layouts.
- ``write_reports``: Markdown reports of both tests from their products.

The adopted catalogues keep Tenv; these are sensitivity scenarios.
"""
from collections.abc import Mapping, Sequence
import hashlib
import itertools
import json
from pathlib import Path
import sys

from astropy import units as u
from astropy.table import Table
import numpy as np
from scipy.stats import pearsonr, spearmanr

from lancet.physical_properties import (
    blackbody_cgs,
    core_optical_depth,
    resolved_core_mask,
    thick_mass,
)

ENVIRONMENT_FIELDS = ["Surface_intercorebkg", "Temp_env_median", "Vsig_env_median",
                      "Surface_env_median", "Interco_env_median"]
COLORS = ["#007C91", "#D77B25", "#76549F"]
MARKERS = ["o", "s", "^"]


def _alignment(qa_scripts: Path | None):
    """Optional panel-alignment audit (audit_panel_alignment.py)."""
    if qa_scripts is None:
        return None
    sys.path.insert(0, str(qa_scripts))
    from audit_panel_alignment import require_matplotlib_panel_alignment

    return require_matplotlib_panel_alignment


def _save(fig, stem: Path, audit) -> None:
    fig.canvas.draw()
    if audit is not None:
        audit(fig, json_out=f"{stem}.alignment.json",
              overlay_svg=f"{stem}.alignment.svg", tolerance_pt=1.5, strict=True)
    fig.savefig(f"{stem}.pdf", bbox_inches="tight", pad_inches=.025)
    fig.savefig(f"{stem}.svg", bbox_inches="tight", pad_inches=.025)
    fig.savefig(f"{stem}.png", dpi=600, bbox_inches="tight", pad_inches=.025)


def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "sans-serif", "font.size": 10,
        "pdf.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": .8,
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False})
    return plt


def run_stress_test(
    inputs: Mapping[str, Path],
    physics: Mapping[str, object],
    caps_k: Sequence[float],
    out: Path,
    *,
    grid_points: int = 101,
    scenario_h: Sequence[float] = (0., .25, .5, .75, 1.),
) -> dict[str, object]:
    """Deterministic core heating; ``inputs`` maps source, physical,
    regional and environment to the catalogues and summary."""
    out.mkdir(parents=True, exist_ok=True)
    source, physical, regional, environment = [Table.read(inputs[k])
        for k in ("source", "physical", "regional", "environment")]
    for table in (source, environment):
        assert len(set(table["ID"])) == len(physical)
        order = {int(row["ID"]): i for i, row in enumerate(table)}
        table[:] = table[[order[int(i)] for i in physical["ID"]]]
        assert np.array_equal(table["ID"], physical["ID"])
    nu = float(physics["frequency_ghz"]) * u.GHz
    temp = np.asarray(physical["Tenv"], dtype=float)
    code = np.asarray(physical["Subregion"])
    assert np.allclose(temp, source["Tenv"])
    assert np.array_equal(code, source["Subregion"])
    caps = dict(zip(("Y", "I", "E"), (float(c) for c in caps_k)))
    ceiling = np.array([caps[r] for r in code])
    assert np.all(ceiling >= temp), "Heating endpoint must not cool any core"
    beam = 1.1331 * float(physics["beam_major_arcsec"]) * float(physics["beam_minor_arcsec"]) * u.arcsec**2
    peak = (np.asarray(source["FXP_BST01"]) * u.Jy / beam).to_value(u.MJy / u.sr)
    background = np.asarray(source["Ibkg2"], dtype=float)
    benv = blackbody_cgs(nu, temp * u.K).to_value(u.MJy / u.sr)
    transmission = np.sqrt(1 - background / benv)
    rear = background / (transmission + 1)

    def evaluate(tcore):
        bc = blackbody_cgs(nu, tcore * u.K).to_value(u.MJy / u.sr)
        contrast = peak / ((bc - rear) * transmission)
        assert np.all((contrast > 0) & (contrast < 1))
        tau = -np.log1p(-contrast)
        correction = tau / -np.expm1(-tau)
        return bc, tau, correction

    _, tau0, f0 = evaluate(temp)
    legacy = core_optical_depth(nu=nu, temp=temp*u.K,
        Isr=peak*u.MJy/u.sr, Ibg=background*u.MJy/u.sr).value
    assert np.allclose(tau0, legacy, rtol=1e-10)
    mass0 = np.asarray(physical["Mass_core"], dtype=float)
    surface0 = np.asarray(physical["Surface_core"], dtype=float)
    density0 = np.asarray(physical["Volden_core"], dtype=float)
    resolved = resolved_core_mask(physical)
    rows, core_rows, correlation_rows = [], [], []
    grid = np.linspace(0, 1, int(grid_points))
    scenarios = {float(h) for h in scenario_h}
    mass_grid = []
    for h in grid:
        tcore = temp + h * (ceiling - temp)
        bc, tau, correction = evaluate(tcore)
        ratio = benv / bc * correction / f0
        mass = mass0 * ratio
        surface = surface0 * tau / tau0
        density = density0 * ratio
        mass_grid.append(mass)
        assert np.all(np.isfinite(mass)) and np.all(mass > 0)
        assert np.all(mass <= mass0 * (1 + 1e-12))
        if h == 0:
            assert np.allclose(mass, mass0, rtol=1e-12)
        # Independent dimensional check using the maintained mass function.
        if h in (0., 1.):
            kwargs = dict(nu=nu, kappa=float(physics["dust_opacity_cm2_per_g"]) /
                          float(physics["gas_to_dust_ratio"]) * u.cm**2/u.g,
                          ftot=np.asarray(source["FXT_BST01"])*u.Jy,
                          dist=float(physics["distance_kpc"])*u.kpc)
            direct = thick_mass(temp=tcore*u.K, tau_c=tau, **kwargs)
            base = thick_mass(temp=temp*u.K, tau_c=tau0, **kwargs)
            assert np.allclose((direct/base).value, ratio, rtol=1e-12)
        for reg in regional:
            r = str(reg["subregion_code"])
            select = code == r
            total = float(mass[select].sum())
            if h == 0:
                assert np.isclose(total, reg["core_mass"], rtol=1e-12)
            rows.append(dict(h=float(h), region=str(reg["label"]), code=r,
                core_count=int(select.sum()), endpoint_k=caps[r],
                core_mass_msun=total, retained_mass_fraction=total/mass0[select].sum(),
                core_cloud_percent=100*total/float(reg["cloud_mass"]),
                core_fragment_percent=100*total/float(reg["fragment_mass"]),
                fragment_cloud_percent=100*float(reg["fragment_mass"])/float(reg["cloud_mass"])))
        if float(h) in scenarios:
            for i in range(len(physical)):
                core_rows.append(dict(h=float(h), ID=int(physical["ID"][i]),
                    Resolvability=str(physical["Resolvability"][i]),
                    code=str(code[i]), Tenv_k=temp[i], Tcore_k=tcore[i],
                    mass_ratio=ratio[i], mass_msun=mass[i], tau=tau[i],
                    tau_correction=correction[i], surface_g_cm2=surface[i],
                    volume_density_catalogue_units=density[i]))
            for prop, values in (("Mass_core", mass), ("Surface_core", surface),
                                 ("Volden_core", density)):
                for env in ENVIRONMENT_FIELDS:
                    x = np.asarray(environment[env], dtype=float)
                    good = np.isfinite(x) & np.isfinite(values) & (x > 0) & (values > 0)
                    if prop == "Volden_core":
                        good &= resolved
                    correlation_rows.append(dict(h=float(h), core_property=prop,
                        environment_property=env, n=int(good.sum()),
                        spearman_rho=float(spearmanr(x[good], values[good]).statistic),
                        pearson_log_r=float(pearsonr(np.log10(x[good]), np.log10(values[good])).statistic)))
            for first, second in itertools.combinations(
                    (("Mass_core", mass), ("Surface_core", surface), ("Volden_core", density)), 2):
                good = resolved if "Volden_core" in {first[0], second[0]} else np.ones(len(mass), dtype=bool)
                correlation_rows.append(dict(h=float(h), core_property=first[0],
                    environment_property=second[0], n=int(good.sum()),
                    spearman_rho=float(spearmanr(first[1][good], second[1][good]).statistic),
                    pearson_log_r=float(pearsonr(np.log10(first[1][good]), np.log10(second[1][good])).statistic)))
    mass_grid = np.asarray(mass_grid)
    assert np.all(np.diff(mass_grid, axis=0) <= 1e-12)
    Table(rows=rows).write(out/"regional_temperature_grid.ecsv", overwrite=True)
    Table(rows=core_rows).write(out/"core_temperature_scenarios.ecsv", overwrite=True)
    Table(rows=correlation_rows).write(out/"correlation_sensitivity.ecsv", overwrite=True)
    corners = []
    for flags in itertools.product((0, 1), repeat=3):
        result = dict(zip(("Y_hot", "I_hot", "E_hot"), flags))
        cloud, fragment = [], []
        for reg, flag in zip(regional, flags):
            r = str(reg["subregion_code"])
            total = float(mass_grid[-1 if flag else 0, code == r].sum())
            result[f"{r}_core_cloud_percent"] = 100*total/float(reg["cloud_mass"])
            result[f"{r}_core_fragment_percent"] = 100*total/float(reg["fragment_mass"])
            cloud.append(result[f"{r}_core_cloud_percent"])
            fragment.append(result[f"{r}_core_fragment_percent"])
        result["core_cloud_Y_lt_I_lt_E"] = bool(cloud[0] < cloud[1] < cloud[2])
        result["core_fragment_Y_gt_I_gt_E"] = bool(fragment[0] > fragment[1] > fragment[2])
        corners.append(result)
    Table(rows=corners).write(out/"independent_regional_endpoints.ecsv", overwrite=True)
    hot = Table(rows=[row for row in core_rows if row["h"] == 1.])
    root = Path.cwd().resolve()
    summary = dict(
        command="PYTHONPATH=src conda run -n py312 python tasks/core_formation_efficiency/run.py --stage stress_test --overwrite",
        input_sha256={name: hashlib.sha256(Path(path).read_bytes()).hexdigest()
                      for name, path in inputs.items()},
        inputs={name: str(Path(path).resolve().relative_to(root)) if Path(path).resolve().is_relative_to(root) else str(path)
                for name, path in inputs.items()},
        caps_k=caps, n_cores=len(physical),
        changed="Core temperature only; every catalogue source retained",
        held_fixed=["Observed peak and total core fluxes", "Source sizes and selection",
                    "Foreground/background temperature and observed intensities",
                    "Dust opacity, gas-to-dust ratio and distance",
                    "Fragment and cloud masses; no free-free correction applied"],
        interpretation="Deterministic sensitivity scenarios, not confidence intervals or measured temperature bounds",
        baseline=[row for row in rows if row["h"] == 0.],
        hot_endpoint=[row for row in rows if row["h"] == 1.],
        temperature_ranges={r: dict(min=float(temp[code==r].min()),
            median=float(np.median(temp[code==r])), max=float(temp[code==r].max())) for r in caps},
        mass_statistics=dict(baseline_min=float(mass0.min()), baseline_median=float(np.median(mass0)),
            baseline_max=float(mass0.max()), hot_min=float(hot["mass_msun"].min()),
            hot_median=float(np.median(hot["mass_msun"])), hot_max=float(hot["mass_msun"].max())),
        tau_statistics=dict(baseline_max=float(tau0.max()), hot_max=float(hot["tau"].max()),
            hot_max_correction=float(hot["tau_correction"].max())),
        coherent_cloud_order_all_grid_points=all(
            rows[i]["core_cloud_percent"] < rows[i+1]["core_cloud_percent"] < rows[i+2]["core_cloud_percent"]
            for i in range(0, len(rows), 3)),
        coherent_fragment_order_all_grid_points=all(
            rows[i]["core_fragment_percent"] > rows[i+1]["core_fragment_percent"] > rows[i+2]["core_fragment_percent"]
            for i in range(0, len(rows), 3)),
        independent_cloud_order_all_corners=all(c["core_cloud_Y_lt_I_lt_E"] for c in corners),
        independent_fragment_order_all_corners=all(c["core_fragment_Y_gt_I_gt_E"] for c in corners),
    )
    (out/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(f"Wrote {out/'summary.json'}")
    return summary


def plot_stress_test(out: Path, *, qa_scripts: Path | None = None) -> None:
    """Regional core/cloud and core/fragment fractions against heating h."""
    plt = _style()
    data = Table.read(out / "regional_temperature_grid.ecsv")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.85))
    fig.subplots_adjust(left=.09, right=.99, bottom=.20, top=.97, wspace=.38)
    metrics = ["core_cloud_percent", "core_fragment_percent"]
    for ax, metric, label in zip(axes, metrics,
            ["Core / cloud mass (%)", "Core / fragment mass (%)"]):
        top = 0.
        for code, color in zip(["Y", "I", "E"], COLORS):
            rows = data[data["code"] == code]
            ax.plot(rows["h"], rows[metric], color=color, linewidth=1.8,
                    label=f'{rows[0]["region"]} ({rows[0]["endpoint_k"]:g} K)')
            ax.scatter([0], [rows[metric][0]], s=22, facecolors="white",
                       edgecolors=color, linewidths=1.1, zorder=3)
            ax.scatter([1], [rows[metric][-1]], s=22, color=color, zorder=3)
            top = max(top, float(np.max(rows[metric])))
        ax.set(xlim=(-.035, 1.035), ylim=(0, 1.15 * top), ylabel=label,
               xlabel="Core-heating fraction, h")
        ax.set_xticks([0, .25, .5, .75, 1])
        ax.tick_params(direction="in", length=3)
        ax.tick_params(axis="x", pad=7)
    for ax, letter in zip(axes, ["(a)", "(b)"]):
        ax.text(.025, .98, letter, transform=ax.transAxes, ha="left", va="top")
    axes[1].legend(loc="upper right", fontsize=8.5, handlelength=1.5)
    _save(fig, out / "temperature_stress_test", _alignment(qa_scripts))
    plt.close(fig)
    print(f"Wrote {out/'temperature_stress_test.png'}")


def plot_sampling(
    sampling: Path,
    regional_summary: Path,
    *,
    qa_scripts: Path | None = None,
) -> None:
    """Prior comparison (independent and shared heating) and the Fig. 9 layouts."""
    from matplotlib.lines import Line2D
    from lancet.core_formation import draw_temperature_efficiency_cascade

    plt = _style()
    audit = _alignment(qa_scripts)
    table = Table.read(sampling / "regional_sampling_summary.ecsv")
    original = Table.read(regional_summary)
    scales = np.array([4., .15, .01])
    for coupling in ("independent", "regional_shared"):
        fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=True)
        fig.subplots_adjust(left=.085, right=.99, bottom=.21, top=.98, wspace=.14)
        for ax, prior, title, letter in zip(axes, ["uniform", "brightness"],
                ["Uniform sampling", "Brightness-dependent"], ["(a)", "(b)"]):
            selected = table[(table["coupling"] == coupling) & (table["prior"] == prior)]
            for row, color, marker in zip(selected, COLORS, MARKERS):
                fragment = row["fragment_mass_msun"]/row["cloud_mass_msun"]
                low, middle, high = [row[f"core_cloud_{p}_percent"]/100 for p in ("p16", "p50", "p84")]
                baseline = row["baseline_mass_msun"]/row["cloud_mass_msun"]
                fractions = np.array([1., fragment, middle])
                assert low > 0 and high >= middle >= low and np.all(fractions > 0)
                ax.plot(scales, fractions, color=color, marker=marker, ms=5, lw=1.5,
                        label=str(row["region"]))
                ax.errorbar([scales[-1]], [middle], yerr=[[middle-low], [high-middle]],
                            fmt="none", color=color, lw=1.3, capsize=4, zorder=4)
                ax.plot(scales[1:], [fragment, baseline], color=color,
                        linestyle="--", alpha=.55, lw=.9)
                ax.scatter([scales[-1]], [baseline], marker=marker, s=32,
                           facecolors="white", edgecolors=color, lw=1, zorder=4)
            ax.set(xscale="log", yscale="log", xlim=(6., .005), ylim=(.0006, 3.))
            ax.set_xticks(scales, ["Cloud\n4 pc", "Fragment\n0.15 pc", "Core\n0.01 pc"])
            ax.set_yticks([.001, .01, .1, 1], ["0.1%", "1%", "10%", "100%"])
            ax.tick_params(direction="in", axis="both", pad=5)
            ax.text(.025, .975, letter, transform=ax.transAxes, ha="left", va="top")
            ax.text(.98, .975, title, fontsize=9, transform=ax.transAxes, ha="right", va="top")
            ax.set_xlabel("Physical scale")
        axes[0].set_ylabel("Mass / cloud mass")
        axes[0].legend(loc="lower left", fontsize=9, handlelength=1.8)
        handles = [Line2D([], [], color="0.3", marker="o", mfc="white", ls="--", lw=.9, label="Adopted Tenv"),
                   Line2D([], [], color="0.3", marker="o", ls="-", lw=1.5, label="Sampled median")]
        axes[1].legend(handles=handles, loc="lower left", fontsize=8.5)
        _save(fig, sampling / f"mass_efficiency_{coupling}", audit)
        plt.close(fig)

    trials = sampling / "figure9_trials"
    trials.mkdir(parents=True, exist_ok=True)
    for layout in ("overlay", "comparison"):
        if layout == "overlay":
            fig, ax = plt.subplots(figsize=(3.9, 3.8))
            fig.subplots_adjust(left=.17, right=.98, bottom=.21, top=.98)
            draw_temperature_efficiency_cascade(ax, table, original)
        else:
            fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6))
            fig.subplots_adjust(left=.085, right=.985, bottom=.22, top=.98, wspace=.36)
            draw_temperature_efficiency_cascade(axes[0], table, original)
            ax = axes[1]
            for j, (base, color, marker) in enumerate(zip(original, COLORS, MARKERS)):
                offset = (j-1)*.13
                y = [100*base["fragment_to_core_efficiency"]]
                lows, highs = [0.], [0.]
                for prior in ("uniform", "brightness"):
                    r = table[(table["coupling"] == "independent") &
                              (table["prior"] == prior) & (table["code"] == base["subregion_code"])][0]
                    p16, p50, p84 = [r[f"core_fragment_{p}_percent"] for p in ("p16", "p50", "p84")]
                    y.append(p50)
                    lows.append(p50-p16)
                    highs.append(p84-p50)
                ax.plot(np.arange(3)+offset, y, color=color, lw=.8, alpha=.45)
                ax.scatter([offset], y[:1], marker=marker, s=30, facecolors="white", edgecolors=color)
                ax.errorbar(np.array([1., 2.])+offset, y[1:], yerr=[lows[1:], highs[1:]],
                    fmt=marker, ms=5, color=color, mfc=color, mec="white", mew=.4,
                    lw=1.2, capsize=3)
            ax.set(xlim=(-.4, 2.4), ylim=(0, 61), ylabel="Core / fragment mass (%)")
            ax.set_xticks([0, 1, 2], ["Adopted\nTenv", "Uniform", "Brightness\ndependent"])
            ax.set_xlabel("Core temperature model")
            ax.tick_params(direction="in", pad=5)
            for a, letter in zip(axes, ["(a)", "(b)"]):
                a.text(.035, .975, letter, transform=a.transAxes, va="top", ha="left")
        _save(fig, trials / f"Figure9_trial_{layout}", audit)
        plt.close(fig)
    print(f"Wrote {sampling/'mass_efficiency_independent.png'} and {trials}")


def _mdtable(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"]*len(headers)) + " |"] +
        ["| " + " | ".join(map(str, row)) + " |" for row in rows])


def _interval(row, prefix, decimals):
    return (f"{row[f'{prefix}_p50_msun']:.{decimals}f} "
            f"[{row[f'{prefix}_p16_msun']:.{decimals}f}, {row[f'{prefix}_p84_msun']:.{decimals}f}]")


def write_reports(stress: Path, sampling: Path, regional_summary: Path) -> None:
    """Markdown reports of the stress test and the temperature sampling."""
    summary = json.loads((stress/"summary.json").read_text())
    meta = json.loads((sampling/"metadata.json").read_text())
    caps = list(summary["caps_k"].values())
    assert caps == list(meta["settings"]["upper_temperatures_k"])
    captext = "/".join(f"{v:g}" for v in caps)
    regional = Table.read(regional_summary)
    fragment_cloud = ", ".join(f'{100*r["fragment_efficiency"]:.3f}%' for r in regional)
    baseline_sums = ", ".join(f'{r["core_mass"]:.2f}' for r in regional)
    base, hot = summary["baseline"], summary["hot_endpoint"]
    endpoints = [[b["region"], f'{h["endpoint_k"]:g}',
        f'{b["core_mass_msun"]:.2f} → {h["core_mass_msun"]:.2f}',
        f'{100*(1-h["retained_mass_fraction"]):.1f}%',
        f'{b["core_cloud_percent"]:.3f}% → {h["core_cloud_percent"]:.3f}%',
        f'{b["core_fragment_percent"]:.2f}% → {h["core_fragment_percent"]:.2f}%'] for b, h in zip(base, hot)]
    bounds = []
    for i, j in itertools.combinations(range(3), 2):
        bounds.append([f'{base[i]["region"]} / {base[j]["region"]}',
            str(base[i]["core_cloud_percent"] < hot[j]["core_cloud_percent"]),
            str(hot[i]["core_fragment_percent"] > base[j]["core_fragment_percent"])])
    corr = Table.read(stress/"correlation_sensitivity.ecsv")
    corrows = []
    for a, b in [("Mass_core", "Surface_core"), ("Surface_core", "Volden_core"),
                 ("Mass_core", "Volden_core"), ("Surface_core", "Surface_intercorebkg")]:
        pair = corr[(corr["core_property"] == a) & (corr["environment_property"] == b)]
        cold, warm = pair[pair["h"] == 0][0], pair[pair["h"] == 1][0]
        corrows.append([a+" / "+b, str(cold["n"]), f'{cold["spearman_rho"]:.3f} → {warm["spearman_rho"]:.3f}'])
    ms, ts = summary["mass_statistics"], summary["tau_statistics"]
    (stress/"report.md").write_text(f"""# Deterministic core-temperature stress test

Upper temperatures: **{captext} K** (Young/Intermediate/Evolved). All
{summary['n_cores']} sources are retained; the adopted catalogue keeps Tenv.

For each core, T(h)=Tenv+h(Tupper−Tenv), with h=0…1. The core temperature
varies; the environment, measured background/foreground brightness, source
flux, size, selection, opacity, gas-to-dust ratio, distance and
fragment/cloud masses stay fixed. There is no free–free correction.

With q=sqrt(1−Ibg/Bnu(Tenv)) and Irear=Ibg/(1+q),
tau(T)=−ln[1−Icore/(q×(Bnu(T)−Irear))] and f(tau)=tau/(1−exp(−tau)); each
adopted mass scales by [Bnu(Tenv)/Bnu(T)]×[f(tau(T))/f(tau(Tenv))], so h=0
reproduces the catalogue exactly. Surface densities scale by the tau ratio
and volume densities by the mass ratio at fixed size.

## Endpoint results

{_mdtable(['Region', 'Tupper (K)', 'Core mass (solar masses)', 'Mass decrease', 'Core/cloud', 'Core/fragment'], endpoints)}

Coherent core/cloud ordering at every h:
**{summary['coherent_cloud_order_all_grid_points']}**; coherent core/fragment
ordering: **{summary['coherent_fragment_order_all_grid_points']}**.
Fragment/cloud fractions stay at {fragment_cloud}.

Median core mass: {ms['baseline_median']:.3f} → {ms['hot_median']:.3f} solar
masses; maximum: {ms['baseline_max']:.2f} → {ms['hot_max']:.2f}. Maximum tau
at the hot endpoint {ts['hot_max']:.5f}, maximum correction
{ts['hot_max_correction']:.5f}.

## Independently varying regional temperatures

The eight endpoint corners are in `independent_regional_endpoints.ecsv`.
Core/cloud ordering at all corners:
**{summary['independent_cloud_order_all_corners']}**; core/fragment ordering:
**{summary['independent_fragment_order_all_corners']}**.

{_mdtable(['Earlier / later region', 'Core/cloud always earlier < later', 'Core/fragment always earlier > later'], bounds)}

## Correlations

Spearman coefficients (baseline → hot endpoint):

{_mdtable(['Pair', 'n', 'rho'], corrows)}

The full table (`correlation_sensitivity.ecsv`) is recomputed at the
scenario points; pairs with core volume density use resolved cores.
""")
    table = Table.read(sampling/"regional_sampling_summary.ecsv")
    ordering = Table.read(sampling/"ordering_draw_fractions.ecsv")
    labels = {"uniform": "Uniform", "brightness": "Brightness-dependent"}
    massrows = []
    for prior in labels:
        for coupling in ("independent", "regional_shared"):
            selected = table[(table["prior"] == prior) & (table["coupling"] == coupling)]
            massrows.append([labels[prior], coupling]+[_interval(r, "mass", 2) for r in selected])
    fractionrows = []
    for r in table[table["coupling"] == "independent"]:
        vals = []
        for prefix in ("core_cloud", "core_fragment"):
            decimals = 3 if prefix == "core_cloud" else 2
            vals.append(f'{r[prefix+"_p50_percent"]:.{decimals}f} [{r[prefix+"_p16_percent"]:.{decimals}f}, {r[prefix+"_p84_percent"]:.{decimals}f}]')
        fractionrows.append([labels[r["prior"]], r["region"]]+vals)
    orderrows = [[labels[r["prior"]], r["coupling"], f'{100*r["core_cloud_Y_lt_I_lt_E"]:.3f}%',
        f'{100*r["core_fragment_Y_gt_I_gt_E"]:.3f}%'] for r in ordering]
    maxsplit = max(table["split_quantile_difference_over_68_width"])*100
    (sampling/"report.md").write_text(f"""# Sampled core temperatures and mass efficiencies

Upper temperatures: **{captext} K** (Young/Intermediate/Evolved). The
adopted catalogue keeps Tenv; Fig. 9 overlays the brightness-prior result.

Tcore_i=Tenv_i+F_i(Tupper_r−Tenv_i). Uniform: F_i~Uniform(0,1). Brightness:
q_i is the global peak-brightness percentile; F_i~Beta(1+2q_i,3−2q_i), mean
0.25+0.5q_i. {meta['settings']['sample_count']} joint realizations per
prior/dependence choice, seed {meta['settings']['random_seed']}. Optical
depths are recomputed with the foreground/background at Tenv; fluxes, sizes,
selection, distance, opacity, gas-to-dust ratio and fragment/cloud masses
are fixed. Independent: every core drawn separately. Regional shared: one
percentile per region, applied through each core's inverse beta CDF.

## Core-mass sums

Median [16th, 84th percentile], solar masses (Tenv baseline: {baseline_sums}):

{_mdtable(['Prior', 'Dependence', 'Young', 'Intermediate', 'Evolved'], massrows)}

## Mass fractions, independent draws

Median [16th, 84th percentile], percent (fragment/cloud fixed at
{fragment_cloud}):

{_mdtable(['Prior', 'Region', 'Core/cloud (%)', 'Core/fragment (%)'], fractionrows)}

## Ordering fractions

Fractions of the sampled realizations, not observational confidence:

{_mdtable(['Prior', 'Dependence', 'Core/cloud Y<I<E', 'Core/fragment Y>I>E'], orderrows)}

## Checks

Baseline regional masses and source IDs match; every sampled mass lies
between its temperature endpoints; sampled means agree with a 160-point
Gauss–Legendre integration within six Monte Carlo standard errors;
split-half 16/50/84 differences are at most {maxsplit:.2f}% of the 16–84
width. Settings and input checksums are in `metadata.json`.
""")
    trials = sampling / "figure9_trials"
    trials.mkdir(parents=True, exist_ok=True)
    (trials/"README.md").write_text(f"""# Figure 9 layouts

Upper temperatures {captext} K. `Figure9_trial_overlay` is the single-panel
overlay drawn by the `efficiency` stage; `Figure9_trial_comparison` adds the
core/fragment fractions under the uniform and brightness-dependent priors.
""")
    print(f"Wrote {stress/'report.md'} and {sampling/'report.md'}")
