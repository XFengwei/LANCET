"""Diagnose core-temperature systematics without changing adopted catalogues.

Keep the foreground/background at Tenv; vary every core from Tenv to its
regional heating endpoint. Ratios are anchored to the adopted physical table.
"""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import tomllib

from astropy import units as u
from astropy.table import Table
import numpy as np
from scipy.stats import pearsonr, spearmanr

from lancet.physical_properties import blackbody_cgs, core_optical_depth, thick_mass, resolved_core_mask


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(
        "results/paper_ii_referee_response/diagnostics/R8d/2026-11-01"))
    parser.add_argument("--caps", nargs=3, type=float,
                        help="Override the upper temperatures in the physical-properties task",
                        metavar=("YOUNG_K", "INTERMEDIATE_K", "EVOLVED_K"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    out = root / args.output
    out.mkdir(parents=True, exist_ok=True)
    names = {
        "source": "results/source_catalogue/tables/source_catalogue.ecsv",
        "physical": "results/physical_properties/tables/core_physical_properties.ecsv",
        "regional": "results/core_formation/tables/core_formation_summary.ecsv",
        "environment": "results/environment_measurements/tables/core_environment_catalogue.ecsv",
        "config": "tasks/physical_properties/task.toml",
    }
    source, physical, regional, environment = [Table.read(root / names[k])
        for k in ("source", "physical", "regional", "environment")]
    for table in (source, environment):
        assert len(set(table["ID"])) == len(physical)
        order = {int(row["ID"]): i for i, row in enumerate(table)}
        table[:] = table[[order[int(i)] for i in physical["ID"]]]
        assert np.array_equal(table["ID"], physical["ID"])
    task_config = tomllib.loads((root / names["config"]).read_text())
    cfg = task_config["physics"]
    nu = cfg["frequency_ghz"] * u.GHz
    temp = np.asarray(physical["Tenv"], dtype=float)
    code = np.asarray(physical["Subregion"])
    assert np.allclose(temp, source["Tenv"])
    assert np.array_equal(code, source["Subregion"])
    caps = dict(zip(("Y", "I", "E"), args.caps or task_config["temperature_sampling"]["upper_temperatures_k"]))
    ceiling = np.array([caps[r] for r in code])
    assert np.all(ceiling >= temp), "Heating endpoint must not cool any core"
    beam = 1.1331 * cfg["beam_major_arcsec"] * cfg["beam_minor_arcsec"] * u.arcsec**2
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
    grid = np.linspace(0, 1, 101)
    fields = ["Surface_intercorebkg", "Temp_env_median", "Vsig_env_median",
              "Surface_env_median", "Interco_env_median"]
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
            kwargs = dict(nu=nu, kappa=cfg["dust_opacity_cm2_per_g"] /
                          cfg["gas_to_dust_ratio"] * u.cm**2/u.g,
                          ftot=np.asarray(source["FXT_BST01"])*u.Jy,
                          dist=cfg["distance_kpc"]*u.kpc)
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
        if h in (0., .25, .5, .75, 1.):
            for i in range(len(physical)):
                core_rows.append(dict(h=float(h), ID=int(physical["ID"][i]),
                    Resolvability=str(physical["Resolvability"][i]),
                    code=str(code[i]), Tenv_k=temp[i], Tcore_k=tcore[i],
                    mass_ratio=ratio[i], mass_msun=mass[i], tau=tau[i],
                    tau_correction=correction[i], surface_g_cm2=surface[i],
                    volume_density_catalogue_units=density[i]))
            for prop, values in (("Mass_core", mass), ("Surface_core", surface),
                                 ("Volden_core", density)):
                for env in fields:
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
    summary = dict(
        command="PYTHONPATH=src conda run -n py312 python scripts/core_temperature_stress_test.py",
        input_sha256={name: hashlib.sha256((root/path).read_bytes()).hexdigest()
                      for name, path in names.items()},
        inputs=names, caps_k=caps, n_cores=len(physical),
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
            for i in range(0,len(rows),3)),
        coherent_fragment_order_all_grid_points=all(
            rows[i]["core_fragment_percent"] > rows[i+1]["core_fragment_percent"] > rows[i+2]["core_fragment_percent"]
            for i in range(0,len(rows),3)),
        independent_cloud_order_all_corners=all(c["core_cloud_Y_lt_I_lt_E"] for c in corners),
        independent_fragment_order_all_corners=all(c["core_fragment_Y_gt_I_gt_E"] for c in corners),
    )
    (out/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
