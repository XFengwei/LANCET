"""Sample explicit temperature priors and propagate regional core mass sums.

This is prior-predictive Monte Carlo, with no fitted temperature likelihood.
The adopted catalogues are read-only; foreground/background stays at Tenv.
"""
from pathlib import Path
import argparse
import hashlib
import json
import tomllib

from astropy import units as u
from astropy.table import Table
import numpy as np
from scipy.special import betaincinv
from scipy.stats import beta, rankdata
from numpy.polynomial.legendre import leggauss

from lancet.physical_properties import blackbody_cgs, core_optical_depth


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("tasks/physical_properties/task.toml"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    cfg = tomllib.loads((root/args.config).read_text())
    settings, physics = cfg["temperature_sampling"], cfg["physics"]
    out = root/settings["output_root"]
    out.mkdir(parents=True, exist_ok=True)
    inputs = dict(source=cfg["inputs"]["source_catalogue"],
        physical=cfg["outputs"]["root"]+"/tables/core_physical_properties.ecsv",
        regional="results/core_formation/tables/core_formation_summary.ecsv")
    source, physical, regional = [Table.read(root/inputs[k]) for k in ("source", "physical", "regional")]
    order = {int(row["ID"]): i for i, row in enumerate(source)}
    assert len(order) == len(source) == len(physical)
    source = source[[order[int(i)] for i in physical["ID"]]]
    assert np.array_equal(source["ID"], physical["ID"])
    tenv = np.asarray(physical["Tenv"], dtype=float)
    codes = np.asarray(physical["Subregion"])
    regcodes = [str(r["subregion_code"]) for r in regional]
    capdict = dict(zip(["Y", "I", "E"], settings["upper_temperatures_k"]))
    caps = np.array([capdict[c] for c in codes])
    assert np.all(caps > tenv) and np.allclose(tenv, source["Tenv"])
    mass0 = np.asarray(physical["Mass_core"], dtype=float)
    masks = [codes == c for c in regcodes]
    for reg, mask in zip(regional, masks):
        assert np.isclose(mass0[mask].sum(), reg["core_mass"], rtol=1e-12)
    nu = physics["frequency_ghz"] * u.GHz
    beam = 1.1331*physics["beam_major_arcsec"]*physics["beam_minor_arcsec"]*u.arcsec**2
    peak_jy = np.asarray(source["FXP_BST01"], dtype=float)
    assert np.all(np.isfinite(peak_jy)) and np.all(peak_jy > 0)
    peak = (peak_jy*u.Jy/beam).to(u.MJy/u.sr)
    background = np.asarray(source["Ibkg2"], dtype=float)*u.MJy/u.sr
    benv = blackbody_cgs(nu, tenv*u.K)

    def evaluate(fraction):
        temp = (tenv + fraction*(caps-tenv))*u.K
        bc = blackbody_cgs(nu, temp)
        tau = core_optical_depth(nu=nu, temp=temp, environment_temp=tenv*u.K,
                                Isr=peak, Ibg=background).value
        correction = tau / -np.expm1(-tau)
        return bc, correction

    _, f0 = evaluate(np.zeros(len(tenv)))

    def masses(fraction):
        bc, correction = evaluate(fraction)
        return mass0*(benv/bc).value*correction/f0

    np.testing.assert_allclose(masses(np.zeros(len(tenv))), mass0, rtol=1e-12)
    minimum_mass = masses(np.ones(len(tenv)))
    # Global peak-brightness percentile; logarithms would preserve this rank.
    percentile = (rankdata(peak_jy, method="average")-1)/(len(tenv)-1)
    strength = float(settings["brightness_strength"])
    priors = {"uniform": (np.ones(len(tenv)), np.ones(len(tenv))),
              "brightness": (1+strength*percentile, 1+strength*(1-percentile))}
    draws = int(settings["sample_count"])
    seeds = np.random.SeedSequence(int(settings["random_seed"])).spawn(4)
    summary_rows, core_rows, order_rows, checks, samples = [], [], [], [], {}
    # Independent quadrature checks the sampled mean without relying on RNG.
    nodes, weights = leggauss(160)
    fractions = (nodes[:, None]+1)/2
    quadrature_mass = masses(fractions)
    seed_index = 0
    for prior, (alpha, beta_shape) in priors.items():
        expected = np.sum(quadrature_mass * beta.pdf(fractions, alpha, beta_shape)
                          * weights[:, None]/2, axis=0)
        for i in range(len(tenv)):
            f_quantiles = betaincinv(alpha[i], beta_shape[i], np.array([.16, .5, .84]))
            # Mass decreases with temperature: reverse temperature quantiles.
            single = np.zeros((3, len(tenv)))
            single[:, i] = f_quantiles[::-1]
            mq = masses(single)[:, i]
            core_rows.append(dict(prior=prior, ID=int(physical["ID"][i]), code=str(codes[i]),
                peak_jy_beam=peak_jy[i], brightness_percentile=percentile[i],
                alpha=alpha[i], beta=beta_shape[i], heating_mean=alpha[i]/(alpha[i]+beta_shape[i]),
                Tenv_k=tenv[i], Tupper_k=caps[i],
                T_p16_k=tenv[i]+f_quantiles[0]*(caps[i]-tenv[i]),
                T_p50_k=tenv[i]+f_quantiles[1]*(caps[i]-tenv[i]),
                T_p84_k=tenv[i]+f_quantiles[2]*(caps[i]-tenv[i]),
                mass_p16_msun=mq[0], mass_p50_msun=mq[1], mass_p84_msun=mq[2]))
        for coupling in ("independent", "regional_shared"):
            rng = np.random.default_rng(seeds[seed_index])
            seed_index += 1
            totals = np.empty((draws, len(regcodes)))
            sum_f = np.zeros(len(tenv))
            for start in range(0, draws, 2500):
                stop = min(start+2500, draws)
                n = stop-start
                if coupling == "independent":
                    f = rng.beta(alpha, beta_shape, size=(n, len(tenv)))
                else:
                    # Same percentile within each region; preserve each core's
                    # marginal prior but test maximally correlated heating.
                    region_u = rng.random((n, len(regcodes)))
                    u_core = np.empty((n, len(tenv)))
                    for j, mask in enumerate(masks):
                        u_core[:, mask] = region_u[:, j, None]
                    f = betaincinv(alpha, beta_shape, u_core)
                mass = masses(f)
                assert np.all(np.isfinite(mass))
                assert np.all(mass >= minimum_mass*(1-1e-12))
                assert np.all(mass <= mass0*(1+1e-12))
                sum_f += f.sum(axis=0)
                for j, mask in enumerate(masks):
                    totals[start:stop, j] = mass[:, mask].sum(axis=1)
            key = f"{prior}_{coupling}"
            samples[key] = totals
            checks.append(dict(model=key,
                max_core_heating_mean_residual=float(np.max(np.abs(sum_f/draws-alpha/(alpha+beta_shape))))))
            for j, (reg, mask) in enumerate(zip(regional, masks)):
                vals = totals[:, j]
                analytic_mean = expected[mask].sum()
                mcse = vals.std(ddof=1)/np.sqrt(draws)
                assert abs(vals.mean()-analytic_mean) < 6*mcse
                quantiles = np.quantile(vals, [.025, .16, .5, .84, .975])
                row = dict(prior=prior, coupling=coupling, code=regcodes[j], region=str(reg["label"]),
                    n_cores=int(mask.sum()), baseline_mass_msun=float(reg["core_mass"]),
                    cloud_mass_msun=float(reg["cloud_mass"]), fragment_mass_msun=float(reg["fragment_mass"]),
                    mean_mass_msun=float(vals.mean()), mean_mcse_msun=float(mcse),
                    quadrature_mean_msun=float(analytic_mean),
                    sd_mass_msun=float(vals.std(ddof=1)))
                for tag, value in zip(["p025", "p16", "p50", "p84", "p975"], quantiles):
                    row[f"mass_{tag}_msun"] = float(value)
                    row[f"core_cloud_{tag}_percent"] = float(100*value/reg["cloud_mass"])
                    row[f"core_fragment_{tag}_percent"] = float(100*value/reg["fragment_mass"])
                first, second = np.quantile(vals[:draws//2], [.16,.5,.84]), np.quantile(vals[draws//2:], [.16,.5,.84])
                row["split_quantile_difference_over_68_width"] = float(np.max(abs(first-second))/(quantiles[3]-quantiles[1]))
                summary_rows.append(row)
            cc = totals / np.array(regional["cloud_mass"])
            cf = totals / np.array(regional["fragment_mass"])
            order_rows.append(dict(prior=prior, coupling=coupling,
                core_cloud_Y_lt_I_lt_E=float(np.mean((cc[:,0]<cc[:,1]) & (cc[:,1]<cc[:,2]))),
                core_fragment_Y_gt_I_gt_E=float(np.mean((cf[:,0]>cf[:,1]) & (cf[:,1]>cf[:,2]))),
                note="Draw fractions conditional on assumed priors, not observational confidence"))
            print(f"Completed {key}: {draws} independent joint realizations", flush=True)
    Table(rows=summary_rows).write(out/"regional_sampling_summary.ecsv", overwrite=True)
    Table(rows=core_rows).write(out/"core_temperature_priors.ecsv", overwrite=True)
    Table(rows=order_rows).write(out/"ordering_draw_fractions.ecsv", overwrite=True)
    np.savez_compressed(out/"regional_mass_draws.npz", **samples)
    metadata = dict(settings=settings, inputs=inputs,
        input_sha256={k: hashlib.sha256((root/p).read_bytes()).hexdigest() for k,p in inputs.items()},
        config_sha256=hashlib.sha256((root/args.config).read_bytes()).hexdigest(),
        method="Direct independent Monte Carlo realizations from explicit assumed priors; no temperature likelihood or MCMC posterior",
        temperature_formula="T_i=Tenv_i+F_i*(Tupper_region-Tenv_i)",
        uniform="F_i~Uniform(0,1)",
        brightness="q_i=global peak-brightness percentile; F_i~Beta(1+2*q_i,3-2*q_i); mean=0.25+0.5*q_i",
        coupling="Independent core draws primary; regional_shared uses one uniform percentile per region per draw and each core's inverse CDF",
        fixed="Measured flux, foreground/background at Tenv, opacity, distance, geometry, fragment/cloud masses; no free-free correction",
        intervals="Equal-tail 16-84% and 2.5-97.5% conditional temperature-prior ranges; not total systematic errors or data-constrained confidence intervals",
        region_column_order=regcodes, verification=checks)
    (out/"metadata.json").write_text(json.dumps(metadata, indent=2)+"\n")
    print(Table(rows=summary_rows)[["prior","coupling","region","mass_p16_msun","mass_p50_msun","mass_p84_msun"]])


if __name__ == "__main__":
    main()
