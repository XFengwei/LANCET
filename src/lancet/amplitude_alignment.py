"""TP and interferometer amplitudes on a common primary beam and resolution.

The Figure 1 amplitude panel compares a TP model that keeps its 18.1-arcsec
restoring beam with interferometer visibilities that have no such beam, and
translates ACA to the 12 m PB with a sky model. This stage instead

* translates ACA to the matched 12 m PB with the 12 m+7 m image (no TP), which
  reproduces the measured ACA visibilities in the `aca-forward-test` stage;
* evaluates the TP model through the same PB at the same uv points; and
* removes the TP restoring beam with the effective transfer of that beam
  through the PB, computed from the same sky image, bracketed by plain
  division by the Gaussian beam response.

No flux factor is fitted or applied.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from .aca_forward_test import sky_model, source_path
from .combination_tests import (
    _broad_ddids, matching_cells, require_products, row_values, sample_image_visibilities, subset_rows,
)
from .full_sampling import _cell_statistics, separate_array_profiles
from .tp_beam_test import smooth_to_circular_beam
from .visibility_diagnostics import C_METRES_PER_SECOND

COLORS = {"configuration_c5": "#7254A3", "configuration_c2": "#007A9E", "aca_7m": "#D87821"}
SIDEBANDS_GHZ = ((0, 216.74), (1, 230.95))


def gaussian_transfer(u_klambda, fwhm_arcsec):
    """Fourier response of a unit-peak circular Gaussian beam."""
    theta = fwhm_arcsec*np.pi/648000
    return np.exp(-(np.pi*theta*np.asarray(u_klambda)*1000)**2/(4*np.log(2)))


def response_limit(fwhm_arcsec, minimum_response):
    """uv distance (kλ) where the Gaussian beam response falls to `minimum_response`."""
    theta = fwhm_arcsec*np.pi/648000
    return float(np.sqrt(np.log(1/minimum_response)*4*np.log(2))/(np.pi*theta)/1000)


def _pb_width(settings, key, frequency_hz):
    return float(settings[key])*float(settings["pb_reference_frequency_ghz"])*1e9/frequency_hz


def _native_and_common(rows, meta, targets, image, header, settings):
    """Model visibilities through the native 7 m PB and the matched 12 m PB."""
    kwargs = dict(pixel_arcsec=float(settings["model_pixel_arcsec"]),
                  patch_arcsec=float(settings["model_patch_arcsec"]),
                  padding=int(settings["model_padding_factor"]))
    native, common = np.full(len(rows["x"]), np.nan+0j), np.full(len(rows["x"]), np.nan+0j)
    for pair in targets:
        field = pair["aca_field"]
        ra, dec = np.rad2deg(np.asarray(meta["subtables"]["FIELD"][field]["PHASE_DIR"]).ravel())
        for band in (False, True):
            use = (rows["field"] == field) & ((rows["frequency"] > 225e9) == band)
            if not use.any():
                continue
            frequency = float(np.median(rows["frequency"][use]))
            uv = rows["uvw"][use]*rows["frequency"][use, None]/C_METRES_PER_SECOND
            native[use] = sample_image_visibilities(image, header, ra % 360, dec, uv[:, 0], uv[:, 1],
                fwhm_arcsec=_pb_width(settings, "pb_7m_fwhm_arcsec", frequency), **kwargs)
            common[use] = sample_image_visibilities(image, header, ra % 360, dec, uv[:, 0], uv[:, 1],
                fwhm_arcsec=_pb_width(settings, "pb_12m_fwhm_arcsec", frequency),
                pb_center=(pair["ra_deg"], pair["dec_deg"]), **kwargs)
    if not (np.all(np.isfinite(native)) and np.all(np.isfinite(common))):
        raise ValueError("Missing model predictions for selected ACA rows")
    return native, common


def continuous_rms(images, centres, keys, n_bins, sectors, radius, settings, pb_key):
    """RMS of each image through one PB type on a continuous radial grid.

    `centres` maps pair index to the (RA, Dec) of the PB centre: the matched
    12 m pointing for the 12 m PB, the ACA phase centre for the native 7 m PB.
    The orientation/sideband population is that of the given ACA cell keys,
    as in the Figure 1 TP reference.
    """
    kwargs = dict(pixel_arcsec=float(settings["model_pixel_arcsec"]),
                  patch_arcsec=float(settings["model_patch_arcsec"]),
                  padding=int(settings["model_padding_factor"]))
    power = {name: np.zeros_like(radius) for name in images}
    number = 0
    for pair, (ra, dec) in centres.items():
        for band, frequency_ghz in SIDEBANDS_GHZ:
            selected = keys[keys//(sectors*n_bins*2) == pair]
            selected = selected[(selected//(sectors*n_bins)) % 2 == band]
            if not len(selected):
                continue
            angle = ((selected % sectors)+0.5)*np.pi/sectors
            u = radius[:, None]*1000*np.cos(angle[None, :])
            v = radius[:, None]*1000*np.sin(angle[None, :])
            width = _pb_width(settings, pb_key, frequency_ghz*1e9)
            for name, (image, header) in images.items():
                value = sample_image_visibilities(image, header, ra, dec,
                    u.ravel(), v.ravel(), fwhm_arcsec=width, **kwargs)
                power[name] += np.sum(np.abs(value.reshape(u.shape))**2, axis=1)
            number += len(selected)
    return {name: 1000*np.sqrt(p/number) for name, p in power.items()}


def _bin_values(stats, use, weights, cell_pair, tp_beam):
    """TP/7 m comparison for one bin with pair weights (native and 12 m PB)."""
    w = weights[cell_pair[use]]
    mean = lambda key: float(np.sum(w*stats[key][use])/np.sum(w))
    power = mean("cross_power")
    uv = float(np.mean(stats["radius"][use]))
    gauss = float(gaussian_transfer(uv, tp_beam))
    values = {"uv_klambda": uv, "gaussian_beam_response": gauss}
    for pb, key in (("native", "native"), ("common", "target")):
        translation = np.sqrt(mean("sky_target")/mean("sky_native")) if pb == "common" else 1.0
        aca = 1000*np.sqrt(power)*translation if power > 0 else np.nan
        tp = 1000*np.sqrt(mean("tp_"+key))
        transfer = np.sqrt(mean("sky_beam_"+key)/mean("sky_"+key))
        values.update({f"aca_{pb}_mjy": aca, f"tp_raw_{pb}_mjy": tp,
            f"tp_beam_transfer_{pb}": transfer, f"tp_corrected_{pb}_mjy": tp/transfer,
            f"tp_gaussian_corrected_{pb}_mjy": tp/gauss,
            f"ratio_raw_{pb}": tp/aca, f"ratio_corrected_{pb}": tp/transfer/aca,
            f"ratio_gaussian_{pb}": tp/gauss/aca})
    return values


def overlap_records(stats, edges, sectors, limit, tp_beam, n_pairs, draws=()):
    """Per-bin TP/7 m comparison on identical ACA cells, with pair bootstrap."""
    size = len(stats["count"])
    n_bins = len(edges)-1
    cell = np.arange(size)
    radial = (cell//sectors) % n_bins
    cell_pair = cell//(2*n_bins*sectors)
    records = []
    for i in range(n_bins):
        use = (radial == i) & (stats["count"] >= 4)
        if not use.any() or edges[i+1] > limit+1e-6:
            continue
        record = {"uv_low_klambda": float(edges[i]), "uv_high_klambda": float(edges[i+1]),
                  "cells": int(use.sum()), "rows": int(stats["count"][use].sum()),
                  "pairs": int(len(set(cell_pair[use])))}
        record.update(_bin_values(stats, use, np.ones(n_pairs), cell_pair, tp_beam))
        if len(draws):
            boot = [_bin_values(stats, use, np.bincount(d, minlength=n_pairs).astype(float),
                                cell_pair, tp_beam) for d in draws]
            for key in [k for k in record if k.startswith("ratio_") or k.startswith("aca_")]:
                values = np.array([b[key] for b in boot], float)
                record[key+"_p16"], record[key+"_p84"] = (float(v) for v in np.nanpercentile(values, [16, 84]))
        records.append(record)
    return records


def _write_csv(path, records):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def run_amplitude_alignment(paths, config, *, save_figure, overwrite=False):
    from astropy.io import fits

    require_products(paths, ["jcomb_continuum", "feather_continuum"])
    root = Path(paths["output_root"])
    out = root/"amplitude_alignment"
    if (out/"overlap.csv").exists() and not overwrite:
        raise FileExistsError("Amplitude alignment outputs exist; use --overwrite")
    out.mkdir(parents=True, exist_ok=True)
    task = config["amplitude_alignment"]
    settings = config["overlap_tests"]
    sectors = int(settings["orientation_sectors"])
    sampling = root/"full_sampling"
    coverage = json.loads((sampling/"summary.json").read_text())
    if coverage["settings"] != settings:
        raise ValueError("Full-sampling settings changed; rerun the full-sampling stage")
    noise = json.loads((root/"overlap_tests/summary.json").read_text())["noise"]
    targets = coverage["target_pointings"]
    edges = np.asarray(coverage["uv_edges_klambda"], float)
    n_bins = len(edges)-1
    size = len(targets)*2*n_bins*sectors

    tp_path = Path(paths["inputs"]["tp_dust_ff"])
    tp, tp_header = fits.getdata(tp_path, header=True)
    if tp_header["BUNIT"].strip().lower() != "jy/beam":
        raise ValueError("TP input must have Jy/beam units")
    tp_beam = float(np.sqrt(tp_header["BMAJ"]*tp_header["BMIN"])*3600)
    tp_area = np.pi/(4*np.log(2))*tp_header["BMAJ"]*tp_header["BMIN"]*(np.pi/180)**2
    tp_image = np.where(np.isfinite(tp), tp, 0)*1e-6/tp_area
    limit = response_limit(tp_beam, float(task["tp_minimum_response"]))

    sky_path = source_path(paths, task["sky_model"])
    print("Amplitude alignment: loading sky model", flush=True)
    sky, sky_header = sky_model(sky_path, settings["model_pixel_arcsec"])
    sky_beam, beam_info = smooth_to_circular_beam(sky, sky_header, tp_beam)

    cache = root/"overlap_tests/aca_7m_rows.npz"
    meta_path = root/"overlap_tests/aca_7m_metadata.json"
    meta = json.loads(meta_path.read_text())
    with np.load(cache) as z:
        rows = {k: z[k] for k in z.files}
    use = np.isin(rows["ddid"], _broad_ddids(meta)) & np.isin(rows["field"], [p["aca_field"] for p in targets])
    _, _, radius = row_values(rows)
    rows = subset_rows(rows, use & (radius >= edges[0]) & (radius < edges[-1]))
    key, good = matching_cells(rows, targets, "aca", edges, sectors)
    kappa = np.array([noise["aca_7m"][str(int(d))]["kappa_time"] for d in rows["ddid"]])
    models = {}
    for label, (image, header) in {"sky": (sky, sky_header), "sky_beam": (sky_beam, sky_header),
                                   "tp": (tp_image, tp_header)}.items():
        print(f"Amplitude alignment: ACA predictions for {label}", flush=True)
        models[label+"_native"], models[label+"_target"] = _native_and_common(
            rows, meta, targets, image, header, settings)
    aca = _cell_statistics(rows, key, good, size, kappa, models)

    statistics = {"aca_7m": aca}
    for array in ("configuration_c2", "configuration_c5"):
        with np.load(sampling/f"{array}_rows.npz") as z:
            data = {k: z[k] for k in z.files}
        k, g = matching_cells(data, targets, "alma", edges, sectors)
        kap = np.array([noise[array][str(int(d))]["kappa_time"] for d in data["ddid"]])
        statistics[array] = _cell_statistics(data, k, g, size, kap)
    profiles = separate_array_profiles(statistics, n_bins, sectors, model="sky")
    _write_csv(out/"array_profiles.csv", profiles)
    fine = np.asarray(task["overlap_uv_edges_klambda"], float)
    fine_key, fine_good = matching_cells(rows, targets, "aca", fine, sectors)
    fine_size = len(targets)*2*(len(fine)-1)*sectors
    fine_stats = _cell_statistics(rows, fine_key, fine_good, fine_size, kappa, models)
    rng = np.random.default_rng(int(task["random_seed"]))
    draws = [rng.integers(0, len(targets), len(targets)) for _ in range(int(task["bootstrap_repeats"]))]
    overlap = overlap_records(fine_stats, fine, sectors, limit, tp_beam, len(targets), draws)
    _write_csv(out/"overlap.csv", overlap)

    keys = np.flatnonzero(aca["count"] >= 4)
    first = int(np.min((keys//sectors) % n_bins))
    keys = keys[(keys//sectors) % n_bins == first]
    grid = np.geomspace(float(task["continuous_min_klambda"]), limit, int(task["continuous_points"]))
    print("Amplitude alignment: continuous TP curves", flush=True)
    images = {"tp": (tp_image, tp_header), "sky": (sky, sky_header), "sky_beam": (sky_beam, sky_header)}
    aca_centres = {}
    for p in targets:
        ra, dec = np.rad2deg(np.asarray(meta["subtables"]["FIELD"][p["aca_field"]]["PHASE_DIR"]).ravel())
        aca_centres[p["pair"]] = (float(ra % 360), float(dec))
    frames = {"common": continuous_rms(images, {p["pair"]: (p["ra_deg"], p["dec_deg"]) for p in targets},
                                       keys, n_bins, sectors, grid, settings, "pb_12m_fwhm_arcsec"),
              "native": continuous_rms(images, aca_centres, keys, n_bins, sectors, grid, settings,
                                       "pb_7m_fwhm_arcsec")}
    gauss = gaussian_transfer(grid, tp_beam)
    continuous = {"uv_klambda": grid, "gaussian_beam_response": gauss}
    for frame, curves in frames.items():
        transfer = curves["sky_beam"]/curves["sky"]
        continuous.update({f"tp_raw_{frame}_mjy": curves["tp"], f"sky_{frame}_mjy": curves["sky"],
            f"tp_beam_transfer_{frame}": transfer, f"tp_corrected_{frame}_mjy": curves["tp"]/transfer,
            f"tp_gaussian_corrected_{frame}_mjy": curves["tp"]/gauss})
    _write_csv(out/"tp_continuous.csv", [dict(zip(continuous, values)) for values in zip(*continuous.values())])

    sources = [cache, meta_path, tp_path, sky_path, sampling/"summary.json"]
    manifest = {"date": "2026-11-13", "tp_beam_arcsec": tp_beam, "tp_response_limit_klambda": limit,
        "sky_model": str(sky_path), "tp_beam_kernel": beam_info, "reference_aca_bin": first,
        "selected_aca_rows": int(len(rows["x"])),
        "method": ("ACA signal RMS is translated from the native 7 m PB to the matched 12 m PB with the "
                   "12m+7m image (no TP). The TP model is Fourier-sampled through the same PB at the same "
                   "ACA uv points. The TP restoring beam is removed by dividing by its effective transfer "
                   "through the PB, i.e. the RMS ratio of the sky image smoothed to the TP beam and unsmoothed, "
                   "on identical cells; the Gaussian beam response gives the alternative bracket. "
                   "No fitted or applied flux factor."),
        "limitations": ("The 12m+7m image lacks scales beyond the ACA's shortest spacing, so its effective "
                        "transfer approaches the Gaussian one; a sky with more large-scale emission leaks more "
                        "power through the PB and implies a smaller correction. Gaussian PBs; 19 matched pointing "
                        "pairs only; no spectral-index scaling. The correction is meaningful only while the "
                        "Gaussian beam response exceeds the configured minimum."),
        "sources_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    (out/"manifest.json").write_text(json.dumps(manifest, indent=2, default=float)+"\n")
    plot_alignment(profiles, overlap, continuous, task, limit, out, save_figure)
    return manifest


def plot_alignment(profiles, overlap, continuous, task, limit, out, save_figure):
    from matplotlib import pyplot as plt
    from matplotlib.ticker import NullFormatter, ScalarFormatter

    plt.rcParams.update({"font.family": "sans-serif", "font.size": 8, "axes.linewidth": 0.7,
                         "pdf.fonttype": 42, "svg.fonttype": "none", "legend.frameon": False})
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(task["width_inches"], task["height_inches"]))
    fig.subplots_adjust(left=.075, right=.99, bottom=.17, top=.97, wspace=.36)
    for axis in (a, b, c):
        axis.tick_params(direction="in", top=True, right=True, which="both", length=3, width=.7)
        axis.tick_params(which="minor", length=1.5)
    plain = lambda axis: [(f.set_major_formatter(ScalarFormatter()), f.set_minor_formatter(NullFormatter()))
                          for f in (axis.xaxis, axis.yaxis)]
    u = continuous["uv_klambda"]
    lo_uv, hi_uv = overlap[0]["uv_low_klambda"], limit
    raw_style = dict(color="0.55", linestyle="-.", linewidth=1)
    fixed_style = dict(color="0.1", linewidth=1.3)

    # (a) full range on the matched 12 m PB
    a.axvspan(lo_uv, hi_uv, color="0.9", linewidth=0, zorder=0)
    labels = {"aca_7m": "7 m (12 m PB)", "configuration_c2": "12 m C2", "configuration_c5": "12 m C5"}
    for array, marker in (("aca_7m", "o"), ("configuration_c2", "s"), ("configuration_c5", "^")):
        p = [r for r in profiles if r["array"] == array]
        x = np.array([r["uv_klambda"] for r in p])
        y = np.array([r["rms_mjy"] if r["detected"] else np.nan for r in p])
        a.plot(x, y, color=COLORS[array], marker=marker, markersize=2.2, linewidth=.9, label=labels[array])
        upper = np.array([(r["uv_klambda"], r["upper_rms_mjy"]) for r in p if not r["detected"]])
        if len(upper):
            a.errorbar(upper[:, 0], upper[:, 1], yerr=upper[:, 1]*.25, uplims=True, fmt="none",
                       color=COLORS[array], elinewidth=.7, capsize=1.3)
    a.plot(u, continuous["tp_raw_common_mjy"], label="TP", **raw_style)
    a.plot(u, continuous["tp_corrected_common_mjy"], label="TP, beam-corrected", **fixed_style)
    a.set(xscale="log", yscale="log", xlim=(1, 1000), ylim=(3, 1e5), ylabel="Amplitude (mJy)",
          xlabel=r"uv distance (k$\lambda$)")
    a.set_xticks([1, 10, 100, 1000], labels=["1", "10", "100", "1000"])
    a.legend(loc="upper right", fontsize=7, handlelength=1.6, labelspacing=.2, borderaxespad=.3)
    a.text(.04, .05, "(a) 12 m PB", transform=a.transAxes)

    # (b) overlap zoom on the native 7 m PB, no translation of the data
    zoom = (u >= lo_uv*.95) & (u <= hi_uv)
    b.fill_between(u[zoom], continuous["tp_corrected_native_mjy"][zoom],
                   continuous["tp_gaussian_corrected_native_mjy"][zoom], color="0.2", alpha=.15, linewidth=0)
    b.plot(u[zoom], continuous["tp_raw_native_mjy"][zoom], label="TP", **raw_style)
    b.plot(u[zoom], continuous["tp_corrected_native_mjy"][zoom], label="TP, beam-corrected", **fixed_style)
    x = np.array([r["uv_klambda"] for r in overlap])
    y = np.array([r["aca_native_mjy"] for r in overlap])
    err = np.array([[r["aca_native_mjy"]-r["aca_native_mjy_p16"], r["aca_native_mjy_p84"]-r["aca_native_mjy"]]
                    for r in overlap]).T
    b.errorbar(x, y, yerr=err, fmt="o", color=COLORS["aca_7m"], markersize=3, elinewidth=.8, capsize=1.5,
               label="7 m measured")
    b.set(xscale="log", yscale="log", xlim=(lo_uv*.97, hi_uv*1.03), ylim=(40, 1500),
          xlabel=r"uv distance (k$\lambda$)", ylabel="Amplitude (mJy)")
    b.set_xticks([5, 6, 7, 8, 9]); b.set_yticks([100, 300, 1000])
    plain(b)
    b.legend(loc="lower left", fontsize=7, handlelength=1.6, labelspacing=.2, borderaxespad=.3)
    b.text(.96, .95, "(b) 7 m PB", transform=b.transAxes, ha="right", va="top")

    # (c) ratio on the native 7 m PB with pointing bootstrap intervals
    for key, style, label in (("ratio_raw_native", raw_style, "TP"),
                              ("ratio_corrected_native", fixed_style, "TP, beam-corrected")):
        y = np.array([r[key] for r in overlap])
        c.fill_between(x, [r[key+"_p16"] for r in overlap], [r[key+"_p84"] for r in overlap],
                       color=style["color"], alpha=.18, linewidth=0)
        c.plot(x, y, marker="o", markersize=2.5, label=label, **style)
    c.plot(x, [r["ratio_gaussian_native"] for r in overlap], color="0.1", linestyle=":", linewidth=1,
           label="TP, Gaussian beam-corrected")
    c.axhline(1, color="k", linewidth=.6)
    c.axhline(1/2.7, color="#B2182B", linewidth=.8, linestyle="--")
    c.text(x.max(), 1/2.7*1.08, "1/2.7", color="#B2182B", ha="right", va="bottom", fontsize=7)
    c.set(xscale="log", yscale="log", xlim=(lo_uv*.97, hi_uv*1.03), ylim=(.1, 3),
          xlabel=r"uv distance (k$\lambda$)", ylabel="TP / 7 m")
    c.set_xticks([5, 6, 7, 8, 9]); c.set_yticks([.1, .3, 1, 3])
    plain(c)
    c.legend(loc="lower left", fontsize=7, handlelength=1.8, labelspacing=.2, borderaxespad=.3)
    c.text(.96, .95, "(c) 7 m PB", transform=c.transAxes, ha="right", va="top")
    save_figure(fig, [out/f"amplitude_alignment.{ext}" for ext in task["output_formats"]],
                dpi=int(task["dpi"]), bbox_inches=None)
    plt.close(fig)


def read_products(out):
    """Saved array profiles, overlap bins, continuous TP curves and manifest."""
    out = Path(out)

    def convert(value):
        if value in ("True", "False"):
            return value == "True"
        try:
            return float(value)
        except ValueError:
            return value

    def read(name):
        with (out/name).open() as handle:
            return [{k: convert(v) for k, v in r.items()} for r in csv.DictReader(handle)]

    continuous = read("tp_continuous.csv")
    continuous = {k: np.array([r[k] for r in continuous]) for k in continuous[0]}
    manifest = json.loads((out/"manifest.json").read_text())
    return read("array_profiles.csv"), read("overlap.csv"), continuous, manifest


def replot_alignment(paths, config, *, save_figure):
    """Redraw the figure from the saved CSV products without recomputing."""
    out = Path(paths["output_root"])/"amplitude_alignment"
    profiles, overlap, continuous, manifest = read_products(out)
    plot_alignment(profiles, overlap, continuous, config["amplitude_alignment"],
                   manifest["tp_response_limit_klambda"], out, save_figure)
