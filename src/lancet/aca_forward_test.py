"""Forward test of combined images against measured ACA visibilities.

Each candidate sky image is multiplied by the native 7 m primary beam at the
actual ACA pointing and sampled at the measured ACA uv coordinates. No image
is used to translate ACA data to another primary beam, and no flux factor is
fitted to the measurements. A combined image that is consistent with the ACA
data reproduces the measured short-baseline signal RMS.
"""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np

from .combination_tests import (
    matching_cells, require_products, row_values, subset_rows, _write_rows_csv,
)
from .full_sampling import FULL_SAMPLING_OUTPUTS, OVERLAP_OUTPUTS, _cell_statistics, require_stage_outputs
from .tp_beam_test import _predict, smooth_to_circular_beam


def source_path(paths, source):
    """Path named by a ``group.name`` source (e.g. ``products.jcomb_continuum``) in the task paths."""
    group, name = source.split(".")
    return Path(paths[group][name])


def sky_model(path, pixel_arcsec):
    """Area-binned MJy/sr sky from an MJy/sr or Jy/beam restored image.

    Images whose pixels already exceed `pixel_arcsec` keep their native grid.
    """
    from astropy.io import fits
    from astropy.nddata import block_reduce

    header = fits.getheader(path)
    unit = header.get("BUNIT", "").strip().lower()
    if unit == "mjy/sr":
        scale = 1.0
    elif unit == "jy/beam":
        scale = 1e-6/(np.pi/(4*np.log(2))*header["BMAJ"]*header["BMIN"]*(np.pi/180)**2)
    else:
        raise ValueError(f"Unsupported sky-model unit {unit!r} in {path}")
    data = np.where(np.isfinite(d := np.squeeze(fits.getdata(path, memmap=True))), d, 0)*scale
    factor = int(round(pixel_arcsec/(abs(header["CDELT2"])*3600)))
    if factor > 1:
        data = block_reduce(data, factor, func=np.mean)
        for i in (1, 2):
            header[f"CRPIX{i}"] = (header[f"CRPIX{i}"]-0.5)/factor+0.5
            header[f"CDELT{i}"] *= factor
            header[f"NAXIS{i}"] = data.shape[2-i]
    header["BUNIT"] = "MJy/sr"
    return np.asarray(data, dtype=float), header


def pair_weights(n_pairs, sample=None):
    """Number of times each pointing pair is drawn (all ones without a sample)."""
    if sample is None:
        return np.ones(n_pairs)
    return np.bincount(sample, minlength=n_pairs).astype(float)


def binned_rms(power, use, cell_pair, weights):
    """RMS from the pair-weighted mean cell power; NaN when it is not positive."""
    w = weights[cell_pair[use]]
    if w.sum() == 0:
        return np.nan
    mean = np.sum(w*power[use])/w.sum()
    return float(np.sqrt(mean)) if mean > 0 else np.nan


def _support(rows, meta, pairs, image, header, settings):
    """PB-weighted fraction of each native ACA beam covered by nonzero model pixels."""
    mask = (np.isfinite(image) & (image != 0)).astype(float)
    result = {}
    for pair in pairs:
        one = subset_rows(rows, rows["field"] == pair["aca_field"])
        zero = {k: v[:1].copy() for k, v in one.items()}
        zero["uvw"] = zero["uvw"]*0
        value = _predict(zero, meta, [pair], mask, header, settings, target_pb=False)[0].real
        frequency = float(zero["frequency"][0])
        fwhm = settings["pb_7m_fwhm_arcsec"]*settings["pb_reference_frequency_ghz"]*1e9/frequency
        result[pair["pair"]] = float(value/(1e6*np.pi/(4*np.log(2))*(fwhm*np.pi/648000)**2))
    return result


def run_aca_forward_test(paths, config, *, save_figure, overwrite=False):
    """Compare sky-model predictions with measured ACA RMS at native ACA PBs."""
    from astropy.io import fits
    from matplotlib import pyplot as plt

    test = config["aca_forward_test"]
    sky_paths = {label: source_path(paths, spec["source"]) for label, spec in test["sky_models"].items()}
    require_products(paths, [spec["source"].split(".")[1] for spec in test["sky_models"].values()
                             if spec["source"].startswith("products.")])
    root = Path(paths["output_root"])
    require_stage_outputs(root, {"overlap": ("overlap_tests", OVERLAP_OUTPUTS),
                                 "full_sampling": ("full_sampling", ("summary.json",))})
    out = root/"aca_forward_test"
    if (out/"profiles.csv").exists() and not overwrite:
        raise FileExistsError("ACA forward test exists; use --overwrite")
    out.mkdir(parents=True, exist_ok=True)
    settings = config["overlap_tests"]
    sectors = settings["orientation_sectors"]
    coverage = json.loads((root/"full_sampling/summary.json").read_text())
    noise = json.loads((root/"overlap_tests/summary.json").read_text())["noise"]["aca_7m"]
    pairs = coverage["target_pointings"]
    n_pairs = len(pairs)
    edges = np.asarray(test["uv_edges_klambda"], float)
    n_bins = len(edges)-1

    tp_path = Path(paths["inputs"]["tp_dust_ff"])
    tp, tp_header = fits.getdata(tp_path, header=True)
    if tp_header["BUNIT"].strip().lower() != "jy/beam":
        raise ValueError("TP input must have Jy/beam units")
    tp_beam = float(np.sqrt(tp_header["BMAJ"]*tp_header["BMIN"])*3600)
    half = 2*np.log(2)/(np.pi*(tp_beam*np.pi/648000))/1000
    tp_limit = half*np.sqrt(np.log(1/config["figure2"]["tp_minimum_response"])/np.log(2))
    tp_area = np.pi/(4*np.log(2))*tp_header["BMAJ"]*tp_header["BMIN"]*(np.pi/180)**2
    tp_image = np.where(np.isfinite(tp), tp, 0)*1e-6/tp_area

    cache = root/"overlap_tests/aca_7m_rows.npz"
    meta_path = root/"overlap_tests/aca_7m_metadata.json"
    meta = json.loads(meta_path.read_text())
    with np.load(cache) as z:
        rows = {k: z[k] for k in z.files}
    _, _, radius = row_values(rows)
    use = np.isin(rows["field"], [p["aca_field"] for p in pairs])
    use &= np.isin(rows["ddid"], test["ddids"]) & (radius >= edges[0]) & (radius < edges[-1])
    rows = subset_rows(rows, use)
    key, good = matching_cells(rows, pairs, "aca", edges, sectors)
    size = n_pairs*2*n_bins*sectors
    kappa = np.array([noise[str(int(d))]["kappa_time"] for d in rows["ddid"]])

    predictions = {"tp": _predict(rows, meta, pairs, tp_image, tp_header, settings, target_pb=False)}
    sources = [cache, meta_path, tp_path, root/"full_sampling/summary.json"]
    support, beam_info = {}, {}
    for label, spec in test["sky_models"].items():
        path = sky_paths[label]
        print(f"ACA forward test: {label}", flush=True)
        sky, header = sky_model(path, settings["model_pixel_arcsec"])
        predictions[label] = _predict(rows, meta, pairs, sky, header, settings, target_pb=False)
        support[label] = _support(rows, meta, pairs, sky, header, settings)
        if spec.get("tp_beam_transfer", False):
            smoothed, beam_info[label] = smooth_to_circular_beam(sky, header, tp_beam)
            predictions[label+"_tp_beam"] = _predict(rows, meta, pairs, smoothed, header, settings,
                                                     target_pb=False)
        sources.append(path)
        del sky
    stats = _cell_statistics(rows, key, good, size, kappa, predictions)

    cell = np.arange(size)
    cell_pair = cell//(2*n_bins*sectors)
    cell_band = (cell//(n_bins*sectors)) % 2
    radial = (cell//sectors) % n_bins
    minimum = test["minimum_support"]
    interior = np.array([min(s[p["pair"]] for s in support.values()) >= minimum for p in pairs])
    rng = np.random.default_rng(test["random_seed"])
    draws = [rng.integers(0, n_pairs, n_pairs) for _ in range(test["bootstrap_repeats"])]
    models = [m for m in predictions if m != "tp" and not m.endswith("_tp_beam")]
    profiles = []
    for selection, pair_use in [("all_pairs", np.ones(n_pairs, bool)), ("interior_pairs", interior)]:
        for band_name, band_use in [("both", np.ones(size, bool)), ("lower", cell_band == 0),
                                    ("upper", cell_band == 1)]:
            for i in range(n_bins):
                use = (radial == i) & (stats["count"] >= 4) & band_use & pair_use[cell_pair]
                if not use.any():
                    continue
                weights = pair_weights(n_pairs)
                measured = binned_rms(stats["cross_power"], use, cell_pair, weights)
                record = {"selection": selection, "sideband": band_name,
                    "uv_low_klambda": float(edges[i]), "uv_high_klambda": float(edges[i+1]),
                    "uv_klambda": float(np.mean(stats["radius"][use])), "cells": int(use.sum()),
                    "rows": int(stats["count"][use].sum()), "pairs": int(len(set(cell_pair[use]))),
                    "aca_rms_mjy": 1000*measured}
                columns = models+(["tp"] if edges[i+1] <= tp_limit+1e-6 else [])
                columns += [m+"_tp_beam" for m in beam_info if edges[i+1] <= tp_limit+1e-6]
                for model in columns:
                    value = binned_rms(stats[model], use, cell_pair, weights)
                    record[model+"_rms_mjy"] = 1000*value
                    record[model+"_to_aca"] = value/measured
                    boot = []
                    for sample in draws:
                        w = pair_weights(n_pairs, sample)
                        boot.append(binned_rms(stats[model], use, cell_pair, w)
                                    / binned_rms(stats["cross_power"], use, cell_pair, w))
                    low, high = np.nanpercentile(boot, [16, 84])
                    record[model+"_to_aca_p16"], record[model+"_to_aca_p84"] = float(low), float(high)
                for model in beam_info:
                    if "tp_rms_mjy" in record:
                        record[f"tp_scale_given_{model}_transfer"] = record["tp_rms_mjy"]/record[model+"_tp_beam_rms_mjy"]
                profiles.append(record)
    _write_rows_csv(out/"profiles.csv", profiles)
    np.savez_compressed(out/"cell_statistics.npz", **stats)
    manifest = {"date": datetime.now(timezone.utc).date().isoformat(), "selected_rows": int(len(rows["x"])),
        "uv_edges_klambda": edges.tolist(), "tp_comparison_limit_klambda": float(tp_limit),
        "tp_beam_arcsec": tp_beam, "pb_7m_fwhm_arcsec_at_reference": settings["pb_7m_fwhm_arcsec"],
        "sky_models": test["sky_models"], "pb_weighted_support": support,
        "minimum_support": minimum, "interior_pairs": [p["pair"] for p, ok in zip(pairs, interior) if ok],
        "tp_transfer_models": beam_info,
        "method": ("Every sky image and the TP model are multiplied by the native 7 m Gaussian PB at the "
                   "ACA phase centre and Fourier-sampled at the measured ACA uv coordinates of the 19 "
                   "matched pointings. Measured RMS is the cross-hand signal power; models and data use "
                   "the same natural weights within cells and equal weights between cells. "
                   "Intervals are 16th-84th percentiles of pointing-pair bootstrap resamples. "
                   "No PB translation of data and no fitted normalization."),
        "limitations": ("Gaussian PBs; restored images include residual noise and deconvolution errors; "
                        "zero-filled outside image footprints (see pb_weighted_support); no spectral-index "
                        "scaling between image and sideband frequencies; pointings are not fully independent."),
        "sources_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    (out/"manifest.json").write_text(json.dumps(manifest, indent=2, default=float)+"\n")
    _plot(profiles, test, out, save_figure, plt)
    return manifest


def _plot(profiles, test, out, save_figure, plt):
    main = [r for r in profiles if r["selection"] == test["figure_selection"] and r["sideband"] == "both"]
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "pdf.fonttype": 42,
                         "svg.fonttype": "none", "legend.frameon": False})
    from matplotlib.ticker import NullFormatter, ScalarFormatter

    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2440944882, 3.6))
    fig.subplots_adjust(left=.1, right=.98, bottom=.3, top=.97, wspace=.3)
    x = np.array([r["uv_klambda"] for r in main])
    colors = plt.cm.viridis(np.linspace(0, .85, len(test["sky_models"])))
    for color, model in zip(colors, test["sky_models"]):
        y = np.array([r[model+"_to_aca"] for r in main])
        lo = np.array([r[model+"_to_aca_p16"] for r in main])
        hi = np.array([r[model+"_to_aca_p84"] for r in main])
        a.fill_between(x, lo, hi, color=color, alpha=.15, linewidth=0)
        a.plot(x, y, "o-", color=color, markersize=3, linewidth=1.2, label=test["sky_models"][model]["label"])
    a.axhline(1, color="0.4", linewidth=.7)
    a.set(xscale="log", yscale="log", xlabel=r"uv distance (k$\lambda$)",
          ylabel="Predicted / measured ACA RMS")
    a.set_xticks([6, 10, 20, 30])
    a.set_yticks([0.2, 0.5, 1, 2])
    for axis in (a.xaxis, a.yaxis):
        axis.set_major_formatter(ScalarFormatter())
        axis.set_minor_formatter(NullFormatter())
    handles, labels = a.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=7.5, bbox_to_anchor=(.5, 0))
    tp = [r for r in main if "tp_to_aca" in r]
    b.plot([r["uv_klambda"] for r in tp], [r["aca_rms_mjy"] for r in tp], "o-", color="k",
           markersize=3, linewidth=1.6, label="ACA measured", zorder=5)
    b.plot([r["uv_klambda"] for r in tp], [r["tp_rms_mjy"] for r in tp], "s--", color="0.45",
           markersize=3, label="TP model, native ACA PB")
    for color, model in zip(colors, test["sky_models"]):
        if test["sky_models"][model].get("signal_panel", False):
            b.plot([r["uv_klambda"] for r in tp], [r[model+"_rms_mjy"] for r in tp], "-",
                   color=color, linewidth=.9)
    b.set(yscale="log", xlabel=r"uv distance (k$\lambda$)", ylabel="Signal RMS (mJy)")
    b.yaxis.set_major_formatter(ScalarFormatter())
    b.yaxis.set_minor_formatter(NullFormatter())
    b.legend(fontsize=7.5, loc="lower left")
    for axis, label in zip((a, b), "ab"):
        axis.tick_params(direction="in", top=True, right=True, which="both")
        axis.text(.97, .95, f"({label})", transform=axis.transAxes, va="top", ha="right")
    save_figure(fig, [out/f"aca_forward_test.{ext}" for ext in test["output_formats"]],
                dpi=300, bbox_inches=None)
    plt.close(fig)
