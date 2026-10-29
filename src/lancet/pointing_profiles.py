"""Pointing-resolved visibility RMS profiles, including sampled TP models."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from .combination_tests import (
    binned_sky_model, matching_cells, row_values, subset_rows, _write_rows_csv,
    sample_image_visibilities,
)
from .full_sampling import _cell_statistics, separate_array_profiles
from .tp_beam_test import _predict


COLORS = {"aca_7m": "#D87821", "configuration_c2": "#007A9E",
          "configuration_c5": "#7254A3", "tp": "#555555"}
LABELS = {"aca_7m": "7 m (12 m PB)", "configuration_c2": "C2",
          "configuration_c5": "C5", "tp": "TP model"}
MARKERS = {"aca_7m": "o", "configuration_c2": "s", "configuration_c5": "^", "tp": "D"}


def per_pointing_profiles(statistics, pairs, edges, sectors):
    """Apply the maintained cell/RMS estimator separately to every pointing."""
    nb = len(edges)-1
    stride = 2*nb*sectors
    records = []
    for pair in pairs:
        p = pair["pair"]
        selected = {name: {k: v[p*stride:(p+1)*stride] for k, v in data.items()}
                    for name, data in statistics.items()}
        for record in separate_array_profiles(selected, nb, sectors):
            i = record["bin"]
            records.append({"pair": p, "aca_field": pair["aca_field"], "alma_field": pair["alma_field"],
                "ra_deg": pair["ra_deg"], "dec_deg": pair["dec_deg"],
                "uv_low_klambda": float(edges[i]), "uv_high_klambda": float(edges[i+1]), **record})
    return records


def _draw_panel(axis, records, pair, ylim):
    for name, color in COLORS.items():
        data = [r for r in records if r["pair"] == pair["pair"] and r["array"] == name]
        radius = np.array([r["uv_klambda"] for r in data])
        amplitude = np.array([r["rms_mjy"] for r in data])
        detected = np.array([r["detected"] for r in data], dtype=bool)
        upper = np.array([r["upper_rms_mjy"] for r in data])
        axis.plot(radius, np.where(detected, amplitude, np.nan), color=color,
                  marker=MARKERS[name], markersize=2.8, linewidth=1.1,
                  linestyle="-." if name == "tp" else "-", label=LABELS[name])
        if (~detected).any():
            axis.errorbar(radius[~detected], upper[~detected], yerr=upper[~detected]*.25,
                uplims=True, fmt="none", color=color, elinewidth=.7, capsize=1.8, markersize=2.5)
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set(xlim=(1, 1000), ylim=ylim)
    axis.set_xticks([1, 10, 100, 1000], labels=["1", "10", "100", "1000"])
    axis.tick_params(direction="in", top=True, right=True, labelsize=10, length=3)
    axis.text(.04, .96, f"P{pair['pair']+1:02d}", transform=axis.transAxes,
              va="top", fontsize=11, fontweight="bold")
    axis.text(.22, .85, f"7m{pair['aca_field']}\n12m{pair['alma_field']}",
              transform=axis.transAxes, ha="right", va="top", fontsize=10, linespacing=1.1)
    axis.legend(loc="upper right", fontsize=10, frameon=True, framealpha=1,
                facecolor="white", edgecolor="none", borderpad=.25,
                handlelength=1.2, handletextpad=.35, labelspacing=.15)


def run_pointing_profiles(paths, config, *, save_figure, overwrite=False):
    from astropy.io import fits
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    import pymupdf

    root = Path(paths["output_root"])
    out = root/"pointing_profiles"
    out.mkdir(parents=True, exist_ok=True)
    if (out/"profiles.csv").exists() and not overwrite:
        raise FileExistsError("Per-pointing profiles exist; use --overwrite")
    coverage_path = root/"full_sampling/summary.json"
    noise_path = root/"overlap_tests/summary.json"
    coverage = json.loads(coverage_path.read_text())
    noise = json.loads(noise_path.read_text())["noise"]
    settings = config["overlap_tests"]
    pairs = coverage["target_pointings"]
    tp_path = Path(paths["reduction_inputs"]["total_power"])
    tp, header = fits.getdata(tp_path, header=True)
    if header["BUNIT"].strip().lower() != "jy/beam":
        raise ValueError("TP input must have Jy/beam units")
    theta = np.sqrt(header["BMAJ"]*header["BMIN"])*np.pi/180
    maximum = np.sqrt(4*np.log(2)*np.log(1/config["figure1"]["tp_minimum_response"]))/np.pi/theta/1000
    edges = np.sort(np.r_[config["figure1"]["amplitude_uv_edges_klambda"], maximum])
    sectors = int(settings["orientation_sectors"])
    size = len(pairs)*2*(len(edges)-1)*sectors
    sources = [coverage_path, noise_path, tp_path]
    statistics, counts = {}, {}
    for name in COLORS:
        if name == "tp":
            continue
        directory = root/("overlap_tests" if name == "aca_7m" else "full_sampling")
        path, meta_path = directory/f"{name}_rows.npz", directory/f"{name}_metadata.json"
        sources.extend([path, meta_path])
        with np.load(path) as z:
            rows = {k: z[k] for k in z.files}
        if name == "aca_7m":
            rows = subset_rows(rows, np.isin(rows["field"], [p["aca_field"] for p in pairs])
                               & np.isin(rows["ddid"], [0, 7]))
        counts[name] = len(rows["x"])
        key, good = matching_cells(rows, pairs, "aca" if name == "aca_7m" else "alma", edges, sectors)
        if not np.all(good):
            raise ValueError("Per-pointing bins omit cached rows")
        kappa = np.array([noise[name][str(int(d))]["kappa_time"] for d in rows["ddid"]])
        models = {}
        if name == "aca_7m":
            meta = json.loads(meta_path.read_text())
            sky_path = Path(paths["images"]["allscale_sky"])
            sources.append(sky_path)
            sky, sky_header = binned_sky_model(sky_path, settings["model_pixel_arcsec"])
            print("Per-pointing ACA PB response models", flush=True)
            models["allscale_native"] = _predict(rows, meta, pairs, sky, sky_header, settings, target_pb=False)
            models["allscale_target"] = _predict(rows, meta, pairs, sky, sky_header, settings, target_pb=True)
            _, _, radius = row_values(rows)
            low = radius < maximum
            area = np.pi/(4*np.log(2))*header["BMAJ"]*header["BMIN"]*(np.pi/180)**2
            image = np.where(np.isfinite(tp), tp, 0)*1e-6/area
            models["tp"] = np.full(len(radius), np.nan+0j)
            print("Per-pointing TP predictions at measured ACA uv samples", flush=True)
            models["tp"][low] = _predict(subset_rows(rows, low), meta, pairs, image, header, settings, target_pb=True)
        statistics[name] = _cell_statistics(rows, key, good, size, kappa, models)
    aca = statistics["aca_7m"]
    tp_valid = np.isfinite(aca["tp"])
    statistics["tp"] = {"count": np.where(tp_valid, aca["count"], 0),
        "radius": aca["radius"], "cross_power": aca["tp"], "cross_variance": np.zeros(size)}
    records = per_pointing_profiles(statistics, pairs, edges, sectors)
    _write_rows_csv(out/"profiles.csv", records)
    _write_rows_csv(out/"pointings.csv", pairs)
    np.savez_compressed(out/"cell_statistics.npz", **{name+"_"+k: v
        for name, data in statistics.items() for k, v in data.items()})
    eligible = {name: int(data["count"][data["count"] >= 4].sum()) for name, data in statistics.items()}
    selected = {name: sum(r["rows"] for r in records if r["array"] == name) for name in COLORS}
    positive = np.array([r["rms_mjy"] if r["detected"] else r["upper_rms_mjy"] for r in records])
    positive = positive[np.isfinite(positive) & (positive > 0)]
    ylim = (positive.min()*.55, positive.max()*3.5)
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 10,
                         "pdf.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": .7})
    page_paths = []
    for page, start in enumerate(range(0, len(pairs), 6), 1):
        fig, axes = plt.subplots(3, 2, figsize=(184/25.4, 240/25.4), sharex=True, sharey=True)
        fig.subplots_adjust(left=.12, right=.96, bottom=.09, top=.97, wspace=.15, hspace=.12)
        subset = pairs[start:start+6]
        for axis, pair in zip(axes.flat, subset):
            _draw_panel(axis, records, pair, ylim)
        for axis in list(axes.flat)[len(subset):]:
            fig.delaxes(axis)
        for axis in fig.axes:
            axis.tick_params(labelbottom=True)
        fig.supxlabel(r"uv distance (k$\lambda$)", y=.025, fontsize=10)
        fig.supylabel("RMS visibility amplitude (mJy)", x=.02, fontsize=10)
        outputs = [out/f"pointing_profiles_page{page}.{ext}" for ext in ("pdf", "svg", "png")]
        save_figure(fig, outputs, dpi=300, bbox_inches=None)
        page_paths.append(outputs[0])
        plt.close(fig)
    combined = pymupdf.open()
    for path in page_paths:
        with pymupdf.open(path) as page:
            combined.insert_pdf(page)
    combined.save(out/"pointing_profiles_atlas.pdf")
    combined.close()
    fig, axes = plt.subplots(5, 4, figsize=(368/25.4, 440/25.4), sharex=True, sharey=True)
    fig.subplots_adjust(left=.07, right=.98, bottom=.06, top=.98, wspace=.12, hspace=.12)
    for axis, pair in zip(axes.flat, pairs):
        _draw_panel(axis, records, pair, ylim)
    fig.delaxes(axes.flat[-1])
    for axis in fig.axes:
        axis.tick_params(labelbottom=True)
    fig.supxlabel(r"uv distance (k$\lambda$)", y=.015, fontsize=10)
    fig.supylabel("RMS visibility amplitude (mJy)", x=.015, fontsize=10)
    outputs = [out/f"pointing_profiles_overview.{ext}" for ext in ("pdf", "svg", "png")]
    save_figure(fig, outputs, dpi=300, bbox_inches=None)
    plt.close(fig)
    manifest = {"date": "2026-10-27", "pointings": len(pairs), "uv_edges_klambda": edges.tolist(),
        "tp_maximum_klambda": float(maximum), "tp_beam_arcsec": float(theta*648000/np.pi),
        "input_rows": counts, "eligible_rows": eligible, "plotted_rows": selected,
        "ylim_mjy": list(ylim), "xlim_klambda": [1, 1000], "base_font_pt": 10,
        "method": "Per-pointing independent-hand cross-power RMS. Natural weights within cells, equal weights between cells; common occupied cells across arrays present where possible, own cells otherwise; minimum four rows per cell. ACA is translated using all-scale model power ratios to the actual matched 12m PB. TP model uses the exact same ACA uv samples, weights and cells below its 10% beam-response limit. No fitted normalization, restoring-beam matching or extrapolation below the ACA sampling support.",
        "limitations": "ACA PB translation is model-dependent. TP retains its 18.1 arcsec beam, so a smooth transition is not guaranteed. Limits omit calibration/model covariance. Row counts reflect cell selection and are not independent sample counts.",
        "sources_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    (out/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return manifest


def run_pointing_range(paths, config, *, save_figure, overwrite=False):
    """Overlay all pointings and evaluate their full low-frequency TP curves."""
    import csv
    from astropy.io import fits
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.lines import Line2D

    root = Path(paths["output_root"])
    out = root/"pointing_profiles"
    display_scale = float(config.get("pointing_range", {}).get("tp_display_scale", 1.0))
    stem = out/("pointing_profiles_overlay_display_scaled" if display_scale != 1 else "pointing_profiles_overlay")
    if stem.with_suffix(".pdf").exists() and not overwrite:
        raise FileExistsError("Single-axis pointing plot exists; use --overwrite")
    with (out/"profiles.csv").open() as handle:
        profiles = list(csv.DictReader(handle))
    coverage_path = root/"full_sampling/summary.json"
    coverage = json.loads(coverage_path.read_text())
    pairs = coverage["target_pointings"]
    summary_path = out/"manifest.json"
    summary = json.loads(summary_path.read_text())
    nb, ns = len(summary["uv_edges_klambda"])-1, int(config["overlap_tests"]["orientation_sectors"])
    settings = config["overlap_tests"]
    with np.load(out/"cell_statistics.npz") as z:
        reference = np.flatnonzero(z["aca_7m_count"] >= 4)
    reference = reference[(reference//ns) % nb == 0]
    meta_path = root/"overlap_tests/aca_7m_metadata.json"
    cache_path = root/"overlap_tests/aca_7m_rows.npz"
    metadata = json.loads(meta_path.read_text())
    with np.load(cache_path) as z:
        fields, frequencies, ddids = z["field"], z["frequency"], z["ddid"]
    tp_path = Path(paths["reduction_inputs"]["total_power"])
    tp, header = fits.getdata(tp_path, header=True)
    if header["BUNIT"].strip().lower() != "jy/beam":
        raise ValueError("TP input must have Jy/beam units")
    area = np.pi/(4*np.log(2))*header["BMAJ"]*header["BMIN"]*(np.pi/180)**2
    image = np.where(np.isfinite(tp), tp, 0)*1e-6/area
    radius = np.r_[0., np.geomspace(1., summary["tp_maximum_klambda"], 96)]
    tp_records = []
    for pair in pairs:
        power, ncells = np.zeros(len(radius)), 0
        ra, dec = np.rad2deg(np.asarray(metadata["subtables"]["FIELD"][pair["aca_field"]]["PHASE_DIR"]).ravel())
        selected = reference[reference//(2*nb*ns) == pair["pair"]]
        for band in (0, 1):
            cells = selected[(selected//(nb*ns)) % 2 == band]
            if not len(cells):
                continue
            use = (fields == pair["aca_field"]) & np.isin(ddids, [0, 7]) & ((frequencies > 225e9) == bool(band))
            frequency = float(np.median(frequencies[use]))
            angles = ((cells % ns)+.5)*np.pi/ns
            u, v = radius[:, None]*1000*np.cos(angles), radius[:, None]*1000*np.sin(angles)
            values = sample_image_visibilities(image, header, ra % 360, dec, u.ravel(), v.ravel(),
                fwhm_arcsec=settings["pb_12m_fwhm_arcsec"]*settings["pb_reference_frequency_ghz"]*1e9/frequency,
                pixel_arcsec=settings["model_pixel_arcsec"], patch_arcsec=settings["model_patch_arcsec"],
                padding=settings["model_padding_factor"], pb_center=(pair["ra_deg"], pair["dec_deg"]))
            power += np.sum(np.abs(values.reshape(u.shape))**2, axis=1)
            ncells += len(cells)
        if not ncells:
            raise ValueError("A pointing has no TP reference cells")
        amplitude = 1000*np.sqrt(power/ncells)
        tp_records.extend({"pair": pair["pair"], "uv_klambda": float(u),
            "rms_mjy": float(a), "reference_cells": ncells} for u, a in zip(radius, amplitude))
    tp_table = out/"tp_continuous_profiles.csv"
    _write_rows_csv(tp_table, tp_records)

    plt.rcParams.update({"font.family": "sans-serif", "font.size": 10,
                         "pdf.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": .7})
    fig, axis = plt.subplots(figsize=(184/25.4, 110/25.4))
    fig.subplots_adjust(left=.105, right=.97, bottom=.135, top=.98)
    positive_values = [r["rms_mjy"]*display_scale for r in tp_records if r["uv_klambda"] >= 1]
    for pair in pairs:
        for name, color in COLORS.items():
            rows = [r for r in profiles if int(r["pair"]) == pair["pair"] and r["array"] == name]
            x = np.array([float(r["uv_klambda"]) for r in rows])
            y = np.array([float(r["rms_mjy"]) for r in rows])
            detected = np.array([r["detected"] == "True" for r in rows])
            upper = np.array([float(r["upper_rms_mjy"]) for r in rows])
            if name == "tp":
                y = y*display_scale
                reference_rows = [r for r in tp_records if r["pair"] == pair["pair"] and r["uv_klambda"] >= 1]
                axis.plot([r["uv_klambda"] for r in reference_rows], [r["rms_mjy"]*display_scale for r in reference_rows],
                          color=color, linewidth=.8, alpha=.45)
                axis.scatter(x, y, color=color, marker="D", s=7, linewidths=0, alpha=.45, zorder=4)
            else:
                axis.plot(x, np.where(detected, y, np.nan), color=color, linewidth=.8, alpha=.45)
                if (~detected).any():
                    axis.errorbar(x[~detected], upper[~detected], yerr=upper[~detected]*.15,
                        uplims=True, fmt="none", color=color, alpha=.22, elinewidth=.6, capsize=1.2, markersize=2)
            positive_values.extend(y[detected].tolist())
            positive_values.extend(upper[~detected].tolist())
    positive_values = np.array(positive_values)
    positive_values = positive_values[np.isfinite(positive_values) & (positive_values > 0)]
    axis.set(xlabel=r"uv distance (k$\lambda$)", ylabel="Amplitude (mJy)",
             xlim=(1, 1000), ylim=(positive_values.min()*.55, positive_values.max()*1.8))
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xticks([1, 10, 100, 1000], labels=["1", "10", "100", "1000"])
    axis.tick_params(direction="in", top=True, right=True, labelsize=10)
    handles = [Line2D([], [], color=COLORS[name], linewidth=1.8, label=LABELS[name]) for name in COLORS]
    axis.legend(handles=handles, loc="upper right", fontsize=10, frameon=False)
    axis.text(.03, .06, "19 pointings per colour\n↓ Upper limits", transform=axis.transAxes,
              fontsize=10, va="bottom", linespacing=1.3)
    outputs = [stem.with_suffix("."+ext) for ext in ("pdf", "svg", "png")]
    save_figure(fig, outputs, dpi=600, bbox_inches=None)
    plt.close(fig)
    sources = [out/"profiles.csv", out/"pointings.csv", out/"cell_statistics.npz", summary_path,
               coverage_path, cache_path, meta_path, tp_path]
    manifest = {"date": "2026-10-28", "pointings": len(pairs), "width_mm": 184,
        "height_mm": 110, "tp_display_scale": display_scale,
        "tp_display_scale_scope": "Only plotted TP curves and diamonds. Native source amplitudes remain unchanged. The factor must be disclosed in the caption; it is not an adopted flux-calibration correction.",
        "tp_display_scale_choice": "2.7 is the rounded median ACA/TP amplitude ratio over 74 detected, matched pointing/bin comparisons below 9.165 klambda (unrounded log-space median 2.710615348473197). Illustrative normalization only.",
        "base_font_pt": 10, "xlim_klambda": [1, 1000], "tp_source_csv": str(tp_table),
        "tp_continuous_reference": "Each pointing uses its first-bin sideband/orientation-cell population, sampled at sector-centre angles. Direct TP Fourier evaluation through actual matched 12m PB from 0 to 9.165 klambda. Positive 1-9.165 klambda shown on log axes; zero retained in CSV. Curves are model references, not new observed visibilities.",
        "tp_matched_points": "Grey diamonds retain exact measured ACA uv sampling and cell weights in each comparison bin.",
        "display": "All 19 curves per array, thin lines at alpha=0.45. TP has the separately recorded display normalization; no averaging or confidence shading. Non-detections remain faint conditional upper-limit arrows. The shortened axis label denotes the same RMS signal-amplitude statistic, defined in the caption.",
        "limitations": summary["limitations"],
        "sources_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    stem.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return manifest
