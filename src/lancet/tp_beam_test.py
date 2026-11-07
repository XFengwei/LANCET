"""Conditional TP/ACA comparison with smoothing applied before the PB."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from .combination_tests import (
    binned_sky_model, matching_cells, row_values, sample_image_visibilities,
    subset_rows, _write_rows_csv,
)
from .full_sampling import _cell_statistics
from .visibility_diagnostics import C_METRES_PER_SECOND


def smooth_to_circular_beam(image, header, target_arcsec):
    """Convolve MJy/sr with a unit-integral, deconvolved Gaussian kernel."""
    from astropy.convolution import Gaussian2DKernel
    from astropy.wcs import WCS
    from scipy.signal import fftconvolve

    sigma = np.array([header["BMAJ"], header["BMIN"]])*3600/np.sqrt(8*np.log(2))
    angle = np.deg2rad(header.get("BPA", 0))
    rotation = np.array([[np.sin(angle), np.cos(angle)], [np.cos(angle), -np.sin(angle)]])
    native = rotation @ np.diag(sigma**2) @ rotation.T
    target = np.eye(2)*(target_arcsec/np.sqrt(8*np.log(2)))**2
    matrix = WCS(header, naxis=2).celestial.pixel_scale_matrix*3600
    inverse = np.linalg.inv(matrix)
    covariance = inverse @ (target-native) @ inverse.T
    eigenvalues, vectors = np.linalg.eigh(covariance)
    if np.any(eigenvalues <= 0):
        raise ValueError("Target beam must exceed the native sky-model beam")
    extent = 2*int(np.ceil(6*np.sqrt(eigenvalues.max())))+1
    kernel = Gaussian2DKernel(*np.sqrt(eigenvalues), theta=np.arctan2(vectors[1, 0], vectors[0, 0]),
                              x_size=extent, y_size=extent)
    result = fftconvolve(image, kernel.array, mode="same")
    info = {"native_beam_arcsec": (np.array([header["BMAJ"], header["BMIN"]])*3600).tolist(),
            "native_pa_deg": float(header.get("BPA", 0)), "target_circular_beam_arcsec": target_arcsec,
            "kernel_sum": float(kernel.array.sum()), "kernel_covariance_pixels2": covariance.tolist(),
            "kernel_support": "six standard deviations of the widest axis"}
    return result, info


def _predict(rows, meta, pairs, image, header, settings, *, target_pb):
    prediction = np.full(len(rows["x"]), np.nan+0j)
    for pair in pairs:
        field = pair["aca_field"]
        ra, dec = np.rad2deg(np.asarray(meta["subtables"]["FIELD"][field]["PHASE_DIR"]).ravel())
        for band in (False, True):
            use = (rows["field"] == field) & ((rows["frequency"] > 225e9) == band)
            if not use.any():
                continue
            frequency = float(np.median(rows["frequency"][use]))
            width = settings["pb_12m_fwhm_arcsec" if target_pb else "pb_7m_fwhm_arcsec"]
            width *= settings["pb_reference_frequency_ghz"]*1e9/frequency
            uv = rows["uvw"][use]*rows["frequency"][use, None]/C_METRES_PER_SECOND
            prediction[use] = sample_image_visibilities(image, header, ra % 360, dec,
                uv[:, 0], uv[:, 1], fwhm_arcsec=width,
                pixel_arcsec=settings["model_pixel_arcsec"], patch_arcsec=settings["model_patch_arcsec"],
                padding=settings["model_padding_factor"],
                pb_center=(pair["ra_deg"], pair["dec_deg"]) if target_pb else None)
    if not np.all(np.isfinite(prediction)):
        raise ValueError("Missing model predictions for selected ACA rows")
    return prediction


def run_tp_beam_test(paths, config, *, save_figure, overwrite=False):
    """Test smoothing without fitting a relative flux scale or changing data."""
    from astropy.io import fits
    from matplotlib import pyplot as plt

    root = Path(paths["output_root"])
    out = root/"tp_beam_test"
    out.mkdir(parents=True, exist_ok=True)
    if (out/"profiles.csv").exists() and not overwrite:
        raise FileExistsError("TP beam test exists; use --overwrite")
    coverage = json.loads((root/"full_sampling/summary.json").read_text())
    noise = json.loads((root/"overlap_tests/summary.json").read_text())["noise"]["aca_7m"]
    settings = config["overlap_tests"]
    tp_path = Path(paths["reduction_inputs"]["total_power"])
    tp, tp_header = fits.getdata(tp_path, header=True)
    if tp_header["BUNIT"].strip().lower() != "jy/beam":
        raise ValueError("TP input must have Jy/beam units")
    beam = np.sqrt(tp_header["BMAJ"]*tp_header["BMIN"])*3600
    half = 2*np.log(2)/(np.pi*(beam*np.pi/648000))/1000
    maximum = half*np.sqrt(np.log(1/config["figure1"]["tp_minimum_response"])/np.log(2))
    edges = np.r_[4.8, 6, 7, 8, maximum]
    tp_area = np.pi/(4*np.log(2))*tp_header["BMAJ"]*tp_header["BMIN"]*(np.pi/180)**2
    tp_image = np.where(np.isfinite(tp), tp, 0)*1e-6/tp_area
    cache = root/"overlap_tests/aca_7m_rows.npz"
    meta_path = root/"overlap_tests/aca_7m_metadata.json"
    meta = json.loads(meta_path.read_text())
    with np.load(cache) as z:
        rows = {k: z[k] for k in z.files}
    _, _, radius = row_values(rows)
    use = np.isin(rows["field"], [p["aca_field"] for p in coverage["target_pointings"]])
    use &= np.isin(rows["ddid"], [0, 7]) & (radius >= edges[0]) & (radius < edges[-1])
    rows = subset_rows(rows, use)
    pairs = coverage["target_pointings"]
    key, good = matching_cells(rows, pairs, "aca", edges, settings["orientation_sectors"])
    size = len(pairs)*2*(len(edges)-1)*settings["orientation_sectors"]
    kappa = np.array([noise[str(int(d))]["kappa_time"] for d in rows["ddid"]])
    predictions = {"tp": _predict(rows, meta, pairs, tp_image, tp_header, settings, target_pb=True)}
    beam_info = {}
    sources = [cache, meta_path, tp_path, root/"full_sampling/summary.json"]
    for label, image_key in [("allscale", "allscale_sky"), ("twelve_m", "twelve_m_pbcor")]:
        path = Path(paths["images"][image_key])
        print(f"TP beam test: {label}", flush=True)
        sky, header = binned_sky_model(path, settings["model_pixel_arcsec"])
        smoothed, beam_info[label] = smooth_to_circular_beam(sky, header, beam)
        predictions[label+"_native"] = _predict(rows, meta, pairs, sky, header, settings, target_pb=False)
        predictions[label+"_pb"] = _predict(rows, meta, pairs, sky, header, settings, target_pb=True)
        predictions[label+"_beam_pb"] = _predict(rows, meta, pairs, smoothed, header, settings, target_pb=True)
        sources.append(path)
    stats = _cell_statistics(rows, key, good, size, kappa, predictions)
    radial = (np.arange(size)//settings["orientation_sectors"]) % (len(edges)-1)
    profiles = []
    for i in range(len(edges)-1):
        use = (radial == i) & (stats["count"] >= 4)
        if not use.any():
            continue
        power = float(np.mean(stats["cross_power"][use]))
        if power <= 0:
            raise ValueError("Short-baseline ACA power is not positive")
        rms = 1000*np.sqrt(power)
        record = {"uv_low_klambda": float(edges[i]), "uv_high_klambda": float(edges[i+1]),
            "uv_klambda": float(np.mean(stats["radius"][use])), "rows": int(stats["count"][use].sum()),
            "cells": int(use.sum()), "aca_native_rms_mjy": float(rms),
            "tp_rms_mjy": float(1000*np.sqrt(np.mean(stats["tp"][use])))}
        for model in beam_info:
            for response in ("pb", "beam_pb"):
                factor = float(np.sqrt(np.mean(stats[model+"_"+response][use])/np.mean(stats[model+"_native"][use])))
                record[model+"_"+response+"_factor"] = factor
                record[model+"_"+response+"_rms_mjy"] = float(rms*factor)
                record[model+"_"+response+"_to_tp"] = float(rms*factor/record["tp_rms_mjy"])
        profiles.append(record)
    _write_rows_csv(out/"profiles.csv", profiles)
    np.savez_compressed(out/"cell_statistics.npz", **stats)
    support = {"aca_7m": {"selected_rows": len(rows["x"]),
                           "eligible_rows": sum(r["rows"] for r in profiles),
                           "eligible_cells": sum(r["cells"] for r in profiles)}}
    for name in ("configuration_c2", "configuration_c5"):
        path = root/"full_sampling"/f"{name}_rows.npz"
        with np.load(path) as z:
            data = {k: z[k] for k in z.files}
        _, _, radius = row_values(data)
        data = subset_rows(data, (radius >= edges[0]) & (radius < edges[-1]))
        keys, valid = matching_cells(data, pairs, "alma", edges, settings["orientation_sectors"])
        counts = np.bincount(keys[valid], minlength=size)
        support[name] = {"selected_rows": len(data["x"]), "eligible_cells": int((counts >= 4).sum()),
                         "eligible_rows": int(counts[counts >= 4].sum())}
        sources.append(path)
    manifest = {"date": "2026-10-27", "selected_rows": len(rows["x"]), "uv_edges_klambda": edges.tolist(),
        "target_beam_arcsec": float(beam), "beam_models": beam_info,
        "sampling_support": support,
        "tp_nonfinite_pixels_zero_filled": int((~np.isfinite(tp)).sum()),
        "method": "Model power ratios translate ACA native 7m PB to the actual paired 12m PB, with and without convolution to the TP restoring beam before PB multiplication. TP is evaluated at the exact same measured ACA uv coordinates and same naturally weighted/equal-cell population in each bin. No fitted normalization or visibility modification.",
        "limitations": "Conditional on restored, data-derived sky models, Gaussian beams and finite/zero-filled image support. The all-scale model contains TP information; the 12m-only model is a sensitivity check, not an external flux standard. No adopted calibration correction. Frequencies change PB widths and uv coordinates; no spectral-index scaling is inferred from the TP header.",
        "sources_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    (out/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 10, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axis = plt.subplots(figsize=(7.2440944882, 4.6))
    fig.subplots_adjust(left=.13, right=.97, bottom=.17, top=.95)
    x = [r["uv_klambda"] for r in profiles]
    for column, color, style, marker, label in [
        ("allscale_pb_rms_mjy", "#D87821", "--", "o", "ACA: PB matched"),
        ("allscale_beam_pb_rms_mjy", "#D87821", "-", "o", "ACA: beam + PB (all-scale model)"),
        ("twelve_m_beam_pb_rms_mjy", "#007A9E", ":", "^", "ACA: beam + PB (12 m model)"),
        ("tp_rms_mjy", "0.3", "-.", "s", "TP: same PB and measured uv samples")]:
        axis.plot(x, [r[column] for r in profiles], color=color, linestyle=style, marker=marker,
                  linewidth=1.5, markersize=4, label=label)
    axis.set(xlabel=r"uv distance (k$\lambda$)", ylabel="Signal RMS (mJy)", xlim=(4.8, maximum))
    axis.set_yscale("log")
    axis.set_ylim(bottom=min(r["tp_rms_mjy"] for r in profiles)*.6,
                  top=max(r["allscale_pb_rms_mjy"] for r in profiles)*4)
    axis.tick_params(direction="in", top=True, right=True)
    axis.legend(loc="upper right", fontsize=10, frameon=False)
    outputs = [out/f"tp_matched_beam_test.{extension}" for extension in ("pdf", "svg", "png")]
    save_figure(fig, outputs, dpi=600, bbox_inches=None)
    plt.close(fig)
    return manifest
