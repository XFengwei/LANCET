"""Continuum pointings, amplitudes and uv-coverage figures."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from .combination_tests import sample_image_visibilities
from .visibility_diagnostics import C_METRES_PER_SECOND


COLORS = {"configuration_c5": "#7254A3", "configuration_c2": "#007A9E", "aca_7m": "#D87821"}
LABELS = {"configuration_c5": "12 m C5", "configuration_c2": "12 m C2", "aca_7m": "ACA 7 m"}


def tp_amplitude_reference(path, statistics_path, summary, minimum_response):
    """Low-frequency TP-model RMS through the same pair-centred 12 m PB.

    The native 18.1-arcsec restoring beam is retained. The reference follows
    the ACA pointing/band/orientation population of its shortest sampled
    bin, with the actual matched 12 m pointing PB applied before the FFT. It is
    a model reference, not sampled TP visibilities or an absolute flux test.
    """
    from astropy.io import fits

    data, header = fits.getdata(path, header=True)
    if header.get("BUNIT", "").strip().lower() != "jy/beam":
        raise ValueError("Expected Jy/beam TP input")
    area = np.pi/(4*np.log(2))*header["BMAJ"]*header["BMIN"]*(np.pi/180)**2
    image = np.where(np.isfinite(data), data, 0) * 1e-6 / area
    theta = np.sqrt(header["BMAJ"] * header["BMIN"]) * np.pi/180
    half_klambda = 2*np.log(2)/(np.pi*theta)/1000
    maximum = half_klambda * np.sqrt(np.log(1/minimum_response)/np.log(2))
    nyquist = 1/(2*max(abs(header["CDELT1"]), abs(header["CDELT2"]))*np.pi/180)/1000
    if maximum >= nyquist:
        raise ValueError("Requested TP reference extends beyond the native pixel Nyquist limit")
    radius = np.linspace(0, maximum, 65)
    settings = summary["settings"]
    nb, ns = len(summary["uv_edges_klambda"])-1, settings["orientation_sectors"]
    with np.load(statistics_path) as z:
        keys = np.flatnonzero(z["aca_count"] >= 4)
    first_bin = int(np.min((keys//ns) % nb))
    keys = keys[(keys//ns) % nb == first_bin]
    total_power = np.zeros_like(radius)
    number = 0
    for pair in summary["target_pointings"]:
        for band, frequency_ghz in [(0, 216.74), (1, 230.95)]:
            selected = keys[keys//(ns*nb*2) == pair["pair"]]
            selected = selected[(selected//(ns*nb)) % 2 == band]
            if not len(selected):
                continue
            angle = ((selected % ns)+0.5)*np.pi/ns
            u = radius[:, None]*1000*np.cos(angle[None, :])
            v = radius[:, None]*1000*np.sin(angle[None, :])
            prediction = sample_image_visibilities(image, header, pair["ra_deg"], pair["dec_deg"],
                u.ravel(), v.ravel(), fwhm_arcsec=settings["pb_12m_fwhm_arcsec"]*settings["pb_reference_frequency_ghz"]/frequency_ghz,
                pixel_arcsec=settings["model_pixel_arcsec"], patch_arcsec=settings["model_patch_arcsec"],
                padding=settings["model_padding_factor"])
            total_power += np.sum(np.abs(prediction.reshape(u.shape))**2, axis=1)
            number += len(selected)
    return radius, np.sqrt(total_power/number)*1000, {
        "beam_fwhm_arcsec": float(np.rad2deg(theta)*3600),
        "half_response_klambda": float(half_klambda), "maximum_plotted_klambda": float(maximum),
        "minimum_beam_response": minimum_response, "reference_aca_cells": number,
        "reference_aca_bin_klambda": summary["uv_edges_klambda"][first_bin:first_bin+2],
        "native_pixel_nyquist_klambda": float(nyquist),
        "nonfinite_pixels_zero_filled": int((~np.isfinite(data)).sum()),
        "description": "Noise-free TP model multiplied by the actual matched 12m pointing PB before Fourier transformation. Native TP restoring beam retained; shortest-ACA-bin pointing/band/sector population; no TP beam deconvolution or fitted normalization.",
    }


def make_combination_figure(paths, config, *, save_figure, overwrite=False, three_panel=False):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Circle, ConnectionPatch, Rectangle
    from matplotlib.ticker import FixedLocator

    root = Path(paths["output_root"])
    settings = config["figure1"]
    tp_display_scale = float(settings.get("tp_display_scale", 1.0))
    stem = "Figure1_combination_diagnostics" + ("_three_panel" if three_panel else "")
    outputs = [root/"figures"/f"{stem}.{suffix}" for suffix in settings["output_formats"]]
    for output in outputs:
        if output.exists() and not overwrite:
            raise FileExistsError(f"Output exists: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
    test = root/"overlap_tests"
    sampling = root/"full_sampling"
    if not (sampling/"amplitudes.csv").is_file():
        from .full_sampling import build_full_sampling_profiles

        build_full_sampling_profiles(paths, config)
    coverage = json.loads((sampling/"summary.json").read_text())
    if coverage["uv_edges_klambda"] != config["figure1"]["amplitude_uv_edges_klambda"] or coverage["settings"] != config["overlap_tests"]:
        raise ValueError("Full-sampling settings changed; rerun the full-sampling stage")
    from .full_sampling import build_separate_array_profiles

    array_profiles = build_separate_array_profiles(paths, coverage)
    # Preserve observation-level density without treating noisy individual
    # scalar amplitudes as PB-translated, noise-corrected measurements.
    observation_sources = {
        "aca_7m": test/"aca_7m_rows.npz",
        "configuration_c2": sampling/"configuration_c2_rows.npz",
        "configuration_c5": sampling/"configuration_c5_rows.npz",
    }
    density_edges = np.geomspace(1, 1000, 46)
    density_centres = np.sqrt(density_edges[:-1]*density_edges[1:])
    density_rows, observation_density = [], []
    for name, source in observation_sources.items():
        with np.load(source) as z:
            use = np.ones(len(z["observation"]), dtype=bool)
            if name == "aca_7m":
                use &= np.isin(z["field"], [p["aca_field"] for p in coverage["target_pointings"]])
                use &= np.isin(z["ddid"], [0, 7])
            radius = np.hypot(z["uvw"][:, 0], z["uvw"][:, 1])*z["frequency"]/C_METRES_PER_SECOND/1000
            for observation in np.unique(z["observation"][use]):
                selected = use & (z["observation"] == observation)
                counts, _ = np.histogram(radius[selected], bins=density_edges)
                if counts.sum() != selected.sum():
                    raise ValueError("Observation samples extend beyond the displayed density range")
                observation_density.append((name, int(observation), counts))
                density_rows.extend((name, int(observation), lo, hi, int(n))
                    for lo, hi, n in zip(density_edges[:-1], density_edges[1:], counts))
    density_table = root/"figures"/(stem+"_observation_sampling.csv" if three_panel else "Figure1_observation_sampling.csv")
    with density_table.open("w") as handle:
        writer = csv.writer(handle)
        writer.writerow(["array", "observation_id", "uv_min_klambda", "uv_max_klambda", "visibility_rows"])
        writer.writerows(density_rows)
    metadata = {name: json.loads((test/f"{name}_metadata.json").read_text()) for name in ("configuration_c5", "aca_7m")}
    coords = {}
    for name in COLORS:
        with np.load(root/"uv_coordinates"/f"{name}.npz") as z:
            if not bool(z["unflagged_only"]):
                raise ValueError("Figure 1 requires the flag-aware UV cache")
            coords[name] = (z["u"], z["v"])
    tp_u, tp_amp, tp_meta = tp_amplitude_reference(
        paths["reduction_inputs"]["total_power"], sampling/"cell_statistics.npz",
        coverage, float(settings["tp_minimum_response"]))
    tp_table = root/"figures"/(stem+"_tp_reference.csv" if three_panel else "Figure1_tp_reference.csv")
    np.savetxt(tp_table, np.column_stack([tp_u, tp_amp]), delimiter=",",
               header="uv_klambda,tp_model_rms_mjy", comments="")
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 10,
        "axes.linewidth": 0.7, "axes.spines.top": True, "axes.spines.right": True,
        "pdf.fonttype": 42, "svg.fonttype": "none", "legend.frameon": False})
    if three_panel:
        # Keep three equal square plot areas at the manuscript text width.
        height = settings["width_inches"]*.885/3.64/.74
        fig, axes = plt.subplots(1, 3, figsize=(settings["width_inches"], height))
        fig.subplots_adjust(left=.10, right=.985, bottom=.20, top=.94, wspace=.32)
        a, c, b = axes
        d = c.inset_axes([.66, .66, .34, .34], zorder=8)
        d.patch.set_facecolor("white")
    else:
        fig, axes = plt.subplots(2, 2, figsize=(settings["width_inches"], settings["height_inches"]))
        fig.subplots_adjust(left=0.105, right=0.96, bottom=0.105, top=0.96, wspace=0.25, hspace=0.25)
        a, b, c, d = axes.flat
    for axis in axes.flat:
        axis.tick_params(direction="in", top=True, right=True, length=3, width=0.7)
    def title(axis, text):
        axis.text(0.04, .96, text, transform=axis.transAxes,
                  ha="left", va="top", fontsize=10 if three_panel else 11,
                  bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.94, "pad": 2}, zorder=20)

    ra0, dec0 = 221.25, -59.8
    for name, marker, label in [("aca_7m", "o", "ACA 7 m"), ("configuration_c5", "+", "12 m C2+C5")]:
        pos = np.array([f["PHASE_DIR"] for f in metadata[name]["subtables"]["FIELD"]]).squeeze()
        ra, dec = np.rad2deg(pos[:, 0]) % 360, np.rad2deg(pos[:, 1])
        a.scatter((ra-ra0)*np.cos(np.deg2rad(dec0))*60, (dec-dec0)*60, s=7 if three_panel else 12,
                  color=COLORS[name] if name == "aca_7m" else COLORS["configuration_c2"], marker=marker,
                  linewidths=0.7, alpha=0.75, label=("7 m" if name == "aca_7m" else "12 m") if three_panel else label)
    a.set(xlabel="RA offset (arcmin)", ylabel="Dec offset (arcmin)")
    a.set_aspect("equal", adjustable="datalim")
    a.invert_xaxis()
    a.margins(0.12)
    if three_panel:
        a.set_xticks([5, 0, -5])
    title(a, "(a)")
    a.legend(loc="upper left", bbox_to_anchor=(0.03, 0.88), fontsize=10, labelspacing=0.35, handletextpad=0.4)

    values = list(tp_amp*tp_display_scale)
    for name, marker, label in [("aca_7m", "o", "7 m (12 m PB)"),
                                ("configuration_c2", "s", "12 m C2"),
                                ("configuration_c5", "^", "12 m C5")]:
        profile = [r for r in array_profiles if r["array"] == name]
        radius = np.array([r["uv_klambda"] for r in profile])
        amplitude = np.array([r["rms_mjy"] for r in profile])
        detected = np.array([r["detected"] for r in profile])
        upper = np.array([r["upper_rms_mjy"] for r in profile])
        b.plot(radius, np.where(detected, amplitude, np.nan), color=COLORS[name],
               linewidth=.9 if three_panel else 1.4, marker=marker, markersize=1.5 if three_panel else 3,
               label={"aca_7m": "7 m", "configuration_c2": "C2", "configuration_c5": "C5"}[name] if three_panel else label)
        if (~detected).any():
            b.errorbar(radius[~detected], upper[~detected], yerr=upper[~detected]*0.25,
                uplims=True, fmt="none", color=COLORS[name], elinewidth=0.8,
                capsize=1.2 if three_panel else 2, markersize=1.5 if three_panel else 3)
        values.extend(amplitude[detected])
        values.extend(upper[~detected])
    positive_tp = tp_u >= 1 if three_panel else tp_u > 0
    b.plot(tp_u[positive_tp], tp_amp[positive_tp]*tp_display_scale, color="0.45", linestyle="-.", linewidth=.9 if three_panel else 1.2, label="TP" if three_panel else "TP (12 m PB)")
    b.set(xlabel=r"uv distance (k$\lambda$)", ylabel="Amplitude (mJy)", xlim=(1, 1000))
    b.set_xscale("log")
    b.set_xticks([1, 10, 100, 1000], labels=["1", "10", "100", "1000"])
    b.set_yscale("log")
    values = np.asarray(values)
    positive = values[np.isfinite(values) & (values > 0)]
    count_top = .40 if three_panel else .35
    b.set_ylim(positive.min()*(.001 if three_panel else .006), positive.max()*(8.0 if three_panel else 4.6))
    sampling_top = b.get_ylim()[0]*(b.get_ylim()[1]/b.get_ylim()[0])**count_top
    major_ticks = [10, 100, 1000, 10000] if three_panel else b.yaxis.get_majorticklocs()
    b.yaxis.set_major_locator(FixedLocator([v for v in major_ticks if v > sampling_top]))
    b.yaxis.set_minor_locator(FixedLocator([v for v in b.yaxis.get_minorticklocs() if v > sampling_top]))
    overlap = coverage["overlap_klambda"]
    maximum_count = max(counts.max() for _, _, counts in observation_density)
    group = np.array([list(observation_sources).index(name) for name, _, _ in observation_density])
    track_y = np.linspace(.24, .04, len(observation_density)) if three_panel else 0.25 - 0.015*np.arange(len(observation_density)) - 0.035*group
    for (name, _, counts), y in zip(observation_density, track_y):
        occupied = counts > 0
        b.scatter(density_centres[occupied], np.full(occupied.sum(), y),
                  s=(4 if three_panel else 10)*counts[occupied]/maximum_count, color=COLORS[name], linewidths=0,
                  transform=b.get_xaxis_transform(), zorder=3)
    for name in (() if three_panel else observation_sources):
        indices = [i for i, (array, _, _) in enumerate(observation_density) if array == name]
        label = {"aca_7m": "7 m", "configuration_c2": "C2", "configuration_c5": "C5"}[name]
        b.text(0.025, track_y[indices].mean(), label, color=COLORS[name], fontsize=10,
               va="center", transform=b.transAxes)
    b.plot([0, 1], [count_top, count_top], color="0.8", linewidth=0.5, transform=b.transAxes)
    b.text(0.025, .30 if three_panel else .29, "Visibility counts", transform=b.transAxes, fontsize=10)
    if three_panel:
        b.get_xticklabels()[-1].set_horizontalalignment("right")
    title(b, "(c)" if three_panel else "(b)")
    b.legend(loc="upper right", bbox_to_anchor=(0.99, 0.98), fontsize=10, labelspacing=0.25,
             ncol=2 if three_panel else 1, columnspacing=.4,
             borderpad=.2 if three_panel else .3, handlelength=1.1 if three_panel else 1.8,
             handletextpad=0.4, frameon=True, facecolor="white", framealpha=1, edgecolor="none")

    zoom = 25.0 if three_panel else float(config["uv_plot"]["zoom_limit_metres"])
    limit = 1.10*max(np.abs(np.r_[u, v]).max() for u, v in coords.values())
    wavelength = C_METRES_PER_SECOND/(config["visibility_profile"]["reference_frequency_ghz"]*1e9)
    tp_half_metres = tp_meta["half_response_klambda"]*1000*wavelength
    tp_outer_metres = tp_meta["maximum_plotted_klambda"]*1000*wavelength
    d.add_patch(Circle((0, 0), tp_outer_metres, facecolor="0.80", edgecolor="0.55", linewidth=0.6, zorder=0))
    d.add_patch(Circle((0, 0), tp_half_metres, facecolor="0.65", edgecolor="0.4", linewidth=0.6, linestyle=":", zorder=1))
    for axis, extent in [(c, limit), (d, zoom)]:
        for name, (u, v) in coords.items():
            if three_panel and axis is d and name != "aca_7m":
                continue
            use = (np.abs(u) <= extent) & (np.abs(v) <= extent)
            axis.scatter(np.r_[u[use], -u[use]], np.r_[v[use], -v[use]], s=(.25 if three_panel else .45) if axis is c else .75,
                         color=COLORS[name], alpha=0.25 if name != "aca_7m" else 0.60,
                         linewidths=0, rasterized=True, zorder=2)
        axis.set(xlim=(-extent, extent), ylim=(-extent, extent), xlabel="u (m)", ylabel="v (m)")
        axis.set_aspect("equal", adjustable="box")
    c.add_patch(Rectangle((-zoom, -zoom), 2*zoom, 2*zoom, facecolor="none", edgecolor="0.2", linewidth=0.7, zorder=3))
    d.scatter([0], [0], marker="+", s=28, linewidths=1.1, color="0.15", zorder=5)
    handles = [Line2D([], [], marker="o", linestyle="none", color=COLORS[k], markersize=4, label=LABELS[k]) for k in COLORS]
    c.legend(handles=handles, loc="lower left", fontsize=10, labelspacing=0.25, handletextpad=0.4,
             frameon=True, facecolor="white", framealpha=1, edgecolor="none")
    tp_handle = Line2D([], [], marker="+", linestyle="none", color="0.3", markersize=6, label="TP model")
    if three_panel:
        c.get_legend().remove()
        short = {"configuration_c5": "C5", "configuration_c2": "C2", "aca_7m": "7 m"}
        handles = [Line2D([], [], marker="o", linestyle="none", color=COLORS[k], markersize=3, label=short[k]) for k in COLORS]
        tp_handle.set_label("TP")
        c.legend(handles=handles+[tp_handle], loc="lower left", ncol=2, fontsize=10,
                 columnspacing=.35, handlelength=.7, handletextpad=.25, labelspacing=.15,
                 frameon=True, facecolor="white", framealpha=1, edgecolor="none")
        c.set(xlabel="u (km)", ylabel="v (km)")
        c.set_xticks([-1000, 0, 1000], labels=["−1", "0", "1"])
        c.set_yticks([-1000, 0, 1000], labels=["−1", "0", "1"])
        d.set(xlabel="", ylabel="", xticks=[], yticks=[])
        d.minorticks_off()
        for origin, destination in [((-zoom, zoom), (0, 0)), ((zoom, -zoom), (1, 0))]:
            c.add_artist(ConnectionPatch(xyA=origin, coordsA=c.transData,
                xyB=destination, coordsB=d.transAxes, color="black", linewidth=.65,
                linestyle="-", zorder=7, clip_on=False))
    else:
        d.legend(handles=handles+[tp_handle], loc="lower left", ncol=2, fontsize=10, columnspacing=0.5,
                 handletextpad=0.3, labelspacing=0.25, frameon=True, facecolor="white", framealpha=1, edgecolor="none")
    title(c, "(b)" if three_panel else "(c)")
    if not three_panel:
        title(d, "(d)")
    save_figure(fig, outputs, dpi=int(settings["dpi"]), bbox_inches=None,
                exclude_axes=[d] if three_panel else ())
    plt.close(fig)
    sources = [sampling/"amplitudes.csv", sampling/"separate_array_amplitudes.csv", sampling/"summary.json", sampling/"cell_statistics.npz", test/"summary.json",
               *observation_sources.values(),
               *[root/"uv_coordinates"/f"{n}.npz" for n in COLORS], Path(paths["reduction_inputs"]["total_power"])]
    amplitude_panel = "c" if three_panel else "b"
    manifest = {"date": "2026-10-28", "status": "Three-panel Figure 1 selected for manuscript; four-panel output retained" if three_panel else "Four-panel output retained; verified adoption provenance is recorded by the referee-response task",
                "tp": tp_meta, "tp_half_response_metres": tp_half_metres,
                "tp_ten_percent_response_metres": tp_outer_metres,
                f"panel_{amplitude_panel}_statistic": "Separate ACA/C2/C5 RMS from independent-hand cross-power in the same 19 pointing pairs; common occupied cells across arrays present in each radial bin where possible, own occupied cells otherwise (CSV matched_arrays flag). ACA model translation targets actual 12m pointing PB. TP has this same PB applied before FFT, with its native restoring beam retained. Non-detections have conditional two-standard-error limits.",
                "figure_width_mm": settings["width_inches"]*25.4, "base_font_pt": 10,
                "interferometer_display_scale": 1.0,
                "tp_display_scale": tp_display_scale,
                "tp_display_scale_scope": f"Panel {amplitude_panel} curve only; native TP source CSV is unscaled. Disclosed as a display normalization in the manuscript caption; no adopted calibration correction.",
                "amplitude_ranges_klambda": coverage["ranges_klambda"], "amplitude_overlap_klambda": overlap,
                f"panel_{amplitude_panel}_xscale": "log; 1-1000 klambda; zero spacing remains in the short-baseline view and the source CSV",
                "panel_order": ["pointings", "uv coverage with short-baseline inset", "amplitudes and visibility counts"] if three_panel else ["pointings", "amplitudes and visibility counts", "full uv coverage", "short baselines"],
                "panel_labels": ["(a)", "(b)", "(c)"] if three_panel else ["(a)", "(b)", "(c)", "(d)"],
                "inset_limit_metres": zoom if three_panel else None,
                "main_uv_display_unit": "km" if three_panel else "m",
                "inset_uv_display_unit": "m" if three_panel else None,
                "inset_arrays": ["aca_7m", "TP model response"] if three_panel else None,
                "inset_background": "white; flush with upper-right parent edges; two solid black zoom connectors" if three_panel else None,
                "inset_ticks": False if three_panel else True,
                "inset_tick_labels": False if three_panel else True,
                f"panel_{amplitude_panel}_marker_size_pt": 1.5 if three_panel else 3,
                "alignment_exclusion": "Upper-right inset in b is intentionally a smaller zoom, not an equal-span main panel" if three_panel else None,
                f"panel_{amplitude_panel}_native_aca_curve": False, f"panel_{amplitude_panel}_shading": False,
                f"panel_{amplitude_panel}_sampling": {
                    "source_csv": str(density_table), "log_bins": len(density_centres),
                    "marker_area": f"{4 if three_panel else 10} pt^2 times row count / global maximum bin count",
                    "maximum_bin_count": int(maximum_count),
                    "selection": "Same 19 pointing pairs and broad windows; all valid cached rows before common-cell selection; no imaging weights or thinning. Vertical track position is an observation index, not amplitude.",
                    "observations": [{"array": name, "observation_id": obs, "rows": int(n.sum())}
                                     for name, obs, n in observation_density]},
                "matched_pair_markers": False,
                "sources_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    outputs[0].with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return outputs
