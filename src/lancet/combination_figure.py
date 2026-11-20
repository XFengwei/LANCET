"""Continuum pointings, amplitudes and uv-coverage figures."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np

from .amplitude_alignment import ALIGNMENT_OUTPUTS, pb_smoothing_note
from .combination_tests import require_products, sample_image_visibilities
from .full_sampling import FULL_SAMPLING_OUTPUTS, OVERLAP_OUTPUTS, require_stage_outputs
from .visibility_diagnostics import C_METRES_PER_SECOND


COLORS = {"configuration_c5": "#7254A3", "configuration_c2": "#007A9E", "aca_7m": "#D87821"}


def figure_dependencies(amplitude_source):
    """Outputs of the earlier diagnostic stages that Fig. 2 reads, by stage, in stage order."""
    stages = {"uv": ("uv_coordinates", tuple(f"{name}.npz" for name in COLORS)),
              "overlap": ("overlap_tests", (*OVERLAP_OUTPUTS, "configuration_c5_metadata.json")),
              "full_sampling": ("full_sampling", (*FULL_SAMPLING_OUTPUTS, "amplitudes.csv", "cell_statistics.npz"))}
    if amplitude_source == "amplitude_alignment":
        stages["amplitude_alignment"] = ("amplitude_alignment", ALIGNMENT_OUTPUTS)
    return stages


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


# Hues outside the array colours and the TP greys (here hue means array). With the array colours,
# every pair involving these is >= 15 OKLab dE x100 in normal vision and >= 6.7 under protan/deutan
# simulation (crimson/green, pink/C2); line style is the secondary encoding.
IMAGE_CURVES = {"interferometer": ("12m+7m", "--", "#B2182B"), "jcomb": ("J-comb", "-", "#006837"),
                "feather": ("Feather", ":", "#CC79A7")}


def image_curve_styles(shown, labels):
    """Style of each shown image curve; the legend labels come from the configuration."""
    return {name: (labels.get(name, IMAGE_CURVES[name][0]), *IMAGE_CURVES[name][1:]) for name in shown}


def check_image_curve_settings(shown, labels, available):
    """Fail early on a curve or label the configuration names but the figure cannot use."""
    unknown = [n for n in shown if n not in IMAGE_CURVES]
    absent = [n for n in shown if n in IMAGE_CURVES and n not in available]
    stray = [n for n in labels if n not in shown]
    if unknown or absent or stray:
        raise ValueError(f"[figure2] image curves: not in IMAGE_CURVES {unknown}; not in the alignment CSV {absent}; "
                         f"image_curve_labels keys not in image_curves_shown {stray}")


def shown_curve_values(curves, shown):
    """Amplitudes of the drawn image curves; only the alignment path has curves."""
    return [v for name in shown for v in curves[name][1]] if curves else []


def _plot_image_curves(axis, curves, shown, labels):
    """Amplitudes of the shown image curves (beam and pixel binning divided out), with their legend."""
    handles = []
    for name, (label, style, color) in image_curve_styles(shown, labels).items():
        u, rms = curves[name]
        handles += axis.plot(u, rms, color=color, linestyle=style, linewidth=1.1, label=label, zorder=4)
    # Below the direct "TP, beam-corrected" label, above the long-baseline profiles.
    axis.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.975, 0.875), fontsize=9,
                labelspacing=0.1, borderpad=.15, borderaxespad=.2, handlelength=1.7, handletextpad=0.35,
                frameon=False)


def make_combination_figure(paths, config, *, save_figure, overwrite=False):
    """Paper II Fig. 2: (a) pointings, (b) uv coverage with the short-baseline inset, (c) amplitudes."""
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Circle, ConnectionPatch, Rectangle
    from matplotlib.ticker import FixedLocator

    require_products(paths, ["jcomb_continuum", "feather_continuum"])
    root = Path(paths["output_root"])
    settings = config["figure2"]
    stem = "figure2_combination_diagnostics" + settings.get("output_suffix", "")
    outputs = [root/"figures"/f"{stem}.{suffix}" for suffix in settings["output_formats"]]
    for output in outputs:
        if output.exists() and not overwrite:
            raise FileExistsError(f"Output exists: {output}")
    amplitude_source = settings.get("amplitude_source", "full_sampling")
    require_stage_outputs(root, figure_dependencies(amplitude_source))
    for output in outputs:
        output.parent.mkdir(parents=True, exist_ok=True)
    test = root/"overlap_tests"
    sampling = root/"full_sampling"
    coverage = json.loads((sampling/"summary.json").read_text())
    if coverage["uv_edges_klambda"] != settings["amplitude_uv_edges_klambda"] or coverage["settings"] != config["overlap_tests"]:
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
    density_table = root/"figures"/f"{stem}_observation_sampling.csv"
    with density_table.open("w") as handle:
        writer = csv.writer(handle)
        writer.writerow(["array", "observation_id", "uv_min_klambda", "uv_max_klambda", "visibility_rows"])
        writer.writerows(density_rows)
    metadata = {name: json.loads((test/f"{name}_metadata.json").read_text()) for name in ("configuration_c5", "aca_7m")}
    coords = {}
    for name in COLORS:
        with np.load(root/"uv_coordinates"/f"{name}.npz") as z:
            if not bool(z["unflagged_only"]):
                raise ValueError("Figure 2 requires the flag-aware UV cache")
            coords[name] = (z["u"], z["v"])
    tp_u, tp_amp, tp_meta = tp_amplitude_reference(
        paths["inputs"]["tp_dust_ff"], sampling/"cell_statistics.npz",
        coverage, float(settings["tp_minimum_response"]))
    shown_curves = list(settings.get("image_curves_shown", IMAGE_CURVES))
    curve_labels = dict(settings.get("image_curve_labels", {}))
    tp_fixed, image_curves, image_meta, alignment_sources = None, {}, None, []
    if amplitude_source == "amplitude_alignment":
        # 7 m translated with the 12m+7m image; TP with and without its beam; image curves.
        from .amplitude_alignment import read_products

        alignment_root = root/"amplitude_alignment"
        array_profiles, _, curves, alignment, image_curves = read_products(alignment_root)
        check_image_curve_settings(shown_curves, curve_labels, image_curves)
        alignment_sources = [alignment_root/name for name in ALIGNMENT_OUTPUTS]
        tp_u, tp_amp, tp_fixed = (curves["uv_klambda"], curves["tp_raw_common_mjy"],
                                  curves["tp_corrected_common_mjy"])
        tp_meta = {**tp_meta, "amplitude_source": str(alignment_root),
                   "translation_sky_model": alignment["sky_model"],
                   "tp_response_limit_klambda": alignment["tp_response_limit_klambda"],
                   "description": alignment["method"]}
        image_table = root/"figures"/f"{stem}_image_curves.csv"
        shutil.copyfile(alignment_root/"image_curves.csv", image_table)
        image_meta = {**alignment["image_curves"], "source_csv": str(image_table),
                      "low_uv_pb_smoothing": pb_smoothing_note(config["overlap_tests"]),
                      "styles": {name: {"label": label, "linestyle": style, "color": color}
                                 for name, (label, style, color) in image_curve_styles(shown_curves, curve_labels).items()},
                      "shown": list(shown_curves)}
    tp_table = root/"figures"/f"{stem}_tp_reference.csv"
    if tp_fixed is None:
        np.savetxt(tp_table, np.column_stack([tp_u, tp_amp]), delimiter=",",
                   header="uv_klambda,tp_model_rms_mjy", comments="")
    else:
        np.savetxt(tp_table, np.column_stack([tp_u, tp_amp, tp_fixed]), delimiter=",",
                   header="uv_klambda,tp_model_rms_mjy,tp_beam_corrected_rms_mjy", comments="")
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 10,
        "axes.linewidth": 0.7, "axes.spines.top": True, "axes.spines.right": True,
        "pdf.fonttype": 42, "svg.fonttype": "none", "legend.frameon": False})
    # Keep three equal square plot areas at the manuscript text width.
    height = settings["width_inches"]*.885/3.64/.74
    fig, axes = plt.subplots(1, 3, figsize=(settings["width_inches"], height))
    fig.subplots_adjust(left=.10, right=.985, bottom=.20, top=.94, wspace=.32)
    a, b, c = axes
    inset = b.inset_axes([.66, .66, .34, .34], zorder=8)
    inset.patch.set_facecolor("white")
    for axis in axes.flat:
        axis.tick_params(direction="in", top=True, right=True, length=3, width=0.7)

    def title(axis, text):
        axis.text(0.04, .96, text, transform=axis.transAxes, ha="left", va="top", fontsize=10,
                  bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.94, "pad": 2}, zorder=20)

    # (a) pointings
    ra0, dec0 = 221.25, -59.8
    for name, marker, label in [("aca_7m", "o", "7 m"), ("configuration_c5", "+", "12 m")]:
        pos = np.array([f["PHASE_DIR"] for f in metadata[name]["subtables"]["FIELD"]]).squeeze()
        ra, dec = np.rad2deg(pos[:, 0]) % 360, np.rad2deg(pos[:, 1])
        a.scatter((ra-ra0)*np.cos(np.deg2rad(dec0))*60, (dec-dec0)*60, s=7,
                  color=COLORS[name] if name == "aca_7m" else COLORS["configuration_c2"], marker=marker,
                  linewidths=0.7, alpha=0.75, label=label)
    a.set(xlabel="RA offset (arcmin)", ylabel="Dec offset (arcmin)")
    a.set_aspect("equal", adjustable="datalim")
    a.invert_xaxis()
    a.margins(0.12)
    a.set_xticks([5, 0, -5])
    title(a, "(a)")
    a.legend(loc="upper left", bbox_to_anchor=(0.03, 0.88), fontsize=10, labelspacing=0.35, handletextpad=0.4)

    # (c) amplitudes and visibility counts
    short = {"aca_7m": "7 m", "configuration_c2": "C2", "configuration_c5": "C5"}
    values = list(tp_amp) + ([] if tp_fixed is None else list(tp_fixed))
    values += shown_curve_values(image_curves, shown_curves)
    for name, marker in [("aca_7m", "o"), ("configuration_c2", "s"), ("configuration_c5", "^")]:
        profile = [r for r in array_profiles if r["array"] == name]
        radius = np.array([r["uv_klambda"] for r in profile])
        amplitude = np.array([r["rms_mjy"] for r in profile])
        detected = np.array([r["detected"] for r in profile])
        upper = np.array([r["upper_rms_mjy"] for r in profile])
        c.plot(radius, np.where(detected, amplitude, np.nan), color=COLORS[name],
               linewidth=.9, marker=marker, markersize=1.5, label=short[name])
        if (~detected).any():
            c.errorbar(radius[~detected], upper[~detected], yerr=upper[~detected]*0.25,
                       uplims=True, fmt="none", color=COLORS[name], elinewidth=0.8, capsize=1.2, markersize=1.5)
        values.extend(amplitude[detected])
        values.extend(upper[~detected])
    positive_tp = tp_u >= 1
    c.plot(tp_u[positive_tp], tp_amp[positive_tp], color="0.45", linestyle="-.", linewidth=.9, label="TP")
    if tp_fixed is not None:
        c.plot(tp_u[positive_tp], tp_fixed[positive_tp], color="0.1", linewidth=1.1, label="TP, beam-corrected")
    if image_curves:
        _plot_image_curves(c, image_curves, shown_curves, curve_labels)
    c.set(xlabel=r"uv distance (k$\lambda$)", ylabel="Amplitude (mJy)", xlim=(1, 1000))
    c.set_xscale("log")
    c.set_xticks([1, 10, 100, 1000], labels=["1", "10", "100", "1000"])
    c.set_yscale("log")
    values = np.asarray(values)
    positive = values[np.isfinite(values) & (values > 0)]
    count_top = .40
    c.set_ylim(positive.min()*.001, positive.max()*8.0)
    sampling_top = c.get_ylim()[0]*(c.get_ylim()[1]/c.get_ylim()[0])**count_top
    c.yaxis.set_major_locator(FixedLocator([v for v in [10, 100, 1000, 10000] if v > sampling_top]))
    c.yaxis.set_minor_locator(FixedLocator([v for v in c.yaxis.get_minorticklocs() if v > sampling_top]))
    overlap = coverage["overlap_klambda"]
    maximum_count = max(counts.max() for _, _, counts in observation_density)
    track_y = np.linspace(.24, .04, len(observation_density))
    for (name, _, counts), y in zip(observation_density, track_y):
        occupied = counts > 0
        c.scatter(density_centres[occupied], np.full(occupied.sum(), y),
                  s=4*counts[occupied]/maximum_count, color=COLORS[name], linewidths=0,
                  transform=c.get_xaxis_transform(), zorder=3)
    # The stack must stay below the "Visibility counts" text at 0.30.
    # One column of labels, level with their tracks but kept a text height apart (bottom to top).
    label_gap, label_y = .085, {}
    for name in reversed(list(observation_sources)):
        own = track_y[[i for i, (array, _, _) in enumerate(observation_density) if array == name]].mean()
        label_y[name] = max(own, max(label_y.values(), default=-1) + label_gap)
    for name in observation_sources if tp_fixed is not None else ():
        c.text(0.025, label_y[name], short[name], color=COLORS[name], fontsize=10,
               va="center", transform=c.transAxes)
    c.plot([0, 1], [count_top, count_top], color="0.8", linewidth=0.5, transform=c.transAxes)
    c.text(0.025, .30, "Visibility counts", transform=c.transAxes, fontsize=10)
    c.get_xticklabels()[-1].set_horizontalalignment("right")
    title(c, "(c)")
    if tp_fixed is not None:
        # Array colours are keyed in panel (b) and labelled on the count tracks;
        # the two TP curves are labelled directly.
        end = np.flatnonzero(positive_tp)[-1]
        c.text(tp_u[end]*1.03, tp_amp[end]*.5, "TP", color="0.45", fontsize=10, ha="left", va="center")
        c.annotate("TP, beam-corrected", xy=(2.5, np.interp(2.5, tp_u, tp_fixed)), xytext=(3.4, tp_fixed.max()*1.7),
                   fontsize=10, ha="left", va="bottom", color="0.1",
                   arrowprops=dict(arrowstyle="-", color="0.1", linewidth=.6, shrinkA=1, shrinkB=1))
    else:
        c.legend(loc="upper right", bbox_to_anchor=(0.99, 0.98), fontsize=10, labelspacing=0.25,
                 ncol=2, columnspacing=.4, borderpad=.2, handlelength=1.1,
                 handletextpad=0.4, frameon=True, facecolor="white", framealpha=1, edgecolor="none")

    # (b) uv coverage with the short-baseline inset (7 m and the TP response)
    zoom = 25.0
    limit = 1.10*max(np.abs(np.r_[u, v]).max() for u, v in coords.values())
    wavelength = C_METRES_PER_SECOND/(config["visibility_profile"]["reference_frequency_ghz"]*1e9)
    tp_half_metres = tp_meta["half_response_klambda"]*1000*wavelength
    tp_outer_metres = tp_meta["maximum_plotted_klambda"]*1000*wavelength
    inset.add_patch(Circle((0, 0), tp_outer_metres, facecolor="0.80", edgecolor="0.55", linewidth=0.6, zorder=0))
    inset.add_patch(Circle((0, 0), tp_half_metres, facecolor="0.65", edgecolor="0.4", linewidth=0.6, linestyle=":", zorder=1))
    for axis, extent in [(b, limit), (inset, zoom)]:
        for name, (u, v) in coords.items():
            if axis is inset and name != "aca_7m":
                continue
            use = (np.abs(u) <= extent) & (np.abs(v) <= extent)
            axis.scatter(np.r_[u[use], -u[use]], np.r_[v[use], -v[use]], s=.25 if axis is b else .75,
                         color=COLORS[name], alpha=0.25 if name != "aca_7m" else 0.60,
                         linewidths=0, rasterized=True, zorder=2)
        axis.set(xlim=(-extent, extent), ylim=(-extent, extent))
        axis.set_aspect("equal", adjustable="box")
    b.add_patch(Rectangle((-zoom, -zoom), 2*zoom, 2*zoom, facecolor="none", edgecolor="0.2", linewidth=0.7, zorder=3))
    inset.scatter([0], [0], marker="+", s=28, linewidths=1.1, color="0.15", zorder=5)
    handles = [Line2D([], [], marker="o", linestyle="none", color=COLORS[k], markersize=3, label=short[k]) for k in COLORS]
    tp_handle = Line2D([], [], marker="+", linestyle="none", color="0.3", markersize=6, label="TP")
    b.legend(handles=handles+[tp_handle], loc="lower left", ncol=2, fontsize=10,
             columnspacing=.35, handlelength=.7, handletextpad=.25, labelspacing=.15,
             frameon=False)
    b.set(xlabel="u (km)", ylabel="v (km)")
    b.set_xticks([-1000, 0, 1000], labels=["−1", "0", "1"])
    b.set_yticks([-1000, 0, 1000], labels=["−1", "0", "1"])
    inset.set(xlabel="", ylabel="", xticks=[], yticks=[])
    inset.minorticks_off()
    for origin, destination in [((-zoom, zoom), (0, 0)), ((zoom, -zoom), (1, 0))]:
        b.add_artist(ConnectionPatch(xyA=origin, coordsA=b.transData,
            xyB=destination, coordsB=inset.transAxes, color="black", linewidth=.65,
            linestyle="-", zorder=7, clip_on=False))
    title(b, "(b)")
    save_figure(fig, outputs, dpi=int(settings["dpi"]), bbox_inches=None, exclude_axes=[inset])
    plt.close(fig)
    sources = [sampling/"amplitudes.csv", sampling/"separate_array_amplitudes.csv", sampling/"summary.json", sampling/"cell_statistics.npz", test/"summary.json",
               *observation_sources.values(),
               *[root/"uv_coordinates"/f"{n}.npz" for n in COLORS], Path(paths["inputs"]["tp_dust_ff"]),
               *alignment_sources]
    manifest = {"date": datetime.now(timezone.utc).date().isoformat(), "status": "Paper II Fig. 2",
                "tp": tp_meta, "tp_half_response_metres": tp_half_metres,
                "tp_ten_percent_response_metres": tp_outer_metres,
                "panel_c_statistic": "Separate ACA/C2/C5 RMS from independent-hand cross-power in the same 19 pointing pairs; common occupied cells across arrays present in each radial bin where possible, own occupied cells otherwise (CSV matched_arrays flag). ACA model translation targets actual 12m pointing PB. TP has this same PB applied before FFT, with its native restoring beam retained. Non-detections have conditional two-standard-error limits.",
                "image_curves": image_meta,
                "figure_width_mm": settings["width_inches"]*25.4, "base_font_pt": 10,
                "amplitude_source": amplitude_source,
                "amplitude_ranges_klambda": coverage["ranges_klambda"], "amplitude_overlap_klambda": overlap,
                "panel_c_xscale": "log; 1-1000 klambda; zero spacing remains in the short-baseline view and the source CSV",
                "panel_order": ["pointings", "uv coverage with short-baseline inset", "amplitudes and visibility counts"],
                "panel_labels": ["(a)", "(b)", "(c)"],
                "inset_limit_metres": zoom, "main_uv_display_unit": "km", "inset_uv_display_unit": "m",
                "inset_arrays": ["aca_7m", "TP model response"],
                "inset_background": "white; flush with upper-right parent edges; two solid black zoom connectors",
                "inset_ticks": False, "inset_tick_labels": False,
                "panel_c_marker_size_pt": 1.5,
                "alignment_exclusion": "Upper-right inset in b is intentionally a smaller zoom, not an equal-span main panel",
                "panel_c_native_aca_curve": False, "panel_c_shading": False,
                "panel_c_sampling": {
                    "source_csv": str(density_table), "log_bins": len(density_centres),
                    "marker_area": "4 pt^2 times row count / global maximum bin count",
                    "maximum_bin_count": int(maximum_count),
                    "selection": "Same 19 pointing pairs and broad windows; all valid cached rows before common-cell selection; no imaging weights or thinning. Vertical track position is an observation index, not amplitude.",
                    "observations": [{"array": name, "observation_id": obs, "rows": int(n.sum())}
                                     for name, obs, n in observation_density]},
                "matched_pair_markers": False,
                "sources_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    outputs[0].with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return outputs
