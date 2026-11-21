"""Figures for the linear continuum combination task."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import numpy as np

COLORS = {"interferometer": "#007A9E", "tp": "#D87821", "combined": "#3B3B3B", "band": "#F3D9BF",
          "guide": "#8C8C8C"}
ARCSEC_PER_KLAMBDA = 206.264806   # angular scale (arcsec) = this / uv distance (klambda)


def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.linewidth": 0.7,
                         "pdf.fonttype": 42, "legend.frameon": False})
    return plt


def _thin(data, max_pixels=900):
    step = max(1, int(np.ceil(max(np.shape(data)) / max_pixels)))
    return np.asarray(data)[::step, ::step]


def _save(fig, stem: Path) -> None:
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(stem.with_suffix(f".{ext}"), dpi=200)


def _image(axis, data, title, coverage=None, **kwargs):
    im = axis.imshow(_thin(data), origin="lower", interpolation="nearest", **kwargs)
    if coverage is not None:
        axis.contour(_thin(coverage).astype(float), [0.5], colors="white", linewidths=0.5)
    axis.set_title(title, fontsize=9)
    axis.set_xticks([])
    axis.set_yticks([])
    return im


def _uv_axis(axis, ticks, scale_ticks):
    """Log uv-distance axis with plain tick labels and the angular scale on top."""
    from matplotlib.ticker import NullFormatter, ScalarFormatter
    def convert(x):
        with np.errstate(divide="ignore"):
            return ARCSEC_PER_KLAMBDA / np.asarray(x, dtype=float)

    axis.set_xscale("log")
    axis.set_xticks(ticks)
    axis.xaxis.set_major_formatter(ScalarFormatter())
    axis.xaxis.set_minor_formatter(NullFormatter())
    axis.set_xlabel("uv distance (kλ)")
    top = axis.secondary_xaxis("top", functions=(convert, convert))
    top.set_xticks(scale_ticks)
    top.xaxis.set_major_formatter(ScalarFormatter())
    top.xaxis.set_minor_formatter(NullFormatter())
    top.set_xlabel("angular scale (″)")


def plot_tp_factor(stem, centres, power_interferometer, power_tp, band_klambda, adopted, measured) -> None:
    """How the TP factor is set: power spectra at the TP resolution before and after the shift.

    Left: 12m+7m smoothed to the TP beam, the TP and the TP times the adopted
    factor. Right: TP / 12m+7m amplitude ratio, before and after.
    """
    plt = _style()
    centres = np.asarray(centres, dtype=float)
    p_int, p_tp = np.asarray(power_interferometer, dtype=float), np.asarray(power_tp, dtype=float)
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.3), layout="constrained")
    ax = axes[0]
    ax.axvspan(*band_klambda, color=COLORS["band"], lw=0, label="matching range")
    ax.loglog(centres, p_int, color=COLORS["interferometer"], lw=1.8, label="12m+7m at the TP resolution")
    ax.loglog(centres, p_tp, color=COLORS["tp"], lw=1.4, ls=(0, (4, 2)), label="TP, original")
    ax.loglog(centres, p_tp * adopted**2, color=COLORS["tp"], lw=1.8,
              label=f"TP × {adopted:g} (power × {adopted**2:.2f})")
    ax.set(ylabel="power (MJy² sr⁻² arcsec²)", title="Power spectra in the common region")
    ax.legend(fontsize=7.5, loc="lower left")
    ax = axes[1]
    ratio = np.sqrt(p_tp / p_int)
    ax.axvspan(*band_klambda, color=COLORS["band"], lw=0)
    ax.plot(centres, ratio, color=COLORS["tp"], lw=1.4, ls=(0, (4, 2)), marker="o", ms=2.5,
            label="TP original / 12m+7m")
    ax.plot(centres, ratio * adopted, color=COLORS["tp"], lw=1.8, marker="o", ms=2.5,
            label=f"TP × {adopted:g} / 12m+7m")
    ax.axhline(1.0, color="black", lw=0.6)
    ax.set_yscale("log")
    ax.set_ylim(0.5, 6.0)
    from matplotlib.ticker import NullFormatter, ScalarFormatter
    ax.set_yticks([0.5, 0.7, 1, 1.5, 2, 3, 5])
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set(ylabel="amplitude ratio  √(P_TP / P_12m+7m)",
           title=f"Lowest point in the range needs × {measured:.2f}; adopted × {adopted:g}")
    ax.legend(fontsize=7.5, loc="lower left")
    for ax in axes:
        _uv_axis(ax, [1, 2, 3, 5, 7, 10, 20], [10, 20, 30, 50, 100, 200])
        ax.set_xlim(centres.min() * 0.95, centres.max() * 1.05)
    _save(fig, stem)
    plt.close(fig)


def plot_comparison(stem, interferometer, tp, combined, coverage, method_label, tp_label) -> None:
    """Before and after: the 12m+7m image, the TP and the combined image (Jy per 12m+7m beam)."""
    plt = _style()
    from matplotlib.colors import AsinhNorm
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.6), layout="constrained")
    maps = {"12m+7m (before)": interferometer, f"{tp_label} (before)": tp, f"{method_label} (after)": combined}
    finite = np.concatenate([np.asarray(m)[np.isfinite(m)].ravel()[::97] for m in maps.values()])
    vmax = float(np.percentile(finite, 99.8))
    norm = AsinhNorm(linear_width=vmax / 50, vmin=0.0, vmax=vmax)
    for axis, (name, data) in zip(axes, maps.items()):
        im = _image(axis, data, name, coverage, cmap="magma", norm=norm)
    bar = fig.colorbar(im, ax=axes.tolist(), shrink=0.9, pad=0.01, label="Jy per 12m+7m beam")
    bar.set_ticks([0.0] + [10.0**e for e in range(-6, 3) if vmax / 10 <= 10.0**e <= vmax])
    _save(fig, stem)
    plt.close(fig)


def plot_zero_spacing(stem, tp, combined_at_tp_beam, mosaic, tp_only, method_label, tp_label, stats: Mapping) -> None:
    """Zero-spacing check on the TP grid: the TP, the combined image at the TP beam, their ratio,
    and the pixel-by-pixel comparison. ``mosaic`` and ``tp_only`` are the two regions."""
    plt = _style()
    fig, axes = plt.subplots(1, 4, figsize=(15.5, 3.9), layout="constrained",
                             gridspec_kw={"width_ratios": [1, 1, 1, 1.05]})
    valid = np.isfinite(combined_at_tp_beam) & np.isfinite(tp)      # includes the mosaic edge (low gain)
    ys, xs = np.nonzero(valid)
    crop = np.s_[ys.min():ys.max() + 1, xs.min():xs.max() + 1]
    vmax = float(np.nanpercentile(tp[valid], 99.8))
    for axis, data, title in ((axes[0], tp, f"{tp_label} (18.1″)"),
                              (axes[1], combined_at_tp_beam, f"{method_label} smoothed to the TP beam")):
        im = axis.imshow(np.where(valid, data, np.nan)[crop], origin="lower", cmap="magma", vmin=0, vmax=vmax,
                         interpolation="nearest")
        axis.set_title(title, fontsize=9)
    fig.colorbar(im, ax=axes[:2].tolist(), shrink=0.85, pad=0.01, label="MJy/sr")
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(valid, combined_at_tp_beam / tp, np.nan)
    im = axes[2].imshow(ratio[crop], origin="lower", cmap="RdBu_r", vmin=0.9, vmax=1.1, interpolation="nearest")
    axes[2].set_title("smoothed combined / TP", fontsize=9)
    fig.colorbar(im, ax=axes[2], shrink=0.85, pad=0.02)
    for axis in axes[:3]:
        axis.contour(mosaic[crop].astype(float), [0.5], colors="black" if axis is axes[2] else "white",
                     linewidths=0.6)
        axis.set_xticks([])
        axis.set_yticks([])
    ax = axes[3]
    for region, color, label in ((tp_only & valid, COLORS["guide"], "outside the mosaic (TP only)"),
                                 (mosaic & valid, COLORS["interferometer"], "inside the mosaic")):
        ax.scatter(tp[region], combined_at_tp_beam[region], s=3, color=color, label=label, rasterized=True)
    limit = float(np.nanmax(tp[valid])) * 1.05
    ax.plot([0, limit], [0, limit], color="black", lw=0.7, label="1 : 1")
    ax.set(xlim=(0, limit), ylim=(0, limit), xlabel=f"{tp_label} (MJy/sr)",
           ylabel=f"{method_label} at the TP beam (MJy/sr)",
           title="flux ratio: " + ", ".join(
               f"{name} {stats[key]['flux_ratio']:.3f}" if stats[key].get("n") else f"{name} none"
               for name, key in (("mosaic", "mosaic"), ("TP only", "tp_only"))))
    ax.title.set_fontsize(9)
    ax.legend(fontsize=7, loc="upper left")
    _save(fig, stem)
    plt.close(fig)


def plot_power_spectrum(stem, centres, powers: Mapping[str, np.ndarray], ratios: Mapping[str, np.ndarray],
                        tp_limit_klambda, guides_klambda, method_label, tp_label) -> None:
    """Power spectra of the inputs and the combined image, and the combined image's amplitude
    relative to each input. ``powers``/``ratios`` have keys ``interferometer``, ``tp`` (and
    ``combined`` for powers). The ratio to the TP is drawn up to ``tp_limit_klambda``; vertical
    guides mark ``guides_klambda`` (the J-comb thresholds)."""
    plt = _style()
    centres = np.asarray(centres, dtype=float)
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.3), layout="constrained")
    ax = axes[0]
    ax.loglog(centres, powers["interferometer"], color=COLORS["interferometer"], lw=1.4, label="12m+7m")
    ax.loglog(centres, powers["tp"], color=COLORS["tp"], lw=1.4, label=f"{tp_label} (18.1″ beam)")
    ax.loglog(centres, powers["combined"], color=COLORS["combined"], lw=1.8, ls=(0, (4, 2)), label=method_label)
    ax.set(ylabel="power (MJy² sr⁻² arcsec²)", title="Power spectra inside the mosaic")
    ax.legend(fontsize=7.5, loc="lower left")
    _uv_axis(ax, [1, 3, 10, 30, 100, 300], [1, 3, 10, 30, 100])
    ax = axes[1]
    low = centres <= tp_limit_klambda
    ax.plot(centres[low], np.array(ratios["tp"], dtype=float)[low], color=COLORS["tp"], lw=1.6, marker="o", ms=3,
            label=f"{method_label} / {tp_label} at the 12m+7m beam")
    ax.plot(centres, ratios["interferometer"], color=COLORS["interferometer"], lw=1.6, marker="o", ms=3,
            label=f"{method_label} / 12m+7m")
    ax.axhline(1.0, color="black", lw=0.6)
    ax.set_ylim(0.0, 2.5)
    ax.set(ylabel="amplitude ratio", title="Combined image relative to each input")
    ax.legend(fontsize=7.5, loc="upper right")
    _uv_axis(ax, [1, 2, 3, 5, 7, 10, 20, 50], [5, 10, 20, 30, 50, 100, 200])
    ax.set_xlim(centres.min() * 0.9, 60.0)
    for ax in axes:
        for k in guides_klambda:
            ax.axvline(k, color=COLORS["guide"], lw=0.7, ls="--")
    _save(fig, stem)
    plt.close(fig)


def plot_planck_maps(stem, planck, model, coverage, coverage_min, peak, radius_pixels, pixel_arcsec, label) -> None:
    """Planck at the model frequency, the model at the Planck beam, and their ratio."""
    plt = _style()
    from matplotlib.patches import Circle
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.9), layout="constrained")
    half = planck.shape[0] * pixel_arcsec / 120.0
    extent = [half, -half, -half, half]            # arcmin offsets, east to the left
    shown = coverage >= 0.5
    vmax = float(np.nanmax(planck))
    for axis, data, title in ((axes[0], planck, "Planck at the TP frequency, CO removed"),
                              (axes[1], np.where(shown, model, np.nan), f"{label} at the Planck beam")):
        im = axis.imshow(data, origin="lower", extent=extent, cmap="magma", vmin=0, vmax=vmax)
        axis.set_title(title, fontsize=9)
    fig.colorbar(im, ax=axes[:2].tolist(), shrink=0.85, pad=0.01, label="MJy/sr")
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(shown, model / planck, np.nan)
    im = axes[2].imshow(ratio, origin="lower", extent=extent, cmap="RdBu_r", vmin=0.5, vmax=1.5)
    axes[2].set_title(f"{label} / Planck (no background removed)", fontsize=9)
    fig.colorbar(im, ax=axes[2], shrink=0.85, pad=0.02)
    n = planck.shape[0]
    cx = half - (peak[1] + 0.5) * pixel_arcsec / 60.0
    cy = -half + (peak[0] + 0.5) * pixel_arcsec / 60.0
    for axis in axes:
        axis.contour(np.linspace(half, -half, n), np.linspace(-half, half, n), coverage, [coverage_min],
                     colors="white" if axis is not axes[2] else "black", linewidths=0.7)
        axis.add_patch(Circle((cx, cy), radius_pixels * pixel_arcsec / 60.0, fill=False, color="#00C2A8", lw=1.3))
        axis.set_xlabel("offset east (arcmin)")
    axes[0].set_ylabel("offset north (arcmin)")
    _save(fig, stem)
    plt.close(fig)


def plot_planck_scatter(stem, planck, model, selected, region, reference, fit, bright, label) -> None:
    """Model against Planck: all selected pixels, the brightest beam and the reference pixels."""
    plt = _style()
    fig, ax = plt.subplots(figsize=(5.6, 4.6), layout="constrained")
    ax.scatter(planck[selected], model[selected], s=6, color=COLORS["guide"], label="coverage-selected pixels")
    ax.scatter(planck[reference], model[reference], s=9, color=COLORS["interferometer"], label="reference (faintest)")
    ax.scatter(planck[region], model[region], s=9, color=COLORS["tp"], label="brightest beam")
    x = np.array([0.0, float(np.nanmax(planck[selected])) * 1.05])
    ax.plot(x, x, color="black", lw=0.7, label="1 : 1")
    ax.plot(x, fit["gain"] * x + fit["offset_mjysr"], color=COLORS["combined"], lw=1.2, ls=(0, (4, 2)),
            label=f"fit: {fit['gain']:.2f} × Planck − {-fit['offset_mjysr']:.1f}")
    direct = bright["ratio_model_over_planck"]["region_mean"]["direct"]
    ax.plot(x, direct * x, color=COLORS["tp"], lw=1.0, ls=":", label=f"brightest beam, through origin: × {direct:.2f}")
    ax.set(xlabel="Planck at the TP frequency (MJy/sr)", ylabel=f"{label} at the Planck beam (MJy/sr)",
           xlim=(0, x[1]), ylim=(0, x[1]))
    ax.legend(fontsize=7, loc="upper left")
    _save(fig, stem)
    plt.close(fig)


def plot_box_spectrum(stem, interferometer, tp, interferometer_tp_res, window, half_size_arcsec, centres,
                      spectra: Mapping, factors: Mapping, band_klambda, tp_limit_klambda) -> None:
    """One box. Top: the 12m+7m image, the TP and the 12m+7m image smoothed to the TP beam (both
    MJy/sr, one scale). Bottom: their power spectra, and the TP factor each uv bin asks for.

    ``spectra``: ``interferometer_tp_res``, ``tp``, ``interferometer``, ``tp_high_res``.
    ``factors``: per-bin ``psd_tp_res``, ``amplitude_high_res``, ``cross_high_res`` and the
    aligned value ``aligned``.
    """
    plt = _style()
    centres = np.asarray(centres, dtype=float)
    half_x, half_y = half_size_arcsec                 # half width (RA) and half height (Dec)
    fig = plt.figure(figsize=(12.0, 8.6), layout="constrained")
    top, bottom = fig.subfigures(2, 1, height_ratios=[0.55 + 0.45 * min(half_y / half_x, 1.0), 1.05])
    image_axes = top.subplots(1, 3)
    axes = bottom.subplots(1, 2)
    extent = [half_x, -half_x, -half_y, half_y]
    contour_x = np.linspace(half_x, -half_x, window.shape[1])
    contour_y = np.linspace(-half_y, half_y, window.shape[0])
    smooth_max = float(np.nanpercentile(np.concatenate([np.ravel(tp), np.ravel(interferometer_tp_res)]), 99.5))
    panels = ((image_axes[0], interferometer, "12m+7m (0.71″ × 0.61″)", float(np.nanpercentile(interferometer, 99.5))),
              (image_axes[1], tp, "TP (18.1″)", smooth_max),
              (image_axes[2], interferometer_tp_res, "12m+7m smoothed to the TP beam", smooth_max))
    for axis, data, title, vmax in panels:
        im = axis.imshow(data, origin="lower", extent=extent, cmap="magma", vmin=0, vmax=vmax)
        axis.contour(contour_x, contour_y, window, [0.5], colors="white", linewidths=0.7)
        axis.set_title(title, fontsize=9)
        axis.set_xlabel("offset east (″)")
        if axis is image_axes[0]:
            top.colorbar(im, ax=axis, shrink=0.8, pad=0.02, label="MJy/sr")
    top.colorbar(im, ax=image_axes[1:].tolist(), shrink=0.8, pad=0.02, label="MJy/sr (TP and smoothed 12m+7m)")
    image_axes[0].set_ylabel("offset north (″)")
    ax = axes[0]
    aligned = factors["aligned"]
    ax.axvspan(*band_klambda, color=COLORS["band"], lw=0, label="matching range")
    ax.loglog(centres, spectra["interferometer"], color=COLORS["interferometer"], lw=1.0, ls=(0, (1, 1.5)),
              label="12m+7m, native")
    ax.loglog(centres, spectra["interferometer_tp_res"], color=COLORS["interferometer"], lw=1.8,
              label="12m+7m at the TP resolution")
    ax.loglog(centres, spectra["tp"], color=COLORS["tp"], lw=1.4, ls=(0, (4, 2)), label="TP")
    ax.loglog(centres, np.asarray(spectra["tp"]) * aligned**2, color=COLORS["tp"], lw=1.8,
              label=f"TP × {aligned:.2f} (aligned)")
    ax.set(ylabel="power (MJy² sr⁻² arcsec²)", title="Power spectra in the box")
    ax.legend(fontsize=7, loc="lower left")
    ax = axes[1]
    low = centres <= tp_limit_klambda
    ax.axvspan(*band_klambda, color=COLORS["band"], lw=0)
    ax.plot(centres, factors["psd_tp_res"], color=COLORS["tp"], lw=1.6, marker="o", ms=3,
            label="power spectra at the TP resolution")
    ax.plot(centres[low], np.array(factors["amplitude_high_res"], dtype=float)[low], color=COLORS["interferometer"],
            lw=1.6, marker="o", ms=3, label="amplitudes at the 12m+7m resolution")
    ax.plot(centres[low], np.array(factors["cross_high_res"], dtype=float)[low], color=COLORS["combined"], lw=1.6,
            marker="s", ms=3, label="cross-spectrum at the 12m+7m resolution")
    ax.axhline(1.0, color="black", lw=0.6)
    ax.set_ylim(0.0, 3.0)
    ax.set(ylabel="TP factor that matches 12m+7m", title="Factor asked for at each scale")
    ax.legend(fontsize=7, loc="upper left")
    for ax in axes:
        _uv_axis(ax, [2, 3, 5, 7, 10, 20], [10, 20, 30, 50, 100])
        ax.set_xlim(centres.min() * 0.9, centres.max() * 1.05)
    _save(fig, stem)
    plt.close(fig)
