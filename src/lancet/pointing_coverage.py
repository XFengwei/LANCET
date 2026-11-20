"""Mosaic pointing centres compared with the MIRIAD gain and the CASA primary-beam response.

The MIRIAD mosaic gain is 1 where full primary-beam correction keeps the noise at or below
that of the noisiest pointing, and falls below 1 outside (``mossen``; ``mostab.for``). The
CASA ``.pb`` is the relative sensitivity of the joint mosaic. This module places the
pointings of each array on both and tabulates the gain against the pb.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

C_M_S = 2.99792458e8


def read_field_centres(measurement_set) -> tuple[np.ndarray, np.ndarray]:
    """RA, Dec (deg) of every row of the FIELD table's PHASE_DIR."""
    from casaconfig import config as casa_config

    casa_config.nologfile = True
    casa_config.logfile = None
    casa_config.data_auto_update = False
    casa_config.measures_auto_update = False
    from casatools import table

    handle = table()
    handle.open(f"{measurement_set}/FIELD", nomodify=True)
    try:
        phase = np.asarray(handle.getcol("PHASE_DIR"))  # (2, 1, nfield), radians
    finally:
        handle.close()
    ra = np.degrees(phase[0, 0]) % 360.0
    return ra, np.degrees(phase[1, 0])


def unique_positions(ra_deg, dec_deg, tolerance_arcsec: float) -> tuple[np.ndarray, np.ndarray]:
    """Positions with repeats closer than ``tolerance_arcsec`` removed (first kept)."""
    ra, dec = np.asarray(ra_deg, float), np.asarray(dec_deg, float)
    keep: list[int] = []
    for k in range(ra.size):
        if keep:
            dx = (ra[keep] - ra[k]) * np.cos(np.radians(dec[k]))
            if np.min(np.hypot(dx, dec[keep] - dec[k])) * 3600.0 <= tolerance_arcsec:
                continue
        keep.append(k)
    return ra[keep], dec[keep]


def primary_beam_fwhm_arcsec(dish_m: float, frequency_ghz: float, factor: float = 1.13) -> float:
    """Primary-beam FWHM = factor x lambda / D, in arcsec."""
    return float(np.degrees(factor * C_M_S / (frequency_ghz * 1e9) / dish_m) * 3600.0)


def gain_against_pb(gain, pb, edges: Sequence[float]) -> list[dict]:
    """Percentiles (5, 50, 95) of the gain in bins of pb; pixels with non-finite values are ignored."""
    gain, pb = np.asarray(gain, float), np.asarray(pb, float)
    ok = np.isfinite(gain) & np.isfinite(pb)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = ok & (pb >= lo) & (pb < hi)
        row = {"pb_low": float(lo), "pb_high": float(hi), "n": int(sel.sum())}
        if sel.any():
            p5, p50, p95 = np.percentile(gain[sel], [5, 50, 95])
            row.update(gain_p5=float(p5), gain_median=float(p50), gain_p95=float(p95))
        rows.append(row)
    return rows


def pb_at_gain_edge(gain, pb, low: float = 0.98, high: float = 0.999) -> dict:
    """Percentiles of the pb where the gain is just below 1, i.e. along the gain = 1 edge."""
    gain, pb = np.asarray(gain, float), np.asarray(pb, float)
    sel = np.isfinite(gain) & np.isfinite(pb) & (gain >= low) & (gain < high)
    p5, p50, p95 = np.percentile(pb[sel], [5, 50, 95])
    return {"gain_low": low, "gain_high": high, "n": int(sel.sum()),
            "pb_p5": float(p5), "pb_median": float(p50), "pb_p95": float(p95)}


def plot_coverage(fig, gain, gain_wcs, pb, pb_wcs, arrays, *, pb_levels, region, title):
    """One panel: gain greyscale, gain = 1 and pb contours, pointings with half-power circles.

    ``arrays``: sequence of dicts with ``label``, ``ra``, ``dec``, ``fwhm_arcsec``, ``marker``,
    ``size``, ``colour``. ``region``: (SkyCoord centre, (width, height) arcsec).
    """
    from matplotlib.lines import Line2D
    from matplotlib.patches import Ellipse

    ax = fig.add_subplot(projection=gain_wcs)
    image = ax.imshow(np.nan_to_num(gain), origin="lower", cmap="Greys", vmin=0, vmax=1, alpha=0.4)
    colourbar = fig.colorbar(image, ax=ax, fraction=0.04, pad=0.02)
    colourbar.set_label("MIRIAD gain")
    ax.contour(np.nan_to_num(gain), levels=[0.999], colors="black", linewidths=1.6)
    pb_colours = ["#1f77b4", "#ff7f0e", "#d62728"][: len(pb_levels)]
    levels = sorted(pb_levels)
    ax.contour(np.nan_to_num(pb), levels=levels, colors=pb_colours, linewidths=1.1, linestyles="--",
               transform=ax.get_transform(pb_wcs))
    sky = ax.get_transform("icrs")
    handles = [Line2D([], [], color="black", lw=1.6, label="MIRIAD gain = 1")]
    handles += [Line2D([], [], color=c, ls="--", label=f"CASA pb = {lv:g}")
                for lv, c in sorted(zip(levels, pb_colours), reverse=True)]
    for spec in arrays:
        ra, dec = np.asarray(spec["ra"]), np.asarray(spec["dec"])
        ax.scatter(ra, dec, transform=sky, marker=spec["marker"], s=spec["size"], c=spec["colour"], lw=1.2)
        width = spec["fwhm_arcsec"] / 3600.0
        for r, d in zip(ra, dec):
            ax.add_patch(Ellipse((r, d), width / np.cos(np.radians(d)), width, transform=sky,
                                 fill=False, ec=spec["colour"], lw=0.4, alpha=0.6))
        handles.append(Line2D([], [], color=spec["colour"], marker=spec["marker"], ls="", ms=8,
                              label=f"{spec['label']} ({ra.size} in mosaic), FWHM {spec['fwhm_arcsec']:.1f}\""))
    centre, (width_arcsec, height_arcsec) = region
    x0, y0 = gain_wcs.world_to_pixel(centre)
    pixel_arcsec = abs(gain_wcs.wcs.cdelt[0]) * 3600.0
    ax.set_xlim(x0 - width_arcsec / 2 / pixel_arcsec, x0 + width_arcsec / 2 / pixel_arcsec)
    ax.set_ylim(y0 - height_arcsec / 2 / pixel_arcsec, y0 + height_arcsec / 2 / pixel_arcsec)
    ax.legend(handles=handles, loc="upper left", fontsize=8, framealpha=0.9)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("RA (ICRS)")
    ax.set_ylabel("Dec (ICRS)")
    return ax
