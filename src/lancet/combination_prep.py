"""Preparation of images for uv-plane combination: I/O, units, gain, regridding,
edge taper and padding."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy import fft, ndimage

from lancet.combination import (
    combine_plane,
    feather_weights,
    jcomb_weights,
    uv_distance_klambda,
)

ARCSEC2_SR = (np.pi / 648000.0) ** 2


def read_plane(path: str | Path):
    """Read a Jy/beam image as a float32 2-D plane with its celestial WCS."""
    with fits.open(path, memmap=False) as hdul:
        header = hdul[0].header.copy()
        data = np.squeeze(np.asarray(hdul[0].data))
    if data.ndim != 2:
        raise ValueError(f"{path}: expected one image plane, got shape {data.shape}")
    bunit = str(header.get("BUNIT", "")).replace(" ", "").lower()
    if bunit != "jy/beam":
        raise ValueError(f"{path}: BUNIT must be Jy/beam, got {header.get('BUNIT')!r}")
    if "BMAJ" not in header or "BMIN" not in header:
        raise ValueError(f"{path}: BMAJ/BMIN missing")
    return data.astype(np.float32), header, WCS(header).celestial


def beam_parameters(header) -> tuple[float, float, float]:
    """(bmaj arcsec, bmin arcsec, bpa deg) from a FITS header."""
    return float(header["BMAJ"]) * 3600.0, float(header["BMIN"]) * 3600.0, float(header.get("BPA", 0.0))


def _beam_sr(bmaj_arcsec, bmin_arcsec):
    return np.pi * bmaj_arcsec * bmin_arcsec / (4.0 * np.log(2.0)) * ARCSEC2_SR


def jybeam_to_mjysr(data, bmaj_arcsec, bmin_arcsec):
    return np.asarray(data) / _beam_sr(bmaj_arcsec, bmin_arcsec) / 1e6


def mjysr_to_jybeam(data, bmaj_arcsec, bmin_arcsec):
    return np.asarray(data) * _beam_sr(bmaj_arcsec, bmin_arcsec) * 1e6


def derive_gain(clean, attenuated, maximum=1.01, weak_fraction=1e-9):
    """Mosaic gain = attenuated / clean; 0 where the clean image is blank.

    The attenuated image is the clean image times the gain, so the ratio is
    exact except where |clean| is essentially zero (below ``weak_fraction`` of
    its peak); those pixels take the gain of the nearest well-defined pixel.
    """
    clean = np.asarray(clean, dtype=float)
    attenuated = np.asarray(attenuated, dtype=float)
    blank = ~np.isfinite(clean) | ~np.isfinite(attenuated)
    weak = ~blank & (np.abs(clean) < weak_fraction * np.nanmax(np.abs(clean)))
    good = ~blank & ~weak
    gain = np.zeros(clean.shape)
    gain[good] = attenuated[good] / clean[good]
    if weak.any():
        _, (iy, ix) = ndimage.distance_transform_edt(~good, return_indices=True)
        gain[weak] = gain[iy[weak], ix[weak]]
    values = gain[~blank]
    if values.size and (values.min() < 0 or values.max() > maximum):
        raise ValueError(f"Gain outside [0, {maximum}]: {values.min():.4g} to {values.max():.4g}")
    return gain.astype(np.float32)


def regrid_to(data, source_wcs, target_wcs, target_shape, order=3, block_rows=512):
    """Spline-interpolate ``data`` onto the target grid; NaN outside its footprint."""
    finite = np.isfinite(data)
    filled = np.where(finite, data, 0.0).astype(float)
    coeffs = ndimage.spline_filter(filled, order=order, mode="grid-constant") if order > 1 else filled
    ny, nx = target_shape
    out = np.full(target_shape, np.nan, dtype=np.float32)
    xs = np.arange(nx)
    for y0 in range(0, ny, block_rows):
        ys = np.arange(y0, min(ny, y0 + block_rows))
        xx, yy = np.meshgrid(xs, ys)
        sx, sy = source_wcs.world_to_pixel_values(*target_wcs.pixel_to_world_values(xx, yy))
        coords = np.array([sy, sx])
        values = ndimage.map_coordinates(coeffs, coords, order=order, mode="grid-constant",
                                         cval=0.0, prefilter=False)
        inside = ndimage.map_coordinates(finite.astype(np.float32), coords, order=0,
                                         mode="constant", cval=0.0) > 0.5
        out[y0:y0 + len(ys)] = np.where(inside, values, np.nan)
    return out


def footprint_margin_arcsec(tp_footprint, mosaic, pixel_arcsec) -> float:
    """Minimum signed distance from mosaic pixels to the TP footprint edge.

    Positive inside the footprint, negative for mosaic pixels outside it.
    """
    footprint = np.asarray(tp_footprint, dtype=bool)
    signed = np.where(footprint, ndimage.distance_transform_edt(footprint),
                      -ndimage.distance_transform_edt(~footprint))
    return float(signed[np.asarray(mosaic, dtype=bool)].min() * pixel_arcsec)


def taper_edge(data, footprint, width_arcsec, pixel_arcsec):
    """Keep ``data`` in the footprint; outside, fall from the nearest edge value
    to zero with a cosine over ``width_arcsec``."""
    footprint = np.asarray(footprint, dtype=bool)
    dist, (iy, ix) = ndimage.distance_transform_edt(~footprint, return_indices=True)
    d = dist * pixel_arcsec
    weight = np.where(d < width_arcsec, 0.5 * (1.0 + np.cos(np.pi * d / width_arcsec)), 0.0)
    edge = np.asarray(data)[iy, ix]
    return np.where(footprint, data, np.nan_to_num(edge) * weight).astype(np.float32)


def _rolloff(n_before, n_data, n_after):
    """1 over the data, cosine to 0 across each padding strip (0 at its far end)."""
    w = np.ones(n_before + n_data + n_after)
    if n_before:
        d = np.arange(n_before, 0, -1)
        w[:n_before] = 0.5 * (1.0 + np.cos(np.pi * d / n_before))
    if n_after:
        d = np.arange(1, n_after + 1)
        w[n_before + n_data:] = 0.5 * (1.0 + np.cos(np.pi * d / n_after))
    return w


def pad(data, pad_pixels, mode="zero"):
    """Pad by at least ``pad_pixels`` on each side to fast FFT sizes.

    ``mode="zero"`` pads with zeros. ``mode="reflect"`` mirrors the data into
    the padding and rolls it off to zero with a cosine, so a plane that is cut
    by the grid edge (the TP) continues smoothly instead of stepping to zero.
    """
    ny, nx = data.shape
    py = fft.next_fast_len(ny + 2 * pad_pixels, real=True)
    px = fft.next_fast_len(nx + 2 * pad_pixels, real=True)
    cut = (slice(pad_pixels, pad_pixels + ny), slice(pad_pixels, pad_pixels + nx))
    widths = ((pad_pixels, py - ny - pad_pixels), (pad_pixels, px - nx - pad_pixels))
    if mode == "zero":
        out = np.zeros((py, px), dtype=data.dtype)
        out[cut] = data
    elif mode == "reflect":
        out = np.pad(data, widths, mode="symmetric")
        out *= np.outer(_rolloff(widths[0][0], ny, widths[0][1]),
                        _rolloff(widths[1][0], nx, widths[1][1])).astype(data.dtype)
    else:
        raise ValueError(f"Unknown padding mode: {mode}")
    return out, cut


def azimuthal_power(data, pixel_arcsec, edges_klambda):
    """Mean |FT|^2 in uv annuli; the FT is scaled by the pixel area so that
    images of the same sky at different pixel sizes agree."""
    _, _, k = uv_distance_klambda(data.shape, pixel_arcsec)
    power = np.abs(fft.rfft2(np.asarray(data, dtype=float)) * pixel_arcsec**2) ** 2
    edges = np.asarray(edges_klambda, dtype=float)
    return np.array([power[(k >= lo) & (k < hi)].mean() for lo, hi in zip(edges[:-1], edges[1:])])


def write_product(path, data, template_header, settings, input_hashes, bunit="Jy/beam", beam=None):
    """Write a float32 2-D image with the template's celestial WCS and beam.

    ``settings`` become ``COMB_<KEY>`` keywords; each input hash a HISTORY line.
    ``beam`` = (bmaj", bmin", bpa) overrides the template beam.
    """
    header = WCS(template_header).celestial.to_header()
    bmaj, bmin, bpa = beam if beam is not None else beam_parameters(template_header)
    header["BUNIT"] = bunit
    header["BMAJ"], header["BMIN"], header["BPA"] = bmaj / 3600.0, bmin / 3600.0, bpa
    for key, value in settings.items():
        header[f"HIERARCH COMB_{key.upper()}"] = value
    for name, digest in input_hashes.items():
        header.add_history(f"input {name} sha256 {digest}")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(np.asarray(data, dtype=np.float32), header).writeto(path, overwrite=True)


def combine_prepared(low, high, pixel_arcsec, high_beam, low_beam, method, settings, blank=None):
    """Combine MJy/sr planes and return Jy per high-resolution beam.

    Methods: ``jcomb`` (Butterworth taper with G_high/G_low) and ``feather``.
    The TP is multiplied by ``settings["tp_factor"]``. NaN inputs count as
    zero; ``blank`` pixels become NaN.
    """
    pad_pixels = int(np.ceil(settings.get("padding_arcsec", 0.0) / pixel_arcsec))
    low_p, cut = pad(np.nan_to_num(np.asarray(low, dtype=np.float32)), pad_pixels, mode="reflect")
    high_p, _ = pad(np.nan_to_num(np.asarray(high, dtype=np.float32)), pad_pixels)
    u, v, k = uv_distance_klambda(low_p.shape, pixel_arcsec)
    if method == "jcomb":
        w, ratio = jcomb_weights(u, v, low_beam, high_beam, settings["k1"], settings["k2"],
                                 settings["order"], settings.get("minimum_tp_response", 0.1))
    elif method == "feather":
        w, ratio = feather_weights(u, v, low_beam)
    else:
        raise ValueError(f"Unknown method: {method}")
    del u, v, k
    out = combine_plane(low_p, high_p, w, ratio, settings.get("tp_factor", 1.0))[cut]
    out = mjysr_to_jybeam(out, high_beam[0], high_beam[1]).astype(np.float32)
    if blank is not None:
        out[np.asarray(blank, dtype=bool)] = np.nan
    return out
