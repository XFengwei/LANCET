"""Measurements that validate a uv-plane combination."""

from __future__ import annotations

import numpy as np
from scipy import fft, ndimage

from lancet.combination import gaussian_transfer, uv_distance_klambda
from lancet.combination_prep import combine_prepared, jybeam_to_mjysr, mjysr_to_jybeam, pad


def smooth_to(data, pixel_arcsec, from_beam, to_fwhm_arcsec):
    """Surface-brightness map smoothed from ``from_beam`` to a circular beam.

    NaNs are excluded by normalised convolution, so mosaic edges are not
    pulled towards zero; pixels where less than 0.2 of the kernel weight falls
    on valid data are NaN. The map is zero-padded by three target beams so the
    FFT does not wrap one grid edge onto the other.
    """
    data = np.asarray(data, dtype=float)
    finite = np.isfinite(data)
    pad_pixels = int(np.ceil(3.0 * to_fwhm_arcsec / pixel_arcsec))
    values, cut = pad(np.where(finite, data, 0.0), pad_pixels)
    weights, _ = pad(finite.astype(float), pad_pixels)
    u, v, _ = uv_distance_klambda(values.shape, pixel_arcsec)
    target = gaussian_transfer(u, v, to_fwhm_arcsec, to_fwhm_arcsec, 0.0)
    num = fft.irfft2(fft.rfft2(values) * (target / gaussian_transfer(u, v, *from_beam)), s=values.shape)[cut]
    den = fft.irfft2(fft.rfft2(weights) * target, s=values.shape)[cut]
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0.2, num / den, np.nan)   # 0.25 at a grid corner


def large_scale_ratio(product_jybeam, tp_mjysr, pixel_arcsec, product_beam, tp_fwhm_arcsec, mask) -> dict:
    """Product smoothed to the TP beam over the TP, in surface brightness."""
    product = jybeam_to_mjysr(product_jybeam, product_beam[0], product_beam[1])
    smooth = smooth_to(product, pixel_arcsec, product_beam, tp_fwhm_arcsec)
    sel = np.asarray(mask, dtype=bool) & np.isfinite(smooth) & (np.asarray(tp_mjysr) > 0)
    ratio = smooth[sel] / np.asarray(tp_mjysr)[sel]
    p16, median, p84 = np.percentile(ratio, [16, 50, 84])
    return {"median": float(median), "p16": float(p16), "p84": float(p84), "n": int(sel.sum())}


def negative_fraction(data, sigma, mask, threshold=-3.0) -> float:
    values = np.asarray(data)[np.asarray(mask, dtype=bool)]
    values = values[np.isfinite(values)]
    return float(np.mean(values < threshold * sigma))


def mosaic_window(gain, min_gain, pixel_arcsec, smooth_fwhm_arcsec) -> np.ndarray:
    """Gain >= ``min_gain`` region smoothed by a Gaussian: a common window for uv comparisons."""
    sigma = smooth_fwhm_arcsec / np.sqrt(8.0 * np.log(2.0)) / pixel_arcsec
    return ndimage.gaussian_filter((np.asarray(gain) >= min_gain).astype(np.float32), sigma)


def tp_at_high_resolution(low, pixel_arcsec, low_beam, high_beam, k_max, padding_arcsec):
    """The TP (MJy/sr) brought to the high-resolution beam: FT(low) x G_high/G_low for k <= k_max.

    This is J-comb's TP term with unit weight; the TP is mirror-padded as in the combination.
    """
    pad_pixels = int(np.ceil(padding_arcsec / pixel_arcsec))
    padded, cut = pad(np.nan_to_num(np.asarray(low, dtype=np.float32)), pad_pixels, mode="reflect")
    u, v, k = uv_distance_klambda(padded.shape, pixel_arcsec)
    keep = k <= k_max
    ratio = np.zeros_like(k)
    ratio[keep] = gaussian_transfer(u[keep], v[keep], *high_beam) / gaussian_transfer(u[keep], v[keep], *low_beam)
    del u, v, k, keep
    return fft.irfft2(fft.rfft2(padded) * ratio, s=padded.shape)[cut].astype(np.float32)


def _windowed_ft(data, window):
    return fft.rfft2(np.nan_to_num(np.asarray(data, dtype=float)) * window)


def _bins(shape, pixel_arcsec, edges_klambda):
    _, _, k = uv_distance_klambda(shape, pixel_arcsec)
    edges = np.asarray(edges_klambda, dtype=float)
    return [(k >= lo) & (k < hi) for lo, hi in zip(edges[:-1], edges[1:])]


def windowed_amplitude_ratio(a, b, window, pixel_arcsec, edges_klambda) -> np.ndarray:
    """Per uv bin, mean |FT(a w)| / mean |FT(b w)| for a common window w."""
    fa, fb = np.abs(_windowed_ft(a, window)), np.abs(_windowed_ft(b, window))
    return np.array([fa[m].mean() / fb[m].mean() for m in _bins(np.shape(a), pixel_arcsec, edges_klambda)])


def cross_spectrum_slopes(a, b, window, pixel_arcsec, edges_klambda):
    """Per uv bin, the scale of ``a`` relative to ``b`` from their windowed cross-spectrum.

    Returns (Re<A B*>/<|B|^2>, <|A|^2>/Re<A B*>). Noise or structure in ``b`` that
    ``a`` lacks lowers the first; the same in ``a`` raises the second, so the two
    bracket the scale when the maps are only partly correlated.
    """
    fa, fb = _windowed_ft(a, window), _windowed_ft(b, window)
    s_ab, s_aa = [], []
    for m in _bins(np.shape(a), pixel_arcsec, edges_klambda):
        cross = np.real(np.sum(fa[m] * np.conj(fb[m])))
        s_ab.append(cross / np.sum(np.abs(fb[m]) ** 2))
        s_aa.append(np.sum(np.abs(fa[m]) ** 2) / cross)
    return np.array(s_ab), np.array(s_aa)


def tapered_window(mask, pixel_arcsec, taper_arcsec) -> np.ndarray:
    """1 inside ``mask``, falling to 0 at its edge with a cosine over ``taper_arcsec``."""
    mask = np.asarray(mask, dtype=bool)
    distance = ndimage.distance_transform_edt(np.pad(mask, 1))[1:-1, 1:-1] * pixel_arcsec
    return (0.5 - 0.5 * np.cos(np.pi * np.clip(distance / taper_arcsec, 0.0, 1.0))).astype(np.float32)


def remove_window_mean(data, window) -> np.ndarray:
    """``data`` minus its window-weighted mean (NaN as 0), so the mean does not leak into other scales."""
    data = np.nan_to_num(np.asarray(data, dtype=float))
    return data - np.sum(data * window) / np.sum(window)


def windowed_power(data, window, pixel_arcsec, edges_klambda) -> np.ndarray:
    """Power spectrum per uv bin of ``data`` inside ``window``, mean removed.

    Mean |FT((data - mean) w)|^2 x pixel area / sum(w^2): for surface-brightness
    input in X the unit is X^2 arcsec^2, independent of the window size.
    """
    ft = fft.rfft2(remove_window_mean(data, window) * window)
    power = np.abs(ft) ** 2 * pixel_arcsec**2 / np.sum(np.asarray(window, dtype=float) ** 2)
    return np.array([power[m].mean() for m in _bins(np.shape(data), pixel_arcsec, edges_klambda)])


def psd_match_factor(power_reference, power_tp, centres_klambda, band_klambda) -> dict:
    """TP factor that lifts the TP power spectrum onto the reference at its weakest point in a band.

    The amplitude ratio sqrt(P_reference / P_TP) is taken per bin; the factor is
    its maximum over the bins whose centre lies in ``band_klambda``, i.e. the TP
    times the factor touches the reference from below there.
    """
    centres = np.asarray(centres_klambda, dtype=float)
    ratio = np.sqrt(np.asarray(power_reference, dtype=float) / np.asarray(power_tp, dtype=float))
    sel = np.flatnonzero((centres >= band_klambda[0]) & (centres <= band_klambda[1]))
    if sel.size == 0:
        raise ValueError(f"No uv bin centre inside {band_klambda} klambda")
    best = sel[np.argmax(ratio[sel])]
    return {"factor": float(ratio[best]), "at_klambda": float(centres[best]),
            "ratio_in_band": ratio[sel].tolist(), "bins_klambda": centres[sel].tolist()}


def tp_only_ratio(product_jybeam, tp_mjysr, high_beam, region) -> dict:
    """Percentiles of product / TP (both as surface brightness) inside ``region``."""
    product = jybeam_to_mjysr(np.asarray(product_jybeam, dtype=float), high_beam[0], high_beam[1])
    tp = np.asarray(tp_mjysr, dtype=float)
    sel = np.asarray(region, dtype=bool) & np.isfinite(product) & (tp > 0)
    if not sel.any():
        return {"p01": None, "p16": None, "median": None, "p84": None, "p99": None, "n": 0}
    p = np.percentile(product[sel] / tp[sel], [1, 16, 50, 84, 99])
    return {"p01": float(p[0]), "p16": float(p[1]), "median": float(p[2]), "p84": float(p[3]),
            "p99": float(p[4]), "n": int(sel.sum())}


def core_contrast_ratio(product, reference, positions_pix, pixel_arcsec, annulus_arcsec=(3.0, 6.0), box=5) -> dict:
    """Background-subtracted peak of ``product`` over that of ``reference`` at each position.

    Peak: maximum in a ``box`` x ``box`` pixel box; background: median in an
    annulus. Smooth emission added by the TP cancels, so a ratio of 1 means the
    combination left the compact structure unchanged.
    """
    half = box // 2
    r_in, r_out = (a / pixel_arcsec for a in annulus_arcsec)
    n = int(np.ceil(r_out))
    yy, xx = np.mgrid[-n:n + 1, -n:n + 1]
    ring = (np.hypot(yy, xx) >= r_in) & (np.hypot(yy, xx) <= r_out)
    ny, nx = np.shape(product)
    ratios = []
    for x, y in positions_pix:
        x, y = int(round(x)), int(round(y))
        if not (n <= y < ny - n and n <= x < nx - n):
            continue
        peak = (slice(y - half, y + half + 1), slice(x - half, x + half + 1))
        stamp = (slice(y - n, y + n + 1), slice(x - n, x + n + 1))
        c_ref = np.nanmax(reference[peak]) - np.nanmedian(reference[stamp][ring])
        c_prod = np.nanmax(product[peak]) - np.nanmedian(product[stamp][ring])
        if np.isfinite(c_ref) and np.isfinite(c_prod) and c_ref > 0:
            ratios.append(c_prod / c_ref)
    ratios = np.asarray(ratios)
    if ratios.size == 0:
        return {"median": float("nan"), "p16": float("nan"), "p84": float("nan"), "n": 0}
    p16, median, p84 = np.percentile(ratios, [16, 50, 84])
    return {"median": float(median), "p16": float(p16), "p84": float(p84), "n": int(ratios.size)}


def beam_volume_outside(method_transfer, high_transfer, shape) -> float:
    """Fraction of the effective PSF volume not covered by the high-resolution Gaussian."""
    eff = fft.irfft2(method_transfer, s=shape)
    high = fft.irfft2(high_transfer, s=shape)
    return float(1.0 - np.sum(np.minimum(eff, high)) / np.sum(eff))


def synthetic_gate(gain, pixel_arcsec, high_beam, low_beam, settings, aca_shortest_klambda=6.6,
                   min_gain=0.5, n_sources=30, rng_seed=3168, flux_tolerance=0.02, rms_tolerance=0.01) -> dict:
    """J-comb on a mock sky with a known answer, on the given grid.

    Mock sky (MJy/sr): a smooth field (60" structure) plus compact Gaussian
    sources. Mock TP = sky through the TP beam; mock interferometer = sky
    through the high beam, high-passed with half power at the ACA's shortest
    spacing. The answer is the sky through the high beam.
    """
    shape = np.shape(gain)
    rng = np.random.default_rng(rng_seed)
    u, v, k = uv_distance_klambda(shape, pixel_arcsec)
    smooth = fft.irfft2(fft.rfft2(rng.random(shape)) * gaussian_transfer(u, v, 60.0, 60.0, 0.0), s=shape)
    smooth = (smooth - smooth.min()) / np.ptp(smooth)
    yy, xx = np.indices(shape)
    inside = np.flatnonzero(np.asarray(gain).ravel() >= min_gain)
    for idx in rng.choice(inside, size=min(n_sources, inside.size), replace=False):
        cy, cx = divmod(idx, shape[1])
        fwhm = rng.uniform(1.0, 3.0) / pixel_arcsec
        smooth += rng.uniform(0.5, 3.0) * np.exp(-4 * np.log(2) * ((xx - cx) ** 2 + (yy - cy) ** 2) / fwhm**2)
    sky_ft = fft.rfft2(smooth)
    low = fft.irfft2(sky_ft * gaussian_transfer(u, v, *low_beam), s=shape)
    truth_ft = sky_ft * gaussian_transfer(u, v, *high_beam)
    highpass = 1.0 - np.exp(-np.log(2.0) * (k / aca_shortest_klambda) ** 2)
    high = fft.irfft2(truth_ft * highpass, s=shape)
    truth = mjysr_to_jybeam(fft.irfft2(truth_ft, s=shape), high_beam[0], high_beam[1])
    del u, v, k, sky_ft, truth_ft
    out = combine_prepared(low, high, pixel_arcsec, high_beam, low_beam, "jcomb", settings)
    sel = np.asarray(gain) >= min_gain
    flux_error = abs(out[sel].sum() - truth[sel].sum()) / truth[sel].sum()
    rms = np.sqrt(np.mean((out[sel] - truth[sel]) ** 2)) / truth[sel].max()
    return {"flux_error": float(flux_error), "rms_fraction_of_peak": float(rms),
            "passed": bool(flux_error < flux_tolerance and rms < rms_tolerance), "n": int(sel.sum())}


def zero_spacing_check(combined_at_tp_beam, tp, regions) -> dict:
    """Is the TP (zero-spacing) information kept by a combination?

    ``combined_at_tp_beam`` is the combined image smoothed to the TP beam and ``tp`` the TP it
    was combined with, on the same grid and in the same surface-brightness unit. Per region
    (boolean masks): the ratio of summed brightness (total flux), percentiles of the pixel
    ratio, and the rms of the difference relative to the mean TP brightness.
    """
    a, b = np.asarray(combined_at_tp_beam, dtype=float), np.asarray(tp, dtype=float)
    out = {}
    for name, region in regions.items():
        sel = np.asarray(region, dtype=bool) & np.isfinite(a) & np.isfinite(b) & (b > 0)
        if not sel.any():
            out[name] = {"n": 0}
            continue
        ratio = a[sel] / b[sel]
        p01, p16, median, p84, p99 = np.percentile(ratio, [1, 16, 50, 84, 99])
        out[name] = {"n": int(sel.sum()), "flux_ratio": float(a[sel].sum() / b[sel].sum()),
                     "pixel_ratio": {"p01": float(p01), "p16": float(p16), "median": float(median),
                                     "p84": float(p84), "p99": float(p99)},
                     "rms_difference_over_mean_tp": float(np.sqrt(np.mean((a[sel] - b[sel]) ** 2)) / b[sel].mean())}
    return out
