"""uv-plane combination of a low-resolution (TP) and a high-resolution image.

J-comb (Jiao et al. 2022, Sci. China PMA 65, 299511) and feather differ only in
the weight given to the low-resolution image, ``w_low``, and in the beam ratio
``ratio`` applied to it:

    C(k) = w_low * ratio * tp_factor * FT(low) + (1 - w_low) * FT(high)

All uv distances are in kilo-wavelengths (k-lambda) on the ``rfft2`` layout.
"""

from __future__ import annotations

import numpy as np
from scipy import fft

ARCSEC_RAD = np.pi / 180.0 / 3600.0
FWHM_TO_SIGMA = 1.0 / np.sqrt(8.0 * np.log(2.0))

Beam = tuple[float, float, float]


def uv_distance_klambda(shape: tuple[int, int], pixel_arcsec: float):
    """Return (u, v, k) in k-lambda for the ``rfft2`` layout of ``shape``."""
    ny, nx = shape
    d = pixel_arcsec * ARCSEC_RAD
    u1 = fft.rfftfreq(nx, d) / 1e3
    v1 = fft.fftfreq(ny, d) / 1e3
    u, v = np.meshgrid(u1, v1)
    return u, v, np.hypot(u, v)


def gaussian_transfer(u, v, bmaj_arcsec, bmin_arcsec, bpa_deg, cdelt1_sign=-1):
    """Fourier transfer of an elliptical Gaussian beam, 1 at k = 0.

    BPA is measured east of north. With ``cdelt1_sign = -1`` east is -x, so the
    major axis points along (-sin BPA, cos BPA) in pixel (x, y).
    """
    theta = np.deg2rad(bpa_deg)
    s = cdelt1_sign
    along = s * u * np.sin(theta) + v * np.cos(theta)
    across = u * np.cos(theta) - s * v * np.sin(theta)
    smaj = bmaj_arcsec * FWHM_TO_SIGMA * ARCSEC_RAD
    smin = bmin_arcsec * FWHM_TO_SIGMA * ARCSEC_RAD
    return np.exp(-2.0 * np.pi**2 * 1e6 * ((smaj * along) ** 2 + (smin * across) ** 2))


def butterworth_taper(k, k1, k2, order):
    """1 for k <= k1, 0 for k >= k2, a rescaled Butterworth in between.

    b(k) = 1 / (1 + (k/k0)^(2n)) with k0 = (k1 + k2)/2, rescaled to run from 1
    at k1 to 0 at k2 (the paper's Eq. 18 normalisation applied to the
    well-behaved Butterworth form; its printed Eq. 19 diverges at k0).
    """
    if not k1 < k2:
        raise ValueError(f"Need k1 < k2, got {k1}, {k2}")
    if order < 1:
        raise ValueError(f"Butterworth order must be >= 1, got {order}")
    k0 = 0.5 * (k1 + k2)

    def b(x):
        return 1.0 / (1.0 + (np.asarray(x, dtype=float) / k0) ** (2 * order))

    w = (b(k) - b(k2)) / (b(k1) - b(k2))
    return np.where(k <= k1, 1.0, np.where(k >= k2, 0.0, w))


def jcomb_weights(u, v, low_beam: Beam, high_beam: Beam, k1, k2, order, minimum_tp_response=0.1):
    """J-comb weights: (w_low, ratio = G_high / G_low where w_low > 0, else 0)."""
    w_low = butterworth_taper(np.hypot(u, v), k1, k2, order)
    g_low = gaussian_transfer(u, v, *low_beam)
    used = w_low > 0
    if np.any(g_low[used] < minimum_tp_response):
        raise ValueError(
            f"TP beam response falls to {g_low[used].min():.3g} inside the TP weight "
            f"(minimum {minimum_tp_response}); lower the upper threshold")
    ratio = np.zeros_like(g_low)
    ratio[used] = gaussian_transfer(u[used], v[used], *high_beam) / g_low[used]
    return w_low, ratio


def feather_weights(u, v, low_beam: Beam, minimum_response=1e-6):
    """Feather weights: the TP's FT is added unchanged (w_low * ratio = 1) and the
    interferometer is down-weighted by 1 - G_low.

    Returned as (w_low = G_low, ratio = 1/G_low); where G_low < ``minimum_response``
    the TP term is dropped (its content there is ~1e-6 of its peak).
    """
    g_low = gaussian_transfer(u, v, *low_beam)
    ratio = np.zeros_like(g_low)
    keep = g_low >= minimum_response
    ratio[keep] = 1.0 / g_low[keep]
    return g_low, ratio


def combine_plane(low, high, w_low, ratio, tp_factor=1.0):
    """Combine two same-grid, same-unit, finite planes in the uv plane."""
    low = np.asarray(low, dtype=np.float32)
    high = np.asarray(high, dtype=np.float32)
    if low.shape != high.shape:
        raise ValueError(f"Shapes differ: {low.shape} vs {high.shape}")
    w = np.asarray(w_low, dtype=np.float32)
    lowpart = (w * np.asarray(ratio, dtype=np.float32) * np.float32(tp_factor))
    total = fft.rfft2(low) * lowpart + fft.rfft2(high) * (np.float32(1.0) - w)
    return fft.irfft2(total, s=low.shape).astype(np.float32)
