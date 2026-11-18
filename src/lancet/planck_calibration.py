"""Comparison of a continuum model with Planck HFI maps (absolute flux scale).

Planck maps are in K_CMB on a HEALPix grid. They are converted to brightness at the
model frequency through the band's transmission for a modified black body, sampled on
a regular grid, and compared with the model smoothed to the Planck beam.
"""

from __future__ import annotations

import numpy as np

H, K, C, TCMB = 6.62607015e-34, 1.380649e-23, 2.99792458e8, 2.7255


def planck_function(nu, temperature):
    """B_nu(T) in SI units."""
    return 2 * H * nu**3 / C**2 / np.expm1(np.minimum(H * nu / (K * temperature), 600))


def dbdt_cmb(nu):
    """dB_nu/dT at the CMB temperature (SI units per K)."""
    x = H * nu / (K * TCMB)
    return 2 * K * nu**2 / C**2 * x**2 / (4 * np.sinh(np.minimum(x, 600) / 2) ** 2)


def read_bandpass(rimo_path, band):
    """(frequency Hz, transmission) of an HFI band from the RIMO file; all positive samples."""
    from astropy.io import fits
    table = fits.getdata(rimo_path, f"BANDPASS_F{band}")
    nu = np.asarray(table["WAVENUMBER"], float) * C * 100.0
    tau = np.clip(np.asarray(table["TRANSMISSION"], float), 0, None)
    keep = (nu > 0) & (tau > 0)
    return nu[keep], tau[keep]


def kcmb_to_mjysr(nu, tau, nu_out, beta, dust_temperature) -> float:
    """MJy/sr at ``nu_out`` per K_CMB in a band, for a modified black body nu^beta B_nu(T)."""
    shape = lambda f: f**beta * planck_function(f, dust_temperature)
    return float(np.trapezoid(tau * dbdt_cmb(nu), nu) * shape(nu_out) / np.trapezoid(tau * shape(nu), nu) * 1e20)


def kcmb_to_mjysr_iras(nu, tau, nu_reference) -> float:
    """The same for the IRAS convention (nu I_nu constant) at the band's reference frequency."""
    return float(np.trapezoid(tau * dbdt_cmb(nu), nu) / np.trapezoid(tau * nu_reference / nu, nu) * 1e20)


def healpix_on_grid(pixels, values, nside, header) -> np.ndarray:
    """Bilinear sample of a partial NESTED Galactic HEALPix map on a celestial grid header."""
    from astropy.coordinates import Galactic
    from astropy.wcs import WCS
    from astropy_healpix import HEALPix
    full = np.full(12 * nside**2, np.nan, dtype=np.float32)
    full[np.asarray(pixels)] = values
    wcs = WCS(header)
    y, x = np.mgrid[:header["NAXIS2"], :header["NAXIS1"]]
    coords = wcs.pixel_to_world(x, y)
    return np.asarray(HEALPix(nside=nside, order="nested", frame=Galactic()).interpolate_bilinear_skycoord(coords, full),
                      dtype=float)


def grid_header(centre_ra_deg, centre_dec_deg, pixel_arcsec, width_arcmin):
    """Square RA/Dec TAN grid."""
    from astropy.io import fits
    n = int(round(width_arcmin * 60.0 / pixel_arcsec))
    return fits.Header({"NAXIS": 2, "NAXIS1": n, "NAXIS2": n, "CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN",
                        "CRVAL1": centre_ra_deg, "CRVAL2": centre_dec_deg, "CRPIX1": n / 2 + 0.5,
                        "CRPIX2": n / 2 + 0.5, "CDELT1": -pixel_arcsec / 3600, "CDELT2": pixel_arcsec / 3600,
                        "CUNIT1": "deg", "CUNIT2": "deg"})


def smooth_with_coverage(data, fwhm_arcsec, pixel_arcsec, kernel_pixels=81):
    """Gaussian-smoothed map with NaN as zero, and the fraction of each beam on finite data."""
    from astropy.convolution import Gaussian2DKernel, convolve_fft
    kernel = Gaussian2DKernel(fwhm_arcsec / pixel_arcsec / np.sqrt(8 * np.log(2)), x_size=kernel_pixels,
                              y_size=kernel_pixels)
    kernel.normalize("integral")
    options = dict(normalize_kernel=False, boundary="fill", fill_value=0)
    return (convolve_fft(np.nan_to_num(data), kernel, **options),
            convolve_fft(np.isfinite(data).astype(float), kernel, **options))


def ratio_stats(a, b, sel) -> dict:
    """Percentiles of a/b and the ratio of sums over the selected pixels."""
    r = (np.asarray(a) / np.asarray(b))[sel]
    return {"n": int(np.sum(sel)), "p16": float(np.percentile(r, 16)), "median": float(np.median(r)),
            "p84": float(np.percentile(r, 84)), "sum_ratio": float(np.sum(a[sel]) / np.sum(b[sel]))}


def linear_fit(model, planck, sel, n_boot=2000, seed=1) -> dict:
    """model = gain x Planck + offset over the selected pixels, with a pixel bootstrap.

    The bootstrap ignores the correlation of pixels inside a Planck beam and
    understates the uncertainty.
    """
    x, y = np.asarray(planck)[sel], np.asarray(model)[sel]
    gain, offset = np.polyfit(x, y, 1)
    rng = np.random.default_rng(seed)
    boot = np.array([np.polyfit(x[i], y[i], 1) for i in rng.integers(0, x.size, (n_boot, x.size))])
    return {"n": int(x.size), "gain": float(gain), "offset_mjysr": float(offset),
            "background_mjysr": float(-offset / gain),
            "gain_boot_p16_p84": np.percentile(boot[:, 0], [16, 84]).tolist(),
            "offset_boot_p16_p84": np.percentile(boot[:, 1], [16, 84]).tolist()}


def bright_region(model_smoothed, allowed, radius_pixels):
    """Peak of the smoothed model among ``allowed`` pixels and the disc of ``radius_pixels`` around it.

    Returns ((y, x) of the peak, boolean disc restricted to ``allowed``).
    """
    masked = np.where(allowed, model_smoothed, -np.inf)
    peak = np.unravel_index(int(np.argmax(masked)), masked.shape)
    yy, xx = np.indices(masked.shape)
    disc = np.hypot(yy - peak[0], xx - peak[1]) <= radius_pixels
    return (int(peak[0]), int(peak[1])), disc & np.asarray(allowed, dtype=bool)


def bright_region_calibration(model, planck, peak, region, reference, background_mjysr) -> dict:
    """Model / Planck in the brightest beam, three ways.

    ``direct``: plain ratio (no background treatment). ``background_removed``:
    Planck minus the fitted uniform background. ``differential``: (bright - reference)
    for both maps, where ``reference`` is the faintest well-covered area; a uniform
    level in either map cancels.
    """
    m_ref, p_ref = float(np.median(model[reference])), float(np.median(planck[reference]))
    m_pk, p_pk = float(model[peak]), float(planck[peak])
    m_reg, p_reg = float(np.mean(model[region])), float(np.mean(planck[region]))
    out = {"peak_pixel": {"model_mjysr": m_pk, "planck_mjysr": p_pk},
           "region_mean": {"model_mjysr": m_reg, "planck_mjysr": p_reg, "n": int(np.sum(region))},
           "reference": {"model_mjysr": m_ref, "planck_mjysr": p_ref, "n": int(np.sum(reference))},
           "ratio_model_over_planck": {}}
    for name, (m, p) in (("peak_pixel", (m_pk, p_pk)), ("region_mean", (m_reg, p_reg))):
        out["ratio_model_over_planck"][name] = {
            "direct": m / p,
            "background_removed": m / (p - background_mjysr),
            "differential": (m - m_ref) / (p - p_ref)}
    return out
