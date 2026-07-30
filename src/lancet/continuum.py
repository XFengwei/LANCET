"""Synthesis of dust and free-free continuum emission maps."""

from __future__ import annotations

from astropy.convolution import Gaussian2DKernel, convolve_fft
import numpy as np


def scale_power_law(
    image: np.ndarray,
    *,
    reference_frequency_ghz: float,
    target_frequency_ghz: float,
    spectral_index: float,
) -> np.ndarray:
    """Scale a flux-density image using ``S_nu proportional to nu**alpha``."""
    data = np.asarray(image)
    reference = float(reference_frequency_ghz)
    target = float(target_frequency_ghz)
    index = float(spectral_index)
    if not np.all(np.isfinite([reference, target, index])):
        raise ValueError("Frequencies and spectral index must be finite")
    if reference <= 0 or target <= 0:
        raise ValueError("Frequencies must be positive")
    return data * (target / reference) ** index


def synthesize_continuum(
    dust_image: np.ndarray,
    freefree_image: np.ndarray,
    *,
    reference_frequency_ghz: float = 230.0,
    target_frequency_ghz: float = 230.0,
    dust_spectral_index: float = 3.8,
    freefree_spectral_index: float = -0.1,
    missing_freefree_value: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Scale and combine dust and free-free flux-density images.

    Missing dust values remain missing. Missing free-free values are replaced
    by ``missing_freefree_value`` so that limited radio coverage does not erase
    otherwise valid dust-model pixels.
    """
    dust = np.asarray(dust_image)
    freefree = np.asarray(freefree_image)
    if dust.shape != freefree.shape:
        raise ValueError(
            f"Dust shape {dust.shape} does not match free-free shape "
            f"{freefree.shape}"
        )
    if dust.ndim != 2:
        raise ValueError("Continuum component images must be two-dimensional")
    if not np.isfinite(missing_freefree_value):
        raise ValueError("missing_freefree_value must be finite")

    freefree_filled = np.where(
        np.isfinite(freefree),
        freefree,
        missing_freefree_value,
    )
    dust_scaled = scale_power_law(
        dust,
        reference_frequency_ghz=reference_frequency_ghz,
        target_frequency_ghz=target_frequency_ghz,
        spectral_index=dust_spectral_index,
    )
    freefree_scaled = scale_power_law(
        freefree_filled,
        reference_frequency_ghz=reference_frequency_ghz,
        target_frequency_ghz=target_frequency_ghz,
        spectral_index=freefree_spectral_index,
    )
    return dust_scaled, freefree_scaled, dust_scaled + freefree_scaled


def smooth_circular_jybeam(
    image: np.ndarray,
    *,
    pixel_scale_arcsec: float,
    original_beam_arcsec: float,
    target_beam_arcsec: float,
    replace_nan_with_zero: bool = True,
) -> np.ndarray:
    """Smooth a Jy/beam map and rescale it to Jy per target beam."""
    data = np.asarray(image)
    if data.ndim != 2:
        raise ValueError("Input image must be two-dimensional")

    pixel_scale = float(pixel_scale_arcsec)
    original_beam = float(original_beam_arcsec)
    target_beam = float(target_beam_arcsec)
    if not np.all(np.isfinite([pixel_scale, original_beam, target_beam])):
        raise ValueError("Pixel scale and beam sizes must be finite")
    if pixel_scale <= 0 or original_beam <= 0:
        raise ValueError("Pixel scale and original beam must be positive")
    if target_beam <= original_beam:
        raise ValueError("Target beam must be larger than original beam")

    working = (
        np.where(np.isfinite(data), data, 0.0)
        if replace_nan_with_zero
        else data
    )
    kernel_fwhm_arcsec = np.sqrt(target_beam**2 - original_beam**2)
    kernel_sigma_pixels = (
        kernel_fwhm_arcsec
        / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        / pixel_scale
    )
    kernel = Gaussian2DKernel(kernel_sigma_pixels)
    kernel.normalize()
    convolved = convolve_fft(
        working,
        kernel,
        boundary="fill",
        fill_value=np.nan,
        nan_treatment="interpolate",
        preserve_nan=True,
        allow_huge=True,
    )
    return convolved * (target_beam / original_beam) ** 2


def calibrate_simba(
    raw_image: np.ndarray,
    *,
    invalid_floor: float = -1e5,
    flux_scale: float = 0.20,
) -> np.ndarray:
    """Apply the adopted SIMBA invalid-pixel mask and flux calibration."""
    data = np.asarray(raw_image)
    floor = float(invalid_floor)
    scale = float(flux_scale)
    if not np.all(np.isfinite([floor, scale])):
        raise ValueError("SIMBA invalid floor and flux scale must be finite")
    if scale <= 0:
        raise ValueError("SIMBA flux scale must be positive")
    return np.where(data > floor, data * scale, np.nan)
