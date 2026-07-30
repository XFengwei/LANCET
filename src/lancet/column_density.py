"""Dust-emission and molecular-column-density calculations."""

from __future__ import annotations

from astropy import constants as const
import numpy as np


_H_CGS = const.h.cgs.value
_C_CGS = const.c.cgs.value
_K_B_CGS = const.k_B.cgs.value
_M_P_G = const.m_p.cgs.value
_CGS_PER_MJY = 1.0e-17


def planck_mjysr(
    frequency_ghz: float,
    temperature_k: np.ndarray | float,
) -> np.ndarray:
    """Evaluate the Planck function in MJy/sr."""
    frequency = float(frequency_ghz)
    temperature = np.asarray(temperature_k, dtype=float)
    if not np.isfinite(frequency) or frequency <= 0:
        raise ValueError("Frequency must be finite and positive")
    if np.any(np.isfinite(temperature) & (temperature <= 0)):
        raise ValueError("Finite temperatures must be positive")

    frequency_hz = frequency * 1.0e9
    exponent = _H_CGS * frequency_hz / (_K_B_CGS * temperature)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        intensity_cgs = (
            2.0
            * _H_CGS
            * frequency_hz**3
            / _C_CGS**2
            / np.expm1(exponent)
        )
    return intensity_cgs / _CGS_PER_MJY


def modified_blackbody_mjysr(
    column_density_cm2: np.ndarray | float,
    temperature_k: np.ndarray | float,
    *,
    frequency_ghz: float,
    opacity_cm2_per_g: float,
    mean_molecular_weight: float = 2.8,
) -> np.ndarray:
    """Calculate monochromatic dust intensity without a thin-limit shortcut."""
    column_density = np.asarray(column_density_cm2, dtype=float)
    opacity = float(opacity_cm2_per_g)
    molecular_weight = float(mean_molecular_weight)
    if opacity <= 0 or molecular_weight <= 0:
        raise ValueError("Opacity and mean molecular weight must be positive")
    optical_depth = (
        molecular_weight * _M_P_G * opacity * column_density
    )
    return planck_mjysr(frequency_ghz, temperature_k) * (
        1.0 - np.exp(-optical_depth)
    )


def brightness_to_column_density(
    intensity_mjysr: np.ndarray | float,
    temperature_k: np.ndarray | float,
    *,
    frequency_ghz: float,
    opacity_cm2_per_g: float,
    mean_molecular_weight: float = 2.8,
) -> np.ndarray:
    """Invert dust brightness to H2 column density in cm^-2."""
    intensity = np.asarray(intensity_mjysr, dtype=float)
    opacity = float(opacity_cm2_per_g)
    molecular_weight = float(mean_molecular_weight)
    if opacity <= 0 or molecular_weight <= 0:
        raise ValueError("Opacity and mean molecular weight must be positive")

    blackbody = planck_mjysr(frequency_ghz, temperature_k)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        optical_depth = -np.log1p(-intensity / blackbody)
    return optical_depth / (molecular_weight * _M_P_G * opacity)


def column_density_to_surface_density(
    column_density_cm2: np.ndarray | float,
    *,
    mean_molecular_weight: float = 2.8,
) -> np.ndarray:
    """Convert H2 column density to gas surface density in g/cm^2."""
    molecular_weight = float(mean_molecular_weight)
    if not np.isfinite(molecular_weight) or molecular_weight <= 0:
        raise ValueError("Mean molecular weight must be finite and positive")
    return (
        np.asarray(column_density_cm2, dtype=float)
        * molecular_weight
        * _M_P_G
    )


def dense_mass_msun(
    surface_density_g_cm2: np.ndarray,
    selection_mask: np.ndarray,
    *,
    pixel_scale_x_arcsec: float,
    pixel_scale_y_arcsec: float,
    distance_pc: float,
    threshold_g_cm2: float = 1.0,
) -> float:
    """Sum mass above a surface-density threshold inside a selection."""
    surface_density = np.asarray(surface_density_g_cm2, dtype=float)
    selection = np.asarray(selection_mask, dtype=bool)
    if surface_density.shape != selection.shape:
        raise ValueError("Surface-density and selection shapes must match")
    values = np.asarray(
        [
            pixel_scale_x_arcsec,
            pixel_scale_y_arcsec,
            distance_pc,
            threshold_g_cm2,
        ],
        dtype=float,
    )
    if np.any(~np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("Pixel scales, distance, and threshold must be positive")

    selected = (
        selection
        & np.isfinite(surface_density)
        & (surface_density > threshold_g_cm2)
    )
    pixel_area_cm2 = (
        pixel_scale_x_arcsec
        * distance_pc
        * const.au.cgs.value
        * pixel_scale_y_arcsec
        * distance_pc
        * const.au.cgs.value
    )
    return float(
        np.sum(surface_density[selected], dtype=np.float64)
        * pixel_area_cm2
        / const.M_sun.cgs.value
    )
