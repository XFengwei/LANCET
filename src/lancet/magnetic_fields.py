"""Planck Stokes-map polarization diagnostics."""

from __future__ import annotations

import numpy as np


def debiased_polarized_intensity(
    stokes_q: np.ndarray,
    stokes_u: np.ndarray,
    rms_qu: float,
) -> np.ndarray:
    """Return the modified-asymptotic debiased polarized intensity."""
    q = np.asarray(stokes_q, dtype=float)
    u = np.asarray(stokes_u, dtype=float)
    if q.shape != u.shape:
        raise ValueError("Stokes Q and U must have identical shapes")
    noise = float(rms_qu)
    if not np.isfinite(noise) or noise <= 0:
        raise ValueError("rms_qu must be finite and positive")
    polarized = np.hypot(q, u)
    correction = np.zeros_like(polarized)
    nonzero = polarized > 0
    correction[nonzero] = (
        noise**2
        * (
            1.0
            - np.exp(
                -(polarized[nonzero] ** 2) / noise**2
            )
        )
        / (2.0 * polarized[nonzero])
    )
    output = np.maximum(polarized - correction, 0.0)
    output[~np.isfinite(q) | ~np.isfinite(u)] = np.nan
    return output


def polarization_angle_degrees(
    stokes_q: np.ndarray,
    stokes_u: np.ndarray,
) -> np.ndarray:
    """Return the IAU polarization angle in degrees."""
    q = np.asarray(stokes_q, dtype=float)
    u = np.asarray(stokes_u, dtype=float)
    if q.shape != u.shape:
        raise ValueError("Stokes Q and U must have identical shapes")
    angle = np.degrees(0.5 * np.arctan2(-u, q))
    angle[~np.isfinite(q) | ~np.isfinite(u)] = np.nan
    return angle


def magnetic_field_angle_degrees(
    stokes_q: np.ndarray,
    stokes_u: np.ndarray,
) -> np.ndarray:
    """Return the plane-of-sky field orientation modulo 180 degrees."""
    return (polarization_angle_degrees(stokes_q, stokes_u) + 90.0) % 180.0


def polarization_angle_error_degrees(
    stokes_q: np.ndarray,
    stokes_u: np.ndarray,
    rms_qu: float,
) -> np.ndarray:
    """Return the high-S/N polarization-angle error approximation."""
    q = np.asarray(stokes_q, dtype=float)
    u = np.asarray(stokes_u, dtype=float)
    if q.shape != u.shape:
        raise ValueError("Stokes Q and U must have identical shapes")
    noise = float(rms_qu)
    if not np.isfinite(noise) or noise <= 0:
        raise ValueError("rms_qu must be finite and positive")
    polarized = np.hypot(q, u)
    error = np.full_like(polarized, np.nan)
    valid = np.isfinite(polarized) & (polarized > 0)
    error[valid] = np.degrees(0.5 * noise / polarized[valid])
    return error


def planck_polarization_products(
    stokes_i: np.ndarray,
    stokes_q: np.ndarray,
    stokes_u: np.ndarray,
    *,
    rms_i: float,
    rms_qu: float,
    intensity_snr_thresholds: tuple[float, ...] = (2.0, 3.0),
    polarized_snr_thresholds: tuple[float, ...] = (2.0, 3.0),
    angle_snr_threshold: float = 2.0,
    fraction_intensity_snr_threshold: float = 3.0,
) -> dict[str, np.ndarray]:
    """Build masks, debiased polarization, fractions, and angle products."""
    intensity = np.asarray(stokes_i, dtype=float)
    q = np.asarray(stokes_q, dtype=float)
    u = np.asarray(stokes_u, dtype=float)
    if intensity.shape != q.shape or intensity.shape != u.shape:
        raise ValueError("Stokes I, Q, and U must have identical shapes")
    intensity_noise = float(rms_i)
    qu_noise = float(rms_qu)
    if (
        not np.isfinite(intensity_noise)
        or intensity_noise <= 0
        or not np.isfinite(qu_noise)
        or qu_noise <= 0
    ):
        raise ValueError("rms_i and rms_qu must be finite and positive")

    polarized = debiased_polarized_intensity(q, u, qu_noise)
    output: dict[str, np.ndarray] = {
        "polarized_intensity_debiased": polarized,
    }
    for threshold in intensity_snr_thresholds:
        output[f"intensity_snr{threshold:g}_mask"] = (
            np.isfinite(intensity)
            & (intensity > float(threshold) * intensity_noise)
        )
    for threshold in polarized_snr_thresholds:
        output[f"polarized_snr{threshold:g}_mask"] = (
            np.isfinite(polarized)
            & (polarized >= float(threshold) * qu_noise)
        )

    angle_mask = (
        np.isfinite(polarized)
        & (polarized >= float(angle_snr_threshold) * qu_noise)
    )
    polarization_angle = polarization_angle_degrees(q, u)
    magnetic_angle = magnetic_field_angle_degrees(q, u)
    angle_error = polarization_angle_error_degrees(q, u, qu_noise)
    for values in (polarization_angle, magnetic_angle, angle_error):
        values[~angle_mask] = np.nan
    fraction = np.full_like(polarized, np.nan)
    fraction_mask = (
        np.isfinite(intensity)
        & (intensity > float(fraction_intensity_snr_threshold) * intensity_noise)
        & np.isfinite(polarized)
    )
    fraction[fraction_mask] = (
        polarized[fraction_mask] / intensity[fraction_mask]
    )
    output.update(
        {
            "polarization_fraction_debiased": fraction,
            "polarization_angle": polarization_angle,
            "polarization_angle_error": angle_error,
            "magnetic_field_angle": magnetic_angle,
        }
    )
    return output
