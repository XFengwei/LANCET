"""Deconvolution of two-dimensional elliptical Gaussian sources."""

from __future__ import annotations

from typing import Literal

import numpy as np


FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
SIGMA_TO_FWHM = 1.0 / FWHM_TO_SIGMA


class UnresolvedSourceError(ValueError):
    """Raised when beam subtraction does not leave a physical Gaussian."""


def fwhm_to_sigma(fwhm: float) -> float:
    """Convert a Gaussian FWHM to sigma in the same units."""
    return float(fwhm) * FWHM_TO_SIGMA


def sigma_to_fwhm(sigma: float) -> float:
    """Convert a Gaussian sigma to FWHM in the same units."""
    return float(sigma) * SIGMA_TO_FWHM


def _validate_ellipse(
    major: float,
    minor: float,
    pa_deg: float,
    *,
    name: str,
) -> tuple[float, float, float]:
    major = float(major)
    minor = float(minor)
    pa_deg = float(pa_deg)

    if not np.all(np.isfinite([major, minor, pa_deg])):
        raise ValueError(f"{name} axes and position angle must be finite")
    if major <= 0 or minor <= 0:
        raise ValueError(f"{name} axes must be positive")
    if minor > major:
        major, minor = minor, major
        pa_deg += 90.0
    return major, minor, pa_deg % 180.0


def _rotation_matrix(theta_rad: float) -> np.ndarray:
    cosine = np.cos(theta_rad)
    sine = np.sin(theta_rad)
    return np.array([[cosine, -sine], [sine, cosine]])


def _ellipse_to_covariance(
    major: float,
    minor: float,
    pa_deg: float,
) -> np.ndarray:
    """Convert FWHM axes and astronomical PA to a covariance matrix.

    The internal image plane has +x toward image right (west) and +y north.
    Astronomical position angle is measured east of north.
    """
    theta_rad = np.deg2rad(90.0 + pa_deg)
    rotation = _rotation_matrix(theta_rad)
    variances = np.diag(
        [fwhm_to_sigma(major) ** 2, fwhm_to_sigma(minor) ** 2]
    )
    return rotation @ variances @ rotation.T


def _covariance_to_ellipse(
    covariance: np.ndarray,
    *,
    tolerance: float,
) -> tuple[float, float, float]:
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    minor_variance, major_variance = eigenvalues

    major = sigma_to_fwhm(np.sqrt(max(major_variance, 0.0)))
    minor = sigma_to_fwhm(np.sqrt(max(minor_variance, 0.0)))

    if np.isclose(major, minor, rtol=tolerance, atol=tolerance):
        return major, minor, 0.0

    major_vector = eigenvectors[:, 1]
    theta_deg = np.rad2deg(np.arctan2(major_vector[1], major_vector[0]))
    pa_deg = (theta_deg - 90.0) % 180.0
    return major, minor, pa_deg


def _unresolved(
    message: str,
    on_unresolved: Literal["nan", "raise"],
) -> tuple[float, float, float]:
    if on_unresolved == "raise":
        raise UnresolvedSourceError(message)
    if on_unresolved != "nan":
        raise ValueError("on_unresolved must be 'nan' or 'raise'")
    return np.nan, np.nan, np.nan


def deconvolve_gaussian(
    observed_major: float,
    observed_minor: float,
    observed_pa_deg: float,
    beam_major: float,
    beam_minor: float,
    beam_pa_deg: float,
    *,
    absolute_tolerance: float = 1e-12,
    relative_tolerance: float = 1e-10,
    on_unresolved: Literal["nan", "raise"] = "nan",
) -> tuple[float, float, float]:
    """Deconvolve an elliptical Gaussian beam from an observed Gaussian.

    Parameters
    ----------
    observed_major, observed_minor
        Observed major and minor FWHM in any common angular or linear unit.
    observed_pa_deg
        Observed position angle in degrees east of north.
    beam_major, beam_minor
        Beam major and minor FWHM in the same unit as the observed axes.
    beam_pa_deg
        Beam position angle in degrees east of north.
    absolute_tolerance, relative_tolerance
        Numerical tolerances applied to the covariance eigenvalues.
    on_unresolved
        Return ``(nan, nan, nan)`` or raise :class:`UnresolvedSourceError`
        when the deconvolution is non-physical.

    Returns
    -------
    source_major, source_minor, source_pa_deg
        Deconvolved FWHM axes in the input unit and PA normalized to
        ``[0, 180)`` degrees.
    """
    observed_major, observed_minor, observed_pa_deg = _validate_ellipse(
        observed_major,
        observed_minor,
        observed_pa_deg,
        name="observed",
    )
    beam_major, beam_minor, beam_pa_deg = _validate_ellipse(
        beam_major,
        beam_minor,
        beam_pa_deg,
        name="beam",
    )
    if absolute_tolerance < 0 or relative_tolerance < 0:
        raise ValueError("tolerances must be non-negative")
    if on_unresolved not in {"nan", "raise"}:
        raise ValueError("on_unresolved must be 'nan' or 'raise'")

    observed_covariance = _ellipse_to_covariance(
        observed_major, observed_minor, observed_pa_deg
    )
    beam_covariance = _ellipse_to_covariance(
        beam_major, beam_minor, beam_pa_deg
    )
    source_covariance = observed_covariance - beam_covariance
    source_covariance = 0.5 * (source_covariance + source_covariance.T)

    eigenvalues = np.linalg.eigvalsh(source_covariance)
    scale = max(
        float(np.max(np.abs(observed_covariance))),
        float(np.max(np.abs(beam_covariance))),
        1.0,
    )
    threshold = absolute_tolerance + relative_tolerance * scale

    if eigenvalues[1] <= threshold:
        return _unresolved(
            "The observed Gaussian is not larger than the beam",
            on_unresolved,
        )
    if eigenvalues[0] < -threshold:
        return _unresolved(
            "Beam subtraction produced a negative source variance",
            on_unresolved,
        )

    if eigenvalues[0] < 0:
        eigenvectors = np.linalg.eigh(source_covariance)[1]
        source_covariance = (
            eigenvectors
            @ np.diag([0.0, eigenvalues[1]])
            @ eigenvectors.T
        )

    return _covariance_to_ellipse(
        source_covariance,
        tolerance=relative_tolerance,
    )


def convolve_gaussians(
    first_major: float,
    first_minor: float,
    first_pa_deg: float,
    second_major: float,
    second_minor: float,
    second_pa_deg: float,
) -> tuple[float, float, float]:
    """Convolve two elliptical Gaussians using covariance addition."""
    first_major, first_minor, first_pa_deg = _validate_ellipse(
        first_major,
        first_minor,
        first_pa_deg,
        name="first",
    )
    second_major, second_minor, second_pa_deg = _validate_ellipse(
        second_major,
        second_minor,
        second_pa_deg,
        name="second",
    )
    covariance = _ellipse_to_covariance(
        first_major, first_minor, first_pa_deg
    ) + _ellipse_to_covariance(second_major, second_minor, second_pa_deg)
    return _covariance_to_ellipse(covariance, tolerance=1e-10)
