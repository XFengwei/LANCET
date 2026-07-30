"""Reusable operations for continuum-image preparation."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np


def gaussian_beam_solid_angle_sr(
    major_fwhm_deg: float,
    minor_fwhm_deg: float,
) -> float:
    """Return the solid angle of a Gaussian beam in steradians."""
    major = float(major_fwhm_deg)
    minor = float(minor_fwhm_deg)
    if not np.all(np.isfinite([major, minor])) or major <= 0 or minor <= 0:
        raise ValueError("Beam FWHM axes must be positive and finite")
    major_radian = np.deg2rad(major)
    minor_radian = np.deg2rad(minor)
    return np.pi * major_radian * minor_radian / (4.0 * np.log(2.0))


def jybeam_to_mjysr(
    image: np.ndarray,
    *,
    major_fwhm_deg: float,
    minor_fwhm_deg: float,
) -> np.ndarray:
    """Convert an image from Jy/beam to MJy/sr."""
    beam_area = gaussian_beam_solid_angle_sr(
        major_fwhm_deg,
        minor_fwhm_deg,
    )
    return np.asarray(image) / beam_area * 1e-6


def jybeam_to_mjysr_from_header(
    image: np.ndarray,
    header: Mapping[str, object],
) -> np.ndarray:
    """Convert Jy/beam data using ``BMAJ`` and ``BMIN`` FITS cards."""
    missing = [key for key in ("BMAJ", "BMIN") if key not in header]
    if missing:
        raise KeyError("Missing FITS beam cards: " + ", ".join(missing))
    return jybeam_to_mjysr(
        image,
        major_fwhm_deg=float(header["BMAJ"]),
        minor_fwhm_deg=float(header["BMIN"]),
    )


def apply_primary_beam_mask(
    image: np.ndarray,
    primary_beam: np.ndarray,
    *,
    minimum_response: float = 0.2,
) -> np.ndarray:
    """Mask image pixels below a primary-beam response threshold."""
    data = np.squeeze(np.asarray(image))
    response = np.squeeze(np.asarray(primary_beam))
    if data.ndim != 2 or response.ndim != 2:
        raise ValueError("Image and primary beam must be two-dimensional")
    if data.shape != response.shape:
        raise ValueError(
            f"Image shape {data.shape} does not match primary-beam shape "
            f"{response.shape}"
        )
    threshold = float(minimum_response)
    if not np.isfinite(threshold) or not 0.0 <= threshold <= 1.0:
        raise ValueError("Primary-beam threshold must be between zero and one")
    return np.where(
        np.isfinite(response) & (response >= threshold),
        data,
        np.nan,
    )
