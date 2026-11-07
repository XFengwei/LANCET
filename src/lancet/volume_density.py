"""Resolution-aware safeguards for volume-density-mapper products."""

from __future__ import annotations

import numpy as np
from astropy import constants as const


_GAUSSIAN_THICKNESS_FACTOR = 2.0 * np.sqrt(np.pi) / np.sqrt(8.0 * np.log(2.0))


def angular_size_cm(angle_arcsec: np.ndarray | float, distance_pc: float) -> np.ndarray:
    """Return the small-angle physical size in centimetres."""
    angle = np.asarray(angle_arcsec, dtype=float)
    distance = float(distance_pc)
    if not np.isfinite(distance) or distance <= 0:
        raise ValueError("Distance must be finite and positive")
    if np.any(~np.isfinite(angle)) or np.any(angle <= 0):
        raise ValueError("Angular sizes must be finite and positive")
    return angle * distance * const.au.cgs.value


def density_from_width(
    surface_density_g_cm2: np.ndarray,
    width_cm: np.ndarray,
) -> np.ndarray:
    """Convert surface density to mass density using the mapper convention.

    Non-positive or non-finite widths remain NaN instead of being converted to
    the largest representable float, as upstream ``numpy.nan_to_num`` does.
    """
    surface = np.asarray(surface_density_g_cm2, dtype=float)
    width = np.asarray(width_cm, dtype=float)
    if surface.shape != width.shape:
        raise ValueError("Surface-density and width arrays must have matching shapes")
    density = np.full(surface.shape, np.nan, dtype=float)
    valid = np.isfinite(surface) & (surface >= 0) & np.isfinite(width) & (width > 0)
    density[valid] = surface[valid] / (width[valid] * _GAUSSIAN_THICKNESS_FACTOR)
    return density


def apply_resolution_floor(
    width_cm: np.ndarray,
    resolution_arcsec: np.ndarray,
    *,
    distance_pc: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Floor characteristic widths at the local angular resolution.

    Returns the effective width, a flag identifying a valid upstream width,
    and a flag identifying pixels where the resolution floor was applied.
    """
    width = np.asarray(width_cm, dtype=float)
    resolution = np.asarray(resolution_arcsec, dtype=float)
    if width.shape != resolution.shape:
        raise ValueError("Width and resolution arrays must have matching shapes")
    beam_cm = angular_size_cm(resolution, distance_pc)
    upstream_defined = np.isfinite(width) & (width > 0)
    effective = np.where(upstream_defined, np.maximum(width, beam_cm), beam_cm)
    floor_applied = ~upstream_defined | (width < beam_cm)
    return effective, upstream_defined, floor_applied


def number_density_h2(
    mass_density_g_cm3: np.ndarray,
    *,
    mean_molecular_weight: float = 2.8,
) -> np.ndarray:
    """Convert mass density to H2 number density in cm^-3."""
    molecular_weight = float(mean_molecular_weight)
    if not np.isfinite(molecular_weight) or molecular_weight <= 0:
        raise ValueError("Mean molecular weight must be finite and positive")
    return np.asarray(mass_density_g_cm3, dtype=float) / (
        molecular_weight * const.m_p.cgs.value
    )


def cube_storage_bytes(shape_yx: tuple[int, int], *, dtype_bytes: int = 4) -> int:
    """Return storage required by the mapper's ``(max(ny,nx), ny, nx)`` cube."""
    ny, nx = (int(value) for value in shape_yx)
    if ny <= 0 or nx <= 0 or dtype_bytes <= 0:
        raise ValueError("Shape and dtype size must be positive")
    return max(ny, nx) * ny * nx * int(dtype_bytes)

