"""Geometry for drawing catalogue sources on continuum maps."""

from __future__ import annotations

from dataclasses import dataclass

from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
import numpy as np


@dataclass(frozen=True)
class EllipseGeometry:
    """A source ellipse in pixel coordinates, ready for matplotlib's Ellipse."""

    x: float
    y: float
    width_px: float
    height_px: float
    angle_deg: float


def _north_up_pixel_arcsec(wcs: WCS) -> float:
    celestial = wcs.celestial
    matrix = celestial.wcs.get_pc()
    if not np.allclose(matrix, np.diag(np.diag(matrix))):
        raise ValueError("Source ellipses need a north-up grid without rotation")
    scales = proj_plane_pixel_scales(celestial) * 3600.0
    if not np.isclose(scales[0], scales[1], rtol=1e-6):
        raise ValueError("Source ellipses need square pixels")
    return float(scales[0])


def ellipse_geometry(
    ra_deg: float,
    dec_deg: float,
    major_arcsec: float,
    minor_arcsec: float,
    position_angle_deg: float,
    wcs: WCS,
    *,
    scale: float = 1.0,
) -> EllipseGeometry:
    """Pixel ellipse of a source with FWHM axes and a position angle measured
    east of north, on a north-up grid (east to the left). matplotlib measures
    the angle counter-clockwise from +x, so north (+y) is 90 degrees."""
    pixel_arcsec = _north_up_pixel_arcsec(wcs)
    x, y = wcs.celestial.world_to_pixel_values(ra_deg, dec_deg)
    return EllipseGeometry(
        x=float(x),
        y=float(y),
        width_px=scale * float(major_arcsec) / pixel_arcsec,
        height_px=scale * float(minor_arcsec) / pixel_arcsec,
        angle_deg=90.0 + float(position_angle_deg),
    )


def box_window(
    wcs: WCS,
    ra_deg: float,
    dec_deg: float,
    size_deg: float,
    *,
    shape: tuple[int, int],
) -> tuple[slice, slice]:
    """Row and column slices of a square box centred on a sky position,
    clipped to an image of the given shape."""
    pixel_arcsec = _north_up_pixel_arcsec(wcs)
    half = int(round(size_deg * 3600.0 / pixel_arcsec / 2.0))
    x, y = wcs.celestial.world_to_pixel_values(ra_deg, dec_deg)
    col, row = int(round(float(x))), int(round(float(y)))
    rows = slice(max(row - half, 0), min(row + half, shape[0]))
    cols = slice(max(col - half, 0), min(col + half, shape[1]))
    return rows, cols
