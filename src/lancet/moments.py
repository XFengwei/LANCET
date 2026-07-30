"""Spectral-cube preparation and moment-map reprojection."""

from __future__ import annotations

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from reproject import reproject_interp


def reproject_image(
    data: np.ndarray,
    source_header: fits.Header,
    target_header: fits.Header,
    *,
    order: str = "bilinear",
) -> tuple[np.ndarray, np.ndarray, fits.Header]:
    """Reproject a 2D image and preserve its physical metadata."""
    image = np.squeeze(np.asarray(data))
    if image.ndim != 2:
        raise ValueError("Input image must be two-dimensional")
    shape_out = (
        int(target_header["NAXIS2"]),
        int(target_header["NAXIS1"]),
    )
    reprojected, footprint = reproject_interp(
        (image, WCS(source_header).celestial),
        output_projection=WCS(target_header).celestial,
        shape_out=shape_out,
        order=order,
        parallel=True,
    )
    reprojected[footprint <= 0] = np.nan
    output_header = target_header.copy()
    for key in ("BUNIT", "BMAJ", "BMIN", "BPA"):
        if key in source_header:
            output_header[key] = source_header[key]
    if "BLANK" in output_header:
        del output_header["BLANK"]
    output_header.add_history(
        "LANCET c18o_moments: reprojected to delivered ALMA continuum grid"
    )
    return reprojected, footprint, output_header
