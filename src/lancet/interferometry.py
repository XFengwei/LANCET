"""Utilities for interferometric coverage and image comparison."""

from __future__ import annotations

from pathlib import Path

from astropy import units as u
from astropy.convolution import convolve_fft
from astropy.wcs import WCS
import numpy as np
from radio_beam import Beam
from reproject import reproject_interp


def deduplicate_uv_coordinates(
    u_metres: np.ndarray,
    v_metres: np.ndarray,
    *,
    decimals: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """Keep one UV coordinate for each rounded ``(u, v)`` pair."""
    u_values = np.asarray(u_metres, dtype=float)
    v_values = np.asarray(v_metres, dtype=float)
    if u_values.ndim != 1 or v_values.ndim != 1:
        raise ValueError("UV coordinates must be one-dimensional")
    if u_values.shape != v_values.shape:
        raise ValueError("U and V coordinate shapes must match")
    if decimals < 0:
        raise ValueError("decimals must be non-negative")
    coordinates = np.column_stack([u_values, v_values])
    rounded = np.round(coordinates, decimals=decimals)
    _, indices = np.unique(rounded, axis=0, return_index=True)
    indices.sort()
    return u_values[indices], v_values[indices]


def read_unique_uv_coordinates(
    measurement_set: str | Path,
    *,
    decimals: int = 1,
    chunk_rows: int = 250_000,
) -> tuple[np.ndarray, np.ndarray]:
    """Read and deduplicate the UVW column without loading all rows at once."""
    path = Path(measurement_set)
    if not path.is_dir():
        raise FileNotFoundError(f"Measurement set not found: {path}")
    if chunk_rows <= 0:
        raise ValueError("chunk_rows must be positive")

    # casatools creates a timestamped CASA log in the current directory at
    # import time unless its configuration is changed first.
    from casaconfig import config as casa_config

    # disable the log file
    casa_config.nologfile = True
    casa_config.logfile = None
    from casatools import table

    handle = table()
    handle.open(str(path), nomodify=True)
    try:
        number_of_rows = int(handle.nrows())
        chunks: list[np.ndarray] = []
        for start in range(0, number_of_rows, chunk_rows):
            row_count = min(chunk_rows, number_of_rows - start)
            uvw = np.asarray(
                handle.getcol(
                    "UVW",
                    startrow=start,
                    nrow=row_count,
                )
            )
            if uvw.ndim != 2 or uvw.shape[0] < 2:
                raise ValueError(
                    f"Unexpected UVW shape in {path}: {uvw.shape}"
                )
            u_unique, v_unique = deduplicate_uv_coordinates(
                uvw[0],
                uvw[1],
                decimals=decimals,
            )
            chunks.append(np.column_stack([u_unique, v_unique]))
    finally:
        handle.close()

    if not chunks:
        return np.array([], dtype=float), np.array([], dtype=float)
    combined = np.concatenate(chunks, axis=0)
    return deduplicate_uv_coordinates(
        combined[:, 0],
        combined[:, 1],
        decimals=decimals,
    )


def beam_match_and_reproject(
    image: np.ndarray,
    image_header,
    reference_image: np.ndarray,
    reference_header,
) -> tuple[np.ndarray, object]:
    """Convolve an image to a reference beam and reproject to its WCS grid."""
    data = np.squeeze(np.asarray(image, dtype=float))
    reference = np.squeeze(np.asarray(reference_image))
    if data.ndim != 2 or reference.ndim != 2:
        raise ValueError("Input and reference images must be two-dimensional")

    current_beam = Beam.from_fits_header(image_header)
    target_beam = Beam.from_fits_header(reference_header)
    try:
        kernel_beam = target_beam.deconvolve(current_beam)
    except ValueError as error:
        raise ValueError(
            "Reference beam must be larger than the input beam"
        ) from error

    pixel_scale = abs(float(image_header["CDELT2"])) * u.deg
    kernel = kernel_beam.as_kernel(pixel_scale)
    convolved = convolve_fft(
        data,
        kernel,
        boundary="fill",
        fill_value=np.nan,
        nan_treatment="interpolate",
        normalize_kernel=True,
        preserve_nan=True,
        allow_huge=True,
    )
    reprojected, footprint = reproject_interp(
        (convolved, WCS(image_header).celestial),
        WCS(reference_header).celestial,
        shape_out=reference.shape,
    )
    reprojected[footprint <= 0] = np.nan

    output_header = reference_header.copy()
    output_header["BMAJ"] = target_beam.major.to_value(u.deg)
    output_header["BMIN"] = target_beam.minor.to_value(u.deg)
    output_header["BPA"] = target_beam.pa.to_value(u.deg)
    output_header["BUNIT"] = image_header.get("BUNIT", "")
    output_header.add_history(
        "LANCET: convolved to reference beam and reprojected to reference WCS"
    )
    return reprojected, output_header
