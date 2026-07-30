"""Source-catalogue selection, sampling, and size measurements."""

from __future__ import annotations

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.wcs import WCS
import numpy as np
from regions import (
    PixCoord,
    RectanglePixelRegion,
    RectangleSkyRegion,
    Regions,
)


def select_science_sources(catalogue: Table) -> Table:
    """Apply the adopted two-scale getsf reliability and size cuts."""
    flux = (catalogue["FXT_BST01"] > 0) | (
        catalogue["FXT_BST02"] > 0
    )
    peak = (catalogue["FXP_BST01"] > 0) | (
        catalogue["FXP_BST02"] > 0
    )
    significance = (catalogue["SIGNM01"] >= 5) | (
        catalogue["SIGNM02"] >= 5
    )
    compact = (
        (catalogue["AFWHM01"] <= 3.0)
        & (catalogue["BFWHM01"] <= 3.0)
    ) | (
        (catalogue["AFWHM02"] <= 3.0)
        & (catalogue["BFWHM02"] <= 3.0)
    )
    return catalogue[flux & peak & significance & compact]


def peak_to_total_flux_ratio(
    peak_flux: np.ndarray,
    peak_error: np.ndarray,
    total_flux: np.ndarray,
    total_error: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return peak/total flux and its independent-error propagation."""
    peak = np.asarray(peak_flux, dtype=float)
    peak_uncertainty = np.asarray(peak_error, dtype=float)
    total = np.asarray(total_flux, dtype=float)
    total_uncertainty = np.asarray(total_error, dtype=float)
    if not (
        peak.shape
        == peak_uncertainty.shape
        == total.shape
        == total_uncertainty.shape
    ):
        raise ValueError("Flux values and errors must have matching shapes")
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = peak / total
        uncertainty = np.sqrt(
            (peak_uncertainty / total) ** 2
            + (total_uncertainty * peak / total**2) ** 2
        )
    return ratio, uncertainty


def relative_flux_difference(
    reference_flux: np.ndarray,
    comparison_flux: np.ndarray,
) -> np.ndarray:
    """Return ``(comparison - reference) / reference``."""
    reference = np.asarray(reference_flux, dtype=float)
    comparison = np.asarray(comparison_flux, dtype=float)
    if reference.shape != comparison.shape:
        raise ValueError("Reference and comparison flux shapes must match")
    with np.errstate(divide="ignore", invalid="ignore"):
        return (comparison - reference) / reference


def uniform_brightness_flux_ratio(
    source_fwhm_arcsec: np.ndarray,
    *,
    beam_fwhm_arcsec: float,
) -> np.ndarray:
    """Return the peak/total ratio of a circular uniform-brightness source."""
    size = np.asarray(source_fwhm_arcsec, dtype=float)
    beam = float(beam_fwhm_arcsec)
    if not np.isfinite(beam) or beam <= 0:
        raise ValueError("Beam FWHM must be finite and positive")
    if np.any(np.isfinite(size) & (size <= 0)):
        raise ValueError("Finite source FWHM values must be positive")
    return (beam / size) ** 2


def getsf_peak_indices(
    x_coordinates: np.ndarray,
    y_coordinates: np.ndarray,
    *,
    shape: tuple[int, int],
) -> tuple[np.ndarray, np.ndarray]:
    """Convert getsf peak coordinates to the adopted NumPy pixel indices."""
    height, width = shape
    x = np.asarray(x_coordinates, dtype=float)
    y = np.asarray(y_coordinates, dtype=float)
    x_index = np.where(x < width / 2.0, x.astype(int), x.astype(int) + 1)
    y_index = np.where(
        y < height / 2.0,
        y.astype(int),
        y.astype(int) + 1,
    )
    if (
        np.any(x_index < 0)
        or np.any(x_index >= width)
        or np.any(y_index < 0)
        or np.any(y_index >= height)
    ):
        raise IndexError("At least one getsf peak falls outside the map")
    return x_index, y_index


def sample_at_getsf_peaks(
    image: np.ndarray,
    x_coordinates: np.ndarray,
    y_coordinates: np.ndarray,
) -> np.ndarray:
    """Sample an image at getsf peak coordinates."""
    data = np.squeeze(np.asarray(image))
    if data.ndim != 2:
        raise ValueError("Sampled image must be two-dimensional")
    x_index, y_index = getsf_peak_indices(
        x_coordinates,
        y_coordinates,
        shape=data.shape,
    )
    return np.asarray(data[y_index, x_index])


def effective_source_sizes(
    observed_major_arcsec: np.ndarray,
    observed_minor_arcsec: np.ndarray,
    *,
    beam_major_arcsec: float,
    beam_minor_arcsec: float,
    distance_pc: float,
    unresolved_threshold_factor: float = 1.25,
    unresolved_size_factor: float = 0.5,
) -> dict[str, np.ndarray]:
    """Calculate the adopted circularized deconvolved source sizes."""
    major = np.asarray(observed_major_arcsec, dtype=float)
    minor = np.asarray(observed_minor_arcsec, dtype=float)
    beam_area = float(beam_major_arcsec) * float(beam_minor_arcsec)
    observed_area = major * minor
    resolved = (
        observed_area - unresolved_threshold_factor * beam_area
    ) >= 0
    size = np.where(
        resolved,
        np.sqrt(np.maximum(observed_area - beam_area, 0.0)),
        unresolved_size_factor * np.sqrt(beam_area),
    )
    return {
        "major_arcsec": size,
        "minor_arcsec": size.copy(),
        "position_angle_deg": np.full(size.shape, np.nan),
        "fwhm_arcsec": size.copy(),
        "fwhm_au": size * float(distance_pc),
        "resolved": np.where(resolved, "Y", "N"),
    }


def rectangle_membership(
    catalogue: Table,
    region_path: str,
    *,
    wcs: WCS,
    ra_column: str = "WCS_ACOOR",
    dec_column: str = "WCS_DCOOR",
) -> np.ndarray:
    """Return membership in the first rectangular DS9 region."""
    coordinates = SkyCoord(
        np.asarray(catalogue[ra_column], dtype=float) * u.deg,
        np.asarray(catalogue[dec_column], dtype=float) * u.deg,
        frame="icrs",
    )
    rectangle = next(
        (
            region
            for region in Regions.read(region_path, format="ds9")
            if isinstance(
                region,
                (RectangleSkyRegion, RectanglePixelRegion),
            )
        ),
        None,
    )
    if rectangle is None:
        raise ValueError(f"No rectangle found in region file: {region_path}")
    if isinstance(rectangle, RectanglePixelRegion):
        x, y = coordinates.to_pixel(wcs)
        membership = rectangle.contains(PixCoord(x=x, y=y))
    else:
        membership = rectangle.contains(coordinates, wcs=wcs)
    return np.asarray(membership, dtype=bool)
