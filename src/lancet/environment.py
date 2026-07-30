"""Aperture measurements of core environments."""

from __future__ import annotations

from collections.abc import Mapping

from astropy.table import Table
from astropy.wcs import WCS
import numpy as np
from photutils.aperture import ApertureStats, SkyCircularAperture

from .continuum_plotting import image_plane


ENVIRONMENT_STATISTICS = (
    "mean",
    "median",
    "std",
    "mad",
    "validfrac",
)


def measure_aperture_statistics(
    data: np.ndarray,
    apertures: SkyCircularAperture,
    wcs: WCS,
    *,
    minimum_valid_value: float = 0.0,
) -> dict[str, np.ndarray]:
    """Measure positive, finite map values inside sky apertures."""
    image = image_plane(data)
    invalid = (~np.isfinite(image)) | (image <= minimum_valid_value)
    statistics = ApertureStats(
        image,
        apertures,
        wcs=wcs.celestial,
        mask=invalid,
    )
    coverage = ApertureStats(
        (~invalid).astype(float),
        apertures,
        wcs=wcs.celestial,
    )
    valid_count = np.atleast_1d(
        np.asarray(coverage.sum, dtype=float)
    )
    total_count = np.atleast_1d(
        np.asarray(coverage.sum_aper_area.value, dtype=float)
    )
    return {
        "mean": np.atleast_1d(
            np.asarray(statistics.mean, dtype=float)
        ),
        "median": np.atleast_1d(
            np.asarray(statistics.median, dtype=float)
        ),
        "std": np.atleast_1d(
            np.asarray(statistics.std, dtype=float)
        ),
        "mad": np.atleast_1d(
            np.asarray(statistics.mad_std, dtype=float)
        ),
        "validfrac": valid_count / total_count,
    }


def add_environment_measurements(
    catalogue: Table,
    measurements: Mapping[str, Mapping[str, np.ndarray]],
    *,
    units: Mapping[str, str | None] | None = None,
) -> Table:
    """Return a catalogue copy with named environment statistics appended."""
    output = catalogue.copy(copy_data=True)
    unit_map = {} if units is None else units
    for prefix, values in measurements.items():
        missing = set(ENVIRONMENT_STATISTICS) - set(values)
        if missing:
            raise KeyError(
                f"{prefix} measurements are missing: "
                + ", ".join(sorted(missing))
            )
        for statistic in ENVIRONMENT_STATISTICS:
            array = np.asarray(values[statistic], dtype=float)
            if array.shape != (len(output),):
                raise ValueError(
                    f"{prefix} {statistic} has shape {array.shape}; "
                    f"expected {(len(output),)}"
                )
            name = f"{prefix}_env_{statistic}"
            output[name] = array
            if statistic != "validfrac" and unit_map.get(prefix):
                output[name].unit = unit_map[prefix]
    return output
