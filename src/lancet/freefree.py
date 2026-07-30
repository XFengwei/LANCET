"""Free-free emission models and pixel-wise continuum fitting."""

from __future__ import annotations

from collections.abc import Callable, Mapping

from astropy import constants as const
import numpy as np
from scipy.optimize import curve_fit


def beam_solid_angle(
    major_fwhm_arcsec: float,
    minor_fwhm_arcsec: float,
) -> float:
    """Return the solid angle of a Gaussian beam in steradians."""
    major = float(major_fwhm_arcsec)
    minor = float(minor_fwhm_arcsec)
    if not np.all(np.isfinite([major, minor])) or major <= 0 or minor <= 0:
        raise ValueError("Beam FWHM axes must be positive and finite")

    arcsec_to_radian = np.pi / (180.0 * 3600.0)
    return (
        np.pi
        * (major * arcsec_to_radian)
        * (minor * arcsec_to_radian)
        / (4.0 * np.log(2.0))
    )


def freefree_optical_depth(
    frequency_ghz: np.ndarray | float,
    emission_measure: np.ndarray | float,
    *,
    electron_temperature_k: float = 8000.0,
) -> np.ndarray:
    """Calculate free-free optical depth.

    Parameters
    ----------
    frequency_ghz
        Frequency in GHz.
    emission_measure
        Emission measure in pc cm^-6.
    electron_temperature_k
        Electron temperature in K.
    """
    frequency = np.asarray(frequency_ghz, dtype=float)
    emission_measure = np.asarray(emission_measure, dtype=float)
    temperature = float(electron_temperature_k)

    if np.any(~np.isfinite(frequency)) or np.any(frequency <= 0):
        raise ValueError("Frequencies must be positive and finite")
    if np.any(~np.isfinite(emission_measure)) or np.any(emission_measure < 0):
        raise ValueError("Emission measure must be non-negative and finite")
    if not np.isfinite(temperature) or temperature <= 0:
        raise ValueError("Electron temperature must be positive and finite")

    return (
        0.08235
        * temperature**-1.35
        * frequency**-2.1
        * emission_measure
    )


def freefree_flux_density(
    frequency_ghz: np.ndarray | float,
    emission_measure: np.ndarray | float,
    *,
    electron_temperature_k: float = 8000.0,
    beam_solid_angle_sr: float,
) -> np.ndarray:
    """Calculate free-free flux density in Jy/beam."""
    frequency = np.asarray(frequency_ghz, dtype=float)
    omega = float(beam_solid_angle_sr)
    if not np.isfinite(omega) or omega <= 0:
        raise ValueError("Beam solid angle must be positive and finite")

    optical_depth = freefree_optical_depth(
        frequency,
        emission_measure,
        electron_temperature_k=electron_temperature_k,
    )
    brightness_temperature = float(electron_temperature_k) * (
        1.0 - np.exp(-optical_depth)
    )
    frequency_hz = frequency * 1e9
    specific_intensity = (
        2.0
        * const.k_B.value
        * brightness_temperature
        * frequency_hz**2
        / const.c.value**2
    )
    return specific_intensity * omega / 1e-26


def fractional_flux_difference(
    comparison_image: np.ndarray,
    reference_image: np.ndarray,
    *,
    absolute: bool = False,
) -> np.ndarray:
    """Return ``(comparison - reference) / reference`` pixel by pixel.

    Pixels with non-finite inputs or a zero reference value are returned as
    NaN. Set ``absolute=True`` for the magnitude of the fractional change.
    """
    comparison = np.asarray(comparison_image, dtype=float)
    reference = np.asarray(reference_image, dtype=float)
    if comparison.shape != reference.shape:
        raise ValueError(
            f"Comparison shape {comparison.shape} does not match "
            f"reference shape {reference.shape}"
        )
    with np.errstate(divide="ignore", invalid="ignore"):
        difference = (comparison - reference) / reference
    valid = (
        np.isfinite(comparison)
        & np.isfinite(reference)
        & (reference != 0)
        & np.isfinite(difference)
    )
    difference = np.where(valid, difference, np.nan)
    return np.abs(difference) if absolute else difference


def fit_freefree_pixel(
    frequency_ghz: np.ndarray,
    flux_jy_per_beam: np.ndarray,
    uncertainty_jy_per_beam: np.ndarray,
    *,
    electron_temperature_k: float = 8000.0,
    beam_solid_angle_sr: float,
) -> tuple[float, float, float, float, np.ndarray]:
    """Fit one spectrum with a fixed-temperature free-free model."""
    frequency = np.asarray(frequency_ghz, dtype=float)
    flux = np.asarray(flux_jy_per_beam, dtype=float)
    uncertainty = np.asarray(uncertainty_jy_per_beam, dtype=float)

    if frequency.ndim != 1:
        raise ValueError("frequency_ghz must be one-dimensional")
    if flux.shape != frequency.shape or uncertainty.shape != frequency.shape:
        raise ValueError("Frequency, flux, and uncertainty shapes must match")
    if len(frequency) < 2:
        raise ValueError("At least two bands are required for a fit")
    if np.any(~np.isfinite(flux)):
        raise ValueError("Flux values must be finite")
    if np.any(~np.isfinite(uncertainty)) or np.any(uncertainty <= 0):
        raise ValueError("Uncertainties must be positive and finite")

    def model(nu_ghz: np.ndarray, emission_measure: float) -> np.ndarray:
        return freefree_flux_density(
            nu_ghz,
            emission_measure,
            electron_temperature_k=electron_temperature_k,
            beam_solid_angle_sr=beam_solid_angle_sr,
        )

    parameters, covariance = curve_fit(
        model,
        frequency,
        flux,
        sigma=uncertainty,
        absolute_sigma=True,
        p0=[1e5],
        bounds=(0.0, np.inf),
        maxfev=20_000,
    )
    emission_measure = float(parameters[0])
    emission_measure_error = (
        float(np.sqrt(covariance[0, 0]))
        if covariance is not None
        and covariance.size
        and covariance[0, 0] >= 0
        else np.nan
    )
    model_flux = model(frequency, emission_measure)
    chi_squared = float(
        np.sum(((flux - model_flux) / uncertainty) ** 2)
    )
    degrees_of_freedom = max(len(flux) - 1, 1)
    reduced_chi_squared = chi_squared / degrees_of_freedom

    return (
        emission_measure,
        emission_measure_error,
        chi_squared,
        reduced_chi_squared,
        model_flux,
    )


def meerkat_subband_frequencies(
    header: Mapping[str, object],
    *,
    number_of_subbands: int = 14,
) -> np.ndarray:
    """Read exact MeerKAT subband frequencies from FITS header cards."""
    if number_of_subbands <= 0:
        raise ValueError("number_of_subbands must be positive")

    keys = [f"FREQ{index:04d}" for index in range(1, number_of_subbands + 1)]
    missing = [key for key in keys if key not in header]
    if missing:
        raise KeyError(
            "Missing MeerKAT subband frequency cards: " + ", ".join(missing)
        )

    frequency_hz = np.asarray([header[key] for key in keys], dtype=float)
    if np.any(~np.isfinite(frequency_hz)) or np.any(frequency_hz <= 0):
        raise ValueError("MeerKAT subband frequencies must be positive and finite")
    if np.any(np.diff(frequency_hz) <= 0):
        raise ValueError("MeerKAT subband frequencies must be increasing")
    return frequency_hz / 1e9


def prepare_sed_cube(
    meerkat_cube: np.ndarray,
    cornish_map: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Combine 14 MeerKAT subbands and a CORNISH map into an SED cube."""
    cube = np.asarray(meerkat_cube)
    cornish = np.squeeze(np.asarray(cornish_map))

    if cube.ndim != 4 or cube.shape[0] != 1 or cube.shape[1] != 16:
        raise ValueError(
            "Expected MeerKAT cube shape (1, 16, ny, nx), "
            f"received {cube.shape}"
        )
    if cornish.ndim != 2:
        raise ValueError(
            f"Expected a two-dimensional CORNISH map, received {cornish.shape}"
        )

    mfs_map = np.asarray(cube[0, 0], dtype=float)
    spectral_index_map = np.asarray(cube[0, 1], dtype=float)
    meerkat_subbands = np.asarray(cube[0, 2:16], dtype=float)
    if cornish.shape != mfs_map.shape:
        raise ValueError(
            f"CORNISH shape {cornish.shape} does not match MeerKAT "
            f"shape {mfs_map.shape}"
        )

    sed_cube = np.concatenate(
        [meerkat_subbands, cornish[np.newaxis, :, :]],
        axis=0,
    )
    return mfs_map, spectral_index_map, sed_cube


def fit_freefree_map(
    sed_cube: np.ndarray,
    frequency_ghz: np.ndarray,
    rms_jy_per_beam: np.ndarray,
    *,
    beam_solid_angle_sr: float,
    electron_temperature_k: float = 8000.0,
    fit_mask: np.ndarray | None = None,
    minimum_bands: int = 6,
    snr_threshold: float | None = None,
    save_model_cube: bool = True,
    progress: Callable[[int, int], None] | None = None,
) -> dict[str, np.ndarray]:
    """Fit a fixed-temperature free-free model independently in each pixel."""
    sed = np.asarray(sed_cube, dtype=float)
    frequency = np.asarray(frequency_ghz, dtype=float)
    rms = np.asarray(rms_jy_per_beam, dtype=float)

    if sed.ndim != 3:
        raise ValueError("sed_cube must have shape (nband, ny, nx)")
    number_of_bands, ny, nx = sed.shape
    if frequency.shape != (number_of_bands,):
        raise ValueError("frequency_ghz length does not match sed_cube")
    if rms.shape != (number_of_bands,):
        raise ValueError("rms_jy_per_beam length does not match sed_cube")
    if np.any(~np.isfinite(rms)) or np.any(rms <= 0):
        raise ValueError("All RMS values must be positive and finite")
    if not 2 <= minimum_bands <= number_of_bands:
        raise ValueError("minimum_bands must be between 2 and nband")
    if snr_threshold is not None and (
        not np.isfinite(snr_threshold) or snr_threshold <= 0
    ):
        raise ValueError("snr_threshold must be positive and finite")

    if fit_mask is None:
        mask = np.ones((ny, nx), dtype=bool)
    else:
        mask = np.asarray(fit_mask, dtype=bool)
        if mask.shape != (ny, nx):
            raise ValueError("fit_mask shape does not match the image")

    outputs: dict[str, np.ndarray] = {
        "EM_map": np.full((ny, nx), np.nan, dtype=np.float32),
        "EMerr_map": np.full((ny, nx), np.nan, dtype=np.float32),
        "ff230_map": np.full((ny, nx), np.nan, dtype=np.float32),
        "chi2_map": np.full((ny, nx), np.nan, dtype=np.float32),
        "rchi2_map": np.full((ny, nx), np.nan, dtype=np.float32),
        "nband_map": np.zeros((ny, nx), dtype=np.int16),
        "fitmask_map": mask.astype(np.int16),
    }
    if save_model_cube:
        outputs["model_cube"] = np.full(
            sed.shape,
            np.nan,
            dtype=np.float32,
        )

    y_pixels, x_pixels = np.where(mask)
    total_pixels = len(y_pixels)
    for index, (iy, ix) in enumerate(
        zip(y_pixels, x_pixels, strict=True),
        start=1,
    ):
        spectrum = sed[:, iy, ix]
        valid = np.isfinite(spectrum) & (spectrum > 0)
        if snr_threshold is not None:
            valid &= spectrum > snr_threshold * rms
        if np.count_nonzero(valid) < minimum_bands:
            continue

        try:
            (
                emission_measure,
                emission_measure_error,
                chi_squared,
                reduced_chi_squared,
                model_flux,
            ) = fit_freefree_pixel(
                frequency[valid],
                spectrum[valid],
                rms[valid],
                electron_temperature_k=electron_temperature_k,
                beam_solid_angle_sr=beam_solid_angle_sr,
            )
        except (RuntimeError, ValueError, FloatingPointError):
            continue

        outputs["EM_map"][iy, ix] = emission_measure
        outputs["EMerr_map"][iy, ix] = emission_measure_error
        outputs["ff230_map"][iy, ix] = freefree_flux_density(
            230.0,
            emission_measure,
            electron_temperature_k=electron_temperature_k,
            beam_solid_angle_sr=beam_solid_angle_sr,
        )
        outputs["chi2_map"][iy, ix] = chi_squared
        outputs["rchi2_map"][iy, ix] = reduced_chi_squared
        outputs["nband_map"][iy, ix] = np.count_nonzero(valid)
        if save_model_cube:
            outputs["model_cube"][valid, iy, ix] = model_flux

        if progress is not None and index % 5000 == 0:
            progress(index, total_pixels)

    return outputs
