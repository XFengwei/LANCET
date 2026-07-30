#!/usr/bin/env python
"""Inspect, prepare, and fit the MeerKAT-CORNISH free-free model."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from pathlib import Path
import tomllib

from astropy.convolution import Gaussian2DKernel, convolve_fft
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from reproject import reproject_interp

from lancet import ProjectPaths
from lancet.continuum import smooth_circular_jybeam
from lancet.freefree import (
    beam_solid_angle,
    fit_freefree_map,
    fractional_flux_difference,
    meerkat_subband_frequencies,
    prepare_sed_cube,
)


DEFAULT_CONFIG = Path("tasks/construct_freefree_model/task.toml")
PREPARED_NATIVE = "CORNISH_5p5GHz_8arcsec_native.fits"
PREPARED_REPROJECTED = "CORNISH_5p5GHz_8arcsec_on_MeerKATgrid.fits"


class MissingRMSConfiguration(ValueError):
    """Raised when a fit is requested without measured per-band RMS values."""


def load_task_config(path: Path) -> dict[str, object]:
    """Read a task TOML file."""
    if not path.is_file():
        raise FileNotFoundError(f"Task configuration not found: {path}")
    with path.open("rb") as handle:
        return tomllib.load(handle)


def _project_path(project: ProjectPaths, value: object) -> Path:
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = project.root / path
    return path.resolve()


def resolve_task_paths(
    project: ProjectPaths,
    config: Mapping[str, object],
) -> dict[str, Path]:
    """Resolve configured inputs and output locations."""
    inputs = config["inputs"]
    outputs = config["outputs"]
    if not isinstance(inputs, Mapping) or not isinstance(outputs, Mapping):
        raise TypeError("Task inputs and outputs must be TOML tables")

    paths = {
        name: _project_path(project, value)
        for name, value in inputs.items()
    }
    paths["output_root"] = project.validate_output(str(outputs["root"]))

    for name in (
        "meerkat_cube",
        "meerkat_reference",
        "cornish_map",
        "dust_reference",
    ):
        path = paths[name]
        if not path.is_relative_to(project.delivered.resolve()):
            raise ValueError(f"{name} must be a delivered input: {path}")
        if not path.is_file():
            raise FileNotFoundError(f"{name} not found: {path}")
    return paths


def _image_description(path: Path) -> tuple[tuple[int, ...], str, str]:
    with fits.open(path, memmap=True) as hdul:
        shape = tuple(hdul[0].shape)
        header = hdul[0].header
        beam = (
            f"{float(header['BMAJ']) * 3600:g} x "
            f"{float(header['BMIN']) * 3600:g} arcsec"
            if "BMAJ" in header and "BMIN" in header
            else "not recorded"
        )
        unit = str(header.get("BUNIT", "not recorded"))
    return shape, beam, unit


def configured_rms(
    config: Mapping[str, object],
) -> np.ndarray:
    """Return the configured 15-band RMS vector or explain what is missing."""
    fit = config["fit"]
    if not isinstance(fit, Mapping):
        raise TypeError("Task fit settings must be a TOML table")

    meerkat = fit.get("rms_meerkat_jy_per_beam")
    cornish = fit.get("rms_cornish_jy_per_beam")
    if meerkat is None or cornish is None:
        raise MissingRMSConfiguration(
            "The fit stage requires measured values for "
            "fit.rms_meerkat_jy_per_beam (14 entries) and "
            "fit.rms_cornish_jy_per_beam in task.toml"
        )

    meerkat_values = np.asarray(meerkat, dtype=float)
    cornish_value = float(cornish)
    if meerkat_values.shape != (14,):
        raise MissingRMSConfiguration(
            "fit.rms_meerkat_jy_per_beam must contain exactly 14 values"
        )
    rms = np.concatenate([meerkat_values, [cornish_value]])
    if np.any(~np.isfinite(rms)) or np.any(rms <= 0):
        raise MissingRMSConfiguration(
            "All configured RMS values must be positive and finite"
        )
    return rms


def configured_frequencies(
    config: Mapping[str, object],
    header: Mapping[str, object],
) -> np.ndarray:
    """Return the 15 fitted frequencies selected by the task profile."""
    fit = config["fit"]
    if not isinstance(fit, Mapping):
        raise TypeError("Task fit settings must be a TOML table")

    source = str(fit.get("frequency_source", "fits_header"))
    if source == "fits_header":
        meerkat = meerkat_subband_frequencies(header)
    elif source == "configured":
        values = fit.get("meerkat_frequencies_ghz")
        if values is None:
            raise ValueError(
                "frequency_source='configured' requires "
                "fit.meerkat_frequencies_ghz"
            )
        meerkat = np.asarray(values, dtype=float)
        if meerkat.shape != (14,):
            raise ValueError(
                "fit.meerkat_frequencies_ghz must contain 14 values"
            )
        if np.any(~np.isfinite(meerkat)) or np.any(meerkat <= 0):
            raise ValueError(
                "Configured MeerKAT frequencies must be positive and finite"
            )
    else:
        raise ValueError(
            "fit.frequency_source must be 'fits_header' or 'configured'"
        )
    return np.concatenate([meerkat, [5.5]])


def inspect_inputs(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
) -> None:
    """Print the input geometry, frequencies, and fit readiness."""
    print("Declared delivered inputs:")
    for name in ("meerkat_cube", "meerkat_reference", "cornish_map"):
        path = paths[name]
        shape, beam, unit = _image_description(path)
        print(f"  {name}: {path}")
        print(f"    shape={shape}, beam={beam}, unit={unit}")

    with fits.open(paths["meerkat_cube"], memmap=True) as hdul:
        header_frequency = meerkat_subband_frequencies(hdul[0].header)
        fit_frequency = configured_frequencies(config, hdul[0].header)
    print("MeerKAT subband frequencies from FITS header [GHz]:")
    print("  " + ", ".join(f"{value:.9f}" for value in header_frequency))
    fit = config["fit"]
    assert isinstance(fit, Mapping)
    print(
        "Frequencies selected by task profile "
        f"({fit.get('frequency_source', 'fits_header')}) [GHz]:"
    )
    print("  " + ", ".join(f"{value:.9f}" for value in fit_frequency))

    try:
        rms = configured_rms(config)
    except MissingRMSConfiguration as error:
        print(f"Fit readiness: BLOCKED — {error}")
    else:
        print(
            "Fit readiness: ready "
            f"({len(rms)} positive per-band RMS values configured)"
        )


def _prepare_output(path: Path, *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists: {path}; pass --overwrite to replace it"
        )
    path.parent.mkdir(parents=True, exist_ok=True)


def prepare_cornish_maps(
    cornish_path: Path,
    reference_path: Path,
    output_directory: Path,
    *,
    original_beam_arcsec: float,
    target_beam_arcsec: float,
    overwrite: bool = False,
) -> tuple[Path, Path]:
    """Smooth CORNISH to a common beam and reproject it to the MeerKAT grid."""
    if target_beam_arcsec <= original_beam_arcsec:
        raise ValueError("Target beam must be larger than the original beam")

    native_output = output_directory / PREPARED_NATIVE
    reprojected_output = output_directory / PREPARED_REPROJECTED
    _prepare_output(native_output, overwrite=overwrite)
    _prepare_output(reprojected_output, overwrite=overwrite)

    with fits.open(cornish_path, memmap=True) as hdul:
        cornish_data = np.squeeze(hdul[0].data).astype(float)
        cornish_header = hdul[0].header.copy()
    if cornish_data.ndim != 2:
        raise ValueError(
            f"Expected a two-dimensional CORNISH map, got {cornish_data.shape}"
        )

    pixel_scale_x = abs(float(cornish_header["CDELT1"])) * 3600.0
    pixel_scale_y = abs(float(cornish_header["CDELT2"])) * 3600.0
    pixel_scale = 0.5 * (pixel_scale_x + pixel_scale_y)
    kernel_fwhm = np.sqrt(
        target_beam_arcsec**2 - original_beam_arcsec**2
    )
    kernel_sigma_pixels = (
        kernel_fwhm
        / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        / pixel_scale
    )
    kernel = Gaussian2DKernel(kernel_sigma_pixels)
    kernel.normalize()
    smoothed = convolve_fft(
        cornish_data,
        kernel,
        boundary="fill",
        fill_value=np.nan,
        nan_treatment="interpolate",
        preserve_nan=True,
        allow_huge=True,
    )
    smoothed *= (target_beam_arcsec / original_beam_arcsec) ** 2

    native_header = cornish_header.copy()
    if "BLANK" in native_header:
        del native_header["BLANK"]
    native_header["BMAJ"] = target_beam_arcsec / 3600.0
    native_header["BMIN"] = target_beam_arcsec / 3600.0
    native_header["BPA"] = 0.0
    native_header["BUNIT"] = "Jy/beam"
    native_header.add_history(
        "LANCET construct_freefree_model: smoothed and beam-rescaled"
    )
    fits.writeto(
        native_output,
        smoothed.astype(np.float32),
        native_header,
        overwrite=overwrite,
    )

    with fits.open(reference_path, memmap=True) as hdul:
        reference_shape = tuple(np.squeeze(hdul[0].data).shape)
        reference_header = hdul[0].header.copy()
    if len(reference_shape) != 2:
        raise ValueError(
            "Expected a two-dimensional MeerKAT reference image, "
            f"got {reference_shape}"
        )

    target_header = WCS(reference_header).celestial.to_header()
    target_header["NAXIS"] = 2
    target_header["NAXIS1"] = reference_shape[1]
    target_header["NAXIS2"] = reference_shape[0]
    target_header["BMAJ"] = target_beam_arcsec / 3600.0
    target_header["BMIN"] = target_beam_arcsec / 3600.0
    target_header["BPA"] = 0.0
    target_header["BUNIT"] = "Jy/beam"
    target_header.add_history(
        "LANCET construct_freefree_model: reprojected from smoothed CORNISH"
    )

    reprojected, footprint = reproject_interp(
        (smoothed, WCS(native_header).celestial),
        target_header,
        shape_out=reference_shape,
        order="bilinear",
    )
    reprojected[footprint <= 0] = np.nan
    fits.writeto(
        reprojected_output,
        reprojected.astype(np.float32),
        target_header,
        overwrite=overwrite,
    )
    return native_output, reprojected_output


def _celestial_header(header: fits.Header, shape: tuple[int, int]) -> fits.Header:
    output = WCS(header).celestial.to_header()
    output["NAXIS"] = 2
    output["NAXIS1"] = shape[1]
    output["NAXIS2"] = shape[0]
    for key in ("BMAJ", "BMIN", "BPA", "BUNIT"):
        if key in header:
            output[key] = header[key]
    return output


def save_fit_results(
    results: Mapping[str, np.ndarray],
    header_2d: fits.Header,
    frequency_ghz: np.ndarray,
    output_directory: Path,
    *,
    electron_temperature_k: float,
    overwrite: bool = False,
) -> list[Path]:
    """Save fitted maps plus explicit metadata for the irregular band axis."""
    map_units = {
        "EM_map": "pc cm-6",
        "EMerr_map": "pc cm-6",
        "ff230_map": "Jy/beam",
        "chi2_map": "",
        "rchi2_map": "",
        "nband_map": "",
        "fitmask_map": "",
    }
    outputs = [
        output_directory / f"{name}.fits"
        for name in map_units
    ]
    if "model_cube" in results:
        outputs.append(output_directory / "model_subbands_cube.fits")
    for output in outputs:
        _prepare_output(output, overwrite=overwrite)

    written: list[Path] = []
    for name, unit in map_units.items():
        output = output_directory / f"{name}.fits"
        header = header_2d.copy()
        header["TE"] = (
            float(electron_temperature_k),
            "Fixed electron temperature [K]",
        )
        if unit:
            header["BUNIT"] = unit
        elif "BUNIT" in header:
            del header["BUNIT"]
        fits.writeto(
            output,
            np.asarray(
                results[name],
                dtype=np.int16
                if name in {"nband_map", "fitmask_map"}
                else np.float32,
            ),
            header,
            overwrite=overwrite,
        )
        written.append(output)

    if "model_cube" in results:
        output = output_directory / "model_subbands_cube.fits"
        cube_header = header_2d.copy()
        cube_header["CTYPE3"] = "BAND"
        cube_header["CRPIX3"] = 1.0
        cube_header["CRVAL3"] = 1.0
        cube_header["CDELT3"] = 1.0
        cube_header["BUNIT"] = "Jy/beam"
        cube_header["TE"] = (
            float(electron_temperature_k),
            "Fixed electron temperature [K]",
        )
        names = [
            *(f"MeerKAT_{index:02d}" for index in range(1, 15)),
            "CORNISH_5.5GHz",
        ]
        bands = fits.BinTableHDU.from_columns(
            [
                fits.Column(
                    name="BAND",
                    format="I",
                    array=np.arange(1, 16, dtype=np.int16),
                ),
                fits.Column(
                    name="NAME",
                    format="20A",
                    array=np.asarray(names),
                ),
                fits.Column(
                    name="FREQUENCY_GHZ",
                    format="D",
                    unit="GHz",
                    array=np.asarray(frequency_ghz, dtype=float),
                ),
            ],
            name="BANDS",
        )
        fits.HDUList(
            [
                fits.PrimaryHDU(
                    data=np.asarray(
                        results["model_cube"],
                        dtype=np.float32,
                    ),
                    header=cube_header,
                ),
                bands,
            ]
        ).writeto(output, overwrite=overwrite)
        written.append(output)
    return written


def run_fit(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> None:
    """Run all configured temperature fits."""
    fit = config["fit"]
    preparation = config["preparation"]
    if not isinstance(fit, Mapping) or not isinstance(preparation, Mapping):
        raise TypeError("Fit and preparation settings must be TOML tables")
    rms = configured_rms(config)

    prepared_cornish = paths["output_root"] / "prepared" / PREPARED_REPROJECTED
    if not prepared_cornish.is_file():
        raise FileNotFoundError(
            f"Prepared CORNISH map not found: {prepared_cornish}; "
            "run --stage prepare first"
        )

    with fits.open(paths["meerkat_cube"], memmap=True) as hdul:
        cube = np.asarray(hdul[0].data)
        cube_header = hdul[0].header.copy()
    with fits.open(prepared_cornish, memmap=True) as hdul:
        cornish = np.asarray(hdul[0].data)

    mfs_map, alpha_map, sed_cube = prepare_sed_cube(cube, cornish)
    frequency = configured_frequencies(config, cube_header)

    fit_mask = (
        np.isfinite(mfs_map)
        & np.isfinite(alpha_map)
        & (mfs_map > float(fit["mfs_threshold_jy_per_beam"]))
        & (alpha_map >= float(fit["alpha_min"]))
        & (alpha_map <= float(fit["alpha_max"]))
    )
    header_2d = _celestial_header(cube_header, mfs_map.shape)
    common_beam = float(preparation["common_beam_arcsec"])
    omega = beam_solid_angle(common_beam, common_beam)

    for temperature in fit["electron_temperatures_k"]:
        temperature = float(temperature)
        print(f"Fitting Te={temperature:g} K")
        results = fit_freefree_map(
            sed_cube,
            frequency,
            rms,
            beam_solid_angle_sr=omega,
            electron_temperature_k=temperature,
            fit_mask=fit_mask,
            minimum_bands=int(fit["minimum_bands"]),
            snr_threshold=float(fit["snr_threshold"]),
            save_model_cube=True,
            progress=lambda done, total: print(
                f"  processed {done}/{total} candidate pixels"
            ),
        )
        output_directory = (
            paths["output_root"] / f"Te_{temperature:g}K"
        )
        written = save_fit_results(
            results,
            header_2d,
            frequency,
            output_directory,
            electron_temperature_k=temperature,
            overwrite=overwrite,
        )
        print(f"  wrote {len(written)} FITS products to {output_directory}")


def postprocess_freefree(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> tuple[Path, Path]:
    """Smooth the reference free-free map and reproject it to the dust grid."""
    preparation = config["preparation"]
    postprocess = config["postprocess"]
    if not isinstance(preparation, Mapping) or not isinstance(
        postprocess, Mapping
    ):
        raise TypeError(
            "Preparation and postprocess settings must be TOML tables"
        )

    temperature = float(postprocess["reference_temperature_k"])
    temperature_directory = (
        paths["output_root"] / f"Te_{temperature:g}K"
    )
    input_map = temperature_directory / "ff230_map.fits"
    if not input_map.is_file():
        raise FileNotFoundError(
            f"Fitted free-free map not found: {input_map}; "
            "run --stage fit first"
        )

    smoothed_output = (
        temperature_directory / "ff230_map_smoothed_18arcsec.fits"
    )
    reprojected_output = (
        temperature_directory
        / "ff230_map_smoothed_18arcsec_regrid2dust.fits"
    )
    _prepare_output(smoothed_output, overwrite=overwrite)
    _prepare_output(reprojected_output, overwrite=overwrite)

    data, header = fits.getdata(input_map, header=True)
    data = np.squeeze(np.asarray(data))
    if data.ndim != 2:
        raise ValueError(f"Expected a two-dimensional map, got {data.shape}")
    pixel_scale = 0.5 * (
        abs(float(header["CDELT1"])) + abs(float(header["CDELT2"]))
    ) * 3600.0
    original_beam = float(preparation["common_beam_arcsec"])
    target_beam = float(postprocess["smoothed_beam_arcsec"])
    smoothed = smooth_circular_jybeam(
        data,
        pixel_scale_arcsec=pixel_scale,
        original_beam_arcsec=original_beam,
        target_beam_arcsec=target_beam,
        replace_nan_with_zero=True,
    )

    smoothed_header = header.copy()
    smoothed_header["BMAJ"] = target_beam / 3600.0
    smoothed_header["BMIN"] = target_beam / 3600.0
    smoothed_header["BPA"] = 0.0
    smoothed_header["BUNIT"] = "Jy/beam"
    smoothed_header.add_history(
        "LANCET construct_freefree_model: smoothed to target beam"
    )
    fits.writeto(
        smoothed_output,
        smoothed.astype(np.float32),
        smoothed_header,
        overwrite=overwrite,
    )

    dust_data, dust_header = fits.getdata(
        paths["dust_reference"],
        header=True,
    )
    dust_shape = np.squeeze(dust_data).shape
    if len(dust_shape) != 2:
        raise ValueError(
            f"Expected a two-dimensional dust reference, got {dust_shape}"
        )
    target_header = WCS(dust_header).celestial.to_header()
    target_header["NAXIS"] = 2
    target_header["NAXIS1"] = dust_shape[1]
    target_header["NAXIS2"] = dust_shape[0]
    for key in ("BMAJ", "BMIN", "BPA"):
        if key in dust_header:
            target_header[key] = dust_header[key]
    target_header["BUNIT"] = "Jy/beam"
    target_header.add_history(
        "LANCET construct_freefree_model: reprojected to dust grid"
    )
    reprojected, footprint = reproject_interp(
        (smoothed, WCS(smoothed_header).celestial),
        target_header,
        shape_out=dust_shape,
        order="bilinear",
    )
    reprojected[footprint <= 0] = np.nan
    fits.writeto(
        reprojected_output,
        reprojected.astype(np.float32),
        target_header,
        overwrite=overwrite,
    )
    return smoothed_output, reprojected_output


def plot_temperature_flux_comparison(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> list[Path]:
    """Plot 230 GHz flux changes caused by the adopted temperatures."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib as mpl
    from matplotlib import pyplot as plt

    settings = config["temperature_comparison_plot"]
    if not isinstance(settings, Mapping):
        raise TypeError(
            "Temperature-comparison plot settings must be a TOML table"
        )
    reference_temperature = float(settings["reference_temperature_k"])
    comparison_temperatures = np.asarray(
        settings["comparison_temperatures_k"],
        dtype=float,
    )
    colors = [str(value) for value in settings["colors"]]
    if comparison_temperatures.ndim != 1 or len(comparison_temperatures) == 0:
        raise ValueError(
            "At least one comparison temperature must be configured"
        )
    if len(colors) != len(comparison_temperatures):
        raise ValueError(
            "A plotting color is required for each comparison temperature"
        )

    reference_path = (
        paths["output_root"]
        / f"Te_{reference_temperature:g}K"
        / "ff230_map.fits"
    )
    if not reference_path.is_file():
        raise FileNotFoundError(
            f"Reference free-free map not found: {reference_path}; "
            "run --stage fit first"
        )
    reference = np.squeeze(np.asarray(fits.getdata(reference_path)))
    if reference.ndim != 2:
        raise ValueError(
            f"Expected a two-dimensional map, got {reference.shape}"
        )

    distributions: list[tuple[float, np.ndarray]] = []
    for temperature in comparison_temperatures:
        comparison_path = (
            paths["output_root"]
            / f"Te_{temperature:g}K"
            / "ff230_map.fits"
        )
        if not comparison_path.is_file():
            raise FileNotFoundError(
                f"Comparison free-free map not found: {comparison_path}; "
                "run --stage fit first"
            )
        comparison = np.squeeze(np.asarray(fits.getdata(comparison_path)))
        difference = fractional_flux_difference(
            comparison,
            reference,
            absolute=True,
        )
        values = difference[np.isfinite(difference) & (difference > 0)]
        if values.size == 0:
            raise ValueError(
                f"No positive finite differences for Te={temperature:g} K"
            )
        distributions.append((float(temperature), values))
        print(
            f"Te={temperature:g} K: {values.size} pixels, "
            f"median difference={np.median(values) * 100:.4g}%"
        )

    all_values = np.concatenate([values for _, values in distributions])
    lower = float(np.min(all_values))
    upper = float(np.max(all_values))
    if not np.isfinite(lower) or not np.isfinite(upper) or lower >= upper:
        raise ValueError("Temperature differences need a finite positive range")
    bins = np.logspace(
        np.log10(lower),
        np.log10(upper),
        int(settings["histogram_bins"]),
    )

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
            "font.size": 11,
            "axes.labelsize": 11,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "axes.linewidth": 0.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    figure, axis = plt.subplots(figsize=(5.0, 4.0))
    for (temperature, values), color in zip(
        distributions,
        colors,
        strict=True,
    ):
        axis.hist(
            values,
            bins=bins,
            density=False,
            histtype="step",
            linewidth=2.0,
            label=rf"$T_e = {temperature:g}\,\mathrm{{K}}$",
            color=color,
        )
    axis.set_xscale("log")
    axis.set_xticks([0.001, 0.01, 0.1])
    axis.set_xticklabels(["0.1", "1", "10"])
    axis.set_xlabel(
        r"$|\Delta S_{230\,\mathrm{GHz}}^{\mathrm{ff}}(T_{\rm e})/"
        rf"S_{{230\,\mathrm{{GHz}}}}^{{\mathrm{{ff}}}}"
        rf"({reference_temperature:g}\,\mathrm{{K}})|$ (%)"
    )
    axis.set_ylabel("Probability density")
    axis.grid(False)
    axis.legend(frameon=False, loc="upper right")
    figure.tight_layout()

    output_formats = [
        str(value).lower() for value in settings["output_formats"]
    ]
    unsupported = set(output_formats) - {"pdf", "png"}
    if unsupported:
        raise ValueError(
            "Unsupported figure formats: " + ", ".join(sorted(unsupported))
        )
    output_stem = (
        paths["output_root"]
        / "figures"
        / "FreeFreeTemplate_Te_varying_flux_difference"
    )
    outputs = [
        output_stem.with_suffix(f".{suffix}") for suffix in output_formats
    ]
    for output in outputs:
        _prepare_output(output, overwrite=overwrite)
    for output in outputs:
        figure.savefig(
            output,
            bbox_inches="tight",
            dpi=int(settings["dpi"]),
        )
    plt.close(figure)
    print(f"Wrote {len(outputs)} temperature-comparison figures")
    return outputs


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="Task TOML file, relative to the project root by default",
    )
    parser.add_argument(
        "--stage",
        choices=("inspect", "prepare", "fit", "postprocess", "plot", "all"),
        default="inspect",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    project = ProjectPaths.discover()
    config_path = args.config
    if not config_path.is_absolute():
        config_path = project.root / config_path
    config = load_task_config(config_path.resolve())
    paths = resolve_task_paths(project, config)
    preparation = config["preparation"]
    if not isinstance(preparation, Mapping):
        raise TypeError("Preparation settings must be a TOML table")

    if args.stage in {"inspect", "all"}:
        inspect_inputs(paths, config)
    if args.stage in {"prepare", "all"}:
        native, reprojected = prepare_cornish_maps(
            paths["cornish_map"],
            paths["meerkat_reference"],
            paths["output_root"] / "prepared",
            original_beam_arcsec=float(
                preparation["cornish_beam_arcsec"]
            ),
            target_beam_arcsec=float(preparation["common_beam_arcsec"]),
            overwrite=args.overwrite,
        )
        print(f"Wrote prepared CORNISH maps:\n  {native}\n  {reprojected}")
    if args.stage in {"fit", "all"}:
        run_fit(paths, config, overwrite=args.overwrite)
    if args.stage in {"postprocess", "all"}:
        smoothed, reprojected = postprocess_freefree(
            paths,
            config,
            overwrite=args.overwrite,
        )
        print(
            "Wrote post-processed free-free maps:"
            f"\n  {smoothed}\n  {reprojected}"
        )
    if args.stage in {"plot", "all"}:
        plot_temperature_flux_comparison(
            paths,
            config,
            overwrite=args.overwrite,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
