#!/usr/bin/env python
"""Construct dust plus free-free continuum models and SIMBA validation maps."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from pathlib import Path
import tomllib

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

from lancet import ProjectPaths
from lancet.continuum import (
    calibrate_simba,
    smooth_circular_jybeam,
    synthesize_continuum,
)


DEFAULT_CONFIG = Path("tasks/synthesize_continuum/task.toml")


def load_task_config(path: Path) -> dict[str, object]:
    """Load task settings from TOML."""
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
    *,
    freefree_override: Path | None = None,
    output_override: Path | None = None,
) -> dict[str, Path]:
    """Resolve the configured inputs and policy-safe output directory."""
    inputs = config["inputs"]
    outputs = config["outputs"]
    if not isinstance(inputs, Mapping) or not isinstance(outputs, Mapping):
        raise TypeError("Task inputs and outputs must be TOML tables")

    paths = {
        name: _project_path(project, value)
        for name, value in inputs.items()
    }
    if freefree_override is not None:
        paths["freefree_230ghz"] = _project_path(
            project,
            freefree_override,
        )
    output_value: object = (
        output_override if output_override is not None else outputs["root"]
    )
    paths["output_root"] = project.validate_output(str(output_value))

    for name in ("dust_230ghz", "simba_250ghz"):
        path = paths[name]
        if not path.is_relative_to(project.delivered.resolve()):
            raise ValueError(f"{name} must be a delivered input: {path}")
        if not path.is_file():
            raise FileNotFoundError(f"{name} not found: {path}")
    return paths


def _describe(path: Path) -> str:
    if not path.is_file():
        return "MISSING"
    with fits.open(path, memmap=True) as hdul:
        header = hdul[0].header
        shape = tuple(hdul[0].shape)
        beam = (
            f"{float(header['BMAJ']) * 3600:g} x "
            f"{float(header['BMIN']) * 3600:g} arcsec"
            if "BMAJ" in header and "BMIN" in header
            else "not recorded"
        )
        unit = header.get("BUNIT", "not recorded")
    return f"shape={shape}, beam={beam}, unit={unit}"


def inspect_inputs(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
) -> None:
    """Print input readiness and configured model assumptions."""
    print("Continuum synthesis inputs:")
    for name in ("dust_230ghz", "freefree_230ghz", "simba_250ghz"):
        print(f"  {name}: {paths[name]}")
        print(f"    {_describe(paths[name])}")

    model = config["model"]
    if not isinstance(model, Mapping):
        raise TypeError("Model settings must be a TOML table")
    print(
        "Configured frequencies [GHz]: "
        + ", ".join(str(value) for value in model["frequencies_ghz"])
    )
    print(
        "Spectral indices: "
        f"dust={model['dust_spectral_index']}, "
        f"free-free={model['freefree_spectral_index']}"
    )
    status = (
        "ready"
        if paths["freefree_230ghz"].is_file()
        else "blocked until the free-free postprocess stage is complete"
    )
    print(f"Model readiness: {status}")


def _prepare_output(path: Path, *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists: {path}; pass --overwrite to replace it"
        )
    path.parent.mkdir(parents=True, exist_ok=True)


def _same_grid(
    first_shape: tuple[int, ...],
    first_header: fits.Header,
    second_shape: tuple[int, ...],
    second_header: fits.Header,
) -> bool:
    if first_shape != second_shape:
        return False
    first = WCS(first_header).celestial
    second = WCS(second_header).celestial
    return bool(
        np.allclose(first.wcs.crpix, second.wcs.crpix, atol=1e-8)
        and np.allclose(first.wcs.crval, second.wcs.crval, atol=1e-10)
        and np.allclose(first.wcs.cdelt, second.wcs.cdelt, atol=1e-12)
        and np.allclose(first.wcs.pc, second.wcs.pc, atol=1e-12)
    )


def _component_header(
    reference: fits.Header,
    *,
    frequency_ghz: float,
    component: str,
) -> fits.Header:
    header = reference.copy()
    if "BLANK" in header and int(header["BITPIX"]) < 0:
        del header["BLANK"]
    header["BUNIT"] = "Jy/beam"
    header["FREQGHZ"] = (float(frequency_ghz), "Model frequency [GHz]")
    header["BTYPE"] = component
    header.add_history(
        "LANCET synthesize_continuum: power-law component synthesis"
    )
    return header


def build_models(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> dict[float, dict[str, Path]]:
    """Build component and combined continuum FITS maps."""
    if not paths["freefree_230ghz"].is_file():
        raise FileNotFoundError(
            f"Free-free input not found: {paths['freefree_230ghz']}; "
            "run the free-free fit and postprocess stages first"
        )

    dust, dust_header = fits.getdata(
        paths["dust_230ghz"],
        header=True,
    )
    freefree, freefree_header = fits.getdata(
        paths["freefree_230ghz"],
        header=True,
    )
    dust = np.squeeze(np.asarray(dust))
    freefree = np.squeeze(np.asarray(freefree))
    if not _same_grid(
        dust.shape,
        dust_header,
        freefree.shape,
        freefree_header,
    ):
        raise ValueError(
            "Dust and free-free inputs must have identical shapes and WCS grids"
        )

    model = config["model"]
    if not isinstance(model, Mapping):
        raise TypeError("Model settings must be a TOML table")
    reference_frequency = float(model["reference_frequency_ghz"])
    model_directory = paths["output_root"] / "model"
    written: dict[float, dict[str, Path]] = {}

    for target_value in model["frequencies_ghz"]:
        target = float(target_value)
        dust_scaled, freefree_scaled, total = synthesize_continuum(
            dust,
            freefree,
            reference_frequency_ghz=reference_frequency,
            target_frequency_ghz=target,
            dust_spectral_index=float(model["dust_spectral_index"]),
            freefree_spectral_index=float(
                model["freefree_spectral_index"]
            ),
            missing_freefree_value=float(
                model["missing_freefree_value_jy_per_beam"]
            ),
        )
        frequency_label = f"{target:g}GHz"
        products = {
            "dust": model_directory / f"dust_{frequency_label}.fits",
            "freefree": (
                model_directory / f"freefree_{frequency_label}.fits"
            ),
            "continuum": (
                model_directory / f"continuum_{frequency_label}.fits"
            ),
        }
        for output in products.values():
            _prepare_output(output, overwrite=overwrite)

        for component, data in (
            ("dust", dust_scaled),
            ("freefree", freefree_scaled),
            ("continuum", total),
        ):
            fits.writeto(
                products[component],
                np.asarray(data, dtype=np.float32),
                _component_header(
                    dust_header,
                    frequency_ghz=target,
                    component=component,
                ),
                overwrite=overwrite,
            )
        written[target] = products
        print(f"Wrote {target:g} GHz model products to {model_directory}")
    return written


def build_simba_validation(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> dict[str, Path]:
    """Calibrate SIMBA, reproject it, and smooth the 250 GHz model."""
    beam = config["beam"]
    simba_settings = config["simba"]
    if not isinstance(beam, Mapping) or not isinstance(
        simba_settings, Mapping
    ):
        raise TypeError("Beam and SIMBA settings must be TOML tables")

    continuum_250 = (
        paths["output_root"] / "model" / "continuum_250GHz.fits"
    )
    if not continuum_250.is_file():
        raise FileNotFoundError(
            f"250 GHz model not found: {continuum_250}; "
            "run --stage model first"
        )

    validation_directory = paths["output_root"] / "validation"
    outputs = {
        "model_smoothed": (
            validation_directory / "continuum_250GHz_smoothed_24arcsec.fits"
        ),
        "simba_calibrated": (
            validation_directory / "simba_250GHz_calibrated.fits"
        ),
        "simba_reprojected": (
            validation_directory / "simba_250GHz_on_model_grid.fits"
        ),
    }
    for output in outputs.values():
        _prepare_output(output, overwrite=overwrite)

    model_data, model_header = fits.getdata(
        continuum_250,
        header=True,
    )
    pixel_scale = 0.5 * (
        abs(float(model_header["CDELT1"]))
        + abs(float(model_header["CDELT2"]))
    ) * 3600.0
    smoothed_model = smooth_circular_jybeam(
        model_data,
        pixel_scale_arcsec=pixel_scale,
        original_beam_arcsec=float(beam["model_fwhm_arcsec"]),
        target_beam_arcsec=float(beam["simba_fwhm_arcsec"]),
    )
    smoothed_header = model_header.copy()
    smoothed_header["BMAJ"] = float(beam["simba_fwhm_arcsec"]) / 3600.0
    smoothed_header["BMIN"] = float(beam["simba_fwhm_arcsec"]) / 3600.0
    smoothed_header["BPA"] = 0.0
    smoothed_header.add_history(
        "LANCET synthesize_continuum: smoothed for SIMBA comparison"
    )
    fits.writeto(
        outputs["model_smoothed"],
        smoothed_model.astype(np.float32),
        smoothed_header,
        overwrite=overwrite,
    )

    simba_raw, simba_header = fits.getdata(
        paths["simba_250ghz"],
        header=True,
    )
    simba_calibrated = calibrate_simba(
        simba_raw,
        invalid_floor=float(simba_settings["invalid_floor"]),
        flux_scale=float(simba_settings["flux_scale"]),
    )
    if "BLANK" in simba_header:
        del simba_header["BLANK"]
    simba_header["BMAJ"] = float(beam["simba_fwhm_arcsec"]) / 3600.0
    simba_header["BMIN"] = float(beam["simba_fwhm_arcsec"]) / 3600.0
    simba_header["BPA"] = 0.0
    simba_header["BUNIT"] = "Jy/beam"
    simba_header["FREQGHZ"] = (250.0, "Representative frequency [GHz]")
    simba_header.add_history(
        "LANCET synthesize_continuum: applied SIMBA calibration"
    )
    fits.writeto(
        outputs["simba_calibrated"],
        simba_calibrated.astype(np.float32),
        simba_header,
        overwrite=overwrite,
    )

    method = str(simba_settings["reprojection_method"])
    if method != "hcongrid":
        raise ValueError("Only the Paper II hcongrid profile is supported")
    from FITS_tools.hcongrid import hcongrid

    simba_reprojected = hcongrid(
        simba_calibrated,
        simba_header,
        model_header,
        preserve_bad_pixels=False,
    )
    reprojected_header = model_header.copy()
    reprojected_header["BMAJ"] = (
        float(beam["simba_fwhm_arcsec"]) / 3600.0
    )
    reprojected_header["BMIN"] = (
        float(beam["simba_fwhm_arcsec"]) / 3600.0
    )
    reprojected_header["BPA"] = 0.0
    reprojected_header["BUNIT"] = "Jy/beam"
    reprojected_header.add_history(
        "LANCET synthesize_continuum: reprojected SIMBA to model grid"
    )
    fits.writeto(
        outputs["simba_reprojected"],
        np.asarray(simba_reprojected, dtype=np.float32),
        reprojected_header,
        overwrite=overwrite,
    )
    print(f"Wrote SIMBA validation products to {validation_directory}")
    return outputs


def _plot_settings(
    config: Mapping[str, object],
) -> tuple[Mapping[str, object], Mapping[str, object], Mapping[str, object]]:
    plot = config["plot"]
    if not isinstance(plot, Mapping):
        raise TypeError("Plot settings must be a TOML table")
    component = plot["component_230ghz"]
    simba = plot["simba_250ghz"]
    if not isinstance(component, Mapping) or not isinstance(simba, Mapping):
        raise TypeError("Plot panel settings must be TOML tables")
    return plot, component, simba


def _configure_plot_style() -> None:
    import matplotlib as mpl

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
            "font.size": 10,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 10,
            "axes.linewidth": 0.8,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.top": True,
            "ytick.right": True,
        }
    )


def _plot_colormap(name: str):
    if name == "mako":
        import seaborn as sns

        return sns.color_palette("mako", n_colors=256, as_cmap=True)

    from matplotlib import colormaps

    try:
        return colormaps[name]
    except KeyError as error:
        raise ValueError(f"Unknown plotting colormap: {name}") from error


def _read_plot_image(path: Path) -> tuple[np.ndarray, fits.Header]:
    data, header = fits.getdata(path, header=True)
    image = np.squeeze(np.asarray(data, dtype=float))
    if image.ndim != 2:
        raise ValueError(f"Expected a two-dimensional image: {path}")
    return np.where(np.isfinite(image), image, 0.0), header


def _draw_map_panel(
    figure,
    subplot_spec,
    path: Path,
    *,
    title: str,
    colormap,
    normalization,
    contour_levels: Sequence[float],
    contour_colormap: str | None,
    contour_smoothing_sigma_pixels: float,
    crop_pixels: Sequence[float],
    ra_tick_spacing_arcmin: float,
    scalebar_length_deg: float,
    scalebar_label: str,
    show_declination_labels: bool,
):
    from astropy import units as u
    from astropy.convolution import Gaussian2DKernel, convolve_fft
    from astropy.visualization.wcsaxes import add_beam, add_scalebar

    data, header = _read_plot_image(path)
    axis = figure.add_subplot(
        subplot_spec,
        projection=WCS(header).celestial,
    )
    image = axis.imshow(
        data,
        origin="lower",
        cmap=colormap,
        norm=normalization,
    )

    kernel = Gaussian2DKernel(contour_smoothing_sigma_pixels)
    kernel.normalize()
    contour_data = convolve_fft(data, kernel)
    contour_options: dict[str, object] = {
        "levels": list(contour_levels),
        "linewidths": 0.4,
    }
    if contour_colormap is None:
        contour_options["colors"] = "black"
    else:
        contour_options["cmap"] = contour_colormap
    axis.contour(contour_data, **contour_options)

    axis.text(
        0.02,
        0.98,
        title,
        transform=axis.transAxes,
        ha="left",
        va="top",
        color="white",
        size=10,
        weight="bold",
    )
    lower, upper = (float(value) for value in crop_pixels)
    axis.set_xlim(lower, upper)
    axis.set_ylim(lower, upper)
    axis.tick_params(axis="both", length=3)
    axis.coords[0].set_ticks(
        spacing=float(ra_tick_spacing_arcmin) * u.arcmin
    )
    axis.set_xlabel("")
    if show_declination_labels:
        axis.set_ylabel("Dec (ICRS)", labelpad=0)
    else:
        axis.set_ylabel("")
        axis.coords[1].set_ticklabel_visible(False)

    add_beam(
        axis,
        header,
        edgecolor="white",
        facecolor="none",
    )
    add_scalebar(
        axis,
        length=float(scalebar_length_deg) * u.deg,
        label=scalebar_label,
        color="white",
    )
    return axis, image


def _save_figure(
    figure,
    output_stem: Path,
    *,
    output_formats: Sequence[str],
    dpi: int,
    overwrite: bool,
) -> list[Path]:
    formats = [str(value).lower() for value in output_formats]
    unsupported = set(formats) - {"pdf", "png"}
    if unsupported:
        raise ValueError(
            "Unsupported figure formats: " + ", ".join(sorted(unsupported))
        )
    outputs = [output_stem.with_suffix(f".{suffix}") for suffix in formats]
    for output in outputs:
        _prepare_output(output, overwrite=overwrite)
    for output in outputs:
        figure.savefig(
            output,
            bbox_inches="tight",
            dpi=int(dpi),
        )
    return outputs


def plot_continuum_components(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> list[Path]:
    """Plot the 230 GHz dust, free-free, and combined continuum maps."""
    import matplotlib

    matplotlib.use("Agg")
    from astropy.visualization import ImageNormalize, SqrtStretch
    from matplotlib import pyplot as plt
    import matplotlib.gridspec as gridspec

    plot, component, _ = _plot_settings(config)
    _configure_plot_style()
    colormap = _plot_colormap(str(plot["colormap"]))
    normalization = ImageNormalize(
        stretch=SqrtStretch(),
        vmin=float(component["vmin_jy_per_beam"]),
        vmax=float(component["vmax_jy_per_beam"]),
    )
    model_directory = paths["output_root"] / "model"
    panels = (
        (
            model_directory / "dust_230GHz.fits",
            "Dust thermal emission",
        ),
        (
            model_directory / "freefree_230GHz.fits",
            "Free-free emission",
        ),
        (
            model_directory / "continuum_230GHz.fits",
            "1.3 mm emission (Dust + Free-free)",
        ),
    )
    for path, _ in panels:
        if not path.is_file():
            raise FileNotFoundError(
                f"Plot input not found: {path}; run --stage model first"
            )

    figure = plt.figure(figsize=(9, 3), dpi=int(plot["dpi"]))
    grid = gridspec.GridSpec(
        nrows=1,
        ncols=4,
        figure=figure,
        left=0.04,
        right=0.95,
        bottom=0.15,
        top=0.95,
        wspace=0.05,
        width_ratios=[1, 1, 1, 0.05],
    )
    last_image = None
    for index, (path, title) in enumerate(panels):
        axis, last_image = _draw_map_panel(
            figure,
            grid[index],
            path,
            title=title,
            colormap=colormap,
            normalization=normalization,
            contour_levels=component[
                "contour_levels_jy_per_beam"
            ],
            contour_colormap=None,
            contour_smoothing_sigma_pixels=float(
                plot["contour_smoothing_sigma_pixels"]
            ),
            crop_pixels=plot["crop_pixels"],
            ra_tick_spacing_arcmin=float(
                plot["ra_tick_spacing_arcmin"]
            ),
            scalebar_length_deg=float(plot["scalebar_length_deg"]),
            scalebar_label=str(plot["scalebar_label"]),
            show_declination_labels=index == 0,
        )
        if index == 1:
            axis.set_xlabel("RA (ICRS)")

    assert last_image is not None
    colorbar_axis = figure.add_subplot(grid[3])
    colorbar = figure.colorbar(last_image, cax=colorbar_axis)
    colorbar.set_label(
        r"Intensity (Jy beam$^{-1}$)",
        fontsize=11,
        rotation=270,
        labelpad=15,
    )
    colorbar.ax.tick_params(
        axis="y",
        labelright=True,
        labelleft=False,
        labelsize=10,
        length=2,
    )
    colorbar.ax.yaxis.set_label_position("right")
    colorbar.ax.yaxis.set_ticks_position("right")

    outputs = _save_figure(
        figure,
        paths["output_root"] / "figures" / "continuum_components_230GHz",
        output_formats=plot["output_formats"],
        dpi=int(plot["dpi"]),
        overwrite=overwrite,
    )
    plt.close(figure)
    return outputs


def plot_simba_comparison(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> list[Path]:
    """Plot the observed and synthesized 250 GHz continuum maps."""
    import matplotlib

    matplotlib.use("Agg")
    from astropy.visualization import ImageNormalize, SqrtStretch
    from matplotlib import pyplot as plt
    import matplotlib.gridspec as gridspec

    plot, _, simba = _plot_settings(config)
    _configure_plot_style()
    colormap = _plot_colormap(str(plot["colormap"]))
    normalization = ImageNormalize(
        stretch=SqrtStretch(),
        vmin=float(simba["vmin_jy_per_beam"]),
        vmax=float(simba["vmax_jy_per_beam"]),
    )
    validation_directory = paths["output_root"] / "validation"
    panels = (
        (
            validation_directory / "simba_250GHz_on_model_grid.fits",
            "1.2 mm emission (SIMBA)",
        ),
        (
            validation_directory
            / "continuum_250GHz_smoothed_24arcsec.fits",
            "1.2 mm emission (Dust + Free-free)",
        ),
    )
    for path, _ in panels:
        if not path.is_file():
            raise FileNotFoundError(
                f"Plot input not found: {path}; run --stage validate first"
            )

    figure = plt.figure(figsize=(6, 3), dpi=int(plot["dpi"]))
    grid = gridspec.GridSpec(
        nrows=1,
        ncols=3,
        figure=figure,
        left=0.04,
        right=0.95,
        bottom=0.15,
        top=0.95,
        wspace=0.05,
        width_ratios=[1, 1, 0.05],
    )
    last_image = None
    for index, (path, title) in enumerate(panels):
        _, last_image = _draw_map_panel(
            figure,
            grid[index],
            path,
            title=title,
            colormap=colormap,
            normalization=normalization,
            contour_levels=simba["contour_levels_jy_per_beam"],
            contour_colormap="viridis_r",
            contour_smoothing_sigma_pixels=float(
                plot["contour_smoothing_sigma_pixels"]
            ),
            crop_pixels=plot["crop_pixels"],
            ra_tick_spacing_arcmin=float(
                plot["ra_tick_spacing_arcmin"]
            ),
            scalebar_length_deg=float(plot["scalebar_length_deg"]),
            scalebar_label=str(plot["scalebar_label"]),
            show_declination_labels=index == 0,
        )

    assert last_image is not None
    colorbar_axis = figure.add_subplot(grid[2])
    colorbar = figure.colorbar(last_image, cax=colorbar_axis)
    colorbar.set_label(
        r"Intensity (Jy beam$^{-1}$)",
        fontsize=11,
        rotation=270,
        labelpad=15,
    )
    colorbar.ax.tick_params(
        axis="y",
        labelright=True,
        labelleft=False,
        labelsize=10,
        length=2,
    )
    colorbar.ax.yaxis.set_label_position("right")
    colorbar.ax.yaxis.set_ticks_position("right")

    outputs = _save_figure(
        figure,
        paths["output_root"] / "figures" / "simba_comparison_250GHz",
        output_formats=plot["output_formats"],
        dpi=int(plot["dpi"]),
        overwrite=overwrite,
    )
    plt.close(figure)
    return outputs


def plot_sed_demonstration(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> list[Path]:
    """Plot an illustrative dust plus free-free spectral energy distribution."""
    import matplotlib

    matplotlib.use("Agg")
    from astropy import constants as const
    from matplotlib import pyplot as plt

    from lancet.freefree import beam_solid_angle, freefree_flux_density

    plot, _, _ = _plot_settings(config)
    sed = plot["sed_demonstration"]
    if not isinstance(sed, Mapping):
        raise TypeError("SED demonstration settings must be a TOML table")
    _configure_plot_style()

    frequency_limits = np.asarray(sed["frequency_range_ghz"], dtype=float)
    if (
        frequency_limits.shape != (2,)
        or np.any(~np.isfinite(frequency_limits))
        or np.any(frequency_limits <= 0)
        or frequency_limits[1] <= frequency_limits[0]
    ):
        raise ValueError(
            "plot.sed_demonstration.frequency_range_ghz must contain "
            "two increasing positive values"
        )
    frequency = np.logspace(
        np.log10(frequency_limits[0]),
        np.log10(frequency_limits[1]),
        int(sed["frequency_samples"]),
    )
    source_fwhm = np.asarray(sed["source_fwhm_arcsec"], dtype=float)
    if source_fwhm.shape != (2,):
        raise ValueError(
            "plot.sed_demonstration.source_fwhm_arcsec must contain "
            "major and minor FWHM values"
        )
    electron_temperature = float(sed["electron_temperature_k"])
    emission_measure = float(sed["emission_measure_pc_cm6"])
    freefree = freefree_flux_density(
        frequency,
        emission_measure,
        electron_temperature_k=electron_temperature,
        beam_solid_angle_sr=beam_solid_angle(*source_fwhm),
    )

    dust_temperature = float(sed["dust_temperature_k"])
    dust_index = float(sed["dust_emissivity_index"])
    dust_reference_frequency = float(
        sed["dust_reference_frequency_ghz"]
    )
    dust_reference_flux = float(sed["dust_reference_flux_jy"])
    frequency_hz = frequency * 1e9
    planck = (
        2.0
        * const.h.value
        * frequency_hz**3
        / const.c.value**2
        / np.expm1(
            const.h.value
            * frequency_hz
            / (const.k_B.value * dust_temperature)
        )
    )
    dust_shape = (
        (frequency / dust_reference_frequency) ** dust_index * planck
    )
    dust = (
        dust_reference_flux
        * dust_shape
        / np.interp(
            dust_reference_frequency,
            frequency,
            dust_shape,
        )
    )
    total = dust + freefree

    radio_frequency = np.asarray(
        sed["radio_frequencies_ghz"],
        dtype=float,
    )
    radio_flux = np.interp(radio_frequency, frequency, freefree)
    prediction_frequency = float(sed["prediction_frequency_ghz"])
    prediction_flux = float(
        np.interp(prediction_frequency, frequency, total)
    )
    wavelength_um = np.asarray(
        sed["far_infrared_wavelengths_um"],
        dtype=float,
    )
    far_infrared_frequency = (
        const.c.value / (wavelength_um * 1e-6) / 1e9
    )
    inside = (
        (far_infrared_frequency >= frequency_limits[0])
        & (far_infrared_frequency <= frequency_limits[1])
    )
    wavelength_um = wavelength_um[inside]
    far_infrared_frequency = far_infrared_frequency[inside]
    far_infrared_flux = np.interp(
        far_infrared_frequency,
        frequency,
        total,
    )

    figure, axis = plt.subplots(figsize=(8.5, 5.4), dpi=int(plot["dpi"]))
    figure.patch.set_facecolor("black")
    axis.set_facecolor("black")
    axis.loglog(
        frequency,
        total,
        color="white",
        linewidth=2.4,
        label="total SED = dust + free-free",
    )
    axis.loglog(
        frequency,
        freefree,
        color="white",
        linewidth=1.9,
        linestyle="--",
        label="free-free physical model",
    )
    axis.loglog(
        frequency,
        dust,
        color="white",
        linewidth=1.9,
        linestyle=":",
        label="dust modified blackbody",
    )
    axis.scatter(
        radio_frequency,
        radio_flux,
        s=65,
        facecolors="black",
        edgecolors="white",
        linewidth=1.5,
        zorder=5,
    )
    axis.scatter(
        far_infrared_frequency,
        far_infrared_flux,
        s=65,
        facecolors="black",
        edgecolors="white",
        linewidth=1.5,
        zorder=5,
    )
    axis.scatter(
        [prediction_frequency],
        [prediction_flux],
        marker="*",
        s=260,
        facecolors="black",
        edgecolors="white",
        linewidth=1.6,
        zorder=6,
    )
    axis.axvline(
        prediction_frequency,
        color="white",
        linewidth=1.0,
        linestyle=":",
        alpha=0.85,
    )

    radio_labels = ("MeerKAT\n1.2 GHz", "ATCA\n6 GHz", "ATCA\n15 GHz")
    alignments = ("left", "center", "center")
    for value, flux, label, alignment in zip(
        radio_frequency,
        radio_flux,
        radio_labels,
        alignments,
        strict=True,
    ):
        axis.text(
            value,
            flux / 1.8,
            label,
            color="white",
            ha=alignment,
            va="top",
            fontsize=9,
        )
    axis.text(
        prediction_frequency,
        prediction_flux * 1.6,
        "1.3 mm\n230 GHz",
        color="white",
        ha="center",
        va="bottom",
        fontsize=10,
    )
    for wavelength, value, flux in zip(
        wavelength_um,
        far_infrared_frequency,
        far_infrared_flux,
        strict=True,
    ):
        axis.text(
            value,
            flux * 1.35,
            f"{wavelength:.0f} $\\mu$m",
            color="white",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    axis.set_xlim(*frequency_limits)
    axis.set_ylim(1e-4, float(np.max(total)) * 2.5)
    axis.text(
        5000,
        axis.get_ylim()[1] / 100.0,
        "Herschel bands",
        color="white",
        ha="right",
        va="top",
        fontsize=12,
    )
    axis.set_xlabel(
        r"Frequency, $\nu$ (GHz)",
        color="white",
        fontsize=12,
    )
    axis.set_ylabel(
        r"Flux density, $S_\nu$ (Jy)",
        color="white",
        fontsize=12,
    )
    axis.tick_params(
        axis="both",
        which="both",
        colors="white",
        direction="in",
        top=True,
        right=True,
        length=5,
    )
    for spine in axis.spines.values():
        spine.set_color("white")

    legend = axis.legend(loc="upper left", frameon=True, fontsize=8.8)
    legend.get_frame().set_facecolor("black")
    legend.get_frame().set_edgecolor("white")
    for text in legend.get_texts():
        text.set_color("white")
    information = (
        rf"$T_e={electron_temperature:.0f}$ K, "
        rf"EM$={emission_measure:.1e}$ pc cm$^{{-6}}$" "\n"
        rf"$S_\nu^{{ff}}=\Omega B_\nu(T_e)(1-e^{{-\tau_\nu}})$" "\n"
        rf"$S_\nu(1.3\,\rm mm)=S_\nu^{{dust}}+S_\nu^{{ff}}$"
    )
    axis.text(
        0.98,
        0.05,
        information,
        transform=axis.transAxes,
        color="white",
        fontsize=8.5,
        ha="right",
        va="bottom",
        bbox={
            "facecolor": "black",
            "edgecolor": "white",
            "boxstyle": "round,pad=0.45",
        },
    )

    outputs = _save_figure(
        figure,
        paths["output_root"] / "figures" / "sed_demonstration_230GHz",
        output_formats=plot["output_formats"],
        dpi=int(plot["dpi"]),
        overwrite=overwrite,
    )
    plt.close(figure)
    return outputs


def plot_products(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> list[Path]:
    """Create all configured continuum-synthesis figures."""
    outputs = plot_continuum_components(
        paths,
        config,
        overwrite=overwrite,
    )
    outputs.extend(
        plot_simba_comparison(
            paths,
            config,
            overwrite=overwrite,
        )
    )
    outputs.extend(
        plot_sed_demonstration(
            paths,
            config,
            overwrite=overwrite,
        )
    )
    print(f"Wrote {len(outputs)} figures to {paths['output_root'] / 'figures'}")
    return outputs


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--stage",
        choices=("inspect", "model", "validate", "plot", "all"),
        default="inspect",
    )
    parser.add_argument(
        "--freefree",
        type=Path,
        help="Override the configured 230 GHz free-free input",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Override the configured task output directory",
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
    paths = resolve_task_paths(
        project,
        config,
        freefree_override=args.freefree,
        output_override=args.output_dir,
    )

    if args.stage in {"inspect", "all"}:
        inspect_inputs(paths, config)
    if args.stage in {"model", "all"}:
        build_models(paths, config, overwrite=args.overwrite)
    if args.stage in {"validate", "all"}:
        build_simba_validation(paths, config, overwrite=args.overwrite)
    if args.stage in {"plot", "all"}:
        plot_products(paths, config, overwrite=args.overwrite)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
