#!/usr/bin/env python
"""Build the project source catalogue from delivered getsf results."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from functools import partial
import json
from pathlib import Path
import tomllib

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Column
from astropy.table import Table, unique, vstack
from astropy.wcs import WCS
import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.ticker import LogLocator

from lancet import ProjectPaths
from lancet.catalogues import (
    effective_source_sizes,
    peak_to_total_flux_ratio,
    rectangle_membership,
    relative_flux_difference,
    sample_at_getsf_peaks,
    select_science_sources,
    uniform_brightness_flux_ratio,
)


DEFAULT_CONFIG = Path("tasks/source_catalogue/task.toml")


def load_task_config(path: Path) -> dict[str, object]:
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
    inputs = config["inputs"]
    outputs = config["outputs"]
    if not isinstance(inputs, Mapping) or not isinstance(outputs, Mapping):
        raise TypeError("Inputs and outputs must be TOML tables")
    paths = {
        name: _project_path(project, value)
        for name, value in inputs.items()
    }
    paths["output_root"] = project.validate_output(str(outputs["root"]))
    for name, path in paths.items():
        if name != "output_root" and not path.is_file():
            raise FileNotFoundError(f"{name} not found: {path}")
    return paths


def inspect_inputs(paths: Mapping[str, Path]) -> None:
    catalogue = Table.read(paths["getsf_catalogue"], format="ascii")
    selected = select_science_sources(catalogue)
    print(f"getsf catalogue: {paths['getsf_catalogue']}")
    print(f"  rows={len(catalogue)}, rows after science cuts={len(selected)}")
    for name in ("intercore_background", "core_background", "temperature"):
        header = fits.getheader(paths[name])
        print(
            f"{name}: {paths[name]} "
            f"shape=({header['NAXIS2']}, {header['NAXIS1']}), "
            f"BITPIX={header['BITPIX']}, BUNIT={header.get('BUNIT')}"
        )
    print(f"Output root: {paths['output_root']}")


def _prepare_output(path: Path, *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists: {path}; pass --overwrite to replace it"
        )
    path.parent.mkdir(parents=True, exist_ok=True)


def output_paths(paths: Mapping[str, Path]) -> dict[str, Path]:
    root = paths["output_root"] / "tables"
    return {
        "catalogue": root / "source_catalogue.ecsv",
        "publication": root / "source_catalogue.tex",
    }


def _roundish(value: float, precision: int = 0) -> str:
    if not np.isfinite(value):
        return "--"
    return f"{int(round(value, precision)):d}"


def build_catalogue(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> Table:
    """Build, tag, and write the analysis source catalogue."""
    catalogue = select_science_sources(
        Table.read(paths["getsf_catalogue"], format="ascii")
    )
    background1 = fits.getdata(paths["intercore_background"], memmap=True)
    background2 = fits.getdata(paths["core_background"], memmap=True)
    temperature, grid_header = fits.getdata(
        paths["temperature"],
        header=True,
        memmap=True,
    )
    x_peak = catalogue["XCO_P"]
    y_peak = catalogue["YCO_P"]
    sampled1 = sample_at_getsf_peaks(background1, x_peak, y_peak)
    sampled2 = sample_at_getsf_peaks(background2, x_peak, y_peak)
    sampled_temperature = sample_at_getsf_peaks(
        temperature,
        x_peak,
        y_peak,
    )

    keys = [
        "NO",
        "XCO_P",
        "YCO_P",
        "WCS_ACOOR",
        "WCS_DCOOR",
        "SIGN",
        "GOOD",
        "SIGNM01",
        "GOODM01",
        "SIGNM02",
        "GOODM02",
        "FXP_BST01",
        "FXP_ERR01",
        "FXT_BST01",
        "FXT_ERR01",
        "FXP_BST02",
        "FXP_ERR02",
        "FXT_BST02",
        "FXT_ERR02",
        "AFWHM01",
        "BFWHM01",
        "THETA01",
        "AFWHM02",
        "BFWHM02",
        "THETA02",
    ]
    output = catalogue[keys]
    output.add_columns(
        [sampled1, sampled2, sampled_temperature],
        names=["Ibkg1", "Ibkg2", "Tenv"],
    )
    observed_major = np.asarray(output["AFWHM01"], dtype=float)
    observed_minor = np.asarray(output["BFWHM01"], dtype=float)
    output["FWHM_OBS_str"] = [
        f"{major:.2f}$\\times${minor:.2f}"
        for major, minor in zip(observed_major, observed_minor)
    ]

    size = config["size"]
    assert isinstance(size, Mapping)
    deconvolved = effective_source_sizes(
        observed_major,
        observed_minor,
        beam_major_arcsec=float(size["beam_major_arcsec"]),
        beam_minor_arcsec=float(size["beam_minor_arcsec"]),
        distance_pc=float(size["distance_pc"]),
        unresolved_threshold_factor=float(
            size["unresolved_threshold_factor"]
        ),
        unresolved_size_factor=float(size["unresolved_size_factor"]),
    )
    output.add_columns(
        [
            deconvolved["major_arcsec"],
            deconvolved["minor_arcsec"],
            deconvolved["position_angle_deg"],
            deconvolved["fwhm_arcsec"],
            deconvolved["fwhm_au"],
            deconvolved["resolved"],
        ],
        names=[
            "Deconv_maj",
            "Deconv_min",
            "Deconv_theta",
            "Deconv_FWHM_angle",
            "Deconv_FWHM_linear",
            "Resolvability",
        ],
    )
    output["FWHM_SRC_str"] = [
        f"{major:.2f}$\\times${minor:.2f}"
        for major, minor in zip(
            deconvolved["major_arcsec"],
            deconvolved["minor_arcsec"],
        )
    ]

    coordinates = SkyCoord(
        output["WCS_ACOOR"],
        output["WCS_DCOOR"],
        unit="deg",
        frame="icrs",
    )
    coordinate_strings = coordinates.fk5.to_string(
        style="hmsdms",
        precision=2,
        sep=":",
    )
    ra_strings, dec_strings = zip(
        *(value.split() for value in coordinate_strings)
    )
    output.add_columns(
        [list(ra_strings), list(dec_strings)],
        indexes=[5, 5],
        names=["RA(hms)", "DEC(dms)"],
    )
    output["FXP_WITH_ERROR"] = [
        f"{peak * 1e3:.1f}({error * 1e3:.1f})"
        for peak, error in zip(output["FXP_BST01"], output["FXP_ERR01"])
    ]
    output["FXT_WITH_ERROR"] = [
        f"{flux * 1e3:.1f}({error * 1e3:.1f})"
        for flux, error in zip(output["FXT_BST01"], output["FXT_ERR01"])
    ]
    output.sort(keys=["SIGN", "GOOD"], reverse=True)
    output.add_column(np.arange(1, len(output) + 1), index=0, name="ID")

    wcs = WCS(grid_header)
    regions = config["regions"]
    assert isinstance(regions, Mapping)
    subsets = []
    for name, label in (
        ("young", "Y"),
        ("intermediate", "I"),
        ("evolved", "E"),
    ):
        mask = rectangle_membership(
            output,
            str(paths[str(regions[name])]),
            wcs=wcs,
        )
        subset = output[mask]
        subset["Subregion"] = np.full(len(subset), label)
        subsets.append(subset)
    merged = unique(vstack(subsets), keys="ID")
    merged.sort("ID")

    products = output_paths(paths)
    for product in products.values():
        _prepare_output(product, overwrite=overwrite)
    merged.write(products["catalogue"], format="ascii.ecsv", overwrite=True)

    publication = merged[
        [
            "ID",
            "RA(hms)",
            "DEC(dms)",
            "Subregion",
            "SIGN",
            "FWHM_OBS_str",
            "THETA01",
            "Deconv_maj",
            "Resolvability",
            "FXP_WITH_ERROR",
            "FXT_WITH_ERROR",
            "Ibkg2",
            "Ibkg1",
        ]
    ]
    publication.write(
        products["publication"],
        format="aastex",
        overwrite=True,
        formats={
            "SIGN": "%d",
            "THETA01": "%d",
            "Deconv_maj": "%.2f",
            "Ibkg2": partial(_roundish, precision=0),
            "Ibkg1": partial(_roundish, precision=0),
        },
    )
    print(
        f"Wrote {products['catalogue']} ({len(merged)} sources; "
        f"{np.count_nonzero(merged['Resolvability'] == 'N')} unresolved)"
    )
    print(f"Wrote {products['publication']}")
    return merged


def _flux_diagnostic_settings(
    config: Mapping[str, object],
) -> Mapping[str, object]:
    settings = config["flux_diagnostics"]
    if not isinstance(settings, Mapping):
        raise TypeError("flux_diagnostics must be a TOML table")
    return settings


def flux_diagnostic_paths(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
) -> dict[str, object]:
    """Return table, summary, and configured diagnostic figure paths."""
    settings = _flux_diagnostic_settings(config)
    formats = [str(value).lower() for value in settings["formats"]]
    unsupported = set(formats) - {"pdf", "png"}
    if unsupported:
        raise ValueError(
            "Unsupported figure formats: " + ", ".join(sorted(unsupported))
        )
    table_root = paths["output_root"] / "tables"
    figure_root = paths["output_root"] / "figures"
    return {
        "table": table_root / "flux_diagnostics.ecsv",
        "summary": table_root / "flux_diagnostics_summary.json",
        "flux_measure": [
            figure_root / f"GetsfResult_FluxMeasure.{extension}"
            for extension in formats
        ],
        "flux_ratio": [
            figure_root / f"GetsfResult_FluxRatio.{extension}"
            for extension in formats
        ],
    }


def _prepare_outputs(
    outputs: Sequence[Path],
    *,
    overwrite: bool,
) -> None:
    for output in outputs:
        if output.exists() and not overwrite:
            raise FileExistsError(
                f"Output already exists: {output}; "
                "pass --overwrite to replace it"
            )
    for output in outputs:
        output.parent.mkdir(parents=True, exist_ok=True)


def calculate_flux_diagnostics(
    catalogue: Table,
    config: Mapping[str, object],
) -> tuple[Table, dict[str, object]]:
    """Calculate the two-array getsf flux diagnostics."""
    settings = _flux_diagnostic_settings(config)
    selected = select_science_sources(catalogue)
    peak_snr_min = float(settings["peak_snr_min"])
    difference_clip = float(settings["relative_difference_clip"])
    if (
        not np.isfinite(peak_snr_min)
        or peak_snr_min <= 0
        or not np.isfinite(difference_clip)
        or difference_clip <= 0
    ):
        raise ValueError(
            "Flux S/N threshold and relative-difference clip must be positive"
        )

    arrays: dict[str, dict[str, np.ndarray]] = {}
    for key, suffix in (("12m7m", "01"), ("12m7mTP", "02")):
        peak = np.asarray(selected[f"FXP_BST{suffix}"], dtype=float)
        peak_error = np.asarray(
            selected[f"FXP_ERR{suffix}"],
            dtype=float,
        )
        total = np.asarray(selected[f"FXT_BST{suffix}"], dtype=float)
        total_error = np.asarray(
            selected[f"FXT_ERR{suffix}"],
            dtype=float,
        )
        major = np.asarray(selected[f"AFWHM{suffix}"], dtype=float)
        minor = np.asarray(selected[f"BFWHM{suffix}"], dtype=float)
        ratio, ratio_error = peak_to_total_flux_ratio(
            peak,
            peak_error,
            total,
            total_error,
        )
        arrays[key] = {
            "peak": peak,
            "peak_error": peak_error,
            "total": total,
            "total_error": total_error,
            "size": np.sqrt(major * minor),
            "ratio": ratio,
            "ratio_error": ratio_error,
        }

    relative_difference = relative_flux_difference(
        arrays["12m7m"]["total"],
        arrays["12m7mTP"]["total"],
    )
    snr_mask = (
        arrays["12m7m"]["peak"]
        >= peak_snr_min * arrays["12m7m"]["peak_error"]
    ) & (
        arrays["12m7mTP"]["peak"]
        >= peak_snr_min * arrays["12m7mTP"]["peak_error"]
    )
    finite = np.isfinite(relative_difference)
    central_mask = (
        snr_mask
        & finite
        & (np.abs(relative_difference) < difference_clip)
    )
    if np.count_nonzero(central_mask) < 2:
        raise ValueError("Too few sources for the relative-flux diagnostic")

    output = Table()
    output["NO"] = selected["NO"]
    for key in ("12m7m", "12m7mTP"):
        values = arrays[key]
        output[f"Peak_{key}"] = Column(
            values["peak"],
            unit=u.Jy / u.beam,
        )
        output[f"err_Peak_{key}"] = Column(
            values["peak_error"],
            unit=u.Jy / u.beam,
        )
        output[f"Integrated_{key}"] = Column(
            values["total"],
            unit=u.Jy,
        )
        output[f"err_Integrated_{key}"] = Column(
            values["total_error"],
            unit=u.Jy,
        )
        output[f"FWHM_{key}"] = Column(
            values["size"],
            unit=u.arcsec,
        )
        output[f"Peak_to_integrated_{key}"] = values["ratio"]
        output[f"err_Peak_to_integrated_{key}"] = values["ratio_error"]
    output["Relative_integrated_flux_difference"] = relative_difference
    output["Peak_SNR_both"] = snr_mask
    output["Central_relative_difference"] = central_mask
    output.meta["reference_measurement"] = "12m7m"
    output.meta["comparison_measurement"] = "12m7mTP"
    output.meta["relative_difference_definition"] = (
        "(Integrated_12m7mTP - Integrated_12m7m) / Integrated_12m7m"
    )
    output.meta["peak_snr_min"] = peak_snr_min
    output.meta["relative_difference_clip"] = difference_clip

    selected_difference = relative_difference[snr_mask & finite]
    central_difference = relative_difference[central_mask]
    central_sigma = float(np.std(central_difference, ddof=0))
    summary: dict[str, object] = {
        "science_source_count": len(selected),
        "relative_difference_source_count": int(
            np.count_nonzero(snr_mask & finite)
        ),
        "central_source_count": int(np.count_nonzero(central_mask)),
        "peak_snr_min": peak_snr_min,
        "relative_difference_clip": difference_clip,
        "relative_difference_definition": output.meta[
            "relative_difference_definition"
        ],
        "relative_difference_median": float(
            np.median(selected_difference)
        ),
        "relative_difference_standard_deviation": float(
            np.std(selected_difference, ddof=0)
        ),
        "central_relative_difference_standard_deviation": central_sigma,
        "display_sigma_percent": int(
            5 * round(100.0 * central_sigma / 5.0)
        ),
    }
    return output, summary


def _configure_flux_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "mathtext.fontset": "dejavusans",
            "font.size": 12,
            "axes.linewidth": 0.8,
            "axes.labelsize": 11,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 4,
            "ytick.major.size": 4,
            "xtick.minor.size": 2,
            "ytick.minor.size": 2,
        }
    )


def _diagnostic_errorbar(
    axis,
    *,
    series: str,
    x: np.ndarray,
    y: np.ndarray,
    xerr: np.ndarray | None,
    yerr: np.ndarray,
    color: str,
    label: str,
    settings: Mapping[str, object],
) -> None:
    filled = bool(settings[f"{series}_filled"])
    lightening = float(settings[f"{series}_error_lightening"])
    if not 0.0 <= lightening <= 1.0:
        raise ValueError("Error-bar lightening must be between zero and one")
    red, green, blue = to_rgb(color)
    error_color = (
        red + (1.0 - red) * lightening,
        green + (1.0 - green) * lightening,
        blue + (1.0 - blue) * lightening,
    )
    edge_color = (
        (0.8 * red, 0.8 * green, 0.8 * blue)
        if filled
        else color
    )
    axis.errorbar(
        x=x,
        y=y,
        xerr=xerr,
        yerr=yerr,
        fmt="o",
        markersize=float(settings[f"{series}_marker_size"]),
        markerfacecolor=color if filled else "white",
        markeredgewidth=0.9,
        markeredgecolor=edge_color,
        ecolor=error_color,
        elinewidth=0.55,
        capsize=1.8,
        alpha=float(settings[f"{series}_marker_alpha"]),
        zorder=int(settings[f"{series}_zorder"]),
        label=label,
    )


def plot_flux_measure(
    table: Table,
    summary: Mapping[str, object],
    outputs: Sequence[Path],
    settings: Mapping[str, object],
) -> None:
    """Plot peak intensity against integrated flux for both arrays."""
    import seaborn as sns

    _configure_flux_plot_style()
    figure, axis = plt.subplots(
        figsize=(
            float(settings["width_inches"]),
            float(settings["height_inches"]),
        ),
        dpi=int(settings["dpi"]),
    )
    lower = float(settings["flux_min_jy"])
    upper = float(settings["flux_max_jy"])
    line = np.logspace(np.log10(lower), np.log10(upper), 100)
    axis.plot(
        line,
        line,
        color="black",
        linestyle="--",
        linewidth=1.5,
        label="y=x: point source",
        zorder=1,
    )
    for key, series, color, label in (
        (
            "12m7mTP",
            "comparison",
            str(settings["comparison_color"]),
            str(settings["comparison_label"]),
        ),
        (
            "12m7m",
            "reference",
            str(settings["reference_color"]),
            str(settings["reference_label"]),
        ),
    ):
        _diagnostic_errorbar(
            axis,
            series=series,
            x=np.asarray(table[f"Peak_{key}"], dtype=float),
            y=np.asarray(table[f"Integrated_{key}"], dtype=float),
            xerr=np.asarray(table[f"err_Peak_{key}"], dtype=float),
            yerr=np.asarray(
                table[f"err_Integrated_{key}"],
                dtype=float,
            ),
            color=color,
            label=label,
            settings=settings,
        )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlim(lower, upper)
    axis.set_ylim(lower, upper)
    axis.xaxis.set_major_locator(LogLocator(base=10))
    axis.yaxis.set_major_locator(LogLocator(base=10))
    axis.set_xlabel(r"Peak intensity (Jy beam$^{-1}$)")
    axis.set_ylabel("Integrated flux (Jy)")
    axis.legend(loc="upper left", frameon=True)

    inset = figure.add_axes([0.60, 0.18, 0.29, 0.35])
    relative = np.asarray(
        table["Relative_integrated_flux_difference"],
        dtype=float,
    )
    mask = np.asarray(table["Peak_SNR_both"], dtype=bool)
    sns.kdeplot(
        relative[mask],
        ax=inset,
        fill=True,
        alpha=0.2,
        color=str(settings["comparison_color"]),
    )
    inset.vlines(
        relative[mask],
        ymin=0,
        ymax=0.2,
        linewidth=0.5,
        alpha=0.5,
        color="gray",
    )
    inset.set_xticks([-0.3, 0.0, 0.3, 1.0])
    inset.set_xticklabels(["-30%", "0", "30%", "100%"])
    inset.set_xlabel("")
    inset.set_ylabel("Density")
    inset.tick_params(axis="x", direction="out")
    inset.text(
        0.47,
        0.96,
        "Relative diff\n"
        + rf"$\sigma\approx{summary['display_sigma_percent']}\%$",
        transform=inset.transAxes,
        ha="left",
        va="top",
    )
    for output in outputs:
        figure.savefig(
            output,
            bbox_inches="tight",
            dpi=int(settings["dpi"]),
        )
        print(f"Wrote {output}")
    plt.close(figure)


def plot_flux_ratio(
    table: Table,
    outputs: Sequence[Path],
    settings: Mapping[str, object],
) -> None:
    """Plot peak/total flux ratio against circularized source size."""
    _configure_flux_plot_style()
    figure, axis = plt.subplots(
        figsize=(
            float(settings["width_inches"]),
            float(settings["height_inches"]),
        ),
        dpi=int(settings["dpi"]),
    )
    size_min = float(settings["size_min_arcsec"])
    size_max = float(settings["size_max_arcsec"])
    theoretical_size = np.linspace(
        float(settings["beam_fwhm_arcsec"]),
        max(3.5, size_max),
        100,
    )
    theoretical_ratio = uniform_brightness_flux_ratio(
        theoretical_size,
        beam_fwhm_arcsec=float(settings["beam_fwhm_arcsec"]),
    )
    axis.plot(
        theoretical_size,
        theoretical_ratio,
        color="black",
        linestyle="-.",
        linewidth=2.0,
        label="uniform brightness",
        zorder=1,
    )
    for key, series, color, label in (
        (
            "12m7mTP",
            "comparison",
            str(settings["comparison_color"]),
            str(settings["comparison_label"]),
        ),
        (
            "12m7m",
            "reference",
            str(settings["reference_color"]),
            str(settings["reference_label"]),
        ),
    ):
        _diagnostic_errorbar(
            axis,
            series=series,
            x=np.asarray(table[f"FWHM_{key}"], dtype=float),
            y=np.asarray(
                table[f"Peak_to_integrated_{key}"],
                dtype=float,
            ),
            xerr=None,
            yerr=np.asarray(
                table[f"err_Peak_to_integrated_{key}"],
                dtype=float,
            ),
            color=color,
            label=label,
            settings=settings,
        )
    axis.set_yscale("log")
    axis.set_xlim(size_min, size_max)
    axis.set_ylim(
        float(settings["ratio_min"]),
        float(settings["ratio_max"]),
    )
    axis.set_xlabel("Source size FWHM (arcsec)")
    axis.set_ylabel("Peak intensity / Integrated flux")
    axis.legend(loc="upper right", frameon=True)
    for output in outputs:
        figure.savefig(
            output,
            bbox_inches="tight",
            dpi=int(settings["dpi"]),
        )
        print(f"Wrote {output}")
    plt.close(figure)


def build_flux_diagnostics(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool,
) -> dict[str, object]:
    """Write the flux-diagnostic data, summary, and figures."""
    catalogue = Table.read(paths["getsf_catalogue"], format="ascii")
    table, summary = calculate_flux_diagnostics(catalogue, config)
    products = flux_diagnostic_paths(paths, config)
    flux_measure = products["flux_measure"]
    flux_ratio = products["flux_ratio"]
    assert isinstance(flux_measure, list)
    assert isinstance(flux_ratio, list)
    outputs = [
        products["table"],
        products["summary"],
        *flux_measure,
        *flux_ratio,
    ]
    _prepare_outputs(outputs, overwrite=overwrite)
    table_output = products["table"]
    summary_output = products["summary"]
    assert isinstance(table_output, Path)
    assert isinstance(summary_output, Path)
    table.write(table_output, format="ascii.ecsv", overwrite=True)
    summary_output.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {table_output}")
    print(f"Wrote {summary_output}")
    settings = _flux_diagnostic_settings(config)
    plot_flux_measure(table, summary, flux_measure, settings)
    plot_flux_ratio(table, flux_ratio, settings)
    return products


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--stage",
        choices=("inspect", "build", "diagnostics", "all"),
        default="inspect",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    project = ProjectPaths.discover()
    config = load_task_config(_project_path(project, args.config))
    paths = resolve_task_paths(project, config)
    if args.stage == "inspect":
        inspect_inputs(paths)
    elif args.stage == "build":
        build_catalogue(paths, config, overwrite=args.overwrite)
    elif args.stage == "diagnostics":
        build_flux_diagnostics(paths, config, overwrite=args.overwrite)
    else:
        build_catalogue(paths, config, overwrite=args.overwrite)
        build_flux_diagnostics(paths, config, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
