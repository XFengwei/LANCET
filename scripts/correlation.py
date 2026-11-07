#!/usr/bin/env python
"""Measure and display correlations between core and environment properties."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
import json
from pathlib import Path
import tomllib

from astropy.table import Table
import matplotlib
import numpy as np
from scipy.stats import gaussian_kde

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter, MaxNLocator

from lancet import ProjectPaths
from lancet.correlation import (
    compare_powerlaw_methods,
    pairwise_correlations,
    prepare_correlation_data,
    prepare_powerlaw_data,
    run_powerlaw_mcmc,
)


DEFAULT_CONFIG = Path("tasks/correlation/task.toml")


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
) -> dict[str, Path]:
    """Resolve the generated input catalogue and policy-safe output root."""
    inputs = config["inputs"]
    outputs = config["outputs"]
    if not isinstance(inputs, Mapping) or not isinstance(outputs, Mapping):
        raise TypeError("Inputs and outputs must be TOML tables")
    paths = {
        name: _project_path(project, value)
        for name, value in inputs.items()
    }
    paths["output_root"] = project.validate_output(str(outputs["root"]))
    for name in inputs:
        if not paths[name].is_file():
            raise FileNotFoundError(f"{name} not found: {paths[name]}")
    return paths


def _mapping(
    config: Mapping[str, object],
    name: str,
) -> Mapping[str, object]:
    value = config[name]
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a TOML table")
    return value


def _sequence(
    value: object,
    *,
    name: str,
) -> list[object]:
    if not isinstance(value, list):
        raise TypeError(f"{name} must be a TOML array")
    return value


def _load_catalogue(path: Path) -> Table:
    return Table.read(path, format="ascii.ecsv")


def _configured_data(
    catalogue: Table,
    config: Mapping[str, object],
) -> tuple[dict[str, np.ndarray], list[str], set[str], dict[str, str]]:
    variables = _mapping(config, "variables")
    columns = [
        str(value)
        for value in _sequence(variables["columns"], name="variables.columns")
    ]
    log_columns = {
        str(value)
        for value in _sequence(
            variables["log_columns"],
            name="variables.log_columns",
        )
    }
    labels_config = _mapping(config, "labels")
    labels = {
        column: str(labels_config.get(column, column))
        for column in columns
    }
    data = prepare_correlation_data(
        catalogue,
        columns,
        log_columns=sorted(log_columns),
        resolved_only_columns=variables.get("resolved_only_columns", []),
    )
    return data, columns, log_columns, labels


def _regression_data(
    catalogue: Table,
    config: Mapping[str, object],
) -> dict[str, np.ndarray | int]:
    regression = _mapping(config, "regression")
    return prepare_powerlaw_data(
        catalogue,
        x_column=str(regression["x_column"]),
        y_column=str(regression["y_column"]),
        y_error_column=str(regression["y_error_column"]),
        x_fractional_error=float(regression["x_fractional_error"]),
        resolved_only_columns=_mapping(config, "variables").get("resolved_only_columns", []),
    )


def inspect_inputs(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
) -> None:
    """Report catalogue, variable, regression, and output readiness."""
    catalogue = _load_catalogue(paths["catalogue"])
    data, columns, log_columns, _ = _configured_data(catalogue, config)
    statistics = pairwise_correlations(data)
    regression = _mapping(config, "regression")
    fit_data = _regression_data(catalogue, config)
    x_name = str(regression["x_column"])
    y_name = str(regression["y_column"])
    x_index = columns.index(x_name)
    y_index = columns.index(y_name)
    print("Correlation input:")
    print(f"  catalogue: {paths['catalogue']} (rows={len(catalogue)})")
    print(f"Variables ({len(columns)}): {', '.join(columns)}")
    print(f"Log10 variables ({len(log_columns)}): {', '.join(sorted(log_columns))}")
    print(
        "Primary regression: "
        f"{y_name} = A * {x_name}^alpha; "
        f"x fractional uncertainty={regression['x_fractional_error']}; "
        f"valid rows={fit_data['n_valid']}/{fit_data['n_total']}"
    )
    print(
        "Configured pair: "
        f"Pearson={statistics['pearson'][y_index, x_index]:.6f}; "
        f"Spearman={statistics['spearman'][y_index, x_index]:.6f}"
    )
    print(
        "MCMC: "
        f"walkers={regression['n_walkers']}, "
        f"steps={regression['n_steps']}, "
        f"burn-in={regression['burn_in']}, thin={regression['thin']}, "
        f"seed={regression['seed']}"
    )
    print(f"Output root: {paths['output_root']}")


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


def correlation_table_path(paths: Mapping[str, Path]) -> Path:
    return paths["output_root"] / "tables" / "pairwise_correlations.ecsv"


def regression_summary_path(paths: Mapping[str, Path]) -> Path:
    return (
        paths["output_root"]
        / "tables"
        / "intercore_surface_density_powerlaw.json"
    )


def posterior_path(paths: Mapping[str, Path]) -> Path:
    return (
        paths["output_root"]
        / "posterior"
        / "intercore_surface_density_powerlaw.npz"
    )


def _figure_paths(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    stem: str,
) -> list[Path]:
    figures = _mapping(config, "figures")
    formats = [
        str(value).lower()
        for value in _sequence(figures["formats"], name="figures.formats")
    ]
    unsupported = set(formats) - {"pdf", "png"}
    if unsupported:
        raise ValueError(
            "Unsupported figure formats: " + ", ".join(sorted(unsupported))
        )
    return [
        paths["output_root"] / "figures" / f"{stem}.{suffix}"
        for suffix in formats
    ]


def build_correlation_table(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> tuple[Table, dict[str, np.ndarray | list[str]], dict[str, np.ndarray]]:
    """Calculate and write all configured pairwise correlations."""
    catalogue = _load_catalogue(paths["catalogue"])
    data, columns, log_columns, _ = _configured_data(catalogue, config)
    statistics = pairwise_correlations(data)
    rows = []
    for row, y_name in enumerate(columns):
        for column, x_name in enumerate(columns):
            rows.append(
                (
                    x_name,
                    y_name,
                    int(statistics["n_points"][row, column]),
                    float(statistics["pearson"][row, column]),
                    float(statistics["pearson_p"][row, column]),
                    float(statistics["spearman"][row, column]),
                    float(statistics["spearman_p"][row, column]),
                )
            )
    output = Table(
        rows=rows,
        names=(
            "x_variable",
            "y_variable",
            "n_points",
            "pearson_r",
            "pearson_p",
            "spearman_rho",
            "spearman_p",
        ),
    )
    output.meta["input_catalogue"] = str(paths["catalogue"])
    output.meta["resolved_only_columns"] = ",".join(_mapping(config, "variables").get("resolved_only_columns", []))
    output.meta["log10_variables"] = ",".join(
        column for column in columns if column in log_columns
    )
    destination = correlation_table_path(paths)
    _prepare_outputs([destination], overwrite=overwrite)
    output.write(destination, format="ascii.ecsv", overwrite=True)
    print(f"Wrote {destination} ({len(output)} ordered pairs)")
    return output, statistics, data


def _log_tick_formatter() -> FuncFormatter:
    def formatter(value: float, _position: float) -> str:
        physical = 10.0**value
        if not np.isfinite(physical) or physical == 0:
            return ""
        if physical >= 1000 or physical < 0.01:
            return f"{physical:.1e}"
        if physical >= 100:
            return f"{physical:.0f}"
        if physical >= 10:
            return f"{physical:.1f}"
        if physical >= 1:
            return f"{physical:.2f}"
        return f"{physical:.2g}"

    return FuncFormatter(formatter)


def _kde_levels(
    density: np.ndarray,
    fractions: Sequence[float],
) -> np.ndarray:
    ordered = np.sort(density.ravel())[::-1]
    cumulative = np.cumsum(ordered)
    if len(cumulative) == 0 or cumulative[-1] <= 0:
        return np.asarray([])
    cumulative /= cumulative[-1]
    values = [
        ordered[np.argmin(np.abs(cumulative - fraction))]
        for fraction in fractions
    ]
    return np.sort(np.unique(values))


def create_corner_figure(
    data: Mapping[str, np.ndarray],
    statistics: Mapping[str, np.ndarray | list[str]],
    labels: Mapping[str, str],
    log_columns: set[str],
    config: Mapping[str, object],
) -> plt.Figure:
    """Create the detailed histogram/scatter/KDE correlation matrix."""
    settings = _mapping(config, "corner_plot")
    keys = [str(value) for value in statistics["keys"]]
    size = float(settings["size_inches"])
    figure, axes = plt.subplots(
        len(keys),
        len(keys),
        figsize=(size, size),
        squeeze=False,
    )
    figure.subplots_adjust(wspace=0.06, hspace=0.06)
    threshold = float(settings["highlight_threshold"])
    padding = float(settings["padding_fraction"])
    contour_fractions = [
        float(value)
        for value in _sequence(
            settings["contour_fractions"],
            name="corner_plot.contour_fractions",
        )
    ]
    formatter = _log_tick_formatter()

    for row, y_name in enumerate(keys):
        for column, x_name in enumerate(keys):
            axis = axes[row, column]
            x = np.asarray(data[x_name], dtype=float)
            y = np.asarray(data[y_name], dtype=float)
            valid = np.isfinite(x) & np.isfinite(y)
            x_valid = x[valid]
            y_valid = y[valid]
            pearson = float(statistics["pearson"][row, column])
            spearman = float(statistics["spearman"][row, column])
            highlighted = (
                np.isfinite(pearson)
                and abs(pearson) > threshold
            ) or (
                np.isfinite(spearman)
                and abs(spearman) > threshold
            )
            for spine in axis.spines.values():
                spine.set_linewidth(0.8)
                spine.set_color("0.5")

            if row == column:
                if len(x_valid):
                    axis.hist(
                        x_valid,
                        bins=int(settings["histogram_bins"]),
                        histtype="stepfilled",
                        color="lightgray",
                        edgecolor="black",
                        alpha=0.8,
                    )
                axis.set_yticks([])
            elif row > column:
                axis.scatter(
                    x_valid,
                    y_valid,
                    s=float(settings["scatter_size"]),
                    alpha=float(settings["scatter_alpha"]),
                    color="black",
                    linewidths=0,
                )
                if len(x_valid) >= 5:
                    x_range = np.ptp(x_valid) or 1.0
                    y_range = np.ptp(y_valid) or 1.0
                    x_limits = (
                        np.min(x_valid) - padding * x_range,
                        np.max(x_valid) + padding * x_range,
                    )
                    y_limits = (
                        np.min(y_valid) - padding * y_range,
                        np.max(y_valid) + padding * y_range,
                    )
                    axis.set_xlim(*x_limits)
                    axis.set_ylim(*y_limits)
                    try:
                        density_estimator = gaussian_kde(
                            np.vstack([x_valid, y_valid])
                        )
                        xx, yy = np.mgrid[
                            x_limits[0]:x_limits[1]:100j,
                            y_limits[0]:y_limits[1]:100j,
                        ]
                        density = density_estimator(
                            np.vstack([xx.ravel(), yy.ravel()])
                        ).reshape(xx.shape)
                        levels = _kde_levels(density, contour_fractions)
                        if len(levels):
                            axis.contour(
                                xx,
                                yy,
                                density,
                                levels=levels,
                                colors=str(settings["kde_color"]),
                                linewidths=1.1,
                            )
                    except (ValueError, np.linalg.LinAlgError):
                        pass
                if highlighted:
                    for spine in axis.spines.values():
                        spine.set_linewidth(2.2)
                        spine.set_color(str(settings["highlight_color"]))
                axis.text(
                    0.04,
                    0.96,
                    f"P = {pearson:.2f}\nS = {spearman:.2f}",
                    transform=axis.transAxes,
                    ha="left",
                    va="top",
                    fontsize=9,
                    color=(
                        str(settings["highlight_color"])
                        if highlighted
                        else "black"
                    ),
                    fontweight="bold" if highlighted else "normal",
                    bbox={
                        "facecolor": "white",
                        "edgecolor": "none",
                        "alpha": 0.65,
                        "pad": 1.5,
                    },
                )
            else:
                axis.axis("off")
                axis.text(
                    0.5,
                    0.55,
                    (
                        f"N = {len(x_valid)}\n"
                        f"Pearson = {pearson:.2f}\n"
                        f"Spearman = {spearman:.2f}"
                    ),
                    ha="center",
                    va="center",
                    fontsize=10,
                    color=(
                        str(settings["highlight_color"])
                        if highlighted
                        else "black"
                    ),
                    fontweight="bold" if highlighted else "normal",
                    transform=axis.transAxes,
                )

            if x_name in log_columns:
                axis.xaxis.set_major_formatter(formatter)
                axis.xaxis.set_major_locator(MaxNLocator(nbins=3))
            if row != column and y_name in log_columns:
                axis.yaxis.set_major_formatter(formatter)
                axis.yaxis.set_major_locator(MaxNLocator(nbins=4))
            axis.tick_params(axis="both", direction="in")
            if row < len(keys) - 1:
                axis.set_xticklabels([])
            else:
                axis.set_xlabel(labels[x_name], fontsize=10)
            if column > 0:
                axis.set_yticklabels([])
            else:
                axis.set_ylabel(labels[y_name], fontsize=10)
    return figure


def create_compact_figure(
    statistics: Mapping[str, np.ndarray | list[str]],
    labels: Mapping[str, str],
    config: Mapping[str, object],
) -> plt.Figure:
    """Create the compact lower-triangle Spearman tile matrix."""
    settings = _mapping(config, "compact_plot")
    keys = [str(value) for value in statistics["keys"]]
    correlation = np.asarray(statistics["spearman"], dtype=float)
    figure, axis = plt.subplots(
        figsize=(
            float(settings["width_inches"]),
            float(settings["height_inches"]),
        )
    )
    normalizer = plt.Normalize(vmin=-1.0, vmax=1.0)
    mapper = plt.cm.ScalarMappable(
        norm=normalizer,
        cmap=str(settings["colormap"]),
    )
    threshold = float(settings["outline_threshold"])
    minimum = float(settings["minimum_box_fraction"])
    maximum = float(settings["maximum_box_fraction"])
    for row in range(len(keys)):
        for column in range(len(keys)):
            if row <= column or not np.isfinite(correlation[row, column]):
                continue
            coefficient = correlation[row, column]
            fraction = minimum + (maximum - minimum) * abs(coefficient)
            strong = abs(coefficient) > threshold
            axis.add_patch(
                Rectangle(
                    (column - fraction / 2, row - fraction / 2),
                    fraction,
                    fraction,
                    facecolor=mapper.to_rgba(coefficient),
                    edgecolor="black" if strong else "none",
                    linewidth=1.3 if strong else 0,
                )
            )
            if bool(settings["show_values"]):
                axis.text(
                    column,
                    row,
                    f"{coefficient:.2f}",
                    ha="center",
                    va="center",
                    fontsize=float(settings["value_fontsize"]),
                )
    count = len(keys)
    axis.set_xlim(-0.5, count - 1.5)
    axis.set_ylim(0.5, count - 0.5)
    axis.invert_yaxis()
    axis.set_xticks(
        np.arange(count - 1),
        [labels[key] for key in keys[:-1]],
        rotation=float(settings["label_rotation_degrees"]),
        ha="right",
        fontsize=float(settings["tick_fontsize"]),
    )
    axis.set_yticks(
        np.arange(1, count),
        [labels[key] for key in keys[1:]],
        fontsize=float(settings["tick_fontsize"]),
    )
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.tick_params(direction="out", length=4)
    axis.set_aspect("equal")
    mapper.set_array([])
    color_axis = axis.inset_axes([0.40, 0.84, 0.55, 0.035])
    colorbar = figure.colorbar(
        mapper,
        cax=color_axis,
        orientation="horizontal",
    )
    colorbar.ax.tick_params(
        labelsize=float(settings["tick_fontsize"]) - 1,
        pad=1,
    )
    colorbar.set_label(
        r"Spearman correlation $\rho$",
        fontsize=float(settings["tick_fontsize"]),
    )
    colorbar.ax.xaxis.set_label_position("top")
    return figure


def _save_figure(
    figure: plt.Figure,
    outputs: Sequence[Path],
    *,
    dpi: int,
) -> None:
    for output in outputs:
        figure.savefig(
            output,
            bbox_inches="tight",
            dpi=dpi,
            facecolor="white",
        )
        print(f"Wrote {output}")
    plt.close(figure)


def build_matrices(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> None:
    """Write the correlation table and the detailed and compact figures."""
    _, statistics, data = build_correlation_table(
        paths,
        config,
        overwrite=overwrite,
    )
    catalogue = _load_catalogue(paths["catalogue"])
    _, _, log_columns, labels = _configured_data(catalogue, config)
    corner_outputs = _figure_paths(
        paths,
        config,
        "core_environment_corner",
    )
    compact_outputs = _figure_paths(
        paths,
        config,
        "core_environment_compact_spearman",
    )
    _prepare_outputs(
        [*corner_outputs, *compact_outputs],
        overwrite=overwrite,
    )
    corner = create_corner_figure(
        data,
        statistics,
        labels,
        log_columns,
        config,
    )
    compact = create_compact_figure(statistics, labels, config)
    figures = _mapping(config, "figures")
    _save_figure(
        corner,
        corner_outputs,
        dpi=int(figures["corner_dpi"]),
    )
    _save_figure(
        compact,
        compact_outputs,
        dpi=int(figures["dpi"]),
    )


def _finite_or_none(values: np.ndarray) -> list[float | None]:
    return [
        float(value) if np.isfinite(value) else None
        for value in np.asarray(values, dtype=float)
    ]


def _regression_summary(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    data: Mapping[str, np.ndarray | int],
    comparison: Mapping[str, object],
    mcmc: Mapping[str, np.ndarray | float | int],
) -> dict[str, object]:
    regression = _mapping(config, "regression")
    quantiles = np.asarray(mcmc["quantiles"], dtype=float)
    scatter = np.asarray(mcmc["scatter_quantiles"], dtype=float)
    ols = comparison["ols"]
    theil = comparison["theil"]
    odr_fit = comparison["odr"]
    fractions = np.asarray(
        comparison["sensitivity_fractions"],
        dtype=float,
    )
    return {
        "input_catalogue": str(paths["catalogue"]),
        "model": "log10(y) = alpha * log10(x) + beta",
        "x_column": str(regression["x_column"]),
        "y_column": str(regression["y_column"]),
        "y_error_column": str(regression["y_error_column"]),
        "x_fractional_error": float(
            regression["x_fractional_error"]
        ),
        "n_total": int(data["n_total"]),
        "n_valid": int(data["n_valid"]),
        "pearson_r": float(comparison["pearson"].statistic),
        "pearson_p": float(comparison["pearson"].pvalue),
        "spearman_rho": float(comparison["spearman"].statistic),
        "spearman_p": float(comparison["spearman"].pvalue),
        "mcmc": {
            "alpha_p2p5_p16_p50_p84_p97p5": quantiles[:, 0].tolist(),
            "beta_p2p5_p16_p50_p84_p97p5": quantiles[:, 1].tolist(),
            "sigma_intrinsic_p2p5_p16_p50_p84_p97p5_dex": (
                scatter.tolist()
            ),
            "acceptance_fraction": float(
                mcmc["acceptance_fraction"]
            ),
            "autocorrelation_time": _finite_or_none(
                np.asarray(mcmc["autocorrelation_time"])
            ),
            "n_walkers": int(regression["n_walkers"]),
            "n_steps": int(regression["n_steps"]),
            "burn_in": int(mcmc["burn_in"]),
            "thin": int(mcmc["thin"]),
            "seed": int(mcmc["seed"]),
        },
        "mle_alpha_beta_lnsigma": np.asarray(
            mcmc["mle"],
            dtype=float,
        ).tolist(),
        "ols": {
            "alpha": float(ols.slope),
            "beta": float(ols.intercept),
            "alpha_stderr": float(ols.stderr),
        },
        "theil_sen": {
            "alpha": float(theil.slope),
            "beta": float(theil.intercept),
            "alpha_95_low": float(theil.low_slope),
            "alpha_95_high": float(theil.high_slope),
        },
        "odr": {
            "alpha": float(odr_fit.beta[0]),
            "beta": float(odr_fit.beta[1]),
            "alpha_stderr": float(odr_fit.sd_beta[0]),
            "beta_stderr": float(odr_fit.sd_beta[1]),
        },
        "sensitivity": [
            {
                "x_fractional_error": float(fraction),
                "likelihood_alpha": float(likelihood),
                "effective_chi2_alpha": float(effective),
            }
            for fraction, likelihood, effective in zip(
                fractions,
                comparison["likelihood_slopes"],
                comparison["effective_chi2_slopes"],
            )
        ],
    }


def create_regression_figure(
    data: Mapping[str, np.ndarray | int],
    comparison: Mapping[str, object],
    mcmc: Mapping[str, np.ndarray | float | int],
    config: Mapping[str, object],
) -> plt.Figure:
    """Create the method-comparison regression diagnostic."""
    regression = _mapping(config, "regression")
    x = np.asarray(data["x"], dtype=float)
    y = np.asarray(data["y"], dtype=float)
    x_error = np.asarray(data["x_error"], dtype=float)
    y_error = np.asarray(data["y_error"], dtype=float)
    samples = np.asarray(mcmc["samples"], dtype=float)
    quantiles = np.asarray(mcmc["quantiles"], dtype=float)
    alpha = quantiles[2, 0]
    beta = quantiles[2, 1]
    x_limits = [
        float(value)
        for value in _sequence(
            regression["x_limits"],
            name="regression.x_limits",
        )
    ]
    y_limits = [
        float(value)
        for value in _sequence(
            regression["y_limits"],
            name="regression.y_limits",
        )
    ]
    x_grid = np.logspace(
        np.log10(x_limits[0]),
        np.log10(x_limits[1]),
        400,
    )
    random = np.random.default_rng(int(regression["plot_seed"]))
    draw_count = min(int(regression["posterior_draws"]), len(samples))
    indices = random.choice(len(samples), draw_count, replace=False)
    curves = (
        10.0 ** samples[indices, 1, None]
        * x_grid[None, :] ** samples[indices, 0, None]
    )
    curve_16, curve_84 = np.percentile(curves, [16, 84], axis=0)
    ols = comparison["ols"]
    theil = comparison["theil"]
    odr_fit = comparison["odr"]
    reference_slope = float(regression["reference_slope"])
    reference_intercept = float(regression["reference_intercept"])

    figure = plt.figure(figsize=(13.0, 6.7), constrained_layout=True)
    grid = figure.add_gridspec(
        2,
        2,
        width_ratios=(1.65, 1.0),
        height_ratios=(1.0, 1.0),
    )
    data_axis = figure.add_subplot(grid[:, 0])
    slope_axis = figure.add_subplot(grid[0, 1])
    sensitivity_axis = figure.add_subplot(grid[1, 1])

    data_axis.errorbar(
        x,
        y,
        xerr=x_error,
        yerr=y_error,
        fmt="none",
        ecolor="0.55",
        elinewidth=0.45,
        alpha=0.17,
    )
    data_axis.scatter(
        x,
        y,
        s=19,
        facecolors="white",
        edgecolors="#666666",
        linewidths=0.6,
        alpha=0.7,
        label=f"Individual cores (N={len(x)})",
    )
    data_axis.fill_between(
        x_grid,
        curve_16,
        curve_84,
        color="#E66101",
        alpha=0.13,
        linewidth=0,
        label="MCMC 68% credible band",
    )
    data_axis.plot(
        x_grid,
        10.0**beta * x_grid**alpha,
        color="#E66101",
        linewidth=2.6,
        label=rf"MCMC: $\alpha={alpha:.3f}$",
    )
    data_axis.plot(
        x_grid,
        10.0**ols.intercept * x_grid**ols.slope,
        color="#1B9E77",
        linewidth=2.0,
        label=rf"OLS $y|x$: $\alpha={ols.slope:.3f}$",
    )
    data_axis.plot(
        x_grid,
        10.0**odr_fit.beta[1] * x_grid**odr_fit.beta[0],
        color="#7570B3",
        linewidth=2.0,
        label=rf"ODR: $\alpha={odr_fit.beta[0]:.3f}$",
    )
    data_axis.plot(
        x_grid,
        10.0**reference_intercept * x_grid**reference_slope,
        color="#D73027",
        linewidth=2.0,
        linestyle="--",
        label=rf"Reference line: $\alpha\simeq{reference_slope:.2f}$",
    )
    sorted_indices = np.argsort(np.asarray(data["log_x"]))
    groups = np.array_split(
        sorted_indices,
        int(regression["median_bins"]),
    )
    median_x = np.asarray(
        [np.median(x[group]) for group in groups if len(group)]
    )
    median_y = np.asarray(
        [np.median(y[group]) for group in groups if len(group)]
    )
    data_axis.scatter(
        median_x,
        median_y,
        marker="D",
        s=45,
        facecolors="black",
        edgecolors="white",
        linewidths=0.7,
        label="Binned medians",
    )
    data_axis.set_xscale("log")
    data_axis.set_yscale("log")
    data_axis.set_xlim(*x_limits)
    data_axis.set_ylim(*y_limits)
    data_axis.set_xlabel(
        r"$\Sigma_{\rm env}^{\rm ic}$ (g cm$^{-2}$)"
    )
    data_axis.set_ylabel(
        r"$\Sigma_{\rm core}^{\rm ic}$ (g cm$^{-2}$)"
    )
    data_axis.grid(True, which="major", color="0.86", linewidth=0.65)
    data_axis.grid(True, which="minor", color="0.93", linewidth=0.4)
    data_axis.tick_params(which="both", direction="in", top=True, right=True)
    data_axis.text(
        0.025,
        0.975,
        (
            rf"Pearson $r={comparison['pearson'].statistic:.3f}$"
            "\n"
            rf"Spearman $\rho={comparison['spearman'].statistic:.3f}$"
        ),
        transform=data_axis.transAxes,
        ha="left",
        va="top",
        bbox={"facecolor": "white", "edgecolor": "0.75", "alpha": 0.9},
    )
    data_axis.legend(
        loc="lower right",
        frameon=True,
        facecolor="white",
        edgecolor="0.8",
        framealpha=0.94,
        fontsize=8.5,
    )
    data_axis.set_title(
        "A   Data and competing power-law fits",
        loc="left",
        fontweight="bold",
    )

    labels = ["MCMC", r"OLS $y|x$", "Theil-Sen", "Reference", "ODR"]
    values = np.asarray(
        [alpha, ols.slope, theil.slope, reference_slope, odr_fit.beta[0]]
    )
    lows = np.asarray(
        [
            quantiles[0, 0],
            ols.slope - 1.96 * ols.stderr,
            theil.low_slope,
            np.nan,
            odr_fit.beta[0] - 1.96 * odr_fit.sd_beta[0],
        ]
    )
    highs = np.asarray(
        [
            quantiles[4, 0],
            ols.slope + 1.96 * ols.stderr,
            theil.high_slope,
            np.nan,
            odr_fit.beta[0] + 1.96 * odr_fit.sd_beta[0],
        ]
    )
    colors = ["#E66101", "#1B9E77", "#4D4D4D", "#D73027", "#7570B3"]
    positions = np.arange(len(labels))[::-1]
    for index, position in enumerate(positions):
        if np.isfinite(lows[index]):
            slope_axis.plot(
                [lows[index], highs[index]],
                [position, position],
                color=colors[index],
                linewidth=3.0,
                solid_capstyle="round",
            )
            marker = "o"
        else:
            marker = "x"
        slope_axis.scatter(
            values[index],
            position,
            color=colors[index],
            marker=marker,
            s=60,
            zorder=3,
        )
    slope_axis.axvline(1.0, color="0.45", linestyle=":", linewidth=1.2)
    slope_axis.set_yticks(positions, labels)
    slope_axis.set_xlabel(r"Power-law slope $\alpha$")
    slope_axis.grid(True, axis="x", color="0.88", linewidth=0.6)
    slope_axis.tick_params(direction="in", top=True)
    slope_axis.set_title(
        "B   Slope comparison (95% intervals)",
        loc="left",
        fontweight="bold",
    )

    fractions = np.asarray(
        comparison["sensitivity_fractions"],
        dtype=float,
    )
    sensitivity_axis.plot(
        fractions,
        comparison["likelihood_slopes"],
        color="#E66101",
        marker="o",
        linewidth=2.2,
        label="Likelihood with intrinsic scatter",
    )
    sensitivity_axis.plot(
        fractions,
        comparison["effective_chi2_slopes"],
        color="#7570B3",
        marker="s",
        linewidth=2.0,
        linestyle="--",
        label=r"Effective $\chi^2$",
    )
    adopted = float(regression["x_fractional_error"])
    sensitivity_axis.axvline(
        adopted,
        color="0.35",
        linestyle=":",
        linewidth=1.2,
        label=f"Adopted fraction = {adopted:.2f}",
    )
    sensitivity_axis.set_xlim(
        float(np.min(fractions)) - 0.02,
        float(np.max(fractions)) + 0.02,
    )
    sensitivity_axis.set_xlabel("Fractional uncertainty assigned to x")
    sensitivity_axis.set_ylabel(r"Fitted slope $\alpha$")
    sensitivity_axis.grid(True, color="0.88", linewidth=0.6)
    sensitivity_axis.tick_params(direction="in", top=True, right=True)
    sensitivity_axis.legend(loc="best", frameon=False, fontsize=8.5)
    sensitivity_axis.set_title(
        "C   Sensitivity to the environmental error model",
        loc="left",
        fontweight="bold",
    )
    figure.suptitle(
        "Intercore surface-density regression",
        fontsize=14,
        fontweight="bold",
    )
    return figure


def create_publication_figure(
    data: Mapping[str, np.ndarray | int],
    comparison: Mapping[str, object],
    mcmc: Mapping[str, np.ndarray | float | int],
    config: Mapping[str, object],
) -> plt.Figure:
    """Create the clean single-panel power-law figure for publication."""
    regression = _mapping(config, "regression")
    settings = _mapping(config, "publication_plot")
    x = np.asarray(data["x"], dtype=float)
    y = np.asarray(data["y"], dtype=float)
    x_error = np.asarray(data["x_error"], dtype=float)
    y_error = np.asarray(data["y_error"], dtype=float)
    samples = np.asarray(mcmc["samples"], dtype=float)
    quantiles = np.asarray(mcmc["quantiles"], dtype=float)
    scatter = np.asarray(mcmc["scatter_quantiles"], dtype=float)
    alpha = quantiles[2, 0]
    beta = quantiles[2, 1]
    x_limits = [
        float(value)
        for value in _sequence(
            regression["x_limits"],
            name="regression.x_limits",
        )
    ]
    y_limits = [
        float(value)
        for value in _sequence(
            regression["y_limits"],
            name="regression.y_limits",
        )
    ]
    x_grid = np.logspace(
        np.log10(np.min(x)) - float(settings["grid_expansion_dex"]),
        np.log10(np.max(x)) + float(settings["grid_expansion_dex"]),
        300,
    )

    plt.rcParams.update(
        {
            "font.size": float(settings["font_size"]),
            "axes.labelsize": float(settings["label_font_size"]),
            "legend.fontsize": float(settings["legend_font_size"]),
            "xtick.labelsize": float(settings["tick_font_size"]),
            "ytick.labelsize": float(settings["tick_font_size"]),
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    figure, axis = plt.subplots(
        figsize=(
            float(settings["width_inches"]),
            float(settings["height_inches"]),
        )
    )
    axis.errorbar(
        x,
        y,
        xerr=x_error,
        yerr=y_error,
        fmt="o",
        markersize=float(settings["marker_size"]),
        markeredgewidth=float(settings["marker_edge_width"]),
        linewidth=float(settings["error_line_width"]),
        alpha=float(settings["point_alpha"]),
        capsize=float(settings["capsize"]),
        markeredgecolor="black",
        markerfacecolor="white",
        ecolor=str(settings["error_color"]),
        elinewidth=float(settings["error_line_width"]),
        zorder=2,
    )

    random = np.random.default_rng(int(settings["posterior_seed"]))
    line_count = min(int(settings["posterior_lines"]), len(samples))
    indices = random.choice(len(samples), line_count, replace=False)
    fit_color = str(settings["fit_color"])
    for index in indices:
        axis.plot(
            x_grid,
            10.0 ** samples[index, 1] * x_grid ** samples[index, 0],
            color=fit_color,
            linewidth=float(settings["posterior_line_width"]),
            alpha=float(settings["posterior_alpha"]),
            zorder=3,
        )
    axis.plot(
        x_grid,
        10.0**beta * x_grid**alpha,
        color=fit_color,
        linewidth=float(settings["best_fit_line_width"]),
        label="Best power-law fit",
        zorder=4,
    )

    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlim(*x_limits)
    axis.set_ylim(*y_limits)
    axis.set_xlabel(
        r"$\Sigma_{\rm env}^{\rm ic}$ (g cm$^{-2}$)"
    )
    axis.set_ylabel(
        r"$\Sigma_{\rm core}^{\rm ic}$ (g cm$^{-2}$)"
    )
    annotation = (
        f"N = {data['n_valid']}\n"
        f"Pearson r = {comparison['pearson'].statistic:.2f}\n"
        f"Spearman $\\rho$ = {comparison['spearman'].statistic:.2f}\n"
        f"$\\alpha$ = {alpha:.3f} "
        f"(+{quantiles[3, 0] - alpha:.3f}/"
        f"-{alpha - quantiles[1, 0]:.3f})\n"
        f"$\\log_{{10}} A$ = {beta:.3f} "
        f"(+{quantiles[3, 1] - beta:.3f}/"
        f"-{beta - quantiles[1, 1]:.3f})\n"
        f"$\\sigma_{{\\rm int}}$ = {scatter[2]:.3f} dex"
    )
    axis.text(
        float(settings["annotation_x"]),
        float(settings["annotation_y"]),
        annotation,
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=float(settings["annotation_font_size"]),
        bbox={
            "facecolor": "white",
            "edgecolor": "0.7",
            "alpha": 0.85,
        },
    )
    axis.legend(
        frameon=False,
        loc=str(settings["legend_location"]),
    )
    axis.tick_params(which="both", direction="in", top=True, right=True)
    return figure


def build_regression(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> dict[str, object]:
    """Run the configured regression, write results, and create its figure."""
    catalogue = _load_catalogue(paths["catalogue"])
    data = _regression_data(catalogue, config)
    regression = _mapping(config, "regression")
    comparison = compare_powerlaw_methods(
        data,
        sensitivity_fractions=[
            float(value)
            for value in _sequence(
                regression["sensitivity_fractions"],
                name="regression.sensitivity_fractions",
            )
        ],
    )
    mcmc = run_powerlaw_mcmc(
        data,
        n_walkers=int(regression["n_walkers"]),
        n_steps=int(regression["n_steps"]),
        burn_in=int(regression["burn_in"]),
        thin=int(regression["thin"]),
        seed=int(regression["seed"]),
    )
    summary = _regression_summary(
        paths,
        config,
        data,
        comparison,
        mcmc,
    )
    summary_output = regression_summary_path(paths)
    samples_output = posterior_path(paths)
    figure_outputs = _figure_paths(
        paths,
        config,
        "intercore_surface_density_regression",
    )
    publication_outputs = _figure_paths(
        paths,
        config,
        "correlation_surface_env-surface_core",
    )
    _prepare_outputs(
        [
            summary_output,
            samples_output,
            *figure_outputs,
            *publication_outputs,
        ],
        overwrite=overwrite,
    )
    with summary_output.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    np.savez_compressed(
        samples_output,
        samples=np.asarray(mcmc["samples"], dtype=np.float32),
        parameter_names=np.asarray(["alpha", "beta", "lnsigma"]),
    )
    figure = create_regression_figure(
        data,
        comparison,
        mcmc,
        config,
    )
    publication_figure = create_publication_figure(
        data,
        comparison,
        mcmc,
        config,
    )
    _save_figure(
        figure,
        figure_outputs,
        dpi=int(_mapping(config, "figures")["dpi"]),
    )
    _save_figure(
        publication_figure,
        publication_outputs,
        dpi=int(_mapping(config, "figures")["publication_dpi"]),
    )
    quantiles = np.asarray(mcmc["quantiles"], dtype=float)
    print(f"Wrote {summary_output}")
    print(f"Wrote {samples_output}")
    print(
        "MCMC alpha="
        f"{quantiles[2, 0]:.6f} "
        f"+{quantiles[3, 0] - quantiles[2, 0]:.6f}"
        f"/-{quantiles[2, 0] - quantiles[1, 0]:.6f}; "
        f"beta={quantiles[2, 1]:.6f}; "
        f"N={data['n_valid']}"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--stage",
        choices=("inspect", "matrices", "regression", "all"),
        default="inspect",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    project = ProjectPaths.discover()
    config = load_task_config(_project_path(project, args.config))
    paths = resolve_task_paths(project, config)
    if args.stage == "inspect":
        inspect_inputs(paths, config)
    if args.stage in {"matrices", "all"}:
        build_matrices(paths, config, overwrite=args.overwrite)
    if args.stage in {"regression", "all"}:
        build_regression(paths, config, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
