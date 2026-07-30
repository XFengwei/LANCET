#!/usr/bin/env python
"""Estimate global and subregion core mass functions."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from pathlib import Path
import tomllib

from astropy import units as u
from astropy.table import Table, vstack
import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt

from lancet import ProjectPaths
from lancet.core_mass_function import (
    bootstrap_powerlaw_alpha,
    diagnose_minimum_mass,
    differential_mass_function,
    powerlaw_pdf,
    sample_alpha_posterior,
    valid_mass_mask,
)


DEFAULT_CONFIG = Path("tasks/core_mass_function/task.toml")


def load_task_config(path: Path) -> dict[str, object]:
    """Load task settings from TOML."""
    if not path.is_file():
        raise FileNotFoundError(f"Task configuration not found: {path}")
    with path.open("rb") as handle:
        return tomllib.load(handle)


def _mapping(
    config: Mapping[str, object],
    name: str,
) -> Mapping[str, object]:
    value = config[name]
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a TOML table")
    return value


def _project_path(project: ProjectPaths, value: object) -> Path:
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = project.root / path
    return path.resolve()


def resolve_task_paths(
    project: ProjectPaths,
    config: Mapping[str, object],
) -> dict[str, Path]:
    """Resolve the generated catalogue and policy-safe output root."""
    inputs = _mapping(config, "inputs")
    outputs = _mapping(config, "outputs")
    paths = {
        name: _project_path(project, value)
        for name, value in inputs.items()
    }
    paths["output_root"] = project.validate_output(str(outputs["root"]))
    for name in inputs:
        if not paths[name].is_file():
            raise FileNotFoundError(f"{name} not found: {paths[name]}")
    return paths


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


def _figure_paths(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    stem: str,
) -> list[Path]:
    figures = _mapping(config, "figures")
    formats = [str(value).lower() for value in figures["formats"]]
    unsupported = set(formats) - {"pdf", "png"}
    if unsupported:
        raise ValueError(
            "Unsupported figure formats: " + ", ".join(sorted(unsupported))
        )
    return [
        paths["output_root"] / "figures" / f"{stem}.{extension}"
        for extension in formats
    ]


def summary_path(paths: Mapping[str, Path]) -> Path:
    return paths["output_root"] / "tables" / "powerlaw_fits.ecsv"


def diagnostics_path(paths: Mapping[str, Path]) -> Path:
    return paths["output_root"] / "tables" / "mmin_diagnostics.ecsv"


def catalogue_export_path(paths: Mapping[str, Path]) -> Path:
    return paths["output_root"] / "tables" / "lancet_cmf_catalogue.csv"


def posterior_path(paths: Mapping[str, Path]) -> Path:
    return paths["output_root"] / "posterior" / "alpha_samples.npz"


def _column_settings(
    config: Mapping[str, object],
) -> Mapping[str, object]:
    return _mapping(config, "columns")


def prepare_catalogue(
    catalogue: Table,
    config: Mapping[str, object],
) -> Table:
    """Select aligned finite-positive mass rows and CMF export columns."""
    columns = _column_settings(config)
    names = [
        str(columns[name])
        for name in (
            "id",
            "mass",
            "mass_error",
            "temperature",
            "size",
            "region",
        )
    ]
    missing = set(names) - set(catalogue.colnames)
    if missing:
        raise KeyError(
            "Physical catalogue is missing columns: "
            + ", ".join(sorted(missing))
        )
    mass = np.asarray(catalogue[str(columns["mass"])], dtype=float)
    selected = catalogue[valid_mass_mask(mass)][names].copy()
    selected.meta["selection"] = (
        f"finite positive {columns['mass']}; row alignment preserved"
    )
    selected.meta["cmf_convention"] = "dN/dM proportional to M^(-alpha)"
    return selected


def _region_definitions(
    config: Mapping[str, object],
) -> list[tuple[str, Mapping[str, object]]]:
    regions = _mapping(config, "regions")
    ordering = _mapping(config, "ordering")
    order = [str(value) for value in ordering["region_order"]]
    definitions: list[tuple[str, Mapping[str, object]]] = []
    for key in order:
        value = regions[key]
        if not isinstance(value, Mapping):
            raise TypeError(f"regions.{key} must be a TOML table")
        definitions.append((key, value))
    return definitions


def _fit_one(
    masses: np.ndarray,
    *,
    key: str,
    label: str,
    code: str,
    minimum_count: int,
    proposal_width: float,
    seed: int,
    config: Mapping[str, object],
) -> dict[str, object]:
    fit = _mapping(config, "fit")
    completeness = float(fit["completeness_mass_msun"])
    diagnosis = diagnose_minimum_mass(
        masses,
        minimum_count=minimum_count,
        completeness_mass=completeness if completeness > 0 else None,
    )
    minimum = float(diagnosis["minimum_mass"])
    posterior = sample_alpha_posterior(
        masses,
        minimum,
        steps=int(fit["mcmc_steps"]),
        burn_in=int(fit["mcmc_burn_in"]),
        proposal_width=proposal_width,
        seed=seed,
    )
    bootstrap = bootstrap_powerlaw_alpha(
        masses,
        minimum,
        samples=int(fit["bootstrap_samples"]),
        seed=seed + int(fit["bootstrap_seed_offset"]),
    )
    return {
        "key": key,
        "label": label,
        "code": code,
        "masses": masses,
        "diagnosis": diagnosis,
        "posterior": posterior,
        "bootstrap": bootstrap,
    }


def analyse_catalogue(
    catalogue: Table,
    config: Mapping[str, object],
) -> dict[str, dict[str, object]]:
    """Fit the global and configured regional CMFs."""
    columns = _column_settings(config)
    fit = _mapping(config, "fit")
    masses = np.asarray(catalogue[str(columns["mass"])], dtype=float)
    region_values = np.asarray(
        catalogue[str(columns["region"])]
    ).astype(str)
    analyses: dict[str, dict[str, object]] = {
        "all": _fit_one(
            masses,
            key="all",
            label="All cores",
            code="all",
            minimum_count=int(fit["global_minimum_count"]),
            proposal_width=float(fit["global_proposal_width"]),
            seed=int(fit["global_seed"]),
            config=config,
        )
    }
    for key, region in _region_definitions(config):
        code = str(region["code"])
        analyses[key] = _fit_one(
            masses[region_values == code],
            key=key,
            label=str(region["label"]),
            code=code,
            minimum_count=int(fit["regional_minimum_count"]),
            proposal_width=float(fit["regional_proposal_width"]),
            seed=int(region["seed"]),
            config=config,
        )
    return analyses


def _summary_table(
    analyses: Mapping[str, Mapping[str, object]],
) -> Table:
    rows: list[dict[str, object]] = []
    for analysis in analyses.values():
        diagnosis = analysis["diagnosis"]
        posterior = analysis["posterior"]
        bootstrap = np.asarray(analysis["bootstrap"], dtype=float)
        masses = np.asarray(analysis["masses"], dtype=float)
        assert isinstance(diagnosis, Mapping)
        assert isinstance(posterior, Mapping)
        rows.append(
            {
                "sample": analysis["key"],
                "label": analysis["label"],
                "subregion_code": analysis["code"],
                "total_count": masses.size,
                "minimum_mass": diagnosis["minimum_mass"],
                "tail_count": diagnosis["count"],
                "ks_distance": diagnosis["ks_distance"],
                "alpha_mle": diagnosis["alpha"],
                "alpha_posterior_median": posterior["median"],
                "alpha_posterior_p16": posterior["p16"],
                "alpha_posterior_p84": posterior["p84"],
                "alpha_bootstrap_median": np.median(bootstrap),
                "alpha_bootstrap_p16": np.percentile(bootstrap, 16),
                "alpha_bootstrap_p84": np.percentile(bootstrap, 84),
                "mcmc_acceptance_fraction": posterior[
                    "acceptance_fraction"
                ],
            }
        )
    table = Table(rows=rows)
    table["minimum_mass"].unit = u.M_sun
    table.meta["cmf_convention"] = "dN/dM proportional to M^(-alpha)"
    table.meta["minimum_mass_selection"] = (
        "minimum KS distance among observed masses satisfying minimum_count"
    )
    table.meta["posterior_prior"] = "flat for alpha > 1"
    table.meta["bootstrap_scope"] = "fixed selected minimum mass"
    return table


def _diagnostics_table(
    analyses: Mapping[str, Mapping[str, object]],
) -> Table:
    tables: list[Table] = []
    for analysis in analyses.values():
        diagnosis = analysis["diagnosis"]
        assert isinstance(diagnosis, Mapping)
        count = len(np.asarray(diagnosis["minimum_masses"]))
        tables.append(
            Table(
                {
                    "sample": np.full(count, analysis["key"]),
                    "label": np.full(count, analysis["label"]),
                    "minimum_mass": diagnosis["minimum_masses"],
                    "tail_count": diagnosis["counts"],
                    "alpha_mle": diagnosis["alphas"],
                    "ks_distance": diagnosis["ks_distances"],
                }
            )
        )
    output = vstack(tables, metadata_conflicts="silent")
    output["minimum_mass"].unit = u.M_sun
    return output


def write_analysis(
    source_catalogue: Table,
    selected_catalogue: Table,
    analyses: Mapping[str, Mapping[str, object]],
    paths: Mapping[str, Path],
    *,
    overwrite: bool,
) -> None:
    """Write fit summaries, diagnostics, samples, and safe catalogue export."""
    del source_catalogue
    outputs = [
        summary_path(paths),
        diagnostics_path(paths),
        catalogue_export_path(paths),
        posterior_path(paths),
    ]
    _prepare_outputs(outputs, overwrite=overwrite)
    _summary_table(analyses).write(
        outputs[0],
        format="ascii.ecsv",
        overwrite=True,
    )
    _diagnostics_table(analyses).write(
        outputs[1],
        format="ascii.ecsv",
        overwrite=True,
    )
    selected_catalogue.write(
        outputs[2],
        format="ascii.csv",
        overwrite=True,
    )
    sample_arrays = {
        str(analysis["key"]): np.asarray(
            analysis["posterior"]["posterior"],
            dtype=float,
        )
        for analysis in analyses.values()
    }
    np.savez_compressed(outputs[3], **sample_arrays)
    for output in outputs:
        print(f"Wrote {output}")


def inspect_inputs(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
) -> None:
    """Report catalogue selection, fitting convention, and dependencies."""
    source = Table.read(paths["catalogue"], format="ascii.ecsv")
    selected = prepare_catalogue(source, config)
    columns = _column_settings(config)
    masses = np.asarray(selected[str(columns["mass"])], dtype=float)
    regions = np.asarray(selected[str(columns["region"])]).astype(str)
    fit = _mapping(config, "fit")
    diagnosis = diagnose_minimum_mass(
        masses,
        minimum_count=int(fit["global_minimum_count"]),
        completeness_mass=(
            float(fit["completeness_mass_msun"])
            if float(fit["completeness_mass_msun"]) > 0
            else None
        ),
    )
    print("Core-mass-function input:")
    print(
        f"  catalogue: {paths['catalogue']} "
        f"(rows={len(source)}, valid positive masses={len(selected)})"
    )
    print(
        f"  mass range: {masses.min():.6g}--{masses.max():.6g} solMass"
    )
    for _, region in _region_definitions(config):
        code = str(region["code"])
        print(
            f"  {region['label']} ({code}): "
            f"{np.count_nonzero(regions == code)} cores"
        )
    print("Convention: dN/dM proportional to M^(-alpha)")
    print(
        "Global KS diagnosis: "
        f"Mmin={diagnosis['minimum_mass']:.6g} solMass, "
        f"Ntail={diagnosis['count']}, "
        f"alpha_MLE={diagnosis['alpha']:.6f}, "
        f"KS={diagnosis['ks_distance']:.6f}"
    )
    if float(fit["completeness_mass_msun"]) <= 0:
        print("Completeness floor: none (Mmin selected statistically)")
    else:
        print(
            "Completeness floor: "
            f"{fit['completeness_mass_msun']} solMass"
        )
    print(f"Output root: {paths['output_root']}")


def _shared_bin_edges(
    selected: Table,
    config: Mapping[str, object],
) -> np.ndarray:
    columns = _column_settings(config)
    figures = _mapping(config, "figures")
    masses = np.asarray(selected[str(columns["mass"])], dtype=float)
    return np.logspace(
        np.log10(masses.min()),
        np.log10(masses.max()),
        int(figures["mass_bins"]) + 1,
    )


def _plot_observed(
    axis: plt.Axes,
    masses: np.ndarray,
    edges: np.ndarray,
    *,
    color: str,
    marker: str,
    label: str,
    alpha: float = 1.0,
    zorder: int = 3,
) -> None:
    cmf = differential_mass_function(masses, edges)
    mask = np.asarray(cmf["counts"]) > 0
    axis.errorbar(
        np.asarray(cmf["centers"])[mask],
        np.asarray(cmf["density"])[mask],
        yerr=np.asarray(cmf["density_error"])[mask],
        fmt=marker,
        markersize=6,
        markerfacecolor=color,
        markeredgecolor="white" if color != "white" else "0.2",
        markeredgewidth=0.8,
        color=color if color != "white" else "0.25",
        ecolor=color if color != "white" else "0.35",
        elinewidth=1.0,
        capsize=2.5,
        alpha=alpha,
        label=label,
        zorder=zorder,
    )


def _model_predictions(
    analysis: Mapping[str, object],
    x_grid: np.ndarray,
    *,
    draws: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    diagnosis = analysis["diagnosis"]
    posterior = analysis["posterior"]
    assert isinstance(diagnosis, Mapping)
    assert isinstance(posterior, Mapping)
    minimum = float(diagnosis["minimum_mass"])
    count = int(diagnosis["count"])
    samples = np.asarray(posterior["posterior"], dtype=float)
    rng = np.random.default_rng(seed)
    indices = rng.choice(
        samples.size,
        size=min(draws, samples.size),
        replace=False,
    )
    prediction = np.asarray(
        [
            count * powerlaw_pdf(x_grid, alpha, minimum)
            for alpha in samples[indices]
        ]
    )
    median_alpha = float(posterior["median"])
    median = count * powerlaw_pdf(x_grid, median_alpha, minimum)
    return (
        median,
        np.percentile(prediction, 16, axis=0),
        np.percentile(prediction, 84, axis=0),
        np.percentile(prediction, 2.5, axis=0),
        np.percentile(prediction, 97.5, axis=0),
    )


def _style_mass_axis(axis: plt.Axes) -> None:
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel(r"Core mass, $M_{\rm core}$ ($M_\odot$)")
    axis.set_ylabel(r"Differential CMF, $\mathrm{d}N/\mathrm{d}M$")
    axis.grid(alpha=0.15, which="both")


def plot_mmin_diagnosis(
    analyses: Mapping[str, Mapping[str, object]],
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool,
) -> list[Path]:
    outputs = _figure_paths(
        paths,
        config,
        "CMF_Mmin_KS_diagnosis",
    )
    _prepare_outputs(outputs, overwrite=overwrite)
    figures = _mapping(config, "figures")
    figure, axis = plt.subplots(
        figsize=(float(figures["width_inches"]), 4.8)
    )
    colors = {"all": "black"}
    for key, region in _region_definitions(config):
        colors[key] = str(region["color"])
    for key, analysis in analyses.items():
        diagnosis = analysis["diagnosis"]
        assert isinstance(diagnosis, Mapping)
        axis.plot(
            diagnosis["minimum_masses"],
            diagnosis["ks_distances"],
            color=colors[key],
            linewidth=1.5,
            label=str(analysis["label"]),
        )
        axis.scatter(
            [diagnosis["minimum_mass"]],
            [diagnosis["ks_distance"]],
            s=35,
            facecolor=colors[key],
            edgecolor="white",
            linewidth=0.6,
            zorder=3,
        )
    axis.set_xscale("log")
    axis.set_xlabel(r"Candidate $M_{\min}$ ($M_\odot$)")
    axis.set_ylabel("KS distance")
    axis.grid(alpha=0.18, which="both")
    axis.legend(frameon=False)
    figure.tight_layout()
    _save_figure(figure, outputs, figures)
    return outputs


def plot_global_fit(
    selected: Table,
    analyses: Mapping[str, Mapping[str, object]],
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool,
) -> list[Path]:
    outputs = _figure_paths(
        paths,
        config,
        "CMF_global_MLE_MCMC_fit",
    )
    _prepare_outputs(outputs, overwrite=overwrite)
    columns = _column_settings(config)
    figures = _mapping(config, "figures")
    analysis = analyses["all"]
    masses = np.asarray(selected[str(columns["mass"])], dtype=float)
    edges = _shared_bin_edges(selected, config)
    diagnosis = analysis["diagnosis"]
    posterior = analysis["posterior"]
    assert isinstance(diagnosis, Mapping)
    assert isinstance(posterior, Mapping)
    minimum = float(diagnosis["minimum_mass"])
    x_grid = np.logspace(
        np.log10(minimum),
        np.log10(masses.max() * 1.15),
        300,
    )
    median, p16, p84, p025, p975 = _model_predictions(
        analysis,
        x_grid,
        draws=int(figures["posterior_draws"]),
        seed=int(figures["plot_seed"]),
    )

    figure, axis = plt.subplots(
        figsize=(
            float(figures["width_inches"]),
            float(figures["height_inches"]),
        )
    )
    _plot_observed(
        axis,
        masses,
        edges,
        color="white",
        marker="o",
        label="Observed CMF",
    )
    fit_color = str(figures["global_fit_color"])
    axis.fill_between(
        x_grid,
        p025,
        p975,
        color=fit_color,
        alpha=0.08,
        linewidth=0,
        label="95% posterior",
    )
    axis.fill_between(
        x_grid,
        p16,
        p84,
        color=fit_color,
        alpha=0.16,
        linewidth=0,
        label="68% posterior",
    )
    upper = float(posterior["p84"]) - float(posterior["median"])
    lower = float(posterior["median"]) - float(posterior["p16"])
    axis.plot(
        x_grid,
        median,
        color=fit_color,
        linewidth=2.4,
        label=(
            rf"Fit: $\alpha={float(posterior['median']):.2f}"
            rf"^{{+{upper:.2f}}}_{{-{lower:.2f}}}$"
        ),
    )
    axis.axvline(
        minimum,
        color="0.45",
        linestyle="--",
        linewidth=1.2,
    )
    axis.text(
        minimum * 1.05,
        0.05,
        rf"$M_{{\min}}={minimum:.2f}\,M_\odot$",
        transform=axis.get_xaxis_transform(),
        rotation=90,
        color="0.4",
        va="bottom",
    )
    axis.text(
        0.97,
        0.96,
        (
            rf"$N_{{\rm core}}={masses.size}$"
            "\n"
            rf"$N(M\geq M_{{\min}})={int(diagnosis['count'])}$"
        ),
        transform=axis.transAxes,
        ha="right",
        va="top",
    )
    _style_mass_axis(axis)
    axis.legend(frameon=False, loc="lower left")
    figure.tight_layout()
    _save_figure(figure, outputs, figures)
    return outputs


def plot_global_posterior(
    analyses: Mapping[str, Mapping[str, object]],
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool,
) -> list[Path]:
    outputs = _figure_paths(
        paths,
        config,
        "CMF_alpha_posterior",
    )
    _prepare_outputs(outputs, overwrite=overwrite)
    figures = _mapping(config, "figures")
    analysis = analyses["all"]
    posterior = analysis["posterior"]
    diagnosis = analysis["diagnosis"]
    assert isinstance(posterior, Mapping)
    assert isinstance(diagnosis, Mapping)
    samples = np.asarray(posterior["posterior"], dtype=float)
    figure, axis = plt.subplots(figsize=(6.0, 4.4))
    axis.hist(
        samples,
        bins=int(figures["posterior_bins"]),
        density=True,
        color=str(figures["global_fit_color"]),
        alpha=0.65,
    )
    for value, style, label in (
        (posterior["median"], "-", "posterior median"),
        (diagnosis["alpha"], "--", "MLE"),
    ):
        axis.axvline(
            float(value),
            color="black",
            linestyle=style,
            linewidth=1.3,
            label=label,
        )
    axis.set_xlabel(r"Power-law exponent $\alpha$")
    axis.set_ylabel("Posterior density")
    axis.legend(frameon=False)
    axis.grid(alpha=0.15)
    figure.tight_layout()
    _save_figure(figure, outputs, figures)
    return outputs


def plot_region_overlay(
    selected: Table,
    analyses: Mapping[str, Mapping[str, object]],
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    with_fits: bool,
    overwrite: bool,
) -> list[Path]:
    stem = (
        "CMF_subregions_MLE_MCMC_fit"
        if with_fits
        else "CMF_global_subregion_overlay"
    )
    outputs = _figure_paths(paths, config, stem)
    _prepare_outputs(outputs, overwrite=overwrite)
    columns = _column_settings(config)
    figures = _mapping(config, "figures")
    mass_column = str(columns["mass"])
    region_column = str(columns["region"])
    masses = np.asarray(selected[mass_column], dtype=float)
    region_values = np.asarray(selected[region_column]).astype(str)
    edges = _shared_bin_edges(selected, config)
    figure, axis = plt.subplots(
        figsize=(
            float(figures["width_inches"]),
            float(figures["height_inches"]),
        )
    )
    _plot_observed(
        axis,
        masses,
        edges,
        color="white",
        marker="o",
        label=f"All cores ($N={masses.size}$)",
        alpha=0.55,
        zorder=1,
    )
    for index, (key, region) in enumerate(
        _region_definitions(config)
    ):
        code = str(region["code"])
        region_masses = masses[region_values == code]
        color = str(region["color"])
        _plot_observed(
            axis,
            region_masses,
            edges,
            color=color,
            marker=str(region["marker"]),
            label=(
                f"{region['label']} ($N={region_masses.size}$)"
                if not with_fits
                else "_nolegend_"
            ),
        )
        if not with_fits:
            continue
        analysis = analyses[key]
        diagnosis = analysis["diagnosis"]
        posterior = analysis["posterior"]
        assert isinstance(diagnosis, Mapping)
        assert isinstance(posterior, Mapping)
        minimum = float(diagnosis["minimum_mass"])
        x_grid = np.logspace(
            np.log10(minimum),
            np.log10(region_masses.max() * 1.15),
            250,
        )
        median, p16, p84, _, _ = _model_predictions(
            analysis,
            x_grid,
            draws=int(figures["posterior_draws"]),
            seed=int(figures["plot_seed"]) + index + 1,
        )
        axis.fill_between(
            x_grid,
            p16,
            p84,
            color=color,
            alpha=0.12,
            linewidth=0,
        )
        upper = float(posterior["p84"]) - float(posterior["median"])
        lower = float(posterior["median"]) - float(posterior["p16"])
        axis.plot(
            x_grid,
            median,
            color=color,
            linewidth=2.2,
            label=(
                rf"{region['label']} fit: "
                rf"$\alpha={float(posterior['median']):.2f}"
                rf"^{{+{upper:.2f}}}_{{-{lower:.2f}}}$"
            ),
        )
        axis.axvline(
            minimum,
            color=color,
            linestyle=":",
            linewidth=1.2,
        )
    _style_mass_axis(axis)
    axis.legend(frameon=False, fontsize=9)
    figure.tight_layout()
    _save_figure(figure, outputs, figures)
    return outputs


def plot_region_posteriors(
    analyses: Mapping[str, Mapping[str, object]],
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool,
) -> list[Path]:
    outputs = _figure_paths(
        paths,
        config,
        "CMF_subregions_alpha_posterior",
    )
    _prepare_outputs(outputs, overwrite=overwrite)
    figures = _mapping(config, "figures")
    figure, axis = plt.subplots(figsize=(6.2, 4.6))
    for key, region in _region_definitions(config):
        posterior = analyses[key]["posterior"]
        assert isinstance(posterior, Mapping)
        axis.hist(
            posterior["posterior"],
            bins=int(figures["posterior_bins"]),
            density=True,
            histtype="step",
            linewidth=1.6,
            color=str(region["color"]),
            label=str(region["label"]),
        )
        axis.axvline(
            float(posterior["median"]),
            color=str(region["color"]),
            linestyle="--",
            linewidth=1.0,
        )
    axis.set_xlabel(r"Power-law exponent $\alpha$")
    axis.set_ylabel("Posterior density")
    axis.legend(frameon=False)
    axis.grid(alpha=0.15)
    figure.tight_layout()
    _save_figure(figure, outputs, figures)
    return outputs


def _save_figure(
    figure: plt.Figure,
    outputs: Sequence[Path],
    settings: Mapping[str, object],
) -> None:
    for output in outputs:
        figure.savefig(
            output,
            dpi=int(settings["dpi"]),
            bbox_inches="tight",
        )
        print(f"Wrote {output}")
    plt.close(figure)


def plot_all(
    selected: Table,
    analyses: Mapping[str, Mapping[str, object]],
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool,
) -> None:
    """Write the six configured CMF diagnostic/publication figure groups."""
    plot_mmin_diagnosis(
        analyses,
        paths,
        config,
        overwrite=overwrite,
    )
    plot_global_fit(
        selected,
        analyses,
        paths,
        config,
        overwrite=overwrite,
    )
    plot_global_posterior(
        analyses,
        paths,
        config,
        overwrite=overwrite,
    )
    plot_region_overlay(
        selected,
        analyses,
        paths,
        config,
        with_fits=False,
        overwrite=overwrite,
    )
    plot_region_overlay(
        selected,
        analyses,
        paths,
        config,
        with_fits=True,
        overwrite=overwrite,
    )
    plot_region_posteriors(
        analyses,
        paths,
        config,
        overwrite=overwrite,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--stage",
        choices=("inspect", "fit", "plot", "all"),
        default="inspect",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    project = ProjectPaths.discover()
    config = load_task_config(_project_path(project, args.config))
    paths = resolve_task_paths(project, config)
    if args.stage == "inspect":
        inspect_inputs(paths, config)
        return

    source = Table.read(paths["catalogue"], format="ascii.ecsv")
    selected = prepare_catalogue(source, config)
    analyses = analyse_catalogue(selected, config)
    if args.stage in {"fit", "all"}:
        write_analysis(
            source,
            selected,
            analyses,
            paths,
            overwrite=args.overwrite,
        )
    if args.stage in {"plot", "all"}:
        plot_all(
            selected,
            analyses,
            paths,
            config,
            overwrite=args.overwrite,
        )


if __name__ == "__main__":
    main()
