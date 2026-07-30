#!/usr/bin/env python
"""Measure and display multiscale core-formation efficiencies."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from pathlib import Path
import tomllib

from astropy.io import fits
from astropy.table import Table
import matplotlib.patheffects as path_effects
import numpy as np

from lancet import ProjectPaths
from lancet.core_formation import (
    create_core_overlay_figure,
    summarize_core_formation,
)


DEFAULT_CONFIG = Path("tasks/core_formation/task.toml")


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
    """Resolve required inputs, an optional reference, and output root."""
    inputs = config["inputs"]
    references = config["references"]
    outputs = config["outputs"]
    if not all(
        isinstance(value, Mapping)
        for value in (inputs, references, outputs)
    ):
        raise TypeError(
            "Inputs, references, and outputs must be TOML tables"
        )

    paths = {
        name: _project_path(project, value)
        for name, value in inputs.items()
    }
    paths.update(
        {
            name: _project_path(project, value)
            for name, value in references.items()
        }
    )
    paths["output_root"] = project.validate_output(str(outputs["root"]))

    for name in inputs:
        if not paths[name].is_file():
            raise FileNotFoundError(f"{name} not found: {paths[name]}")
    reference = paths["paper1_fragment_catalogue"]
    if reference.is_file() and not reference.is_relative_to(
        project.delivered.resolve()
    ):
        raise ValueError(
            "paper1_fragment_catalogue must be a delivered reference"
        )
    return paths


def _region_settings(
    config: Mapping[str, object],
) -> dict[str, Mapping[str, object]]:
    regions = config["regions"]
    if not isinstance(regions, Mapping):
        raise TypeError("Regions must be a TOML table")
    settings: dict[str, Mapping[str, object]] = {}
    for key, value in regions.items():
        if not isinstance(value, Mapping):
            raise TypeError(f"Region {key} must be a TOML table")
        settings[str(key)] = value
    return settings


def _opacity_settings(
    config: Mapping[str, object],
) -> tuple[float, float]:
    opacity = config["dust_opacity"]
    if not isinstance(opacity, Mapping):
        raise TypeError("Dust opacity must be a TOML table")
    return (
        float(opacity["paper1_cm2_per_g"]),
        float(opacity["adopted_cm2_per_g"]),
    )


def _summarize(
    catalogue: Table,
    config: Mapping[str, object],
) -> Table:
    reference_opacity, adopted_opacity = _opacity_settings(config)
    return summarize_core_formation(
        catalogue,
        _region_settings(config),
        reference_opacity_cm2_per_g=reference_opacity,
        adopted_opacity_cm2_per_g=adopted_opacity,
    )


def _fits_description(path: Path) -> str:
    with fits.open(path, memmap=True) as hdul:
        header = hdul[0].header
        return (
            f"shape={tuple(hdul[0].shape)}, "
            f"BITPIX={int(header['BITPIX'])}, "
            f"BUNIT={header.get('BUNIT', 'not recorded')}"
        )


def inspect_inputs(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
) -> None:
    """Report dependencies, adopted values, and output location."""
    physical = Table.read(
        paths["physical_catalogue"],
        format="ascii.ecsv",
    )
    sources = Table.read(
        paths["source_catalogue"],
        format="ascii.ecsv",
    )
    fragments = Table.read(
        paths["fragment_ellipse_catalogue"],
        format="ascii",
    )
    print("Core-formation inputs:")
    print(
        f"  physical_catalogue: {paths['physical_catalogue']} "
        f"(rows={len(physical)})"
    )
    print(
        f"  source_catalogue: {paths['source_catalogue']} "
        f"(rows={len(sources)})"
    )
    print(
        "  fragment_ellipse_catalogue: "
        f"{paths['fragment_ellipse_catalogue']} "
        f"(rows={len(fragments)})"
    )
    for name in ("high_resolution_continuum", "low_resolution_continuum"):
        print(f"  {name}: {paths[name]}")
        print(f"    {_fits_description(paths[name])}")

    reference = paths["paper1_fragment_catalogue"]
    if reference.is_file():
        paper1 = Table.read(reference)
        print(
            "Paper I reference catalogue: "
            f"{reference} (rows={len(paper1)}, not used in calculations)"
        )
    else:
        print(
            "Paper I reference catalogue: "
            f"{reference} (optional reference is unavailable)"
        )

    reference_opacity, adopted_opacity = _opacity_settings(config)
    print(
        "Fragment opacity conversion: "
        f"{reference_opacity:g} -> {adopted_opacity:g} cm2/g; "
        f"mass scale={reference_opacity / adopted_opacity:.3f}"
    )
    summary = _summarize(physical, config)
    for row in summary:
        print(
            f"  {row['label']}: cores={row['core_count']}, "
            f"corrected fragment mass={row['fragment_mass']:.3f} solMass, "
            f"fragment efficiency="
            f"{100 * row['fragment_efficiency']:.3f}%, "
            f"core mass={row['core_mass']:.3f} solMass, "
            f"core efficiency={100 * row['core_efficiency']:.3f}%, "
            f"fragment-to-core efficiency="
            f"{100 * row['fragment_to_core_efficiency']:.3f}%"
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
        output.parent.mkdir(parents=True, exist_ok=True)


def summary_path(paths: Mapping[str, Path]) -> Path:
    """Return the machine-readable summary path."""
    return (
        paths["output_root"]
        / "tables"
        / "core_formation_summary.ecsv"
    )


def build_summary(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> Table:
    """Calculate and write regional masses and efficiencies."""
    catalogue = Table.read(
        paths["physical_catalogue"],
        format="ascii.ecsv",
    )
    summary = _summarize(catalogue, config)
    output = summary_path(paths)
    _prepare_outputs([output], overwrite=overwrite)
    summary.write(output, format="ascii.ecsv", overwrite=True)
    print(f"Wrote {output}")
    for row in summary:
        print(
            f"{row['label']}: M_core={row['core_mass']:.3f} solMass, "
            f"epsilon_fragment={100 * row['fragment_efficiency']:.3f}%, "
            f"epsilon_core={100 * row['core_efficiency']:.3f}%, "
            f"epsilon_fragment_to_core="
            f"{100 * row['fragment_to_core_efficiency']:.3f}%"
        )
    return summary


def _figure_outputs(
    paths: Mapping[str, Path],
    stem: str,
    config: Mapping[str, object],
) -> list[Path]:
    figures = config["figures"]
    if not isinstance(figures, Mapping):
        raise TypeError("Figures must be a TOML table")
    formats = [str(value).lower() for value in figures["formats"]]
    unsupported = set(formats) - {"pdf", "png"}
    if unsupported:
        raise ValueError(
            "Unsupported figure formats: "
            + ", ".join(sorted(unsupported))
        )
    root = paths["output_root"] / "figures"
    return [root / f"{stem}.{suffix}" for suffix in formats]


def _save_figure(
    figure,
    outputs: Sequence[Path],
    *,
    dpi: int,
) -> None:
    from matplotlib import pyplot as plt

    try:
        for output in outputs:
            figure.savefig(output, bbox_inches="tight", dpi=dpi)
            print(f"Wrote {output}")
    finally:
        plt.close(figure)


def _load_or_build_summary(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
) -> Table:
    path = summary_path(paths)
    if path.is_file():
        return Table.read(path, format="ascii.ecsv")
    catalogue = Table.read(
        paths["physical_catalogue"],
        format="ascii.ecsv",
    )
    return _summarize(catalogue, config)


def plot_efficiency(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> list[Path]:
    """Plot the adopted cloud-to-fragment-to-core efficiencies."""
    from matplotlib import pyplot as plt

    summary = _load_or_build_summary(paths, config)
    plot = config["efficiency_plot"]
    figures = config["figures"]
    if not isinstance(plot, Mapping) or not isinstance(figures, Mapping):
        raise TypeError("Plot and figure settings must be TOML tables")

    outputs = _figure_outputs(
        paths,
        "multiscale_mass_efficiency",
        config,
    )
    _prepare_outputs(outputs, overwrite=overwrite)
    scales = np.asarray(plot["scales_pc"], dtype=float)
    if scales.shape != (3,) or np.any(scales <= 0):
        raise ValueError("scales_pc must contain three positive values")

    style = {
        "Young": {"marker": "o", "color": "blue"},
        "Intermediate": {"marker": "s", "color": "orange"},
        "Evolved": {"marker": "^", "color": "red"},
    }
    with plt.rc_context(
        {
            "font.size": 9,
            "axes.linewidth": 1.0,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "xtick.minor.width": 0.8,
            "ytick.minor.width": 0.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    ):
        figure, axis = plt.subplots(figsize=(4.0, 3.4))
        for row in summary:
            label = str(row["label"])
            masses = np.asarray(
                [
                    row["cloud_mass"],
                    row["fragment_mass"],
                    row["core_mass"],
                ],
                dtype=float,
            )
            fractions = masses / masses[0]
            counts = np.asarray(
                [np.inf, row["fragment_count"], row["core_count"]],
                dtype=float,
            )
            errors = fractions * 0.5 / np.sqrt(counts - 1.0)
            region_style = style[label]
            axis.errorbar(
                scales,
                fractions,
                xerr=scales * float(plot["scale_fractional_error"]),
                yerr=errors,
                capsize=2,
                label=label,
                marker=region_style["marker"],
                color=region_style["color"],
                markeredgecolor="white",
                linewidth=1.8,
                markersize=7,
            )
            for scale, fraction in zip(scales[1:], fractions[1:]):
                horizontal_factor = (
                    2.8 if np.isclose(scale, scales[1]) else 2.0
                )
                axis.text(
                    scale * horizontal_factor,
                    fraction * 1.5,
                    f"{fraction * 100:.1f}%",
                    color=region_style["color"],
                    fontsize=12,
                    fontweight="heavy",
                    path_effects=[
                        path_effects.withStroke(
                            linewidth=0.8,
                            foreground="white",
                        )
                    ],
                )

        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.invert_xaxis()
        axis.set_xlabel("Physical scale (pc)")
        axis.set_ylabel("Mass efficiency")
        axis.set_xticks(scales)
        axis.set_xticklabels(
            ["Cloud\n4 pc", "Fragment\n0.15 pc", "Core\n0.01 pc"]
        )
        axis.set_yticks([0.001, 0.01, 0.1, 1.0])
        axis.set_yticklabels(["0.1%", "1%", "10%", "100%"])
        axis.scatter(
            scales[0],
            1.0,
            s=120,
            facecolors="none",
            edgecolors="black",
            linewidths=1.5,
            zorder=5,
        )
        axis.text(
            0.28,
            0.95,
            "Cloud-scale gas reservoir",
            transform=axis.transAxes,
            fontsize=9,
            horizontalalignment="left",
            verticalalignment="center",
            zorder=15,
            fontweight="bold",
            path_effects=[
                path_effects.withStroke(
                    linewidth=0.5,
                    foreground="white",
                )
            ],
        )
        axis.annotate(
            "",
            xy=(0.25, 0.95),
            xytext=(0.07, 0.95),
            xycoords="axes fraction",
            textcoords="axes fraction",
            arrowprops={
                "arrowstyle": "->",
                "linewidth": 1.5,
                "color": "black",
            },
            zorder=15,
        )
        axis.text(
            scales[1],
            0.001,
            r"$\epsilon_{\rm fragment}$",
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=12,
        )
        axis.text(
            scales[2],
            0.001,
            r"$\epsilon_{\rm core}$",
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=12,
        )
        axis.set_xlim(6.0, 0.005)
        axis.set_ylim(0.0006, 1.5)
        axis.legend(frameon=False, loc="lower left", handlelength=2.2)
        _save_figure(
            figure,
            outputs,
            dpi=int(figures["dpi"]),
        )
    return outputs


def plot_overlays(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> list[Path]:
    """Create the three core/fragment catalogue overlay figures."""
    sources = Table.read(
        paths["source_catalogue"],
        format="ascii.ecsv",
    )
    fragments = Table.read(
        paths["fragment_ellipse_catalogue"],
        format="ascii",
    )
    high_data, high_header = fits.getdata(
        paths["high_resolution_continuum"],
        header=True,
        memmap=True,
    )
    low_data, low_header = fits.getdata(
        paths["low_resolution_continuum"],
        header=True,
        memmap=True,
    )
    overlay = config["overlay"]
    panels = config["panels"]
    figures = config["figures"]
    if (
        not isinstance(overlay, Mapping)
        or not isinstance(panels, list)
        or not isinstance(figures, Mapping)
    ):
        raise TypeError(
            "Overlay, panels, and figures must have the expected TOML types"
        )

    all_outputs: list[Path] = []
    for panel in panels:
        if not isinstance(panel, Mapping):
            raise TypeError("Each panel must be a TOML table")
        key = str(panel["key"])
        outputs = _figure_outputs(
            paths,
            f"core_catalogue_lowres_contours_{key}",
            config,
        )
        _prepare_outputs(outputs, overwrite=overwrite)
        figure = create_core_overlay_figure(
            high_data,
            high_header,
            low_data,
            low_header,
            sources,
            fragments,
            center_ra_deg=float(panel["center_ra_deg"]),
            center_dec_deg=float(panel["center_dec_deg"]),
            width_arcsec=float(overlay["width_arcsec"]),
            height_arcsec=float(overlay["height_arcsec"]),
            label=str(panel["label"]),
            image_vmin=float(overlay["image_vmin"]),
            image_vmax=float(panel["image_vmax"]),
            contour_levels=tuple(
                float(value) for value in overlay["contour_levels_mjysr"]
            ),
        )
        _save_figure(
            figure,
            outputs,
            dpi=int(figures["dpi"]),
        )
        all_outputs.extend(outputs)
    return all_outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--stage",
        choices=("inspect", "summary", "efficiency", "overlays", "all"),
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
    if args.stage in {"summary", "all"}:
        build_summary(paths, config, overwrite=args.overwrite)
    if args.stage in {"efficiency", "all"}:
        plot_efficiency(paths, config, overwrite=args.overwrite)
    if args.stage in {"overlays", "all"}:
        plot_overlays(paths, config, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
