#!/usr/bin/env python
"""Calculate and plot foreground/background optical-depth effects."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
import json
from pathlib import Path
import tomllib

from astropy import units as u
from astropy.table import Column, Table
import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt

from lancet import ProjectPaths
from lancet.imaging import jybeam_to_mjysr
from lancet.physical_properties import (
    blackbody_cgs,
    core_optical_depth,
    symmetric_environment_optical_depth,
)


DEFAULT_CONFIG = Path("tasks/background_effects/task.toml")


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


def resolve_output_root(
    project: ProjectPaths,
    config: Mapping[str, object],
) -> Path:
    """Resolve the policy-safe task output root."""
    outputs = _mapping(config, "outputs")
    return project.validate_output(str(outputs["root"]))


def table_path(output_root: Path) -> Path:
    return output_root / "tables" / "background_effects.ecsv"


def summary_path(output_root: Path) -> Path:
    return output_root / "tables" / "background_effects_summary.json"


def figure_paths(
    output_root: Path,
    config: Mapping[str, object],
) -> list[Path]:
    figures = _mapping(config, "figures")
    formats = [str(value).lower() for value in figures["formats"]]
    unsupported = set(formats) - {"pdf", "png"}
    if unsupported:
        raise ValueError(
            "Unsupported figure formats: " + ", ".join(sorted(unsupported))
        )
    return [
        output_root / "figures" / f"Background_effects.{extension}"
        for extension in formats
    ]


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


def _model_values(
    config: Mapping[str, object],
) -> dict[str, object]:
    model = _mapping(config, "model")
    scalar_names = (
        "frequency_ghz",
        "temperature_k",
        "beam_major_arcsec",
        "beam_minor_arcsec",
        "background_min_mjysr",
        "background_max_mjysr",
    )
    values = {name: float(model[name]) for name in scalar_names}
    if np.any(~np.isfinite(list(values.values()))):
        raise ValueError("Model scalar settings must be finite")
    if np.any(np.asarray(list(values.values())) <= 0):
        raise ValueError("Model scalar settings must be positive")
    samples = int(model["background_samples"])
    if samples < 2:
        raise ValueError("background_samples must be at least two")
    brightness = np.asarray(
        model["core_brightness_mjy_beam"],
        dtype=float,
    )
    if (
        brightness.ndim != 1
        or brightness.size == 0
        or np.any(~np.isfinite(brightness))
        or np.any(brightness <= 0)
    ):
        raise ValueError(
            "core_brightness_mjy_beam must contain positive values"
        )
    values["background_samples"] = samples
    values["core_brightness_mjy_beam"] = brightness
    return values


def _brightness_label(value: float) -> str:
    text = f"{value:g}".replace("-", "m").replace(".", "p")
    return f"{text}mJyBeam"


def calculate_table(
    config: Mapping[str, object],
) -> tuple[Table, dict[str, object]]:
    """Calculate opacity corrections across the configured background grid."""
    model = _model_values(config)
    frequency = float(model["frequency_ghz"]) * u.GHz
    temperature = float(model["temperature_k"]) * u.K
    background = np.geomspace(
        float(model["background_min_mjysr"]),
        float(model["background_max_mjysr"]),
        int(model["background_samples"]),
    ) * u.MJy / u.sr
    blackbody = blackbody_cgs(frequency, temperature).to(u.MJy / u.sr)
    if np.max(background) >= blackbody:
        raise ValueError(
            "Every background intensity must be below the Planck function"
        )

    brightness_mjy_beam = np.asarray(
        model["core_brightness_mjy_beam"],
        dtype=float,
    )
    brightness_mjysr = jybeam_to_mjysr(
        brightness_mjy_beam * 1e-3,
        major_fwhm_deg=(
            float(model["beam_major_arcsec"]) / 3600.0
        ),
        minor_fwhm_deg=(
            float(model["beam_minor_arcsec"]) / 3600.0
        ),
    )
    table = Table()
    table["Background_intensity"] = Column(
        background.to_value(u.MJy / u.sr),
        unit=u.MJy / u.sr,
    )
    environmental_tau = symmetric_environment_optical_depth(
        nu=frequency,
        temp=temperature,
        Ibg=background,
    )
    table["Environment_optical_depth"] = np.asarray(
        environmental_tau.to_value(u.dimensionless_unscaled),
        dtype=float,
    )

    zero_background = 0.0 * u.MJy / u.sr
    for brightness_beam, brightness_sr in zip(
        brightness_mjy_beam,
        brightness_mjysr,
        strict=True,
    ):
        core_intensity = brightness_sr * u.MJy / u.sr
        optical_depth = core_optical_depth(
            nu=frequency,
            temp=temperature,
            Isr=core_intensity,
            Ibg=background,
        )
        baseline = core_optical_depth(
            nu=frequency,
            temp=temperature,
            Isr=core_intensity,
            Ibg=zero_background,
        )
        label = _brightness_label(float(brightness_beam))
        table[f"Core_optical_depth_{label}"] = np.asarray(
            optical_depth.to_value(u.dimensionless_unscaled),
            dtype=float,
        )
        table[f"Relative_core_optical_depth_{label}"] = np.asarray(
            (optical_depth / baseline).to_value(
                u.dimensionless_unscaled
            ),
            dtype=float,
        )

    table.meta["frequency_ghz"] = float(model["frequency_ghz"])
    table.meta["temperature_k"] = float(model["temperature_k"])
    table.meta["beam_major_arcsec"] = float(
        model["beam_major_arcsec"]
    )
    table.meta["beam_minor_arcsec"] = float(
        model["beam_minor_arcsec"]
    )
    table.meta["normalization"] = "zero-background core optical depth"
    table.meta["environment_geometry"] = (
        "equal-temperature symmetric foreground and background"
    )
    summary: dict[str, object] = {
        "frequency_ghz": float(model["frequency_ghz"]),
        "temperature_k": float(model["temperature_k"]),
        "planck_function_mjysr": blackbody.to_value(u.MJy / u.sr),
        "beam_major_arcsec": float(model["beam_major_arcsec"]),
        "beam_minor_arcsec": float(model["beam_minor_arcsec"]),
        "core_brightness_mjy_beam": brightness_mjy_beam.tolist(),
        "core_brightness_mjysr": brightness_mjysr.tolist(),
        "background_range_mjysr": [
            background[0].to_value(u.MJy / u.sr),
            background[-1].to_value(u.MJy / u.sr),
        ],
        "background_samples": len(background),
        "environment_geometry": table.meta["environment_geometry"],
        "normalization": table.meta["normalization"],
    }
    return table, summary


def inspect_task(
    output_root: Path,
    config: Mapping[str, object],
) -> None:
    """Report configured assumptions and calculation readiness."""
    _, summary = calculate_table(config)
    print("Background-effects model:")
    print(
        f"  frequency={summary['frequency_ghz']} GHz, "
        f"temperature={summary['temperature_k']} K"
    )
    print(
        "  Planck function="
        f"{float(summary['planck_function_mjysr']):.6g} MJy/sr"
    )
    print(
        "  core brightness [mJy/beam]: "
        + ", ".join(
            f"{float(value):g}"
            for value in summary["core_brightness_mjy_beam"]
        )
    )
    print(
        "  core brightness [MJy/sr]: "
        + ", ".join(
            f"{float(value):.6g}"
            for value in summary["core_brightness_mjysr"]
        )
    )
    print(f"Output root: {output_root}")


def build_table(
    output_root: Path,
    config: Mapping[str, object],
    *,
    overwrite: bool,
) -> tuple[Path, Path]:
    """Write the opacity grid and its calculation summary."""
    table, summary = calculate_table(config)
    products = [table_path(output_root), summary_path(output_root)]
    _prepare_outputs(products, overwrite=overwrite)
    table.write(products[0], format="ascii.ecsv", overwrite=True)
    products[1].write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {products[0]}")
    print(f"Wrote {products[1]}")
    return products[0], products[1]


def plot_table(
    output_root: Path,
    config: Mapping[str, object],
    *,
    overwrite: bool,
) -> list[Path]:
    """Create the background-opacity sensitivity figure."""
    source = table_path(output_root)
    if not source.is_file():
        raise FileNotFoundError(
            f"Calculation table not found: {source}; run --stage table"
        )
    table = Table.read(source, format="ascii.ecsv")
    outputs = figure_paths(output_root, config)
    _prepare_outputs(outputs, overwrite=overwrite)
    model = _model_values(config)
    figures = _mapping(config, "figures")
    colors = [str(value) for value in figures["colors"]]
    brightness = np.asarray(
        model["core_brightness_mjy_beam"],
        dtype=float,
    )
    if len(colors) < len(brightness):
        raise ValueError("figures.colors must cover every core brightness")

    figure, axis = plt.subplots(
        figsize=(
            float(figures["width_inches"]),
            float(figures["height_inches"]),
        ),
        dpi=int(figures["dpi"]),
    )
    background = np.asarray(table["Background_intensity"], dtype=float)
    for value, color in zip(brightness, colors, strict=False):
        label = _brightness_label(float(value))
        axis.plot(
            background,
            np.asarray(
                table[f"Relative_core_optical_depth_{label}"],
                dtype=float,
            ),
            color=color,
            linewidth=2.0,
            marker="o",
            markevery=int(figures["marker_every"]),
            markersize=float(figures["marker_size"]),
            markeredgecolor="black",
            markeredgewidth=0.5,
            label=rf"$I_{{\rm core}}={value:g}$ mJy beam$^{{-1}}$",
        )
    axis.set_xscale(str(figures["xscale"]))
    axis.set_xlabel(
        r"Environmental emission $I_{\nu}^{\rm env}$ "
        r"(MJy sr$^{-1}$)"
    )
    axis.set_ylabel(
        r"$\tau_{\nu}^{\rm core}(I_{\nu}^{\rm env})/"
        r"\tau_{\nu}^{\rm core}(0)$"
    )
    axis.grid(alpha=0.18, which="both")
    axis.legend(loc="upper left", frameon=False, fontsize=8.5)

    right = axis.twinx()
    right.plot(
        background,
        np.asarray(table["Environment_optical_depth"], dtype=float),
        color="0.35",
        linewidth=2.0,
        linestyle="--",
        label=r"$\tau_{\nu}^{\rm env}$",
    )
    right.set_yscale("log")
    right.set_ylabel(r"Environmental optical depth $\tau_{\nu}^{\rm env}$")
    right.legend(loc="lower right", frameon=False, fontsize=8.5)
    figure.tight_layout()
    for output in outputs:
        figure.savefig(
            output,
            dpi=int(figures["dpi"]),
            bbox_inches="tight",
        )
        print(f"Wrote {output}")
    plt.close(figure)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--stage",
        choices=("inspect", "table", "plot", "all"),
        default="inspect",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    project = ProjectPaths.discover()
    config_path = _project_path(project, args.config)
    config = load_task_config(config_path)
    output_root = resolve_output_root(project, config)
    if args.stage == "inspect":
        inspect_task(output_root, config)
    elif args.stage == "table":
        build_table(output_root, config, overwrite=args.overwrite)
    elif args.stage == "plot":
        plot_table(output_root, config, overwrite=args.overwrite)
    else:
        build_table(output_root, config, overwrite=args.overwrite)
        plot_table(output_root, config, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
