#!/usr/bin/env python
"""Plot the delivered all-scale continuum images."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from pathlib import Path
import tomllib

from astropy.io import fits

from lancet import ProjectPaths
from lancet.continuum_plotting import create_dual_colour_figure


DEFAULT_CONFIG = Path("tasks/plot_allscale_continuum/task.toml")


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
    """Resolve all read-only inputs and the policy-safe output root."""
    inputs = config["inputs"]
    outputs = config["outputs"]
    if not isinstance(inputs, Mapping) or not isinstance(outputs, Mapping):
        raise TypeError("Task inputs and outputs must be TOML tables")

    paths = {
        name: _project_path(project, value)
        for name, value in inputs.items()
    }
    paths["output_root"] = project.validate_output(str(outputs["root"]))

    for name in ("combined_continuum", "dust_continuum"):
        path = paths[name]
        if not path.is_relative_to(project.delivered.resolve()):
            raise ValueError(f"{name} must be a delivered input: {path}")
    synthetic = paths["synthetic_continuum"]
    if synthetic.is_relative_to(project.data.resolve()):
        raise ValueError(
            "synthetic_continuum must come from a task result, not data/"
        )

    for name in (
        "combined_continuum",
        "dust_continuum",
        "synthetic_continuum",
        "primary_beam_regions",
    ):
        if not paths[name].is_file():
            raise FileNotFoundError(f"{name} not found: {paths[name]}")
    return paths


def _describe(path: Path) -> str:
    with fits.open(path, memmap=True) as hdul:
        header = hdul[0].header
        shape = tuple(hdul[0].shape)
        bitpix = int(header["BITPIX"])
        unit = header.get("BUNIT", "not recorded")
    return f"shape={shape}, BITPIX={bitpix}, BUNIT={unit}"


def inspect_inputs(paths: Mapping[str, Path]) -> None:
    """Report task readiness without reading full image arrays."""
    print("All-scale continuum inputs:")
    for name in (
        "combined_continuum",
        "dust_continuum",
        "synthetic_continuum",
    ):
        print(f"  {name}: {paths[name]}")
        print(f"    {_describe(paths[name])}")
    print(f"  primary_beam_regions: {paths['primary_beam_regions']}")
    print(f"Output root: {paths['output_root']}")


def _prepare_outputs(
    stem: Path,
    formats: Sequence[object],
    *,
    overwrite: bool,
) -> list[Path]:
    suffixes = [str(value).lower() for value in formats]
    unsupported = set(suffixes) - {"pdf", "png"}
    if unsupported:
        raise ValueError(
            "Unsupported figure formats: " + ", ".join(sorted(unsupported))
        )
    outputs = [stem.with_suffix(f".{suffix}") for suffix in suffixes]
    for output in outputs:
        if output.exists() and not overwrite:
            raise FileExistsError(
                f"Output already exists: {output}; "
                "pass --overwrite to replace it"
            )
        output.parent.mkdir(parents=True, exist_ok=True)
    return outputs


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


def _plot_settings(
    config: Mapping[str, object],
    *,
    style_name: str,
) -> tuple[Mapping[str, object], Mapping[str, object]]:
    view = config["view"]
    styles = config["styles"]
    if not isinstance(view, Mapping) or not isinstance(styles, Mapping):
        raise TypeError("View and styles settings must be TOML tables")
    style = styles[style_name]
    if not isinstance(style, Mapping):
        raise TypeError(f"Style {style_name} must be a TOML table")
    return view, style


def _create_figure(
    path: Path,
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    style_name: str,
    annotated: bool,
    title: str | None = None,
):
    view, style = _plot_settings(config, style_name=style_name)
    with fits.open(path, memmap=True) as hdul:
        data = hdul[0].data
        header = hdul[0].header.copy()
        figure = create_dual_colour_figure(
            data,
            header,
            center_ra_deg=float(view["center_ra_deg"]),
            center_dec_deg=float(view["center_dec_deg"]),
            width_deg=float(view["width_deg"]),
            height_deg=float(view["height_deg"]),
            background_vmin=float(style["background_vmin"]),
            transition=float(style["transition"]),
            highlight_vmax=float(style["highlight_vmax"]),
            figsize=(
                float(view["figure_width_inches"]),
                float(view["figure_height_inches"]),
            ),
            annotated=annotated,
            region_path=(
                paths["primary_beam_regions"] if annotated else None
            ),
            title=title,
            resolution_label=(
                str(view["resolution_label"]) if annotated else None
            ),
            scalebar_length_deg=float(view["scalebar_length_deg"]),
            scalebar_label=str(view["scalebar_label"]),
        )
    return figure


def plot_combined(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> list[Path]:
    """Plot the delivered dust-plus-free-free all-scale map."""
    figures = config["figures"]
    if not isinstance(figures, Mapping):
        raise TypeError("Figure settings must be a TOML table")
    output_root = paths["output_root"] / "figures"
    written: list[Path] = []

    annotated_outputs = _prepare_outputs(
        output_root / "allscale_continuum_dust_freefree",
        figures["annotated_formats"],
        overwrite=overwrite,
    )
    annotated = _create_figure(
        paths["combined_continuum"],
        paths,
        config,
        style_name="allscale",
        annotated=True,
        title="12m7mTP combined 1.3 mm continuum map",
    )
    _save_figure(
        annotated,
        annotated_outputs,
        dpi=int(figures["combined_annotated_dpi"]),
    )
    written.extend(annotated_outputs)

    clean_outputs = _prepare_outputs(
        output_root / "allscale_continuum_dust_freefree_clean",
        figures["clean_formats"],
        overwrite=overwrite,
    )
    clean = _create_figure(
        paths["combined_continuum"],
        paths,
        config,
        style_name="allscale",
        annotated=False,
    )
    _save_figure(clean, clean_outputs, dpi=int(figures["clean_dpi"]))
    written.extend(clean_outputs)
    return written


def plot_synthetic(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> list[Path]:
    """Plot the current synthesized-continuum model."""
    figures = config["figures"]
    if not isinstance(figures, Mapping):
        raise TypeError("Figure settings must be a TOML table")
    outputs = _prepare_outputs(
        paths["output_root"] / "figures" / "synthetic_continuum_clean",
        figures["clean_formats"],
        overwrite=overwrite,
    )
    figure = _create_figure(
        paths["synthetic_continuum"],
        paths,
        config,
        style_name="synthetic",
        annotated=False,
    )
    _save_figure(figure, outputs, dpi=int(figures["clean_dpi"]))
    return outputs


def plot_dust(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> list[Path]:
    """Plot the delivered dust-only all-scale map."""
    figures = config["figures"]
    if not isinstance(figures, Mapping):
        raise TypeError("Figure settings must be a TOML table")
    output_root = paths["output_root"] / "figures"
    written: list[Path] = []

    annotated_outputs = _prepare_outputs(
        output_root / "allscale_continuum_dust",
        figures["annotated_formats"],
        overwrite=overwrite,
    )
    annotated = _create_figure(
        paths["dust_continuum"],
        paths,
        config,
        style_name="allscale",
        annotated=True,
        title="All-scale 1.3 mm dust continuum map",
    )
    _save_figure(
        annotated,
        annotated_outputs,
        dpi=int(figures["dust_annotated_dpi"]),
    )
    written.extend(annotated_outputs)

    clean_outputs = _prepare_outputs(
        output_root / "allscale_continuum_dust_clean",
        figures["clean_formats"],
        overwrite=overwrite,
    )
    clean = _create_figure(
        paths["dust_continuum"],
        paths,
        config,
        style_name="allscale",
        annotated=False,
    )
    _save_figure(clean, clean_outputs, dpi=int(figures["clean_dpi"]))
    written.extend(clean_outputs)
    return written


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="Task TOML file",
    )
    parser.add_argument(
        "--stage",
        choices=("inspect", "combined", "synthetic", "dust", "all"),
        default="inspect",
        help="Task stage to run",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing task outputs",
    )
    args = parser.parse_args()

    project = ProjectPaths.discover()
    config_path = _project_path(project, args.config)
    config = load_task_config(config_path)
    paths = resolve_task_paths(project, config)

    if args.stage == "inspect":
        inspect_inputs(paths)
        return

    import matplotlib

    matplotlib.use("Agg")
    if args.stage in ("combined", "all"):
        plot_combined(paths, config, overwrite=args.overwrite)
    if args.stage in ("synthetic", "all"):
        plot_synthetic(paths, config, overwrite=args.overwrite)
    if args.stage in ("dust", "all"):
        plot_dust(paths, config, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
