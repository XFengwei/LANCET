#!/usr/bin/env python
"""Build provisional Planck polarization and magnetic-field diagnostics."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
import json
from pathlib import Path
import tomllib

from astropy.io import fits
from astropy.visualization import ImageNormalize, SqrtStretch
from astropy.wcs import WCS
import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.collections import LineCollection

from lancet import ProjectPaths
from lancet.magnetic_fields import planck_polarization_products


DEFAULT_CONFIG = Path("tasks/magnetic_fields/task.toml")


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
    """Resolve delivered Stokes inputs and the policy-safe result root."""
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
        if not paths[name].is_relative_to(project.delivered.resolve()):
            raise ValueError(f"{name} must be a delivered input")
    return paths


def _read_stokes(
    paths: Mapping[str, Path],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, fits.Header]:
    intensity, header = fits.getdata(
        paths["stokes_i"],
        header=True,
        memmap=True,
    )
    q, q_header = fits.getdata(
        paths["stokes_q"],
        header=True,
        memmap=True,
    )
    u, u_header = fits.getdata(
        paths["stokes_u"],
        header=True,
        memmap=True,
    )
    arrays = [
        np.squeeze(np.asarray(value, dtype=float))
        for value in (intensity, q, u)
    ]
    if any(value.ndim != 2 for value in arrays):
        raise ValueError("Stokes I, Q, and U must reduce to 2D images")
    if len({value.shape for value in arrays}) != 1:
        raise ValueError("Stokes I, Q, and U image shapes do not match")

    reference_wcs = WCS(header).celestial
    sample_x = np.array([0.0, arrays[0].shape[1] - 1.0])
    sample_y = np.array([0.0, arrays[0].shape[0] - 1.0])
    reference_world = reference_wcs.pixel_to_world_values(sample_x, sample_y)
    for name, comparison_header in (
        ("stokes_q", q_header),
        ("stokes_u", u_header),
    ):
        comparison_world = WCS(
            comparison_header
        ).celestial.pixel_to_world_values(sample_x, sample_y)
        if not all(
            np.allclose(first, second, rtol=0, atol=1e-8)
            for first, second in zip(reference_world, comparison_world)
        ):
            raise ValueError(f"{name} WCS does not match stokes_i")
    return arrays[0], arrays[1], arrays[2], header


def _product_paths(paths: Mapping[str, Path]) -> dict[str, Path]:
    root = paths["output_root"] / "maps"
    return {
        "intensity_snr2_mask": root / "stokes_i_snr2_mask.fits",
        "intensity_snr3_mask": root / "stokes_i_snr3_mask.fits",
        "polarized_intensity_debiased": (
            root / "polarized_intensity_debiased.fits"
        ),
        "polarized_snr2_mask": (
            root / "polarized_intensity_snr2_mask.fits"
        ),
        "polarized_snr3_mask": (
            root / "polarized_intensity_snr3_mask.fits"
        ),
        "polarization_fraction_debiased": (
            root / "polarization_fraction_debiased.fits"
        ),
        "polarization_angle": root / "polarization_angle.fits",
        "polarization_angle_error": (
            root / "polarization_angle_error.fits"
        ),
        "magnetic_field_angle": root / "magnetic_field_angle.fits",
    }


def _summary_path(paths: Mapping[str, Path]) -> Path:
    return paths["output_root"] / "tables" / "polarization_summary.json"


def _figure_paths(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
) -> list[Path]:
    plot = _mapping(config, "plot")
    formats = [str(value).lower() for value in plot["formats"]]
    unsupported = set(formats) - {"pdf", "png"}
    if unsupported:
        raise ValueError(
            "Unsupported figure formats: " + ", ".join(sorted(unsupported))
        )
    return [
        paths["output_root"]
        / "figures"
        / f"planck_magnetic_field_diagnostic.{extension}"
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


def _task_status(config: Mapping[str, object]) -> str:
    task = _mapping(config, "task")
    return str(task["status"])


def inspect_inputs(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
) -> None:
    """Report input compatibility and provisional assumptions."""
    intensity, q, u, header = _read_stokes(paths)
    noise = _mapping(config, "noise")
    thresholds = _mapping(config, "thresholds")
    print("Magnetic-field workflow:")
    print(f"  status: {_task_status(config).upper()} / NOT PAPER-READY")
    for name in ("stokes_i", "stokes_q", "stokes_u"):
        print(f"  {name}: {paths[name]}")
    print(
        f"  common shape: {intensity.shape}; "
        f"BUNIT={header.get('BUNIT', 'not recorded')}; "
        f"finite I/Q/U={np.isfinite(intensity).sum()}/"
        f"{np.isfinite(q).sum()}/{np.isfinite(u).sum()}"
    )
    print(
        "Adopted provisional RMS: "
        f"I={noise['rms_i_kcmb']} K_CMB, "
        f"Q/U={noise['rms_qu_kcmb']} K_CMB"
    )
    print(
        "Angle mask: debiased PI S/N >= "
        f"{thresholds['angle_polarized_snr']}; "
        "no field-strength estimate is performed"
    )
    print(
        "Continuity note: keep this task ongoing until RMS estimation, "
        "angle convention, beam treatment, and a strength method are fixed."
    )
    print(f"Output root: {paths['output_root']}")


def _calculate(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
) -> tuple[
    dict[str, np.ndarray],
    np.ndarray,
    fits.Header,
]:
    intensity, q, u, header = _read_stokes(paths)
    noise = _mapping(config, "noise")
    thresholds = _mapping(config, "thresholds")
    products = planck_polarization_products(
        intensity,
        q,
        u,
        rms_i=float(noise["rms_i_kcmb"]),
        rms_qu=float(noise["rms_qu_kcmb"]),
        intensity_snr_thresholds=tuple(
            float(value)
            for value in thresholds["intensity_snr"]
        ),
        polarized_snr_thresholds=tuple(
            float(value)
            for value in thresholds["polarized_snr"]
        ),
        angle_snr_threshold=float(
            thresholds["angle_polarized_snr"]
        ),
        fraction_intensity_snr_threshold=float(
            thresholds["fraction_intensity_snr"]
        ),
    )
    return products, intensity, header


def _output_header(
    source: fits.Header,
    *,
    btype: str,
    bunit: str,
) -> fits.Header:
    header = source.copy()
    header["BTYPE"] = btype
    header["BUNIT"] = bunit
    header.add_history(
        "LANCET magnetic_fields ongoing workflow; assumptions are provisional"
    )
    return header


def build_products(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool,
) -> dict[str, np.ndarray]:
    """Write masks and polarization diagnostics outside delivered data."""
    products, intensity, source_header = _calculate(paths, config)
    outputs = _product_paths(paths)
    summary_output = _summary_path(paths)
    _prepare_outputs(
        [*outputs.values(), summary_output],
        overwrite=overwrite,
    )
    descriptions = {
        "intensity_snr2_mask": ("Stokes I SNR > 2 mask", "1"),
        "intensity_snr3_mask": ("Stokes I SNR > 3 mask", "1"),
        "polarized_intensity_debiased": (
            "Debiased polarized intensity",
            str(source_header.get("BUNIT", "Kcmb")),
        ),
        "polarized_snr2_mask": ("Debiased PI SNR >= 2 mask", "1"),
        "polarized_snr3_mask": ("Debiased PI SNR >= 3 mask", "1"),
        "polarization_fraction_debiased": (
            "Debiased polarization fraction",
            "1",
        ),
        "polarization_angle": ("Polarization position angle", "deg"),
        "polarization_angle_error": (
            "Polarization angle uncertainty",
            "deg",
        ),
        "magnetic_field_angle": (
            "Plane-of-sky magnetic-field orientation",
            "deg",
        ),
    }
    for name, output in outputs.items():
        values = products[name]
        description, unit = descriptions[name]
        if values.dtype == bool:
            written = values.astype(np.uint8)
        else:
            written = np.asarray(values, dtype=np.float32)
        fits.writeto(
            output,
            written,
            _output_header(
                source_header,
                btype=description,
                bunit=unit,
            ),
            overwrite=True,
        )
        print(f"Wrote {output}")

    finite_angle = np.isfinite(products["magnetic_field_angle"])
    finite_fraction = np.isfinite(
        products["polarization_fraction_debiased"]
    )
    summary = {
        "status": _task_status(config),
        "paper_ready": False,
        "shape": list(intensity.shape),
        "magnetic_angle_pixel_count": int(np.count_nonzero(finite_angle)),
        "polarization_fraction_pixel_count": int(
            np.count_nonzero(finite_fraction)
        ),
        "median_polarization_fraction": float(
            np.nanmedian(products["polarization_fraction_debiased"])
        ),
        "continuity_note": (
            "Revisit RMS estimation, coordinate/angle convention, beam "
            "treatment, and field-strength methodology before scientific use."
        ),
    }
    summary_output.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {summary_output}")
    return products


def plot_diagnostic(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool,
) -> list[Path]:
    """Plot Planck intensity with the provisional field orientations."""
    angle_path = _product_paths(paths)["magnetic_field_angle"]
    intensity_mask_path = _product_paths(paths)["intensity_snr3_mask"]
    if not angle_path.is_file() or not intensity_mask_path.is_file():
        raise FileNotFoundError(
            "Magnetic-field maps are unavailable; run --stage build first"
        )
    intensity, header = fits.getdata(
        paths["stokes_i"],
        header=True,
        memmap=True,
    )
    angle = fits.getdata(angle_path, memmap=True)
    intensity_mask = fits.getdata(intensity_mask_path, memmap=True).astype(bool)
    intensity = np.squeeze(np.asarray(intensity, dtype=float))
    angle = np.squeeze(np.asarray(angle, dtype=float))
    intensity_mask = np.squeeze(intensity_mask)
    outputs = _figure_paths(paths, config)
    _prepare_outputs(outputs, overwrite=overwrite)
    plot = _mapping(config, "plot")

    figure = plt.figure(
        figsize=(
            float(plot["width_inches"]),
            float(plot["height_inches"]),
        )
    )
    axis = figure.add_subplot(projection=WCS(header).celestial)
    normalization = ImageNormalize(
        vmin=float(plot["vmin_kcmb"]),
        vmax=float(plot["vmax_kcmb"]),
        stretch=SqrtStretch(),
        clip=True,
    )
    image = axis.imshow(
        intensity,
        origin="lower",
        cmap=str(plot["colormap"]),
        norm=normalization,
    )
    step = int(plot["vector_step_pixels"])
    length = float(plot["vector_length_pixels"])
    y_grid, x_grid = np.mgrid[
        0 : angle.shape[0] : step,
        0 : angle.shape[1] : step,
    ]
    sampled_angle = angle[::step, ::step]
    valid = (
        np.isfinite(sampled_angle)
        & intensity_mask[::step, ::step]
    )
    radians = np.radians(sampled_angle[valid])
    longitude_sign = np.sign(float(header.get("CDELT1", -1.0)))
    latitude_sign = np.sign(float(header.get("CDELT2", 1.0)))
    dx = length * np.sin(radians) / longitude_sign
    dy = length * np.cos(radians) / latitude_sign
    x_centers = x_grid[valid]
    y_centers = y_grid[valid]
    segments = np.stack(
        [
            np.column_stack(
                [x_centers - 0.5 * dx, y_centers - 0.5 * dy]
            ),
            np.column_stack(
                [x_centers + 0.5 * dx, y_centers + 0.5 * dy]
            ),
        ],
        axis=1,
    )
    vectors = LineCollection(
        segments,
        colors=str(plot["vector_color"]),
        linewidths=float(plot["vector_line_width"]),
        capstyle="round",
        zorder=3,
    )
    axis.add_collection(vectors)
    vectors.set_path_effects(
        [
            path_effects.Stroke(
                linewidth=float(plot["vector_outline_width"]),
                foreground=str(plot["vector_outline_color"]),
            ),
            path_effects.Normal(),
        ]
    )
    axis.coords[0].set_axislabel("Galactic longitude")
    axis.coords[1].set_axislabel("Galactic latitude")
    colorbar = figure.colorbar(image, ax=axis, pad=0.02)
    colorbar.set_label(r"Stokes $I$ (K$_{\rm CMB}$)")
    axis.text(
        0.02,
        0.98,
        "Ongoing diagnostic - assumptions provisional",
        transform=axis.transAxes,
        ha="left",
        va="top",
        color="white",
        fontsize=9,
        bbox={
            "facecolor": "black",
            "edgecolor": "none",
            "alpha": 0.55,
            "pad": 2.0,
        },
        path_effects=[
            path_effects.withStroke(linewidth=2.5, foreground="black")
        ],
    )
    figure.tight_layout()
    for output in outputs:
        figure.savefig(
            output,
            dpi=int(plot["dpi"]),
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
        choices=("inspect", "build", "plot", "all"),
        default="inspect",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    project = ProjectPaths.discover()
    config = load_task_config(_project_path(project, args.config))
    paths = resolve_task_paths(project, config)
    if args.stage == "inspect":
        inspect_inputs(paths, config)
    elif args.stage == "build":
        build_products(paths, config, overwrite=args.overwrite)
    elif args.stage == "plot":
        plot_diagnostic(paths, config, overwrite=args.overwrite)
    else:
        build_products(paths, config, overwrite=args.overwrite)
        plot_diagnostic(paths, config, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
