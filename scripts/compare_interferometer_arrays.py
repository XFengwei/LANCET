#!/usr/bin/env python
"""Compare ALMA 12 m and ACA 7 m spatial coverage and continuum images."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
import json
from pathlib import Path
import tomllib

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

from lancet import ProjectPaths
from lancet.interferometry import (
    beam_match_and_reproject,
    read_unique_uv_coordinates,
)


DEFAULT_CONFIG = Path("tasks/compare_interferometer_arrays/task.toml")


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
) -> dict[str, object]:
    """Resolve declared read-only inputs and the policy-safe output root."""
    measurement_sets = config["measurement_sets"]
    images = config["images"]
    regions = config["regions"]
    outputs = config["outputs"]
    if not all(
        isinstance(value, Mapping)
        for value in (measurement_sets, images, regions, outputs)
    ):
        raise TypeError("Task path settings must be TOML tables")

    resolved_sets: dict[str, Path] = {}
    for name, value in measurement_sets.items():
        path = _project_path(project, value)
        if not path.is_relative_to(project.reduction.resolve()):
            raise ValueError(
                f"Measurement set must be a declared reduction input: {path}"
            )
        if not path.is_dir():
            raise FileNotFoundError(f"Measurement set not found: {path}")
        resolved_sets[name] = path

    resolved_images: dict[str, Path] = {}
    for name, value in images.items():
        path = _project_path(project, value)
        if path.is_relative_to(project.data.resolve()):
            raise ValueError(
                f"Prepared image must come from a task output: {path}"
            )
        if not path.is_file():
            raise FileNotFoundError(
                f"Prepared image not found: {path}; "
                "run prepare_continuum_images first"
            )
        resolved_images[name] = path

    resolved_regions: dict[str, Path] = {}
    for name, value in regions.items():
        path = _project_path(project, value)
        if not path.is_file():
            raise FileNotFoundError(f"Region file not found: {path}")
        resolved_regions[name] = path

    return {
        "measurement_sets": resolved_sets,
        "images": resolved_images,
        "regions": resolved_regions,
        "output_root": project.validate_output(str(outputs["root"])),
    }


def _prepare_output(path: Path, *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists: {path}; pass --overwrite to replace it"
        )
    path.parent.mkdir(parents=True, exist_ok=True)


def inspect_task(paths: Mapping[str, object]) -> None:
    """Report input readiness without reading image arrays or UVW columns."""
    print("Declared read-only measurement sets:")
    measurement_sets = paths["measurement_sets"]
    assert isinstance(measurement_sets, Mapping)
    for name, path in measurement_sets.items():
        print(f"  {name}: {path}")

    print("Prepared continuum images:")
    images = paths["images"]
    assert isinstance(images, Mapping)
    for name, path in images.items():
        header = fits.getheader(path)
        shape = tuple(
            int(header[f"NAXIS{index}"])
            for index in range(int(header["NAXIS"]), 0, -1)
        )
        print(
            f"  {name}: shape={shape}, "
            f"BITPIX={int(header['BITPIX'])}, "
            f"BUNIT={header.get('BUNIT', 'not recorded')}"
        )
    print(f"Output root: {paths['output_root']}")


def _figure_outputs(
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
        _prepare_output(output, overwrite=overwrite)
    return outputs


def _save_figure(figure, outputs: Sequence[Path], *, dpi: int) -> None:
    for output in outputs:
        figure.savefig(output, bbox_inches="tight", dpi=dpi)


def _uv_cache(
    path: Path,
    measurement_set: Path,
    *,
    decimals: int,
    chunk_rows: int,
    refresh: bool,
) -> tuple[np.ndarray, np.ndarray]:
    if path.is_file() and not refresh:
        with np.load(path) as stored:
            return np.asarray(stored["u"]), np.asarray(stored["v"])

    path.parent.mkdir(parents=True, exist_ok=True)
    u_values, v_values = read_unique_uv_coordinates(
        measurement_set,
        decimals=decimals,
        chunk_rows=chunk_rows,
    )
    np.savez_compressed(
        path,
        u=np.asarray(u_values, dtype=np.float32),
        v=np.asarray(v_values, dtype=np.float32),
    )
    return u_values, v_values


def plot_uv_coverage(
    paths: Mapping[str, object],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
    refresh_cache: bool = False,
) -> list[Path]:
    """Create the composed UV-coverage manuscript figure."""
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Circle, ConnectionPatch, Rectangle

    settings = config["uv_plot"]
    if not isinstance(settings, Mapping):
        raise TypeError("UV plot settings must be a TOML table")
    output_root = Path(paths["output_root"])
    outputs = _figure_outputs(
        output_root / "figures" / "UVcoverage_sets",
        settings["output_formats"],
        overwrite=overwrite,
    )

    measurement_sets = paths["measurement_sets"]
    assert isinstance(measurement_sets, Mapping)
    coordinates: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for name, measurement_set in measurement_sets.items():
        cache = output_root / "uv_coordinates" / f"{name}.npz"
        u_values, v_values = _uv_cache(
            cache,
            Path(measurement_set),
            decimals=int(settings["round_decimals"]),
            chunk_rows=int(settings["chunk_rows"]),
            refresh=refresh_cache,
        )
        coordinates[str(name)] = (u_values, v_values)
        print(f"{name}: retained {len(u_values)} unique UV coordinates")

    c5_u, c5_v = coordinates["configuration_c5"]
    c2_u, c2_v = coordinates["configuration_c2"]
    aca_u, aca_v = coordinates["aca_7m"]
    combined_u, combined_v = coordinates["combined_12m"]
    zoom = float(settings["zoom_limit_metres"])

    figure = plt.figure(figsize=(7.2, 13.8))
    grid = figure.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.13)
    overview = figure.add_subplot(grid[0])
    zoomed = figure.add_subplot(grid[1])

    overview.scatter(
        c5_u,
        c5_v,
        s=1.2,
        color="deepskyblue",
        alpha=0.12,
        rasterized=True,
    )
    overview.scatter(
        -c5_u,
        -c5_v,
        s=1.2,
        color="deepskyblue",
        alpha=0.12,
        rasterized=True,
    )
    overview.scatter(
        c2_u,
        c2_v,
        s=1.2,
        color="green",
        alpha=0.14,
        rasterized=True,
    )
    overview.scatter(
        -c2_u,
        -c2_v,
        s=1.2,
        color="green",
        alpha=0.14,
        rasterized=True,
    )
    overview.scatter(
        aca_u,
        aca_v,
        s=1.2,
        color="tomato",
        alpha=0.18,
        rasterized=True,
    )
    overview.scatter(
        -aca_u,
        -aca_v,
        s=1.2,
        color="tomato",
        alpha=0.18,
        rasterized=True,
    )
    zoom_box = Rectangle(
        (-zoom, -zoom),
        2.0 * zoom,
        2.0 * zoom,
        edgecolor="black",
        facecolor="none",
        linewidth=1.2,
        label="Zoom-in UV space",
    )
    overview.add_patch(zoom_box)
    extent = 1.08 * max(
        np.max(np.abs(c5_u)),
        np.max(np.abs(c5_v)),
        np.max(np.abs(c2_u)),
        np.max(np.abs(c2_v)),
        np.max(np.abs(aca_u)),
        np.max(np.abs(aca_v)),
    )
    overview.set_xlim(-extent, extent)
    overview.set_ylim(-extent, extent)
    overview.set_xlabel("u (m)")
    overview.set_ylabel("v (m)")
    overview.set_aspect("equal", adjustable="box")
    overview.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor="deepskyblue",
                markeredgecolor="none",
                markersize=6,
                label="ALMA 12m C5",
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor="green",
                markeredgecolor="none",
                markersize=6,
                label="ALMA 12m C2",
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor="tomato",
                markeredgecolor="none",
                markersize=6,
                label="ACA 7m",
            ),
            zoom_box,
        ],
        frameon=False,
        loc="upper right",
    )

    overlap = Circle(
        (0, 0),
        float(settings["overlap_radius_metres"]),
        edgecolor="black",
        facecolor="gray",
        linewidth=1.5,
        hatch="//",
        label="Overlap",
        zorder=1,
    )
    total_power = Circle(
        (0, 0),
        float(settings["total_power_radius_metres"]),
        edgecolor="black",
        facecolor="0.6",
        linewidth=1.5,
        label="Total power",
        zorder=2,
    )
    zoomed.add_patch(overlap)
    zoomed.add_patch(total_power)
    zoomed.scatter(
        combined_u,
        combined_v,
        s=1.3,
        color="deepskyblue",
        alpha=0.30,
        zorder=3,
        rasterized=True,
    )
    zoomed.scatter(
        -combined_u,
        -combined_v,
        s=1.3,
        color="deepskyblue",
        alpha=0.30,
        zorder=3,
        rasterized=True,
    )
    zoomed.scatter(
        aca_u,
        aca_v,
        s=1.3,
        color="tomato",
        alpha=0.45,
        zorder=4,
        rasterized=True,
    )
    zoomed.scatter(
        -aca_u,
        -aca_v,
        s=1.3,
        color="tomato",
        alpha=0.45,
        zorder=4,
        rasterized=True,
    )
    zoomed.set_xlim(-zoom, zoom)
    zoomed.set_ylim(-zoom, zoom)
    zoomed.set_xlabel("u (m)")
    zoomed.set_ylabel("v (m)")
    zoomed.set_aspect("equal", adjustable="box")
    zoomed.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor="deepskyblue",
                markeredgecolor="none",
                markersize=6,
                label="ALMA 12m (C5+C2)",
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor="tomato",
                markeredgecolor="none",
                markersize=6,
                label="ACA 7m",
            ),
            overlap,
            total_power,
        ],
        frameon=True,
        loc="upper right",
    )
    for x_value in (-zoom, zoom):
        connector = ConnectionPatch(
            xyA=(x_value, -zoom),
            coordsA=overview.transData,
            xyB=(x_value, zoom),
            coordsB=zoomed.transData,
            color="black",
            linewidth=1.0,
            clip_on=False,
        )
        figure.add_artist(connector)

    _save_figure(figure, outputs, dpi=int(settings["dpi"]))
    plt.close(figure)
    print(f"Wrote UV-coverage figures to {outputs[0].parent}")
    return outputs


def _read_2d(path: Path) -> tuple[np.ndarray, fits.Header]:
    data, header = fits.getdata(path, header=True)
    image = np.squeeze(np.asarray(data, dtype=float))
    if image.ndim != 2:
        raise ValueError(f"Expected a two-dimensional image: {path}")
    return image, header


def _recenter_wcs(
    axis,
    header: fits.Header,
    *,
    ra_deg: float,
    dec_deg: float,
    width_deg: float,
    height_deg: float,
) -> None:
    wcs = WCS(header).celestial
    x_center, y_center = wcs.world_to_pixel_values(ra_deg, dec_deg)
    pixel_scale_x = abs(float(wcs.wcs.cdelt[0]))
    pixel_scale_y = abs(float(wcs.wcs.cdelt[1]))
    width_pixels = width_deg / pixel_scale_x
    height_pixels = height_deg / pixel_scale_y
    axis.set_xlim(
        x_center - width_pixels / 2.0,
        x_center + width_pixels / 2.0,
    )
    axis.set_ylim(
        y_center - height_pixels / 2.0,
        y_center + height_pixels / 2.0,
    )


def _add_regions(axis, header: fits.Header, path: Path) -> None:
    from regions import Regions

    for region in Regions.read(path, format="ds9"):
        patch = region.to_pixel(WCS(header).celestial).as_artist()
        patch.set_edgecolor("yellow")
        patch.set_facecolor("none")
        patch.set_linewidth(1.0)
        axis.add_patch(patch)


def plot_image_comparison(
    paths: Mapping[str, object],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> list[Path]:
    """Create the vertical ACA, 12 m, and combined image comparison."""
    import matplotlib

    matplotlib.use("Agg")
    from astropy import units as u
    from astropy.visualization import ImageNormalize, SqrtStretch
    from astropy.visualization.wcsaxes import add_beam, add_scalebar
    from matplotlib import pyplot as plt
    import matplotlib as mpl
    import seaborn as sns

    settings = config["image_comparison_plot"]
    if not isinstance(settings, Mapping):
        raise TypeError("Image-comparison settings must be a TOML table")
    output_root = Path(paths["output_root"])
    outputs = _figure_outputs(
        output_root
        / "figures"
        / "Interferometry_combine_comparison_vertical",
        settings["output_formats"],
        overwrite=overwrite,
    )
    images = paths["images"]
    regions = paths["regions"]
    assert isinstance(images, Mapping) and isinstance(regions, Mapping)
    panels = (
        (
            Path(images["aca_pbcor"]),
            "ACA 7m",
            r"$6^{\prime\prime}\!.2$=0.08 pc",
        ),
        (
            Path(images["twelve_m_pbcor"]),
            "ALMA 12m (C5 + C2)",
            r"$0^{\prime\prime}\!.44$=1240 au",
        ),
        (
            Path(images["combined_pbcut"]),
            "ALMA 12m (C5+C2) + ACA 7m",
            r"$0^{\prime\prime}\!.65$=1820 au",
        ),
    )
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
            "font.size": 10,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "axes.linewidth": 0.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    colormap = sns.color_palette(
        str(settings["colormap"]),
        n_colors=256,
        as_cmap=True,
    )
    normalization = ImageNormalize(
        stretch=SqrtStretch(),
        vmin=float(settings["vmin_mjysr"]),
        vmax=float(settings["vmax_mjysr"]),
    )
    figure = plt.figure(figsize=(5, 9), dpi=int(settings["dpi"]))
    grid = figure.add_gridspec(
        4,
        1,
        height_ratios=[0.05, 1, 1, 1],
        hspace=0.02,
    )
    last_image = None
    for index, (path, title, resolution) in enumerate(panels):
        data, header = _read_2d(path)
        axis = figure.add_subplot(
            grid[index + 1, 0],
            projection=WCS(header).celestial,
        )
        last_image = axis.imshow(
            data,
            origin="lower",
            cmap=colormap,
            norm=normalization,
        )
        axis.text(
            0.02,
            0.98,
            title,
            transform=axis.transAxes,
            ha="left",
            va="top",
            color="black",
            size=10,
            weight="bold",
        )
        axis.text(
            0.05,
            0.08,
            resolution,
            transform=axis.transAxes,
            ha="left",
            va="top",
            color="black" if index else "white",
            size=8,
            weight="bold",
        )
        add_beam(
            axis,
            header,
            edgecolor="black",
            facecolor="white",
        )
        add_scalebar(
            axis,
            length=float(settings["scalebar_length_deg"]) * u.deg,
            label=str(settings["scalebar_label"]),
            color="black",
        )
        _add_regions(
            axis,
            header,
            Path(regions["successful_combinations"]),
        )
        _recenter_wcs(
            axis,
            header,
            ra_deg=float(settings["center_ra_deg"]),
            dec_deg=float(settings["center_dec_deg"]),
            width_deg=float(settings["width_deg"]),
            height_deg=float(settings["height_deg"]),
        )
        if index < 2:
            axis.set_xlabel("")
            axis.set_ylabel("")
            axis.coords[0].set_ticklabel_visible(False)
            axis.coords[1].set_ticklabel_visible(False)
        else:
            axis.set_xlabel("RA (ICRS)")
            axis.set_ylabel("Dec (ICRS)", labelpad=-1)

    assert last_image is not None
    colorbar_axis = figure.add_subplot(grid[0, 0])
    colorbar = figure.colorbar(
        last_image,
        cax=colorbar_axis,
        orientation="horizontal",
    )
    colorbar.set_label(r"Intensity (MJy sr$^{-1}$)", fontsize=12)
    colorbar.ax.tick_params(
        axis="x",
        labelbottom=False,
        bottom=False,
        top=True,
        labeltop=True,
    )
    colorbar.ax.xaxis.set_label_position("top")
    colorbar.ax.xaxis.set_ticks_position("top")

    _save_figure(figure, outputs, dpi=int(settings["dpi"]))
    plt.close(figure)
    print(f"Wrote image-comparison figures to {outputs[0].parent}")
    return outputs


def run_beam_match(
    paths: Mapping[str, object],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> dict[str, Path]:
    """Beam-match the combined image and quantify its ACA residual."""
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    settings = config["beam_match"]
    if not isinstance(settings, Mapping):
        raise TypeError("Beam-match settings must be a TOML table")
    images = paths["images"]
    assert isinstance(images, Mapping)
    output_root = Path(paths["output_root"]) / "beam_match"
    outputs = {
        "matched": output_root / "combined_convolved_to_ACA.fits",
        "residual": output_root / "combined_minus_ACA_residual.fits",
        "summary": output_root / "residual_summary.json",
        "figure": (
            output_root
            / f"combined_minus_ACA_residual.{settings['output_format']}"
        ),
    }
    for output in outputs.values():
        _prepare_output(output, overwrite=overwrite)

    combined, combined_header = _read_2d(Path(images["combined_pbcut"]))
    aca_pbcor, aca_header = _read_2d(Path(images["aca_pbcor"]))
    matched, matched_header = beam_match_and_reproject(
        combined,
        combined_header,
        aca_pbcor,
        aca_header,
    )
    matched_header["BTYPE"] = "Intensity"
    fits.writeto(
        outputs["matched"],
        np.asarray(matched, dtype=np.float32),
        matched_header,
        overwrite=overwrite,
    )

    reference_name = str(settings["residual_reference"])
    reference, _ = _read_2d(Path(images[reference_name]))
    if reference.shape != matched.shape:
        raise ValueError("Residual reference does not match the ACA grid")
    residual = matched - reference
    residual_header = matched_header.copy()
    residual_header["BTYPE"] = "Residual"
    residual_header.add_history(
        f"LANCET residual: combined beam match minus {reference_name}"
    )
    fits.writeto(
        outputs["residual"],
        np.asarray(residual, dtype=np.float32),
        residual_header,
        overwrite=overwrite,
    )

    finite = residual[np.isfinite(residual)]
    if finite.size == 0:
        raise ValueError("Beam-match residual has no finite pixels")
    summary = {
        "reference": reference_name,
        "finite_pixels": int(finite.size),
        "mean_mjysr": float(np.mean(finite)),
        "median_mjysr": float(np.median(finite)),
        "standard_deviation_mjysr": float(np.std(finite)),
        "rms_mjysr": float(np.sqrt(np.mean(finite**2))),
        "median_absolute_deviation_mjysr": float(
            np.median(np.abs(finite - np.median(finite)))
        ),
    }
    outputs["summary"].write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    limit = float(np.nanpercentile(np.abs(finite), 99.0))
    figure, axis = plt.subplots(figsize=(6, 4.5))
    image = axis.imshow(
        residual,
        origin="lower",
        cmap="RdBu_r",
        vmin=-limit,
        vmax=limit,
    )
    axis.set_title("Combined interferometer image - ACA")
    axis.set_xlabel("Pixel x")
    axis.set_ylabel("Pixel y")
    colorbar = figure.colorbar(image, ax=axis)
    colorbar.set_label(r"Residual (MJy sr$^{-1}$)")
    figure.tight_layout()
    figure.savefig(
        outputs["figure"],
        bbox_inches="tight",
        dpi=int(settings["dpi"]),
    )
    plt.close(figure)
    print(json.dumps(summary, indent=2))
    return outputs


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--stage",
        choices=("inspect", "uv", "compare-images", "beam-match", "all"),
        default="inspect",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--refresh-uv",
        action="store_true",
        help="Reread measurement sets instead of reusing cached coordinates",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    project = ProjectPaths.discover()
    config_path = args.config
    if not config_path.is_absolute():
        config_path = project.root / config_path
    config = load_task_config(config_path.resolve())
    paths = resolve_task_paths(project, config)

    if args.stage in {"inspect", "all"}:
        inspect_task(paths)
    if args.stage in {"uv", "all"}:
        plot_uv_coverage(
            paths,
            config,
            overwrite=args.overwrite,
            refresh_cache=args.refresh_uv,
        )
    if args.stage in {"compare-images", "all"}:
        plot_image_comparison(paths, config, overwrite=args.overwrite)
    if args.stage in {"beam-match", "all"}:
        run_beam_match(paths, config, overwrite=args.overwrite)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
