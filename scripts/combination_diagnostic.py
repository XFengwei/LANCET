#!/usr/bin/env python
"""Diagnose continuum data combination: coverage, visibilities and images."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
import csv
import hashlib
import json
import os
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
from lancet.visibility_diagnostics import (
    profile_columns,
    radial_visibility_sums,
    read_visibility_profiles,
)


DEFAULT_CONFIG = Path("tasks/combination_diagnostic/task.toml")


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
        if path.is_relative_to(project.data.resolve()) and name != "allscale_sky":
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

    reduction_inputs = {}
    for name, value in config.get("reduction_inputs", {}).items():
        path = _project_path(project, value)
        if not path.is_relative_to(project.reduction.resolve()) or not path.is_file():
            raise ValueError(f"Expected a declared, read-only reduction file: {path}")
        reduction_inputs[name] = path

    return {
        "measurement_sets": resolved_sets,
        "images": resolved_images,
        "regions": resolved_regions,
        "reduction_inputs": reduction_inputs,
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
    for name, path in paths.get("reduction_inputs", {}).items():
        print(f"Read-only reduction input {name}: {path}")


def _figure_outputs(
    stem: Path,
    formats: Sequence[object],
    *,
    overwrite: bool,
) -> list[Path]:
    suffixes = [str(value).lower() for value in formats]
    unsupported = set(suffixes) - {"pdf", "svg", "png"}
    if unsupported:
        raise ValueError(
            "Unsupported figure formats: " + ", ".join(sorted(unsupported))
        )
    outputs = [stem.with_suffix(f".{suffix}") for suffix in suffixes]
    for output in outputs:
        _prepare_output(output, overwrite=overwrite)
    return outputs


def _save_figure(figure, outputs: Sequence[Path], *, dpi: int, bbox_inches="tight", exclude_axes=()) -> None:
    if os.environ.get("LANCET_FIGURE_QA") == "1":
        from audit_panel_alignment import require_matplotlib_panel_alignment

        require_matplotlib_panel_alignment(
            figure, json_out=outputs[0].with_suffix(".alignment.json"), strict=True,
            exclude_axes=exclude_axes,
        )
    for output in outputs:
        figure.savefig(output, bbox_inches=bbox_inches, dpi=dpi)


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
            if "unflagged_only" in stored and bool(stored["unflagged_only"]):
                return np.asarray(stored["u"]), np.asarray(stored["v"])

    path.parent.mkdir(parents=True, exist_ok=True)
    u_values, v_values = read_unique_uv_coordinates(
        measurement_set,
        decimals=decimals,
        chunk_rows=chunk_rows,
        unflagged_only=True,
    )
    np.savez_compressed(
        path,
        u=np.asarray(u_values, dtype=np.float32),
        v=np.asarray(v_values, dtype=np.float32),
        unflagged_only=True,
    )
    return u_values, v_values


ARRAY_STYLE = {
    "configuration_c5": ("12 m C5", "#7254A3", "-"),
    "configuration_c2": ("12 m C2", "#007A9E", "--"),
    "aca_7m": ("ACA 7 m", "#D87821", "-."),
}


def _plot_style():
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    plt.rcParams.update({
        "font.family": "sans-serif", "font.size": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "svg.fonttype": "none",
        "legend.frameon": False, "axes.linewidth": 0.7,
    })
    return plt


def plot_uv_coverage(paths, config, *, overwrite=False, refresh_cache=False):
    """Show measured full and short-baseline geometry with equal plot areas."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle

    plt = _plot_style()
    settings = config["uv_plot"]
    root = Path(paths["output_root"])
    outputs = _figure_outputs(root / "figures" / "UVcoverage_sets",
                              settings["output_formats"], overwrite=overwrite)
    coordinates = {}
    for name, measurement_set in paths["measurement_sets"].items():
        coordinates[name] = _uv_cache(
            root / "uv_coordinates" / f"{name}.npz", Path(measurement_set),
            decimals=int(settings["round_decimals"]),
            chunk_rows=int(settings["chunk_rows"]), refresh=refresh_cache,
        )
        print(f"{name}: {len(coordinates[name][0])} distinct unflagged UV positions", flush=True)

    figure, axes = plt.subplots(1, 2, figsize=(7.5, 4.0))
    figure.subplots_adjust(left=0.09, right=0.98, bottom=0.17, top=0.84, wspace=0.32)
    zoom = float(settings["zoom_limit_metres"])
    extent = 1.04 * max(np.abs(np.concatenate(coordinates[name])).max()
                        for name in ARRAY_STYLE)
    for index, axis in enumerate(axes):
        limit = extent if index == 0 else zoom
        for name, (label, color, _) in ARRAY_STYLE.items():
            u, v = coordinates[name]
            keep = (np.abs(u) <= limit) & (np.abs(v) <= limit)
            # Plot every distinct position plus its conjugate; no sample thinning.
            axis.scatter(np.r_[u[keep], -u[keep]], np.r_[v[keep], -v[keep]],
                         s=0.65 if index == 0 else 1.0, color=color,
                         alpha=0.24 if name != "aca_7m" else 0.50,
                         linewidths=0, rasterized=True)
        axis.set(xlim=(-limit, limit), ylim=(-limit, limit),
                 xlabel="u (m)", ylabel="v (m)")
        axis.set_aspect("equal", adjustable="box")
        axis.tick_params(direction="out", length=3, width=0.7)
        axis.set_title(["a  Full coverage", "b  Short baselines"][index],
                       loc="left", fontsize=10, pad=10)
    axes[0].add_patch(Rectangle((-zoom, -zoom), 2 * zoom, 2 * zoom,
                               facecolor="none", edgecolor="0.2", linewidth=0.8))
    handles = [Line2D([], [], marker="o", linestyle="none", color=color,
                      markersize=5, label=label) for label, color, _ in ARRAY_STYLE.values()]
    figure.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.54, 1.0),
                  ncol=3, columnspacing=2.0, handletextpad=0.5)
    _save_figure(figure, outputs, dpi=int(settings["dpi"]))
    plt.close(figure)
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


def _write_profile_csv(path, profiles):
    """Write profiles as long-form source data, including empty radial bins."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = None
        for name, columns in profiles.items():
            if writer is None:
                writer = csv.writer(handle)
                writer.writerow(["dataset", *columns])
            writer.writerows([name, *row] for row in zip(*columns.values()))


def run_visibility_diagnostics(paths, config, *, overwrite=False):
    """Generate amplitude/weight profiles from all eligible continuum rows."""
    plt = _plot_style()
    settings = config["visibility_profile"]
    root = Path(paths["output_root"])
    edges = np.geomspace(float(settings["uv_min_klambda"]),
                         float(settings["uv_max_klambda"]),
                         int(settings["radial_bins"]) + 1)
    table_path = root / "visibility_profiles" / "radial_profiles.csv"
    _prepare_output(table_path, overwrite=overwrite)
    profiles, metadata = {}, {}
    for name, ms in paths["measurement_sets"].items():
        print(f"Reading visibility profiles: {name}", flush=True)
        sums, fields, ddids, meta = read_visibility_profiles(
            Path(ms), edges, chunk_rows=int(settings["chunk_rows"]))
        if meta["counts"]["rows_outside_bins"]:
            raise ValueError(f"Radial bins exclude usable data in {name}; extend the range")
        profiles[name] = profile_columns(sums, edges)
        metadata[name] = meta
        np.savez_compressed(table_path.parent / f"{name}.npz", edges=edges,
                            sums=sums, field_sums=fields, ddid_sums=ddids)
        _write_profile_csv(table_path.parent / f"{name}_per_field.csv", {
            str(i): profile_columns(values, edges) for i, values in enumerate(fields)
        })
        _write_profile_csv(table_path.parent / f"{name}_per_ddid.csv", {
            str(i): profile_columns(values, edges) for i, values in enumerate(ddids)
        })
    _write_profile_csv(table_path, profiles)
    metadata_path = table_path.parent / "measurement_sets.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")

    separate = [profiles[k] for k in ("configuration_c5", "configuration_c2")]
    combined = profiles["combined_12m"]
    total_weight = lambda profile: float(profile["sum_weight"].sum())
    total_amp = lambda profile: float(np.nansum(
        profile["mean_amplitude_jy"] * profile["sum_weight"]))
    w_separate = sum(total_weight(p) for p in separate)
    a_separate = sum(total_amp(p) for p in separate) / w_separate
    control = {
        "combined_to_separate_total_weight_ratio": total_weight(combined) / w_separate,
        "combined_to_separate_global_scalar_amplitude_ratio": (
            total_amp(combined) / total_weight(combined)) / a_separate,
        "combined_usable_rows": int(combined["averaged_row_count"].sum()),
        "separate_usable_rows": int(sum(p["averaged_row_count"].sum() for p in separate)),
        "interpretation": "Checks channel averaging and concatenation accounting; not an independent flux calibration test.",
        "radial_support_intersection_klambda": [
            max(metadata[k]["ranges"]["uv_klambda"][0] for k in ARRAY_STYLE),
            min(metadata[k]["ranges"]["uv_klambda"][1] for k in ARRAY_STYLE),
        ],
    }
    (table_path.parent / "consistency_checks.json").write_text(json.dumps(control, indent=2) + "\n")

    # Keep the concatenated 12 m set as a control; never double-count it with C2/C5.
    figures = []
    for zoomed in (False, True):
        stem = "visibility_profiles_short_baselines" if zoomed else "visibility_profiles"
        outputs = _figure_outputs(root / "figures" / stem,
                                  settings["output_formats"], overwrite=overwrite)
        figure, axes = plt.subplots(2, 2, figsize=(7.5, 6.0), sharex=True)
        figure.subplots_adjust(left=0.12, right=0.97, bottom=0.13, top=0.84,
                               hspace=0.40, wspace=0.40)
        metrics = [("mean_amplitude_jy", "Mean amplitude (Jy)"),
                   ("mean_weight", "Mean row weight (MS units)"),
                   ("sum_weight", "Summed weight (MS units)"),
                   ("averaged_row_count", "Averaged visibility rows")]
        for letter, axis, (metric, ylabel) in zip("abcd", axes.flat, metrics):
            for name, (label, color, linestyle) in ARRAY_STYLE.items():
                values = profiles[name]
                use = ((values["averaged_row_count"] > 0)
                       & np.isfinite(values[metric]) & (values[metric] > 0))
                if zoomed:
                    use &= values["uv_mean_klambda"] <= float(settings["zoom_max_klambda"])
                axis.plot(values["uv_mean_klambda"][use], values[metric][use],
                          color=color, linestyle=linestyle, linewidth=1.4, label=label)
            axis.set_yscale("log")
            axis.set_ylabel(ylabel)
            axis.set_title(letter, loc="left", fontweight="bold", pad=8)
            axis.grid(axis="y", alpha=0.12, linewidth=0.5)
            if zoomed:
                axis.set_xlim(0, float(settings["zoom_max_klambda"]))
            else:
                axis.set_xscale("log")
                axis.set_xlim(3, 1200)
        for axis in axes[1]:
            axis.set_xlabel(r"uv distance (k$\lambda$)")
        figure.legend(*axes[0, 0].get_legend_handles_labels(), loc="upper center",
                      bbox_to_anchor=(0.55, 0.98), ncol=3, columnspacing=1.7)
        _save_figure(figure, outputs, dpi=int(settings["dpi"]))
        plt.close(figure)
        figures.extend(outputs)
    return figures


def run_total_power_diagnostic(paths, config, *, overwrite=False):
    """Show the finite-map Fourier amplitudes separately from measured MS data."""
    from astropy.wcs.utils import proj_plane_pixel_area
    from scipy.ndimage import distance_transform_edt

    plt = _plot_style()
    settings = config["visibility_profile"]
    path = Path(paths["reduction_inputs"]["total_power"])
    data, header = _read_2d(path)
    if header.get("BUNIT", "").strip().lower() != "jy/beam":
        raise ValueError("TP Fourier conversion requires a Jy/beam input")
    pixel_area = float(proj_plane_pixel_area(WCS(header).celestial)) * (np.pi / 180)**2
    beam_area = np.pi / (4 * np.log(2)) * float(header["BMAJ"]) * float(header["BMIN"]) * (np.pi / 180)**2
    window = np.outer(np.hanning(data.shape[0]), np.hanning(data.shape[1]))
    # Also taper the internal finite-data mask, whose edge lies inside the frame.
    beam_pixels = np.sqrt(float(header["BMAJ"]) * float(header["BMIN"])) / np.sqrt(pixel_area) * np.pi / 180
    taper_pixels = 2 * beam_pixels
    mask_distance = distance_transform_edt(np.isfinite(data))
    window *= np.sin(0.5 * np.pi * np.minimum(mask_distance / taper_pixels, 1))**2
    transformed = np.fft.fftshift(np.fft.fft2(
        np.where(np.isfinite(data), data, 0) * pixel_area / beam_area * window))
    u = np.fft.fftshift(np.fft.fftfreq(data.shape[1], abs(header["CDELT1"]) * np.pi / 180))
    v = np.fft.fftshift(np.fft.fftfreq(data.shape[0], abs(header["CDELT2"]) * np.pi / 180))
    radius = np.hypot(u[None, :], v[:, None]) / 1000
    # Use only complete Fourier annuli within the rectangular grid's inscribed circle.
    nyquist = min(np.max(np.abs(u)), np.max(np.abs(v))) / 1000
    edges = np.linspace(0, nyquist, 61)
    sums = radial_visibility_sums(radius.ravel(), transformed.ravel(),
                                 np.ones(data.size), edges)
    columns = profile_columns(sums, edges)
    root = Path(paths["output_root"])
    _write_profile_csv(root / "total_power" / "fourier_profile.csv", {"TP_map": columns})
    outputs = _figure_outputs(root / "figures" / "total_power_fourier_diagnostic",
                              settings["output_formats"], overwrite=overwrite)
    figure, axes = plt.subplots(1, 2, figsize=(7.5, 3.5))
    figure.subplots_adjust(left=0.11, right=0.97, bottom=0.2, top=0.84, wspace=0.4)
    axes[0].semilogy(columns["uv_mean_klambda"], columns["mean_amplitude_jy"], color="0.25")
    axes[0].set_ylabel("Mean Fourier amplitude (Jy)")
    axes[0].set_title("a  Windowed TP image", loc="left", fontsize=10)
    rho = np.linspace(0, nyquist * 1000, 300)
    theta = np.sqrt(header["BMAJ"] * header["BMIN"]) * np.pi / 180
    transfer = np.exp(-np.pi**2 * theta**2 * rho**2 / (4 * np.log(2)))
    axes[1].plot(rho / 1000, transfer, color="0.25")
    axes[1].set_ylabel("Gaussian beam response")
    axes[1].set_title("b  TP beam transfer", loc="left", fontsize=10)
    for axis in axes:
        axis.set(xlabel=r"uv distance (k$\lambda$)", xlim=(0, nyquist))
    _save_figure(figure, outputs, dpi=int(settings["dpi"]))
    plt.close(figure)
    script = Path(paths["reduction_inputs"]["combination_script"])
    summary = {
        "input": str(path), "shape": list(data.shape),
        "nonfinite_pixels_zero_filled": int((~np.isfinite(data)).sum()),
        "apodization": "Full-frame Hann times cosine-squared finite-mask edge taper over two beam FWHM; no beam deconvolution",
        "mask_edge_taper_pixels": float(taper_pixels),
        "beam_fwhm_arcsec": float(np.sqrt(header["BMAJ"] * header["BMIN"]) * 3600),
        "reference_frequency_ghz": settings["reference_frequency_ghz"],
        "half_response_klambda": float(2 * np.log(2) / (np.pi * theta) / 1000),
        "windowed_zero_spacing_jy": float(np.abs(transformed[data.shape[0] // 2, data.shape[1] // 2])),
        "complete_annulus_max_klambda": nyquist,
        "combination_script_sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
        "interpretation": "Finite, windowed, beam-convolved model-image transform; not measured TP visibilities, and not comparable to raw mosaic MS amplitudes without a common response. No visibility weights exist for this FITS image.",
    }
    (root / "total_power" / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return outputs


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--stage",
        choices=("inspect", "uv", "visibilities", "overlap-tests", "pb-curve", "full-sampling", "figure1", "figure1-three-panel", "tp-beam-test", "pointing-profiles", "pointing-range", "joint-images", "dirty-beam", "total-power", "compare-images", "beam-match", "all"),
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
    if args.stage in {"visibilities", "all"}:
        run_visibility_diagnostics(paths, config, overwrite=args.overwrite)
    if args.stage in {"total-power", "all"}:
        run_total_power_diagnostic(paths, config, overwrite=args.overwrite)
    if args.stage == "overlap-tests":
        from lancet.combination_tests import run_overlap_tests

        run_overlap_tests(paths, config, save_figure=_save_figure, overwrite=args.overwrite)
    if args.stage == "pb-curve":
        from lancet.combination_tests import run_pb_amplitude_curve

        run_pb_amplitude_curve(paths, config, save_figure=_save_figure, overwrite=args.overwrite)
    if args.stage == "figure1":
        from lancet.combination_figure import make_combination_figure

        make_combination_figure(paths, config, save_figure=_save_figure, overwrite=args.overwrite)
    if args.stage == "figure1-three-panel":
        from lancet.combination_figure import make_combination_figure

        make_combination_figure(paths, config, save_figure=_save_figure, overwrite=args.overwrite, three_panel=True)
    if args.stage == "full-sampling":
        from lancet.full_sampling import build_full_sampling_profiles

        build_full_sampling_profiles(paths, config, overwrite=args.overwrite)
    if args.stage == "tp-beam-test":
        from lancet.tp_beam_test import run_tp_beam_test

        run_tp_beam_test(paths, config, save_figure=_save_figure, overwrite=args.overwrite)
    if args.stage == "pointing-profiles":
        from lancet.pointing_profiles import run_pointing_profiles

        run_pointing_profiles(paths, config, save_figure=_save_figure, overwrite=args.overwrite)
    if args.stage == "pointing-range":
        from lancet.pointing_profiles import run_pointing_range

        run_pointing_range(paths, config, save_figure=_save_figure, overwrite=args.overwrite)
    if args.stage == "joint-images":
        from lancet.joint_image_comparison import run_joint_image_comparison
        run_joint_image_comparison(project.root, config, save_figure=_save_figure, overwrite=args.overwrite)
        return 0

    if args.stage == "dirty-beam":
        from lancet.dirty_beam import run_dirty_beam

        run_dirty_beam(project.root, config, save_figure=_save_figure, overwrite=args.overwrite)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
