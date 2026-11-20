"""Diagnostic (phase 2) runners of the combine_continuum task.

uv coverage, the ACA / 12 m / 12m+7m image comparison and the pointing coverage, with the
figure helpers that every phase 2 stage shares. Each runner takes the task paths with
``output_root`` set to the diagnostics folder.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import csv
import json
import os
from pathlib import Path

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

from .imaging import read_frequency_ghz, read_mjysr
from .interferometry import read_unique_uv_coordinates


def _prepare_output(path: Path, *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists: {path}; pass --overwrite to replace it"
        )
    path.parent.mkdir(parents=True, exist_ok=True)


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
    inputs = paths["inputs"]
    regions = paths["regions"]
    assert isinstance(inputs, Mapping) and isinstance(regions, Mapping)
    panels = (
        (
            Path(inputs["aca_pbcor"]),
            "ACA 7m",
            r"$6^{\prime\prime}\!.2$=0.08 pc",
        ),
        (
            Path(inputs["twelve_m_pbcor"]),
            "ALMA 12m (C5 + C2)",
            r"$0^{\prime\prime}\!.44$=1240 au",
        ),
        (
            Path(inputs["interferometer"]),
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
        data, header = read_mjysr(path)
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


def run_amplitude_alignment(paths, config, *, save_figure, overwrite=False):
    """Compute the TP / interferometer amplitude alignment, then redraw it from the saved products."""
    from .amplitude_alignment import replot_alignment
    from .amplitude_alignment import run_amplitude_alignment as compute

    manifest = compute(paths, config, save_figure=save_figure, overwrite=overwrite)
    replot_alignment(paths, config, save_figure=save_figure)
    return manifest


def run_pointing_coverage(paths, config, *, overwrite=False):
    """12m and 7m pointings on the MIRIAD gain and CASA pb, one figure per region."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from .pointing_coverage import (
        gain_against_pb, pb_at_gain_edge, plot_coverage, primary_beam_fwhm_arcsec,
        read_field_centres, unique_positions,
    )

    settings = config["pointing_coverage"]
    root = Path(paths["output_root"]) / "pointing_coverage"
    stems = [root / f"pointing_coverage_{r['name']}" for r in settings["regions"]]
    outputs = [s.with_suffix(e) for s in stems for e in (".pdf", ".png")]
    outputs += [root / "pointings.csv", root / "gain_vs_pb.json"]
    for path in outputs:
        _prepare_output(path, overwrite=overwrite)
    images = {"gain": Path(paths["inputs"]["gain"]), "pb": Path(paths["inputs"]["casa_pb"])}

    def centres(names):
        ra, dec = [], []
        for name in names:
            r, d = read_field_centres(paths["measurement_sets"][name])
            ra.append(r); dec.append(d)
        return unique_positions(np.concatenate(ra), np.concatenate(dec),
                                settings["duplicate_tolerance_arcsec"])

    # Primary beams at the TP model frequency, the frequency of the combined products.
    freq, factor = read_frequency_ghz(paths["inputs"]["tp_dust_ff"]), settings["beam_factor"]
    arrays = [
        {"label": "12m pointings", "ra": None, "dec": None, "sets": settings["twelve_m_sets"],
         "fwhm_arcsec": primary_beam_fwhm_arcsec(12.0, freq, factor), "marker": ".", "size": 40,
         "colour": "#9467bd"},
        {"label": "7m pointings", "ra": None, "dec": None, "sets": settings["seven_m_sets"],
         "fwhm_arcsec": primary_beam_fwhm_arcsec(7.0, freq, factor), "marker": "+", "size": 90,
         "colour": "#2ca02c"},
    ]
    for spec in arrays:
        spec["ra"], spec["dec"] = centres(spec["sets"])
        print(f"{spec['label']}: {spec['ra'].size} unique centres from {spec['sets']}")

    def load(key, step):
        with fits.open(images[key]) as handle:
            data = np.squeeze(handle[0].data)[::step, ::step].astype(float)
            wcs = WCS(handle[0].header).celestial[::step, ::step]
        return data, wcs

    step = settings["display_step"]
    gain, gain_wcs = load("gain", step)
    pb, pb_wcs = load("pb", step)
    for spec_region, stem in zip(settings["regions"], stems):
        fig = plt.figure(figsize=(9, 8) if spec_region["size_arcsec"][0] < 400 else (12, 7.5))
        centre = SkyCoord(*spec_region["centre"].split(), unit=(u.hourangle, u.deg))
        plot_coverage(fig, gain, gain_wcs, pb, pb_wcs, arrays, pb_levels=settings["pb_levels"],
                      region=(centre, tuple(spec_region["size_arcsec"])),
                      title=f"{spec_region['title']}: 12m and 7m pointings, MIRIAD gain and CASA 12m+7m pb")
        _save_figure(fig, [stem.with_suffix(".pdf"), stem.with_suffix(".png")], dpi=160)
        plt.close(fig)

    with (root / "pointings.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["array", "ra_deg", "dec_deg"])
        for spec in arrays:
            for r, d in zip(spec["ra"], spec["dec"]):
                writer.writerow([spec["label"].split()[0], f"{r:.7f}", f"{d:.7f}"])

    from reproject import reproject_interp

    step = settings["statistics_step"]
    gain, gain_wcs = load("gain", step)
    pb, pb_wcs = load("pb", step)
    pb_on_gain, _ = reproject_interp((pb, pb_wcs), gain_wcs, shape_out=gain.shape)
    summary = {
        "gain": str(images["gain"]), "pb": str(images["pb"]), "frequency_ghz": freq,
        "statistics_step_pixels": step,
        "primary_beam_fwhm_arcsec": {s["label"]: s["fwhm_arcsec"] for s in arrays},
        "pointings": {s["label"]: int(s["ra"].size) for s in arrays},
        "gain_by_pb_bin": gain_against_pb(gain, pb_on_gain, settings["pb_bin_edges"]),
        "pb_along_gain_edge": pb_at_gain_edge(gain, pb_on_gain),
    }
    (root / "gain_vs_pb.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary["pb_along_gain_edge"]))
    return outputs
