#!/usr/bin/env python
"""Stack multi-wavelength views of the LANCET filament on one rotated strip grid."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from pathlib import Path
import tomllib
import warnings

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import FITSFixedWarning, WCS
import astropy.units as u
import numpy as np
from reproject import reproject_adaptive
from scipy import ndimage

from lancet import ProjectPaths
from lancet.continuum_plotting import truncated_mako_colormap
from lancet.infrared_mosaic import (
    boost_saturation,
    coadd_tiles,
    load_tile,
    stretch_channel,
)


DEFAULT_CONFIG = Path("tasks/multiwavelength_strip/task.toml")


def load_toml(path: Path) -> dict[str, object]:
    """Load settings from TOML."""
    if not path.is_file():
        raise FileNotFoundError(f"Configuration not found: {path}")
    with path.open("rb") as handle:
        return tomllib.load(handle)


def strip_wcs(grid: Mapping[str, object]) -> tuple[WCS, tuple[int, int]]:
    """TAN grid whose +x axis points along *long_axis_pa_deg* (east of north)."""
    scale = float(grid["pixel_scale_arcsec"])
    width = int(round(float(grid["width_arcmin"]) * 60.0 / scale))
    height = int(round(float(grid["height_arcmin"]) * 60.0 / scale))
    # In a north-up, east-left image +x points to PA 270; rotate from there.
    angle = np.deg2rad(float(grid["long_axis_pa_deg"]) - 270.0)
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [float(grid["center_ra_deg"]), float(grid["center_dec_deg"])]
    wcs.wcs.crpix = [(width + 1) / 2.0, (height + 1) / 2.0]
    wcs.wcs.cdelt = [-scale / 3600.0, scale / 3600.0]
    wcs.wcs.pc = [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
    wcs.wcs.radesys = "ICRS"
    return wcs, (height, width)


def read_image(path: Path) -> tuple[np.ndarray, WCS]:
    """Return the first 2-D celestial image in *path*."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        with fits.open(path) as hdul:
            hdu = next(h for h in hdul if h.data is not None)
            data = np.squeeze(hdu.data).astype(np.float32)
            wcs = WCS(hdu.header).celestial
    if data.ndim != 2:
        raise ValueError(f"Expected a 2-D image in {path}, got shape {data.shape}")
    return data, wcs


def fill_interior_holes(data: np.ndarray) -> np.ndarray:
    """Give enclosed NaN holes (e.g. saturated peaks) their nearest valid value."""
    finite = np.isfinite(data)
    holes = ndimage.binary_fill_holes(finite) & ~finite
    if not holes.any():
        return data
    nearest = ndimage.distance_transform_edt(
        ~finite, return_distances=False, return_indices=True
    )
    filled = data.copy()
    filled[holes] = data[tuple(index[holes] for index in nearest)]
    return filled


def reproject_image(
    path: Path,
    wcs_out: WCS,
    shape_out: tuple[int, int],
    fill_holes: bool = False,
) -> np.ndarray:
    """Resample a surface-brightness image onto the strip grid (anti-aliased)."""
    data, wcs = read_image(path)
    if fill_holes:
        data = fill_interior_holes(data)
    image, footprint = reproject_adaptive(
        (data, wcs),
        wcs_out,
        shape_out=shape_out,
        conserve_flux=False,
        kernel="gaussian",
        boundary_mode="strict",
    )
    image = image.astype(np.float32)
    image[footprint < 0.5] = np.nan
    return image


def vvv_band_images(
    project: ProjectPaths,
    vvv_config: Mapping[str, object],
    wcs_out: WCS,
    shape_out: tuple[int, int],
    pixel_arcsec: float,
) -> dict[str, np.ndarray]:
    """Seamless VVV mosaics of every band on the strip grid."""
    inputs = vvv_config["inputs"]
    blending = vvv_config["blending"]
    zero_point = float(vvv_config["photometry"]["reference_zero_point"])
    directory = project.root / str(inputs["tile_directory"])
    edge_trim = int(blending["edge_trim_pixels"])
    images = {}
    for band in inputs["bands"]:
        tiles = [
            load_tile(directory / f"{tile}_{band}.fits", zero_point, edge_trim)
            for tile in inputs["tiles"]
        ]
        images[band] = coadd_tiles(
            tiles,
            wcs_out,
            shape_out,
            block=int(blending["block_pixels"]),
            smoothness=float(blending["smoothness"]),
            decay_pixels=float(blending["decay_arcsec"]) / pixel_arcsec,
            feather_pixels=float(blending["feather_arcsec"]) / pixel_arcsec,
        )
    return images


def grid_path(output_root: Path, panel: str, key: str) -> Path:
    return output_root / "grid" / f"{panel}_{key}.fits"


def run_reproject(project: ProjectPaths, config: Mapping[str, object]) -> None:
    """Write every panel input on the common strip grid."""
    output_root = project.validate_output(str(config["outputs"]["root"]))
    (output_root / "grid").mkdir(parents=True, exist_ok=True)
    wcs_out, shape_out = strip_wcs(config["grid"])
    header = wcs_out.to_header()
    pixel_arcsec = float(config["grid"]["pixel_scale_arcsec"])

    for panel in config["panels"]:
        name = str(panel["name"])
        if panel["kind"] == "vvv":
            vvv_config = load_toml(project.root / str(panel["vvv_config"]))
            images = vvv_band_images(
                project, vvv_config, wcs_out, shape_out, pixel_arcsec
            )
        elif panel["kind"] == "rgb":
            images = {
                channel: reproject_image(project.root / path, wcs_out, shape_out)
                for channel, path in zip(("r", "g", "b"), panel["inputs"])
            }
        else:
            image = reproject_image(
                project.root / str(panel["input"]),
                wcs_out,
                shape_out,
                fill_holes=bool(panel.get("fill_holes", False)),
            )
            if "mask_input" in panel:
                response = reproject_image(
                    project.root / str(panel["mask_input"]), wcs_out, shape_out
                )
                with np.errstate(invalid="ignore"):
                    image[~(response >= float(panel["mask_minimum"]))] = np.nan
            images = {"i": image}
        for key, image in images.items():
            path = grid_path(output_root, name, key)
            fits.PrimaryHDU(image, header=header).writeto(path, overwrite=True)
            coverage = np.isfinite(image).mean()
            print(f"Wrote {path.name} coverage={coverage:.3f}")


def asinh_scale(
    image: np.ndarray, low_percentile: float, high_percentile: float, softening: float
) -> np.ndarray:
    """Percentile-clipped asinh stretch onto [0, 1]; NaN stays NaN."""
    low, high = np.nanpercentile(image, [low_percentile, high_percentile])
    scaled = np.clip((image - low) / (high - low), 0.0, 1.0)
    return np.arcsinh(scaled / softening) / np.arcsinh(1.0 / softening)


def dual_colour_rgb(
    image: np.ndarray,
    background_vmin: float,
    transition: float,
    highlight_vmax: float,
) -> np.ndarray:
    """Paper dual-colour scaling: linear greys up to *transition*, then a
    square-root truncated mako map up to *highlight_vmax*."""
    import matplotlib

    finite = np.isfinite(image)
    value = np.where(finite, image, background_vmin)
    grey = np.clip((value - background_vmin) / (transition - background_vmin), 0, 1)
    rgb = matplotlib.colormaps["Greys"](grey)[..., :3]
    highlight = finite & (value > transition)
    level = np.sqrt(
        np.clip((value[highlight] - transition) / (highlight_vmax - transition), 0, 1)
    )
    rgb[highlight] = truncated_mako_colormap()(level)[..., :3]
    return rgb


def panel_rgb(
    project: ProjectPaths,
    output_root: Path,
    panel: Mapping[str, object],
) -> np.ndarray:
    """Render one panel as an RGB array with rows running bottom-up."""
    import matplotlib

    name = str(panel["name"])
    kind = panel["kind"]
    if kind == "vvv":
        rgb_config = load_toml(project.root / str(panel["vvv_config"]))["rgb"]
        channels = []
        for colour, weight in zip(("red", "green", "blue"), rgb_config["weights"]):
            data = fits.getdata(grid_path(output_root, name, rgb_config[colour]))
            channel = stretch_channel(
                data,
                float(rgb_config["high_percentile"]),
                float(rgb_config["asinh_softening"]),
                float(rgb_config["black_sigma"]),
            )
            channels.append(np.clip(channel * float(weight), 0.0, 1.0))
        return boost_saturation(
            np.stack(channels, axis=-1), float(rgb_config["saturation"])
        )
    if kind == "rgb":
        channels = []
        for index, key in enumerate(("r", "g", "b")):
            data = fits.getdata(grid_path(output_root, name, key))
            channel = asinh_scale(
                data,
                float(panel["low_percentile"][index]),
                float(panel["high_percentile"][index]),
                float(panel["asinh_softening"][index]),
            )
            channels.append(np.clip(channel * float(panel["weights"][index]), 0, 1))
        return np.nan_to_num(np.stack(channels, axis=-1), nan=0.0)

    data = fits.getdata(grid_path(output_root, name, "i"))
    if panel.get("style") == "dual_colour":
        return dual_colour_rgb(
            data,
            float(panel["background_vmin"]),
            float(panel["transition"]),
            float(panel["highlight_vmax"]),
        )
    scaled = asinh_scale(
        data,
        float(panel["low_percentile"]),
        float(panel["high_percentile"]),
        float(panel["asinh_softening"]),
    )
    cmap = matplotlib.colormaps[str(panel["cmap"])].with_extremes(bad="black")
    return cmap(np.ma.masked_invalid(scaled))[..., :3]


def compass_vectors(wcs: WCS, shape: tuple[int, int]) -> dict[str, np.ndarray]:
    """Unit pixel-space vectors pointing north and east at the grid centre."""
    cy, cx = (shape[0] - 1) / 2.0, (shape[1] - 1) / 2.0
    centre = wcs.pixel_to_world(cx, cy)
    vectors = {}
    for label, pa in (("N", 0.0), ("E", 90.0)):
        tip = centre.directional_offset_by(pa * u.deg, 10 * u.arcsec)
        x, y = wcs.world_to_pixel(tip)
        vector = np.array([x - cx, y - cy])
        vectors[label] = vector / np.hypot(*vector)
    return vectors


def run_plot(project: ProjectPaths, config: Mapping[str, object]) -> None:
    """Stack the panels into one figure without coordinate axes."""
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as effects

    output_root = project.validate_output(str(config["outputs"]["root"]))
    figure_config = config["figure"]
    wcs, shape = strip_wcs(config["grid"])
    panels = config["panels"]
    pixel_arcsec = float(config["grid"]["pixel_scale_arcsec"])

    width = float(figure_config["width_inches"])
    panel_height = width * shape[0] / shape[1]
    gap = float(figure_config["panel_gap_inches"])
    height = len(panels) * panel_height + (len(panels) - 1) * gap
    fig = plt.figure(figsize=(width, height), facecolor="black")
    outline = [effects.withStroke(linewidth=2.5, foreground="black")]
    fontsize = float(figure_config["label_fontsize"])

    for index, panel in enumerate(panels):
        bottom = 1.0 - ((index + 1) * panel_height + index * gap) / height
        ax = fig.add_axes((0.0, bottom, 1.0, panel_height / height))
        ax.imshow(
            panel_rgb(project, output_root, panel),
            origin="lower",
            interpolation="none",
        )
        ax.set_axis_off()
        label_colour = str(panel.get("label_colour", "white"))
        label_outline = "white" if label_colour == "black" else "black"
        ax.text(
            0.008,
            0.93,
            str(panel["label"]),
            transform=ax.transAxes,
            color=label_colour,
            fontsize=fontsize,
            va="top",
            path_effects=[
                effects.withStroke(linewidth=2.5, foreground=label_outline)
            ],
        )

        if index == 0 and figure_config["compass"]:
            origin = np.array([shape[1] - 140.0, 120.0])
            for label, vector in compass_vectors(wcs, shape).items():
                tip = origin + 80.0 * vector
                ax.annotate(
                    "",
                    xy=tip,
                    xytext=origin,
                    arrowprops=dict(arrowstyle="-|>", color="white", lw=1.2),
                )
                ax.text(
                    *(origin + 105.0 * vector),
                    label,
                    color="white",
                    fontsize=fontsize * 0.85,
                    ha="center",
                    va="center",
                    path_effects=outline,
                )

        if index == len(panels) - 1:
            length = float(figure_config["scalebar_length_deg"]) * 3600 / pixel_arcsec
            x1 = shape[1] - 60.0
            y = 70.0
            ax.plot([x1 - length, x1], [y, y], color="white", lw=2.5)
            ax.text(
                x1 - length / 2,
                y + 25,
                str(figure_config["scalebar_label"]),
                color="white",
                fontsize=fontsize * 0.9,
                ha="center",
                va="bottom",
                path_effects=outline,
            )

    for fmt in figure_config["formats"]:
        output = output_root / f"lancet_multiwavelength_strip.{fmt}"
        fig.savefig(output, dpi=int(figure_config["dpi"]), facecolor="black")
        print(f"Wrote {output}")
    plt.close(fig)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--stage",
        choices=("reproject", "plot", "all"),
        default="all",
        help="reproject: cache inputs on the strip grid; plot: render the figure",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    project = ProjectPaths.discover()
    config_path = args.config if args.config.is_absolute() else project.root / args.config
    config = load_toml(config_path)
    if args.stage in ("reproject", "all"):
        run_reproject(project, config)
    if args.stage in ("plot", "all"):
        run_plot(project, config)


if __name__ == "__main__":
    main()
