"""Plotting utilities for all-scale continuum images."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from astropy import units as u
from astropy.io import fits
from astropy.visualization import (
    ImageNormalize,
    LinearStretch,
    SqrtStretch,
)
from astropy.visualization.wcsaxes import add_beam, add_scalebar
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
import matplotlib.patheffects as path_effects
from matplotlib.colors import ListedColormap
import numpy as np
from regions import Regions


def image_plane(data: np.ndarray) -> np.ndarray:
    """Return a two-dimensional image plane without changing its precision."""
    image = np.squeeze(np.asanyarray(data))
    if image.ndim != 2:
        raise ValueError(
            f"Continuum image must be two-dimensional after squeezing: "
            f"{image.shape}"
        )
    return image


def world_crop_limits(
    header: fits.Header,
    *,
    center_ra_deg: float,
    center_dec_deg: float,
    width_deg: float,
    height_deg: float,
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Convert a sky-centred rectangular view to pixel-axis limits."""
    if width_deg <= 0 or height_deg <= 0:
        raise ValueError("Plot width and height must be positive")

    wcs = WCS(header).celestial
    center_x, center_y = wcs.world_to_pixel_values(
        center_ra_deg,
        center_dec_deg,
    )
    scales = np.asarray(proj_plane_pixel_scales(wcs), dtype=float)
    if scales.shape != (2,) or np.any(~np.isfinite(scales)):
        raise ValueError("Could not determine two-dimensional WCS pixel scales")

    width_pixels = width_deg / abs(scales[0])
    height_pixels = height_deg / abs(scales[1])
    return (
        (center_x - width_pixels / 2.0, center_x + width_pixels / 2.0),
        (center_y - height_pixels / 2.0, center_y + height_pixels / 2.0),
    )


def highlighted_emission(
    data: np.ndarray,
    *,
    transition: float,
) -> np.ma.MaskedArray:
    """Mask emission at and below the dual-colour transition value."""
    image = image_plane(data)
    return np.ma.masked_where(
        (~np.isfinite(image)) | (image <= transition),
        image,
        copy=False,
    )


def truncated_mako_colormap(
    *,
    sample_count: int = 256,
    remove_brightest: int = 20,
) -> ListedColormap:
    """Return the paper colour map with its brightest end removed."""
    if sample_count < 2:
        raise ValueError("sample_count must be at least two")
    if not 0 <= remove_brightest < sample_count:
        raise ValueError(
            "remove_brightest must be non-negative and below sample_count"
        )

    import seaborn as sns

    base = sns.color_palette(
        "mako",
        n_colors=sample_count,
        as_cmap=True,
    )
    colours = base(np.linspace(0.0, 1.0, sample_count))
    if remove_brightest:
        colours = colours[:-remove_brightest]
    return ListedColormap(colours, name="mako_cut")


def _apply_crop(
    axis: Any,
    limits: tuple[tuple[float, float], tuple[float, float]],
) -> None:
    axis.set_xlim(limits[0])
    axis.set_ylim(limits[1])


def _add_regions(axis: Any, path: Path) -> None:
    for region in Regions.read(path, format="ds9"):
        pixel_region = region.to_pixel(axis.wcs)
        axis.add_patch(
            pixel_region.as_artist(
                facecolor="none",
                edgecolor="yellow",
                linestyle="-",
                linewidth=0.5,
            )
        )


def create_dual_colour_figure(
    data: np.ndarray,
    header: fits.Header,
    *,
    center_ra_deg: float,
    center_dec_deg: float,
    width_deg: float,
    height_deg: float,
    background_vmin: float,
    transition: float,
    highlight_vmax: float,
    figsize: tuple[float, float] = (10.0, 6.0),
    annotated: bool = False,
    region_path: str | Path | None = None,
    title: str | None = None,
    resolution_label: str | None = None,
    scalebar_length_deg: float = 0.02,
    scalebar_label: str = "1 pc",
    colourbar_label: str = r"Intensity (MJy sr$^{-1}$)",
    background_ticks: tuple[float, ...] | None = None,
    image_axes_rect: tuple[float, float, float, float] = (
        0.05,
        0.05,
        0.80,
        0.90,
    ),
) -> Any:
    """Create the dual-colour all-scale continuum view."""
    if not background_vmin < transition < highlight_vmax:
        raise ValueError(
            "Expected background_vmin < transition < highlight_vmax"
        )

    from matplotlib import pyplot as plt

    image = image_plane(data)
    wcs = WCS(header).celestial
    limits = world_crop_limits(
        header,
        center_ra_deg=center_ra_deg,
        center_dec_deg=center_dec_deg,
        width_deg=width_deg,
        height_deg=height_deg,
    )
    background_norm = ImageNormalize(
        stretch=LinearStretch(),
        vmin=background_vmin,
        vmax=transition,
    )
    highlight_norm = ImageNormalize(
        stretch=SqrtStretch(),
        vmin=transition,
        vmax=highlight_vmax,
    )
    highlight = highlighted_emission(image, transition=transition)
    colour_map = truncated_mako_colormap()

    figure = plt.figure(figsize=figsize, dpi=100)
    background_axis = figure.add_axes(
        image_axes_rect,
        projection=wcs,
    )
    background_artist = background_axis.imshow(
        image,
        origin="lower",
        cmap="Greys",
        norm=background_norm,
    )
    _apply_crop(background_axis, limits)

    overlay_axis = figure.add_axes(
        image_axes_rect,
        projection=wcs,
    )
    overlay_axis.patch.set_alpha(0.0)
    highlight_artist = overlay_axis.imshow(
        highlight,
        origin="lower",
        cmap=colour_map,
        norm=highlight_norm,
    )
    _apply_crop(overlay_axis, limits)
    overlay_axis.set_xlabel("")
    overlay_axis.set_ylabel("")

    if not annotated:
        background_axis.set_xlabel("")
        background_axis.set_ylabel("")
        for axis in (background_axis, overlay_axis):
            axis.tick_params(
                axis="x",
                bottom=False,
                labelbottom=False,
                top=False,
            )
            axis.tick_params(
                axis="y",
                left=False,
                labelleft=False,
                right=False,
            )
        return figure

    if region_path is not None:
        _add_regions(overlay_axis, Path(region_path))
    background_axis.set_xlabel("RA (ICRS)")
    background_axis.set_ylabel("Dec (ICRS)", labelpad=0)
    background_axis.tick_params(
        axis="x",
        bottom=True,
        labelbottom=True,
        top=False,
    )
    background_axis.tick_params(
        axis="y",
        left=True,
        labelleft=True,
        right=False,
    )

    if all(key in header for key in ("BMAJ", "BMIN", "BPA")):
        add_beam(
            overlay_axis,
            header,
            edgecolor="black",
            facecolor="white",
            frame=True,
        )
    add_scalebar(
        overlay_axis,
        length=scalebar_length_deg * u.deg,
        label=scalebar_label,
        color="black",
    )
    if resolution_label:
        overlay_axis.text(
            0.04,
            0.035,
            resolution_label,
            transform=overlay_axis.transAxes,
            horizontalalignment="left",
            verticalalignment="top",
            color="black",
            size=10,
        )
    if title:
        overlay_axis.text(
            0.02,
            0.98,
            title,
            transform=overlay_axis.transAxes,
            horizontalalignment="left",
            verticalalignment="top",
            color="black",
            size=15,
            weight="bold",
            path_effects=[
                path_effects.withStroke(
                    linewidth=2.5,
                    foreground="white",
                )
            ],
        )

    background_colourbar_axis = figure.add_axes(
        [0.86, 0.05, 0.03, 0.90]
    )
    background_colourbar = figure.colorbar(
        background_artist,
        cax=background_colourbar_axis,
    )
    if background_ticks is not None:
        background_colourbar.set_ticks(background_ticks)
    background_colourbar.ax.tick_params(
        axis="y",
        right=True,
        labelright=True,
        labelleft=False,
        left=False,
        labelsize=8,
        length=2,
        rotation=0,
    )

    highlight_colourbar_axis = figure.add_axes([0.93, 0.05, 0.03, 0.90])
    highlight_colourbar = figure.colorbar(
        highlight_artist,
        cax=highlight_colourbar_axis,
    )
    highlight_colourbar.set_label(
        colourbar_label,
        fontsize=12,
    )
    highlight_colourbar.ax.tick_params(
        axis="y",
        labelright=True,
        labelleft=False,
        left=False,
        right=True,
        labelsize=8,
        length=2,
        rotation=0,
    )
    return figure
