"""Core-formation summaries and catalogue-overlay plotting utilities."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import Table
from astropy.visualization import ImageNormalize, SqrtStretch
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
import matplotlib.patheffects as path_effects
from matplotlib.patches import Ellipse
import numpy as np

from .continuum_plotting import image_plane


def counting_fractional_uncertainty(count: float) -> float:
    """Return the adopted half-Poisson fractional uncertainty."""
    if np.isinf(count):
        return 0.0
    if not np.isfinite(count) or count <= 1:
        raise ValueError("Finite source counts must be greater than one")
    return 0.5 / np.sqrt(count - 1.0)


def dust_mass_opacity_scale(
    reference_opacity_cm2_per_g: float,
    adopted_opacity_cm2_per_g: float,
) -> float:
    """Return the mass scale for changing an optically thin dust opacity."""
    reference = float(reference_opacity_cm2_per_g)
    adopted = float(adopted_opacity_cm2_per_g)
    if reference <= 0 or adopted <= 0:
        raise ValueError("Dust opacities must be positive")
    return reference / adopted


def summarize_core_formation(
    catalogue: Table,
    region_settings: Mapping[str, Mapping[str, object]],
    *,
    reference_opacity_cm2_per_g: float,
    adopted_opacity_cm2_per_g: float,
) -> Table:
    """Summarize cloud, fragment, and core masses by evolutionary region."""
    required = {"Subregion", "Mass_core"}
    missing = required - set(catalogue.colnames)
    if missing:
        raise KeyError(
            "Physical catalogue is missing columns: "
            + ", ".join(sorted(missing))
        )

    subregions = np.asarray(catalogue["Subregion"]).astype(str)
    opacity_scale = dust_mass_opacity_scale(
        reference_opacity_cm2_per_g,
        adopted_opacity_cm2_per_g,
    )
    rows: list[dict[str, object]] = []
    for key, settings in region_settings.items():
        code = str(settings["code"])
        mask = subregions == code
        core_count = int(np.count_nonzero(mask))
        if core_count == 0:
            raise ValueError(f"No cores found for subregion code {code!r}")

        core_mass = float(np.nansum(catalogue["Mass_core"][mask]))
        if "err_Mass_core" in catalogue.colnames:
            core_mass_error = float(
                np.sqrt(
                    np.nansum(
                        np.asarray(
                            catalogue["err_Mass_core"][mask],
                            dtype=float,
                        )
                        ** 2
                    )
                )
            )
        else:
            core_mass_error = np.nan

        cloud_mass = float(settings["cloud_mass_msun"])
        paper1_fragment_mass = float(
            settings["paper1_fragment_mass_msun"]
        )
        fragment_mass = paper1_fragment_mass * opacity_scale
        fragment_count = int(settings["fragment_count"])
        if cloud_mass <= 0 or fragment_mass <= 0:
            raise ValueError("Cloud and fragment masses must be positive")

        fragment_efficiency = fragment_mass / cloud_mass
        core_efficiency = core_mass / cloud_mass
        fragment_to_core_efficiency = core_mass / fragment_mass
        fragment_to_core_efficiency_error = (
            core_mass_error / fragment_mass
            if np.isfinite(core_mass_error)
            else np.nan
        )
        rows.append(
            {
                "region": key,
                "label": str(settings["label"]),
                "subregion_code": code,
                "cloud_mass": cloud_mass,
                "paper1_fragment_mass": paper1_fragment_mass,
                "fragment_mass": fragment_mass,
                "core_mass": core_mass,
                "core_mass_error": core_mass_error,
                "fragment_count": fragment_count,
                "core_count": core_count,
                "fragment_efficiency": fragment_efficiency,
                "core_efficiency": core_efficiency,
                "fragment_to_core_efficiency": (
                    fragment_to_core_efficiency
                ),
                "fragment_to_core_efficiency_core_mass_error": (
                    fragment_to_core_efficiency_error
                ),
                "fragment_efficiency_counting_error": (
                    fragment_efficiency
                    * counting_fractional_uncertainty(fragment_count)
                ),
                "core_efficiency_counting_error": (
                    core_efficiency
                    * counting_fractional_uncertainty(core_count)
                ),
            }
        )

    output = Table(rows=rows)
    for name in (
        "cloud_mass",
        "paper1_fragment_mass",
        "fragment_mass",
        "core_mass",
        "core_mass_error",
    ):
        output[name].unit = u.M_sun
    output.meta["mass_efficiency_definition"] = (
        "mass at each smaller scale divided by the cloud mass"
    )
    output.meta["fragment_to_core_efficiency_definition"] = (
        "core mass divided by opacity-corrected fragment mass"
    )
    output.meta["fragment_to_core_efficiency_error"] = (
        "quadrature core-mass error divided by corrected fragment mass; "
        "fragment-mass uncertainty is not included"
    )
    output.meta["counting_uncertainty"] = (
        "efficiency * 0.5 / sqrt(number - 1)"
    )
    output.meta["paper1_fragment_opacity_cm2_per_g"] = float(
        reference_opacity_cm2_per_g
    )
    output.meta["adopted_fragment_opacity_cm2_per_g"] = float(
        adopted_opacity_cm2_per_g
    )
    output.meta["fragment_mass_opacity_scale"] = opacity_scale
    return output


def add_catalogue_ellipses(
    axis: Any,
    wcs: WCS,
    catalogue: Table,
    *,
    ra_column: str,
    dec_column: str,
    major_column: str,
    minor_column: str,
    position_angle_column: str,
    position_angle_type: str = "astronomical",
    size_unit: str = "arcsec",
    expansion_factor: float = 1.0,
    edgecolor: str = "white",
    linewidth: float = 1.2,
    alpha: float = 0.9,
    zorder: int = 10,
) -> list[Ellipse]:
    """Overlay catalogue ellipses on an image axis and return the artists."""
    required = {
        ra_column,
        dec_column,
        major_column,
        minor_column,
        position_angle_column,
    }
    missing = required - set(catalogue.colnames)
    if missing:
        raise KeyError(
            "Ellipse catalogue is missing columns: "
            + ", ".join(sorted(missing))
        )
    if expansion_factor <= 0:
        raise ValueError("expansion_factor must be positive")
    if position_angle_type not in {"astronomical", "mathematical"}:
        raise ValueError(
            "position_angle_type must be astronomical or mathematical"
        )

    celestial = wcs.celestial
    scales_deg = np.asarray(
        proj_plane_pixel_scales(celestial),
        dtype=float,
    )
    pixel_scale_arcsec = (
        np.sqrt(abs(scales_deg[0] * scales_deg[1])) * u.deg
    ).to_value(u.arcsec)
    artists: list[Ellipse] = []
    for row in catalogue:
        coordinates = np.asarray(
            [
                row[ra_column],
                row[dec_column],
                row[major_column],
                row[minor_column],
                row[position_angle_column],
            ],
            dtype=float,
        )
        if not np.all(np.isfinite(coordinates)):
            continue

        sky = SkyCoord(
            coordinates[0] * u.deg,
            coordinates[1] * u.deg,
            frame="icrs",
        )
        x_pixel, y_pixel = sky.to_pixel(celestial)
        major = coordinates[2] * expansion_factor
        minor = coordinates[3] * expansion_factor
        if size_unit == "arcsec":
            width = major / pixel_scale_arcsec
            height = minor / pixel_scale_arcsec
        elif size_unit == "deg":
            width = major * 3600.0 / pixel_scale_arcsec
            height = minor * 3600.0 / pixel_scale_arcsec
        elif size_unit == "pixel":
            width = major
            height = minor
        else:
            raise ValueError("size_unit must be arcsec, deg, or pixel")

        angle = coordinates[4]
        if position_angle_type == "astronomical":
            angle += 90.0
        artist = Ellipse(
            (x_pixel, y_pixel),
            width=width,
            height=height,
            angle=angle,
            edgecolor=edgecolor,
            facecolor="none",
            linewidth=linewidth,
            alpha=alpha,
            zorder=zorder,
        )
        axis.add_patch(artist)
        artists.append(artist)
    return artists


def create_core_overlay_figure(
    high_resolution_data: np.ndarray,
    high_resolution_header: fits.Header,
    low_resolution_data: np.ndarray,
    low_resolution_header: fits.Header,
    source_catalogue: Table,
    fragment_catalogue: Table,
    *,
    center_ra_deg: float,
    center_dec_deg: float,
    width_arcsec: float,
    height_arcsec: float,
    label: str,
    image_vmin: float,
    image_vmax: float,
    contour_levels: tuple[float, ...],
    excluded_fragment_catalogue: Table | None = None,
    figsize: tuple[float, float] = (4.0, 4.0),
) -> Any:
    """Create one high-resolution core view with low-resolution contours."""
    if width_arcsec <= 0 or height_arcsec <= 0:
        raise ValueError("Panel width and height must be positive")
    if image_vmin >= image_vmax:
        raise ValueError("image_vmin must be below image_vmax")

    from matplotlib import pyplot as plt

    high_image = image_plane(high_resolution_data)
    low_image = image_plane(low_resolution_data)
    center = SkyCoord(
        center_ra_deg * u.deg,
        center_dec_deg * u.deg,
        frame="icrs",
    )
    cutout = Cutout2D(
        high_image,
        center,
        (height_arcsec * u.arcsec, width_arcsec * u.arcsec),
        wcs=WCS(high_resolution_header).celestial,
        mode="partial",
        fill_value=0.0,
    )
    image = np.nan_to_num(
        cutout.data,
        copy=True,
        nan=0.0,
        posinf=0.0,
        neginf=0.0,
    )
    norm = ImageNormalize(
        stretch=SqrtStretch(),
        vmin=image_vmin,
        vmax=image_vmax,
    )

    figure = plt.figure(figsize=figsize, dpi=100)
    axis = figure.add_axes(
        [0.05, 0.05, 0.9, 0.9],
        projection=cutout.wcs,
    )
    axis.imshow(
        image,
        origin="lower",
        cmap="gray",
        norm=norm,
        zorder=10,
    )
    add_catalogue_ellipses(
        axis,
        cutout.wcs,
        source_catalogue,
        ra_column="WCS_ACOOR",
        dec_column="WCS_DCOOR",
        major_column="AFWHM01",
        minor_column="BFWHM01",
        position_angle_column="THETA01",
        position_angle_type="astronomical",
        expansion_factor=1.5,
        edgecolor="cyan",
        linewidth=0.6,
        zorder=30,
    )
    low_wcs = WCS(low_resolution_header).celestial
    axis.contour(
        low_image,
        levels=contour_levels,
        transform=axis.get_transform(low_wcs),
        colors="yellow",
        linewidths=0.5,
        zorder=15,
    )
    add_catalogue_ellipses(
        axis,
        cutout.wcs,
        fragment_catalogue,
        ra_column="x_cen",
        dec_column="y_cen",
        major_column="major_sigma",
        minor_column="minor_sigma",
        position_angle_column="position_angle",
        position_angle_type="mathematical",
        expansion_factor=2.35,
        edgecolor="red",
        linewidth=1.5,
        zorder=20,
    )
    if excluded_fragment_catalogue is not None:
        add_catalogue_ellipses(
            axis, cutout.wcs, excluded_fragment_catalogue,
            ra_column="x_cen", dec_column="y_cen",
            major_column="major_sigma", minor_column="minor_sigma",
            position_angle_column="position_angle",
            position_angle_type="mathematical", expansion_factor=2.35,
            edgecolor="orange", linewidth=1.5, alpha=0.9, zorder=20,
        )
    axis.set_xlim(-0.5, image.shape[1] - 0.5)
    axis.set_ylim(-0.5, image.shape[0] - 0.5)

    axis.text(
        0.02,
        0.98,
        label,
        transform=axis.transAxes,
        horizontalalignment="left",
        verticalalignment="top",
        color="white",
        fontsize=15,
        fontweight="bold",
        zorder=40,
        path_effects=[
            path_effects.withStroke(linewidth=2.5, foreground="black")
        ],
    )
    axis.set_xlabel("")
    axis.set_ylabel("")
    axis.tick_params(
        axis="x",
        which="both",
        bottom=False,
        top=False,
        labelbottom=False,
    )
    axis.tick_params(
        axis="y",
        which="both",
        left=False,
        right=False,
        labelleft=False,
    )
    axis.coords[0].set_ticks_visible(False)
    axis.coords[1].set_ticks_visible(False)
    axis.coords[0].set_ticklabel_visible(False)
    axis.coords[1].set_ticklabel_visible(False)
    try:
        axis.coords.frame.set_visible(False)
    except AttributeError:
        pass
    for spine in axis.spines.values():
        spine.set_visible(False)
    return figure


def draw_temperature_efficiency_cascade(ax, table, original):
    """Retain the adopted hierarchy and overlay the preferred heating model."""
    from matplotlib.lines import Line2D

    COLORS = ["#246FA8", "#D4851E", "#BD3654"]
    MARKERS = ["o", "s", "^"]
    scales = np.array([4., .15, .01])
    assert np.all(scales > 0)
    selected = table[(table["coupling"]=="independent") & (table["prior"]=="brightness")]
    for row, base, color, marker in zip(selected, original, COLORS, MARKERS):
        assert str(row["code"]) == str(base["subregion_code"])
        for sampled, adopted in [("baseline_mass_msun", "core_mass"),
                                 ("fragment_mass_msun", "fragment_mass"),
                                 ("cloud_mass_msun", "cloud_mass")]:
            assert np.isclose(row[sampled], base[adopted]), "Temperature scenarios require the current regional masses"
        fragment = row["fragment_mass_msun"]/row["cloud_mass_msun"]
        baseline = row["baseline_mass_msun"]/row["cloud_mass_msun"]
        assert fragment > 0 and baseline > 0
        low, middle, high = [row[f"core_cloud_{p}_percent"]/100 for p in ("p16", "p50", "p84")]
        assert 0 < low <= middle <= high <= baseline
        # Interpolation is a guide between measured reservoir scales only.
        x = np.geomspace(scales[1], scales[2], 80)
        weight = np.linspace(0, 1, len(x))
        lower = np.exp((1-weight)*np.log(fragment)+weight*np.log(low))
        upper = np.exp((1-weight)*np.log(fragment)+weight*np.log(high))
        ax.fill_between(x, lower, upper, color=color, alpha=.18, linewidth=0)
        ax.plot(scales, [1., fragment, baseline], color=color, lw=1.6,
                marker=marker, ms=5, mfc="white", mew=1.1)
        # Keep the original counting-error convention separate from temperature.
        ax.errorbar(scales[1:], [fragment, baseline],
            yerr=[base["fragment_efficiency_counting_error"], base["core_efficiency_counting_error"]],
            fmt="none", color=color, alpha=.65, lw=.75, capsize=2)
        ax.plot(scales[1:], [fragment, middle], color=color, ls="--", lw=1.25)
        ax.errorbar([scales[2]], [middle], yerr=[[middle-low], [high-middle]],
            fmt=marker, ms=4.5, mfc=color, mec="white", mew=.5,
            color=color, lw=1.2, capsize=3, zorder=5)
    ax.set(xscale="log", yscale="log", xlim=(6., .0045), ylim=(.0007, 2.4))
    ax.set_xticks(scales, ["Cloud\n4 pc", "Fragment\n0.15 pc", "Core\n0.01 pc"])
    ax.set_yticks([.001, .01, .1, 1], ["0.1%", "1%", "10%", "100%"])
    ax.set_xlabel("Physical scale")
    ax.set_ylabel("Mass / cloud mass")
    ax.tick_params(direction="in", pad=5)
    regions = [Line2D([], [], color=c, marker=m, lw=1.5, label=str(r["label"]))
               for r,c,m in zip(original, COLORS, MARKERS)]
    legend = ax.legend(handles=regions, loc="lower left", fontsize=9, handlelength=1.6)
    ax.add_artist(legend)
    models = [Line2D([], [], color=".3", marker="o", mfc="white", lw=1.5, label="Adopted Tenv"),
              Line2D([], [], color=".3", marker="o", ls="--", lw=1.25, label="Brightness prior")]
    ax.legend(handles=models, loc="upper right", bbox_to_anchor=(1., .95),
              fontsize=8.5, handlelength=1.8)

