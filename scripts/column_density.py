#!/usr/bin/env python
"""Build and display dust-derived column-density products."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
import json
from pathlib import Path
import tomllib

from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
import numpy as np

from lancet import ProjectPaths
from lancet.column_density import (
    brightness_to_column_density,
    column_density_to_surface_density,
    dense_mass_msun,
)
from lancet.continuum_plotting import create_dual_colour_figure


DEFAULT_CONFIG = Path("tasks/column_density/task.toml")


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
    """Resolve declared read-only inputs and the policy-safe output root."""
    inputs = config["inputs"]
    outputs = config["outputs"]
    if not isinstance(inputs, Mapping) or not isinstance(outputs, Mapping):
        raise TypeError("Task inputs and outputs must be TOML tables")
    paths = {
        name: _project_path(project, value)
        for name, value in inputs.items()
    }
    paths["output_root"] = project.validate_output(str(outputs["root"]))

    for name in (
        "temperature_lowres",
        "column_density_lowres",
        "column_density_multiscale",
        "dust_continuum_allscale",
    ):
        if not paths[name].is_relative_to(project.delivered.resolve()):
            raise ValueError(f"{name} must be delivered: {paths[name]}")
    for name, path in paths.items():
        if name != "output_root" and not path.is_file():
            raise FileNotFoundError(f"{name} not found: {path}")
    return paths


def _describe(path: Path) -> str:
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
    """Report dependency readiness without reading full image arrays."""
    print("Column-density inputs:")
    for name in (
        "temperature_lowres",
        "column_density_lowres",
        "column_density_multiscale",
        "dust_continuum_allscale",
        "intercore_background",
        "freefree_smoothed",
    ):
        print(f"  {name}: {paths[name]}")
        print(f"    {_describe(paths[name])}")
    model = config["model"]
    assert isinstance(model, Mapping)
    print(
        "Model: "
        f"frequency={model['frequency_ghz']} GHz, "
        f"opacity={model['opacity_cm2_per_g']} cm2/g, "
        f"mu_H2={model['mean_molecular_weight']}"
    )
    print(f"Output root: {paths['output_root']}")


def _prepare_output(path: Path, *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists: {path}; pass --overwrite to replace it"
        )
    path.parent.mkdir(parents=True, exist_ok=True)


def _float_header(
    header: fits.Header,
    *,
    unit: str,
    history: str,
) -> fits.Header:
    output = header.copy()
    output["BUNIT"] = unit
    if "BLANK" in output:
        del output["BLANK"]
    output.add_history(history)
    return output


def _write_float32(
    path: Path,
    data: np.ndarray,
    header: fits.Header,
    *,
    overwrite: bool,
) -> None:
    _prepare_output(path, overwrite=overwrite)
    fits.writeto(
        path,
        np.asarray(data, dtype=np.float32),
        header,
        overwrite=overwrite,
    )
    print(f"Wrote {path}")


def prepared_paths(paths: Mapping[str, Path]) -> dict[str, Path]:
    root = paths["output_root"] / "prepared"
    return {
        "temperature_allscale": root / "temperature_allscale.fits",
        "surface_density_lowres": root / "surface_density_lowres.fits",
        "surface_density_multiscale": (
            root / "surface_density_multiscale.fits"
        ),
    }


def map_paths(paths: Mapping[str, Path]) -> dict[str, Path]:
    root = paths["output_root"] / "maps"
    return {
        "column_density_allscale": root / "column_density_allscale.fits",
        "column_density_intercore": root / "column_density_intercore.fits",
        "surface_density_allscale": root / "surface_density_allscale.fits",
        "surface_density_intercore": root / "surface_density_intercore.fits",
    }


def prepare_products(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> dict[str, Path]:
    """Regrid temperature and convert retained lower-resolution maps."""
    from FITS_tools.hcongrid import hcongrid

    products = prepared_paths(paths)
    temperature, temperature_header = fits.getdata(
        paths["temperature_lowres"],
        header=True,
        memmap=True,
    )
    target_header = fits.getheader(paths["dust_continuum_allscale"])
    regridded = hcongrid(temperature, temperature_header, target_header)
    output_header = target_header.copy()
    for key in ("BUNIT", "BMAJ", "BMIN", "BPA"):
        if key in temperature_header:
            output_header[key] = temperature_header[key]
    _write_float32(
        products["temperature_allscale"],
        regridded,
        _float_header(
            output_header,
            unit="K",
            history="LANCET column_density: low-resolution temperature regrid",
        ),
        overwrite=overwrite,
    )
    del regridded

    model = config["model"]
    assert isinstance(model, Mapping)
    molecular_weight = float(model["mean_molecular_weight"])
    for input_name, output_name in (
        ("column_density_lowres", "surface_density_lowres"),
        ("column_density_multiscale", "surface_density_multiscale"),
    ):
        column, header = fits.getdata(
            paths[input_name],
            header=True,
            memmap=True,
        )
        surface = column_density_to_surface_density(
            column,
            mean_molecular_weight=molecular_weight,
        )
        _write_float32(
            products[output_name],
            surface,
            _float_header(
                header,
                unit="g cm^-2",
                history=(
                    "LANCET column_density: H2 column to gas surface density"
                ),
            ),
            overwrite=overwrite,
        )
        del surface
    return products


def build_maps(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> dict[str, Path]:
    """Create all-scale and intercore column/surface-density maps."""
    temperature_path = prepared_paths(paths)["temperature_allscale"]
    if not temperature_path.is_file():
        raise FileNotFoundError(
            f"Prepared temperature not found: {temperature_path}; "
            "run --stage prepare first"
        )
    model = config["model"]
    assert isinstance(model, Mapping)
    products = map_paths(paths)

    with fits.open(temperature_path, memmap=True) as temperature_hdul:
        temperature = temperature_hdul[0].data
        for input_name, column_name, surface_name in (
            (
                "dust_continuum_allscale",
                "column_density_allscale",
                "surface_density_allscale",
            ),
            (
                "intercore_background",
                "column_density_intercore",
                "surface_density_intercore",
            ),
        ):
            intensity, header = fits.getdata(
                paths[input_name],
                header=True,
                memmap=True,
            )
            if intensity.shape != temperature.shape:
                raise ValueError(
                    f"{input_name} shape {intensity.shape} does not match "
                    f"temperature shape {temperature.shape}"
                )
            column = brightness_to_column_density(
                intensity,
                temperature,
                frequency_ghz=float(model["frequency_ghz"]),
                opacity_cm2_per_g=float(model["opacity_cm2_per_g"]),
                mean_molecular_weight=float(
                    model["mean_molecular_weight"]
                ),
            )
            _write_float32(
                products[column_name],
                column,
                _float_header(
                    header,
                    unit="cm^-2",
                    history=(
                        "LANCET column_density: finite-opacity brightness "
                        "inversion"
                    ),
                ),
                overwrite=overwrite,
            )
            surface = column_density_to_surface_density(
                column,
                mean_molecular_weight=float(
                    model["mean_molecular_weight"]
                ),
            )
            _write_float32(
                products[surface_name],
                surface,
                _float_header(
                    header,
                    unit="g cm^-2",
                    history=(
                        "LANCET column_density: H2 column to gas surface "
                        "density"
                    ),
                ),
                overwrite=overwrite,
            )
            del intensity, column, surface
    return products


def _figure_outputs(
    stem: Path,
    formats: Sequence[object],
    *,
    overwrite: bool,
) -> list[Path]:
    outputs = [stem.with_suffix(f".{str(value).lower()}") for value in formats]
    for output in outputs:
        if output.suffix not in {".pdf", ".png"}:
            raise ValueError(f"Unsupported figure format: {output.suffix}")
        _prepare_output(output, overwrite=overwrite)
    return outputs


def _save(figure, outputs: Sequence[Path], *, dpi: int) -> None:
    for output in outputs:
        figure.savefig(output, bbox_inches="tight", dpi=dpi)
        print(f"Wrote {output}")


def plot_products(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> list[Path]:
    """Create the full-field and subregion surface-density figures."""
    import matplotlib

    matplotlib.use("Agg")
    from astropy import units as u
    from astropy.visualization.wcsaxes import add_beam, add_scalebar
    from matplotlib import pyplot as plt
    import matplotlib.patheffects as path_effects
    from regions import Regions

    surface_path = map_paths(paths)["surface_density_allscale"]
    if not surface_path.is_file():
        raise FileNotFoundError(
            f"Surface-density map not found: {surface_path}; "
            "run --stage maps first"
        )
    view = config["view"]
    style = config["plot_style"]
    figures = config["figures"]
    subregions = config["subregions"]
    if not all(
        isinstance(item, Mapping)
        for item in (view, style, figures, subregions)
    ):
        raise TypeError("Plot settings must be TOML tables")
    data, header = fits.getdata(surface_path, header=True, memmap=True)
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
        annotated=True,
        title="H$_2$ surface density map",
        resolution_label=str(view["resolution_label"]),
        scalebar_length_deg=float(view["scalebar_length_deg"]),
        scalebar_label=str(view["scalebar_label"]),
        colourbar_label=r"Surface density (g cm$^{-2}$)",
        background_ticks=tuple(
            float(value) for value in style["background_ticks"]
        ),
    )
    output_root = paths["output_root"] / "figures"
    output_formats = figures["formats"]
    dpi = int(figures["full_field_dpi"])
    written: list[Path] = []

    base_outputs = _figure_outputs(
        output_root / "surface_density_allscale",
        output_formats,
        overwrite=overwrite,
    )
    _save(figure, base_outputs, dpi=dpi)
    written.extend(base_outputs)

    overlay_axis = figure.axes[1]
    freefree, freefree_header = fits.getdata(
        paths["freefree_smoothed"],
        header=True,
        memmap=True,
    )
    contours = overlay_axis.contour(
        np.squeeze(freefree),
        levels=[float(value) for value in style["freefree_levels"]],
        colors="orange",
        linewidths=1.0,
        origin="lower",
        transform=overlay_axis.get_transform(
            WCS(freefree_header).celestial
        ),
    )
    overlay_outputs = _figure_outputs(
        output_root / "surface_density_allscale_freefree",
        output_formats,
        overwrite=overwrite,
    )
    _save(figure, overlay_outputs, dpi=dpi)
    written.extend(overlay_outputs)
    contours.remove()

    for region in Regions.read(paths["box_regions"], format="ds9"):
        overlay_axis.add_patch(
            region.to_pixel(overlay_axis.wcs).as_artist(
                facecolor="none",
                edgecolor="yellow",
                linewidth=1.5,
            )
        )
    box_outputs = _figure_outputs(
        output_root / "surface_density_allscale_subregion_boxes",
        output_formats,
        overwrite=overwrite,
    )
    _save(figure, box_outputs, dpi=dpi)
    written.extend(box_outputs)
    plt.close(figure)

    distance_pc = float(config["model"]["distance_pc"])
    scale_length_pc = float(view["subregion_scalebar_pc"])
    scale_length_deg = np.degrees(scale_length_pc / distance_pc)
    for name, settings in subregions.items():
        assert isinstance(settings, Mapping)
        subfigure = create_dual_colour_figure(
            data,
            header,
            center_ra_deg=float(settings["center_ra_deg"]),
            center_dec_deg=float(settings["center_dec_deg"]),
            width_deg=float(view["subregion_width_deg"]),
            height_deg=float(view["subregion_height_deg"]),
            background_vmin=float(style["background_vmin"]),
            transition=float(style["transition"]),
            highlight_vmax=float(style["highlight_vmax"]),
            figsize=(5.0, 5.0),
            annotated=False,
            image_axes_rect=(0.05, 0.05, 0.90, 0.90),
        )
        axis = subfigure.axes[1]
        add_beam(
            axis,
            header,
            edgecolor="black",
            facecolor="white",
            frame=True,
        )
        add_scalebar(
            axis,
            length=scale_length_deg * u.deg,
            label=f"{scale_length_pc:g} pc",
            color="black",
        )
        axis.text(
            0.07,
            0.05,
            str(view["resolution_label"]),
            transform=axis.transAxes,
            color="black",
            size=12,
            path_effects=[
                path_effects.withStroke(
                    linewidth=2.5,
                    foreground="white",
                )
            ],
        )
        axis.text(
            0.02,
            0.98,
            str(settings["label"]),
            transform=axis.transAxes,
            horizontalalignment="left",
            verticalalignment="top",
            color="black",
            size=20,
            weight="bold",
            path_effects=[
                path_effects.withStroke(
                    linewidth=2.5,
                    foreground="white",
                )
            ],
        )
        outputs = _figure_outputs(
            output_root / f"surface_density_{name}",
            output_formats,
            overwrite=overwrite,
        )
        _save(subfigure, outputs, dpi=int(figures["subregion_dpi"]))
        written.extend(outputs)
        plt.close(subfigure)
    return written


def measure_dense_gas(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> Path:
    """Measure dense mass and dense-gas fraction by evolutionary subregion."""
    import pyregion

    surface_path = map_paths(paths)["surface_density_allscale"]
    surface, header = fits.getdata(
        surface_path,
        header=True,
        memmap=True,
    )
    scales = np.abs(proj_plane_pixel_scales(WCS(header).celestial)) * 3600.0
    dense = config["dense_gas"]
    regions = config["dense_regions"]
    if not isinstance(dense, Mapping) or not isinstance(regions, Mapping):
        raise TypeError("Dense-gas settings must be TOML tables")

    summary: dict[str, object] = {
        "threshold_g_cm2": float(dense["threshold_g_cm2"]),
        "distance_pc": float(config["model"]["distance_pc"]),
        "pixel_scale_arcsec": [float(scales[0]), float(scales[1])],
        "regions": {},
    }
    region_results = summary["regions"]
    assert isinstance(region_results, dict)
    shape = surface.shape
    for name, region_path in regions.items():
        mask = pyregion.open(str(paths[str(region_path)])).as_imagecoord(
            header=header
        ).get_mask(shape=shape)
        cloud_mass = float(dense[f"cloud_mass_{name}_msun"])
        mass = dense_mass_msun(
            surface,
            mask,
            pixel_scale_x_arcsec=float(scales[0]),
            pixel_scale_y_arcsec=float(scales[1]),
            distance_pc=float(config["model"]["distance_pc"]),
            threshold_g_cm2=float(dense["threshold_g_cm2"]),
        )
        region_results[str(name)] = {
            "dense_mass_msun": mass,
            "cloud_mass_msun": cloud_mass,
            "dense_gas_fraction": mass / cloud_mass,
        }
        print(
            f"{name}: dense mass={mass:.3f} Msun, "
            f"fraction={mass / cloud_mass:.6f}"
        )

    output = paths["output_root"] / "measurements" / "dense_gas.json"
    _prepare_output(output, overwrite=overwrite)
    output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {output}")
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--stage",
        choices=("inspect", "prepare", "maps", "plot", "dense", "all"),
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
    if args.stage in ("prepare", "all"):
        prepare_products(paths, config, overwrite=args.overwrite)
    if args.stage in ("maps", "all"):
        build_maps(paths, config, overwrite=args.overwrite)
    if args.stage in ("plot", "all"):
        plot_products(paths, config, overwrite=args.overwrite)
    if args.stage in ("dense", "all"):
        measure_dense_gas(paths, config, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
