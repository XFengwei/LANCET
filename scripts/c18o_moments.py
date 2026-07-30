#!/usr/bin/env python
"""Reconstruct C18O cubes, masked moments, and ALMA-grid products."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from pathlib import Path
import tomllib

from astropy import units as u
from astropy.io import fits
import numpy as np
from spectral_cube import SpectralCube

from lancet import ProjectPaths
from lancet.moments import reproject_image


DEFAULT_CONFIG = Path("tasks/c18o_moments/task.toml")
MOMENT_SUFFIXES = (
    "ecube.fits.gz",
    "mask.fits.gz",
    "flux.csv",
    "mom0.fits.gz",
    "mom1.fits.gz",
    "mom2.fits.gz",
    "emom0.fits.gz",
    "emom1.fits.gz",
    "emom2.fits.gz",
)


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
    """Resolve the single line-data input, WCS template, and output root."""
    inputs = config["inputs"]
    outputs = config["outputs"]
    if not isinstance(inputs, Mapping) or not isinstance(outputs, Mapping):
        raise TypeError("Task inputs and outputs must be TOML tables")
    paths = {
        name: _project_path(project, value)
        for name, value in inputs.items()
    }
    paths["output_root"] = project.validate_output(str(outputs["root"]))
    for name in ("c18o_cube", "alma_grid"):
        path = paths[name]
        if not path.is_relative_to(project.delivered.resolve()):
            raise ValueError(f"{name} must be a delivered input: {path}")
        if not path.is_file():
            raise FileNotFoundError(f"{name} not found: {path}")
    return paths


def cube_paths(paths: Mapping[str, Path]) -> dict[str, Path]:
    root = paths["output_root"] / "cubes"
    return {
        "velocity_kms": root / "c18o_jybeam_velocity_kms.fits",
        "brightness_k": root / "c18o_brightness_temperature.fits",
    }


def moment_paths(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
) -> dict[str, Path]:
    settings = config["moments"]
    assert isinstance(settings, Mapping)
    basename = str(settings["basename"])
    root = paths["output_root"] / "moments"
    return {
        suffix: root / f"{basename}.{suffix}"
        for suffix in MOMENT_SUFFIXES
    }


def _describe(path: Path) -> str:
    with fits.open(path, memmap=True) as hdul:
        header = hdul[0].header
        return (
            f"shape={tuple(hdul[0].shape)}, "
            f"BITPIX={int(header['BITPIX'])}, "
            f"BUNIT={header.get('BUNIT', 'not recorded')}, "
            f"CUNIT3={header.get('CUNIT3', 'n/a')}"
        )


def inspect_inputs(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
) -> None:
    """Report input metadata and declared mask settings."""
    print(f"C18O cube: {paths['c18o_cube']}")
    print(f"  {_describe(paths['c18o_cube'])}")
    print(f"ALMA grid: {paths['alma_grid']}")
    print(f"  {_describe(paths['alma_grid'])}")
    settings = config["moments"]
    assert isinstance(settings, Mapping)
    print(
        "Mask settings: "
        f"snr_hi={settings['snr_hi']}, "
        f"snr_lo={settings['snr_lo']}, "
        f"min_tot_ch={settings['min_tot_ch']}, "
        f"edgech={settings['edge_channels']}"
    )
    print(f"Output root: {paths['output_root']}")


def _prepare_output(path: Path, *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists: {path}; pass --overwrite to replace it"
        )
    path.parent.mkdir(parents=True, exist_ok=True)


def _write_cube_float32(
    cube: SpectralCube,
    path: Path,
    *,
    overwrite: bool,
    history: str,
) -> None:
    _prepare_output(path, overwrite=overwrite)
    header = cube.header.copy()
    if "BLANK" in header:
        del header["BLANK"]
    header.add_history(history)
    fits.writeto(
        path,
        np.asarray(cube.unmasked_data[:].value, dtype=np.float32),
        header,
        overwrite=overwrite,
    )
    print(f"Wrote {path}")


def prepare_cubes(
    paths: Mapping[str, Path],
    *,
    overwrite: bool = False,
) -> dict[str, Path]:
    """Create the velocity-axis and brightness-temperature derivatives."""
    products = cube_paths(paths)
    cube = SpectralCube.read(paths["c18o_cube"])
    velocity_cube = cube.with_spectral_unit(u.km / u.s)
    _write_cube_float32(
        velocity_cube,
        products["velocity_kms"],
        overwrite=overwrite,
        history="LANCET c18o_moments: spectral axis converted to km/s",
    )
    temperature_cube = cube.to(u.K)
    _write_cube_float32(
        temperature_cube,
        products["brightness_k"],
        overwrite=overwrite,
        history="LANCET c18o_moments: Jy/beam converted to brightness K",
    )
    return products


def build_moments(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> dict[str, Path]:
    """Run the declared dilated-mask moment calculation."""
    from maskmoment import maskmoment

    input_cube = cube_paths(paths)["brightness_k"]
    if not input_cube.is_file():
        raise FileNotFoundError(
            f"Brightness-temperature cube not found: {input_cube}; "
            "run --stage prepare first"
        )
    products = moment_paths(paths, config)
    for product in products.values():
        _prepare_output(product, overwrite=overwrite)

    settings = config["moments"]
    assert isinstance(settings, Mapping)
    maskmoment(
        img_fits=str(input_cube),
        outdir=str(paths["output_root"] / "moments"),
        outname=str(settings["basename"]),
        snr_hi=float(settings["snr_hi"]),
        snr_lo=float(settings["snr_lo"]),
        minbeam=float(settings["min_beam"]),
        snr_hi_minch=int(settings["snr_hi_min_channels"]),
        snr_lo_minch=int(settings["snr_lo_min_channels"]),
        min_tot_ch=int(settings["min_tot_ch"]),
        edgech=int(settings["edge_channels"]),
        mom1_minch=int(settings["moment1_min_channels"]),
        mom2_minch=int(settings["moment2_min_channels"]),
    )
    return products


def reproject_moment2(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> dict[str, Path]:
    """Reproject velocity dispersion and its footprint to the ALMA grid."""
    source = moment_paths(paths, config)["mom2.fits.gz"]
    if not source.is_file():
        raise FileNotFoundError(
            f"Moment-2 map not found: {source}; run --stage moments first"
        )
    data, source_header = fits.getdata(source, header=True, memmap=True)
    target_header = fits.getheader(paths["alma_grid"])
    reprojected, footprint, output_header = reproject_image(
        data,
        source_header,
        target_header,
        order=str(config["reprojection"]["order"]),
    )
    outputs = {
        "moment2": (
            paths["output_root"]
            / "reprojected"
            / "c18o_moment2_on_alma_grid.fits"
        ),
        "footprint": (
            paths["output_root"]
            / "reprojected"
            / "c18o_moment2_reprojection_footprint.fits"
        ),
    }
    for name, array in (
        ("moment2", reprojected),
        ("footprint", footprint),
    ):
        _prepare_output(outputs[name], overwrite=overwrite)
        header = output_header.copy()
        header["BUNIT"] = (
            source_header.get("BUNIT", "km s-1")
            if name == "moment2"
            else "1"
        )
        header["BTYPE"] = (
            "Velocity dispersion"
            if name == "moment2"
            else "Reprojection footprint"
        )
        fits.writeto(
            outputs[name],
            np.asarray(array, dtype=np.float32),
            header,
            overwrite=overwrite,
        )
        print(f"Wrote {outputs[name]}")
    return outputs


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


def plot_moments(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> list[Path]:
    """Create the stacked moment-0, moment-1, and moment-2 figure."""
    import matplotlib

    matplotlib.use("Agg")
    from astropy.visualization.wcsaxes import add_beam
    from astropy.wcs import WCS
    from matplotlib import pyplot as plt
    import seaborn as sns

    products = moment_paths(paths, config)
    for suffix in ("mom0.fits.gz", "mom1.fits.gz", "mom2.fits.gz"):
        if not products[suffix].is_file():
            raise FileNotFoundError(
                f"Moment product not found: {products[suffix]}"
            )
    figure = plt.figure(figsize=(5.5, 12.0), dpi=100)
    grid = figure.add_gridspec(
        3,
        2,
        width_ratios=[1.0, 0.05],
        hspace=0.02,
        wspace=0.0,
    )
    panels = (
        (
            "mom0.fits.gz",
            "Moment-0",
            sns.color_palette("mako", as_cmap=True),
            "Line brightness (K km s$^{-1}$)",
            {},
        ),
        (
            "mom1.fits.gz",
            "Moment-1",
            "RdBu_r",
            "Velocity (km s$^{-1}$)",
            {},
        ),
        (
            "mom2.fits.gz",
            "Moment-2",
            "inferno",
            "Velocity dispersion (km s$^{-1}$)",
            {"vmin": 0.5, "vmax": 3.0},
        ),
    )
    for row, (suffix, label, colormap, colourbar_label, image_options) in (
        enumerate(panels)
    ):
        data, header = fits.getdata(
            products[suffix],
            header=True,
            memmap=True,
        )
        axis = figure.add_subplot(
            grid[row, 0],
            projection=WCS(header).celestial,
        )
        artist = axis.imshow(data, cmap=colormap, **image_options)
        axis.tick_params(
            axis="x",
            bottom=True,
            top=True,
            labelbottom=row == 2,
        )
        axis.tick_params(
            axis="y",
            left=True,
            right=True,
            labelleft=row == 2,
        )
        axis.set_xlabel("RA (ICRS)")
        axis.set_ylabel("DEC (ICRS)")
        axis.text(
            0.02,
            0.98,
            label,
            transform=axis.transAxes,
            horizontalalignment="left",
            verticalalignment="top",
            size=10,
            weight="bold",
        )
        axis.set_facecolor("0.8")
        add_beam(
            axis,
            header,
            borderpad=0.2,
            color="black",
            frame=True,
            corner="bottom right",
        )
        if row == 2:
            axis.text(
                0.04,
                0.03,
                '30" = 0.4 pc',
                transform=axis.transAxes,
                horizontalalignment="left",
                verticalalignment="bottom",
                size=10,
                weight="bold",
                bbox={
                    "facecolor": "white",
                    "edgecolor": "none",
                    "alpha": 0.75,
                    "pad": 1.5,
                },
            )
        colourbar_axis = figure.add_subplot(grid[row, 1])
        colourbar = figure.colorbar(artist, cax=colourbar_axis)
        colourbar.set_label(colourbar_label, fontsize=10)
        colourbar.ax.tick_params(labelsize=10, length=2)

    figure_settings = config["figures"]
    assert isinstance(figure_settings, Mapping)
    outputs = _figure_outputs(
        paths["output_root"] / "figures" / "c18o_moment_maps",
        figure_settings["formats"],
        overwrite=overwrite,
    )
    try:
        for output in outputs:
            figure.savefig(
                output,
                bbox_inches="tight",
                dpi=int(figure_settings["dpi"]),
            )
            print(f"Wrote {output}")
    finally:
        plt.close(figure)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--stage",
        choices=(
            "inspect",
            "prepare",
            "moments",
            "reproject",
            "plot",
            "all",
        ),
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
        prepare_cubes(paths, overwrite=args.overwrite)
    if args.stage in ("moments", "all"):
        build_moments(paths, config, overwrite=args.overwrite)
    if args.stage in ("reproject", "all"):
        reproject_moment2(paths, config, overwrite=args.overwrite)
    if args.stage in ("plot", "all"):
        plot_moments(paths, config, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
