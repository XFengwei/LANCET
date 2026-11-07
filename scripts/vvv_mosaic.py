#!/usr/bin/env python
"""Mosaic the delivered VVV J/H/Ks tiles of G316.8 and render a colour image."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from pathlib import Path
import tomllib

from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import numpy as np
from reproject.mosaicking import find_optimal_celestial_wcs

from lancet import ProjectPaths
from lancet.infrared_mosaic import (
    boost_saturation,
    coadd_tiles,
    load_tile,
    stretch_channel,
)


DEFAULT_CONFIG = Path("tasks/vvv_mosaic/task.toml")


def load_task_config(path: Path) -> dict[str, object]:
    """Load task settings from TOML."""
    if not path.is_file():
        raise FileNotFoundError(f"Task configuration not found: {path}")
    with path.open("rb") as handle:
        return tomllib.load(handle)


def tile_paths(
    project: ProjectPaths,
    config: Mapping[str, object],
) -> dict[str, list[Path]]:
    """Return the delivered tile files grouped by band."""
    inputs = config["inputs"]
    directory = (project.root / str(inputs["tile_directory"])).resolve()
    if not directory.is_relative_to(project.delivered.resolve()):
        raise ValueError(f"tile_directory must be a delivered input: {directory}")
    paths: dict[str, list[Path]] = {}
    for band in inputs["bands"]:
        paths[band] = [directory / f"{tile}_{band}.fits" for tile in inputs["tiles"]]
        for path in paths[band]:
            if not path.is_file():
                raise FileNotFoundError(f"VVV tile not found: {path}")
    return paths



def mosaic_header(
    tiles: Sequence[tuple[np.ndarray, WCS]],
    grid: Mapping[str, object],
) -> tuple[WCS, tuple[int, int]]:
    """Return the configured north-up grid, or one covering every tile."""
    frame = str(grid["frame"])
    pixel_scale = float(grid["pixel_scale_arcsec"]) / 3600.0
    if "center_ra_deg" not in grid:
        return find_optimal_celestial_wcs(
            list(tiles),
            frame=frame,
            projection="TAN",
            resolution=pixel_scale * u.deg,
        )
    if frame != "icrs":
        raise ValueError("An explicit grid centre requires frame = 'icrs'")
    width = int(grid["width_pixels"])
    height = int(grid["height_pixels"])
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [float(grid["center_ra_deg"]), float(grid["center_dec_deg"])]
    wcs.wcs.crpix = [(width + 1) / 2.0, (height + 1) / 2.0]
    wcs.wcs.cdelt = [-pixel_scale, pixel_scale]
    wcs.wcs.radesys = "ICRS"
    return wcs, (height, width)


def build_band_mosaic(
    paths: Sequence[Path],
    config: Mapping[str, object],
    target: tuple[WCS, tuple[int, int]] | None = None,
) -> tuple[np.ndarray, WCS, tuple[int, int]]:
    """Reproject the tiles of one band, flatten their seams and feather them."""
    zero_point = float(config["photometry"]["reference_zero_point"])
    grid = config["grid"]
    blending = config["blending"]
    edge_trim = int(blending["edge_trim_pixels"])
    tiles = [load_tile(path, zero_point, edge_trim) for path in paths]
    if target is None:
        target = mosaic_header(tiles, grid)
    wcs_out, shape_out = target
    pixel_arcsec = float(grid["pixel_scale_arcsec"])
    mosaic = coadd_tiles(
        tiles,
        wcs_out,
        shape_out,
        block=int(blending["block_pixels"]),
        smoothness=float(blending["smoothness"]),
        decay_pixels=float(blending["decay_arcsec"]) / pixel_arcsec,
        feather_pixels=float(blending["feather_arcsec"]) / pixel_arcsec,
    )
    return mosaic, wcs_out, shape_out


def run_mosaic(project: ProjectPaths, config: Mapping[str, object]) -> None:
    """Write one coadded FITS mosaic per band."""
    output_root = project.validate_output(str(config["outputs"]["root"]))
    output_root.mkdir(parents=True, exist_ok=True)
    target = None
    for band, paths in tile_paths(project, config).items():
        mosaic, wcs_out, shape_out = build_band_mosaic(paths, config, target)
        target = (wcs_out, shape_out)
        header = wcs_out.to_header()
        header["BUNIT"] = ("ADU", "sky-subtracted, rescaled counts")
        header["MAGZPT"] = (
            float(config["photometry"]["reference_zero_point"]),
            "common zero point of the mosaic",
        )
        header["FILTER"] = band.upper()
        output = output_root / f"vvv_{band}_mosaic.fits"
        fits.PrimaryHDU(mosaic, header=header).writeto(output, overwrite=True)
        print(f"Wrote {output} shape={mosaic.shape}")


def run_rgb(project: ProjectPaths, config: Mapping[str, object]) -> None:
    """Combine the band mosaics into a Ks/H/J colour PNG."""
    from PIL import Image

    output_root = project.validate_output(str(config["outputs"]["root"]))
    rgb_config = config["rgb"]
    channels = []
    for name, weight in zip(("red", "green", "blue"), rgb_config["weights"]):
        band = rgb_config[name]
        data = fits.getdata(output_root / f"vvv_{band}_mosaic.fits")
        channel = stretch_channel(
            data,
            float(rgb_config["high_percentile"]),
            float(rgb_config["asinh_softening"]),
            float(rgb_config["black_sigma"]),
        )
        channels.append(np.clip(channel * float(weight), 0.0, 1.0))

    # FITS rows run bottom-up; images run top-down.
    rgb = np.flipud(np.stack(channels, axis=-1))
    rgb = boost_saturation(rgb, float(rgb_config["saturation"]))
    image = Image.fromarray((rgb * 255.0 + 0.5).astype(np.uint8))

    figures = config["figures"]
    output = output_root / str(figures["png"])
    image.save(output)
    print(f"Wrote {output} size={image.size}")

    preview = image.copy()
    preview.thumbnail((int(figures["preview_max_pixels"]),) * 2)
    preview_path = output.with_name(output.stem + "_preview.png")
    preview.save(preview_path)
    print(f"Wrote {preview_path} size={preview.size}")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--stage",
        choices=("mosaic", "rgb", "all"),
        default="all",
        help="mosaic: coadd FITS per band; rgb: render the colour PNG",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    project = ProjectPaths.discover()
    config_path = args.config if args.config.is_absolute() else project.root / args.config
    config = load_task_config(config_path)
    if args.stage in ("mosaic", "all"):
        run_mosaic(project, config)
    if args.stage in ("rgb", "all"):
        run_rgb(project, config)


if __name__ == "__main__":
    main()
