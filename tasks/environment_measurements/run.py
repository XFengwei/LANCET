#!/usr/bin/env python
"""Measure temperature, kinematic, and surface-density core environments."""

from __future__ import annotations

import argparse
from collections.abc import Mapping
from pathlib import Path
import tomllib

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from photutils.aperture import SkyCircularAperture
import numpy as np

from lancet import ProjectPaths
from lancet.environment import (
    add_environment_measurements,
    measure_aperture_statistics,
)


DEFAULT_CONFIG = Path("tasks/environment_measurements/task.toml")


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
    """Resolve declared inputs and the policy-safe output root."""
    inputs = config["inputs"]
    outputs = config["outputs"]
    if not isinstance(inputs, Mapping) or not isinstance(outputs, Mapping):
        raise TypeError("Inputs and outputs must be TOML tables")
    paths = {
        name: _project_path(project, value)
        for name, value in inputs.items()
    }
    paths["output_root"] = project.validate_output(str(outputs["root"]))
    if not paths["temperature"].is_relative_to(
        project.delivered.resolve()
    ):
        raise ValueError(
            f"temperature must be delivered: {paths['temperature']}"
        )
    for name in inputs:
        if not paths[name].is_file():
            raise FileNotFoundError(f"{name} not found: {paths[name]}")
    return paths


def _fits_description(path: Path) -> str:
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
    """Report catalogue, map, and aperture readiness."""
    catalogue = Table.read(
        paths["physical_catalogue"],
        format="ascii.ecsv",
    )
    aperture = config["aperture"]
    if not isinstance(aperture, Mapping):
        raise TypeError("Aperture must be a TOML table")
    print("Environment-measurement inputs:")
    print(
        f"  physical_catalogue: {paths['physical_catalogue']} "
        f"(rows={len(catalogue)})"
    )
    for name in (
        "temperature",
        "velocity_dispersion",
        "intercore_surface_density",
        "surface_density",
    ):
        print(f"  {name}: {paths[name]}")
        print(f"    {_fits_description(paths[name])}")
    print(
        f"Aperture: radius={aperture['radius_arcsec']} arcsec; "
        f"valid pixels > {aperture['minimum_valid_value']}"
    )
    print(f"Output root: {paths['output_root']}")


def output_path(paths: Mapping[str, Path]) -> Path:
    """Return the generated environment-catalogue path."""
    return (
        paths["output_root"]
        / "tables"
        / "core_environment_catalogue.ecsv"
    )


def _prepare_output(path: Path, *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists: {path}; "
            "pass --overwrite to replace it"
        )
    path.parent.mkdir(parents=True, exist_ok=True)


def build_catalogue(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> Table:
    """Measure every configured map and write the augmented catalogue."""
    catalogue = Table.read(
        paths["physical_catalogue"],
        format="ascii.ecsv",
    )
    required = {"WCS_ACOOR", "WCS_DCOOR"}
    missing = required - set(catalogue.colnames)
    if missing:
        raise KeyError(
            "Physical catalogue is missing coordinates: "
            + ", ".join(sorted(missing))
        )
    aperture_config = config["aperture"]
    maps = config["maps"]
    if not isinstance(aperture_config, Mapping) or not isinstance(
        maps,
        Mapping,
    ):
        raise TypeError("Aperture and maps must be TOML tables")

    coordinates = SkyCoord(
        np.asarray(catalogue["WCS_ACOOR"], dtype=float) * u.deg,
        np.asarray(catalogue["WCS_DCOOR"], dtype=float) * u.deg,
        frame="icrs",
    )
    apertures = SkyCircularAperture(
        coordinates,
        r=float(aperture_config["radius_arcsec"]) * u.arcsec,
    )
    minimum = float(aperture_config["minimum_valid_value"])
    measurements: dict[str, dict[str, np.ndarray]] = {}
    units: dict[str, str | None] = {}
    for prefix, input_name in maps.items():
        name = str(input_name)
        if name not in paths:
            raise KeyError(f"Map input {name!r} is not declared")
        data, header = fits.getdata(
            paths[name],
            header=True,
            memmap=True,
        )
        measurements[str(prefix)] = measure_aperture_statistics(
            data,
            apertures,
            WCS(header).celestial,
            minimum_valid_value=minimum,
        )
        units[str(prefix)] = header.get("BUNIT")

    output = add_environment_measurements(
        catalogue,
        measurements,
        units=units,
    )
    output.meta["aperture_radius_arcsec"] = float(
        aperture_config["radius_arcsec"]
    )
    output.meta["minimum_valid_value"] = minimum
    destination = output_path(paths)
    _prepare_output(destination, overwrite=overwrite)
    output.write(destination, format="ascii.ecsv", overwrite=True)
    print(f"Wrote {destination} ({len(output)} cores)")
    for prefix in maps:
        valid = np.asarray(
            output[f"{prefix}_env_validfrac"],
            dtype=float,
        )
        print(
            f"  {prefix}: median valid fraction="
            f"{np.nanmedian(valid):.3f}, minimum={np.nanmin(valid):.3f}"
        )
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--stage",
        choices=("inspect", "build"),
        default="inspect",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    project = ProjectPaths.discover()
    config = load_task_config(_project_path(project, args.config))
    paths = resolve_task_paths(project, config)
    if args.stage == "inspect":
        inspect_inputs(paths, config)
    else:
        build_catalogue(paths, config, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
