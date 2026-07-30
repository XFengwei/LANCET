#!/usr/bin/env python
"""Calculate physical properties for the project source catalogue."""

from __future__ import annotations

import argparse
from collections.abc import Mapping
from pathlib import Path
import tomllib

from astropy import constants as const
from astropy import units as u
from astropy.table import Table
import numpy as np

from lancet import ProjectPaths
from lancet.physical_properties import (
    brightness_to_column_density,
    core_optical_depth,
    core_surface_density,
    core_volume_density,
    propagate_uncertainty,
    thick_mass,
)


DEFAULT_CONFIG = Path("tasks/physical_properties/task.toml")


def load_task_config(path: Path) -> dict[str, object]:
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
    inputs = config["inputs"]
    outputs = config["outputs"]
    if not isinstance(inputs, Mapping) or not isinstance(outputs, Mapping):
        raise TypeError("Inputs and outputs must be TOML tables")
    catalogue = _project_path(project, inputs["source_catalogue"])
    if not catalogue.is_file():
        raise FileNotFoundError(
            f"Source catalogue not found: {catalogue}; "
            "run source_catalogue --stage build first"
        )
    return {
        "source_catalogue": catalogue,
        "output_root": project.validate_output(str(outputs["root"])),
    }


def inspect_inputs(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
) -> None:
    catalogue = Table.read(paths["source_catalogue"], format="ascii.ecsv")
    physics = config["physics"]
    uncertainty = config["uncertainty"]
    assert isinstance(physics, Mapping)
    assert isinstance(uncertainty, Mapping)
    print(f"Source catalogue: {paths['source_catalogue']}")
    print(f"  rows={len(catalogue)}")
    print(
        "Physics: "
        f"frequency={physics['frequency_ghz']} GHz, "
        f"kappa_dust={physics['dust_opacity_cm2_per_g']} cm2/g, "
        f"gas_to_dust={physics['gas_to_dust_ratio']}, "
        f"distance={physics['distance_kpc']} kpc"
    )
    print(
        f"Uncertainty: samples={uncertainty['sample_count']}, "
        f"seed={uncertainty['random_seed']}"
    )
    print(f"Output root: {paths['output_root']}")


def _prepare_output(path: Path, *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists: {path}; pass --overwrite to replace it"
        )
    path.parent.mkdir(parents=True, exist_ok=True)


def output_paths(paths: Mapping[str, Path]) -> dict[str, Path]:
    root = paths["output_root"] / "tables"
    return {
        "catalogue": root / "core_physical_properties.ecsv",
        "publication": root / "core_physical_properties.tex",
    }


def _monte_carlo_series(
    function,
    parameter_factory,
    count: int,
    *,
    sample_count: int,
    seed: int,
    transform=None,
) -> tuple[np.ndarray, np.ndarray]:
    means = []
    errors = []
    for index in range(count):
        result = propagate_uncertainty(
            function,
            parameter_factory(index),
            sample_count=sample_count,
            random_seed=seed,
        )
        mean = result["mean"]
        error = result["std"]
        if transform is not None:
            mean = transform(mean)
            error = transform(error)
        means.append(mean.value if isinstance(mean, u.Quantity) else mean)
        errors.append(error.value if isinstance(error, u.Quantity) else error)
    return np.asarray(means), np.asarray(errors)


def _round_uncertainty(value: float, error: float) -> str:
    if not np.isfinite(value) or not np.isfinite(error) or error <= 0:
        return f"{value:.2g}"
    exponent = int(np.floor(np.log10(error)))
    significant_digits = 2 if error / 10**exponent < 3 else 1
    round_to = exponent - (significant_digits - 1)
    rounded_error = round(error, -round_to)
    rounded_value = round(value, -round_to)
    if round_to >= 0:
        template = "%.0f"
    else:
        template = f"%.{-round_to}f"
    return (
        template % rounded_value
        + r"$\pm$"
        + template % rounded_error
    )


def build_properties(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> Table:
    """Calculate and write the physical-property catalogue."""
    catalogue = Table.read(paths["source_catalogue"], format="ascii.ecsv")
    physics = config["physics"]
    uncertainty = config["uncertainty"]
    assert isinstance(physics, Mapping)
    assert isinstance(uncertainty, Mapping)
    count = len(catalogue)
    samples = int(uncertainty["sample_count"])
    seed = int(uncertainty["random_seed"])

    gas_to_dust = float(physics["gas_to_dust_ratio"])
    gas_to_dust_error = float(physics["gas_to_dust_error"])
    dust_opacity = float(physics["dust_opacity_cm2_per_g"]) * (
        u.cm**2 / u.g
    )
    kappa = dust_opacity / gas_to_dust
    kappa_error = (
        dust_opacity / gas_to_dust**2 * gas_to_dust_error
    )
    distance = float(physics["distance_kpc"]) * u.kpc
    distance_error = float(physics["distance_error_kpc"]) * u.kpc
    frequency = float(physics["frequency_ghz"]) * u.GHz
    background_error = float(
        uncertainty["background_rms_mjysr"]
    ) * u.MJy / u.sr
    temperature = np.asarray(catalogue["Tenv"]) * u.K
    temperature_error = (
        temperature * float(uncertainty["temperature_fraction"])
    )
    intercore = np.asarray(catalogue["Ibkg1"]) * u.MJy / u.sr
    background = np.asarray(catalogue["Ibkg2"]) * u.MJy / u.sr
    total_flux = np.asarray(catalogue["FXT_BST01"]) * u.Jy
    total_flux_error = np.asarray(catalogue["FXT_ERR01"]) * u.Jy
    beam_area = (
        1.1331
        * float(physics["beam_major_arcsec"])
        * float(physics["beam_minor_arcsec"])
        * u.arcsec**2
    )
    peak = (
        np.asarray(catalogue["FXP_BST01"]) * u.Jy / u.beam
    ).to(u.MJy / u.sr, equivalencies=u.beam_angular_area(beam_area))
    peak_error = (
        np.asarray(catalogue["FXP_ERR01"]) * u.Jy / u.beam
    ).to(u.MJy / u.sr, equivalencies=u.beam_angular_area(beam_area))
    angle = np.asarray(catalogue["Deconv_FWHM_angle"]) * u.arcsec

    def background_parameters(values, index):
        return {
            "nu": {"dist": "fixed", "value": frequency},
            "kappa": {
                "dist": "normal",
                "mean": kappa,
                "std": kappa_error,
            },
            "Iv": {
                "dist": "normal",
                "mean": values[index],
                "std": background_error,
            },
            "temp": {
                "dist": "normal",
                "mean": temperature[index],
                "std": temperature_error[index],
            },
        }

    to_surface = lambda value: (
        value * float(physics["mean_molecular_weight"]) * const.m_p
    ).to(u.g / u.cm**2)
    surface_intercore, error_surface_intercore = _monte_carlo_series(
        brightness_to_column_density,
        lambda index: background_parameters(intercore, index),
        count,
        sample_count=samples,
        seed=seed,
        transform=to_surface,
    )
    surface_background, error_surface_background = _monte_carlo_series(
        brightness_to_column_density,
        lambda index: background_parameters(background, index),
        count,
        sample_count=samples,
        seed=seed,
        transform=to_surface,
    )

    def tau_parameters(index):
        return {
            "nu": {"dist": "fixed", "value": frequency},
            "Isr": {
                "dist": "normal",
                "mean": peak[index],
                "std": peak_error[index],
            },
            "temp": {
                "dist": "normal",
                "mean": temperature[index],
                "std": temperature_error[index],
            },
            "Ibg": {
                "dist": "normal",
                "mean": background[index],
                "std": background_error,
            },
        }

    tau, error_tau = _monte_carlo_series(
        core_optical_depth,
        tau_parameters,
        count,
        sample_count=samples,
        seed=seed,
    )

    def mass_parameters(index):
        return {
            "nu": {"dist": "fixed", "value": frequency},
            "kappa": {
                "dist": "normal",
                "mean": kappa,
                "std": kappa_error,
            },
            "ftot": {
                "dist": "normal",
                "mean": total_flux[index],
                "std": total_flux_error[index],
            },
            "temp": {
                "dist": "normal",
                "mean": temperature[index],
                "std": temperature_error[index],
            },
            "tau_c": {
                "dist": "normal",
                "mean": tau[index],
                "std": error_tau[index],
            },
            "dist": {
                "dist": "normal",
                "mean": distance,
                "std": distance_error,
            },
        }

    mass, error_mass = _monte_carlo_series(
        thick_mass,
        mass_parameters,
        count,
        sample_count=samples,
        seed=seed,
    )

    def surface_parameters(index):
        parameters = tau_parameters(index)
        return {
            "nu": parameters["nu"],
            "kappa": {
                "dist": "normal",
                "mean": kappa,
                "std": kappa_error,
            },
            "Isr": parameters["Isr"],
            "temp": parameters["temp"],
            "Ibg": parameters["Ibg"],
        }

    core_surface, error_core_surface = _monte_carlo_series(
        core_surface_density,
        surface_parameters,
        count,
        sample_count=samples,
        seed=seed,
    )

    def volume_parameters(index):
        parameters = mass_parameters(index)
        parameters["angle"] = {
            "dist": "fixed",
            "value": angle[index],
        }
        return parameters

    volume, error_volume = _monte_carlo_series(
        core_volume_density,
        volume_parameters,
        count,
        sample_count=samples,
        seed=seed,
    )
    volume_quantity = volume * 1e-16 * u.g / u.cm**3
    error_volume_quantity = error_volume * 1e-16 * u.g / u.cm**3
    number_density = (
        volume_quantity
        / float(physics["mean_molecular_weight"])
        / const.m_p
    ).to(1e8 / u.cm**3)
    error_number_density = (
        error_volume_quantity
        / float(physics["mean_molecular_weight"])
        / const.m_p
    ).to(1e8 / u.cm**3)

    major = np.asarray(catalogue["Deconv_maj"], dtype=float)
    minor = np.asarray(catalogue["Deconv_min"], dtype=float)
    ellipticity = np.sqrt(1 - minor**2 / major**2)
    output = catalogue[
        ["ID", "WCS_ACOOR", "WCS_DCOOR", "Tenv", "Subregion"]
    ]
    output.add_columns(
        [
            angle,
            np.asarray(catalogue["Deconv_FWHM_linear"]) * u.au,
            ellipticity,
            surface_intercore * u.g / u.cm**2,
            error_surface_intercore * u.g / u.cm**2,
            surface_background * u.g / u.cm**2,
            error_surface_background * u.g / u.cm**2,
            tau,
            error_tau,
            core_surface * u.g / u.cm**2,
            error_core_surface * u.g / u.cm**2,
            mass * u.M_sun,
            error_mass * u.M_sun,
            volume_quantity.to(1e-16 * u.g / u.cm**3),
            error_volume_quantity.to(1e-16 * u.g / u.cm**3),
            number_density,
            error_number_density,
        ],
        names=[
            "FWHM_angle",
            "FWHM_linear",
            "Ellipticity",
            "Surface_intercorebkg",
            "err_Surface_intercorebkg",
            "Surface_corebkg",
            "err_Surface_corebkg",
            "Tau_core",
            "err_Tau_core",
            "Surface_core",
            "err_Surface_core",
            "Mass_core",
            "err_Mass_core",
            "Volden_core",
            "err_Volden_core",
            "nVolden_core",
            "err_nVolden_core",
        ],
    )

    products = output_paths(paths)
    for path in products.values():
        _prepare_output(path, overwrite=overwrite)
    output.write(products["catalogue"], format="ascii.ecsv", overwrite=True)

    publication = output[["ID", "Tenv", "FWHM_linear"]]
    for name, values, errors in (
        (
            "Surf_intercorebkg",
            surface_intercore,
            error_surface_intercore,
        ),
        (
            "Surf_corebkg",
            surface_background,
            error_surface_background,
        ),
        ("Tau_c", tau, error_tau),
        ("Surf_c", core_surface, error_core_surface),
        ("M_c", mass, error_mass),
        ("Volden_c", volume, error_volume),
        (
            "nVolden_c",
            number_density.value,
            error_number_density.value,
        ),
    ):
        publication[name] = [
            _round_uncertainty(value, error)
            for value, error in zip(values, errors)
        ]
    publication.write(
        products["publication"],
        format="aastex",
        overwrite=True,
        formats={
            "Tenv": "%.1f",
            "FWHM_linear": lambda value: (
                "--" if not np.isfinite(value) else f"{int(round(value, -2)):d}"
            ),
        },
    )
    print(f"Wrote {products['catalogue']} ({len(output)} sources)")
    print(f"Wrote {products['publication']}")
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
        build_properties(paths, config, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
