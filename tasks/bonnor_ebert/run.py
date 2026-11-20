#!/usr/bin/env python
"""Calculate configurable Bonnor--Ebert stability diagnostics."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
import json
from pathlib import Path
import tomllib

from astropy import constants as const
from astropy import units as u
from astropy.table import Column, Table
import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt

from lancet import ProjectPaths
from lancet.physical_properties import resolved_core_mask
from lancet.stability import (
    bonnor_ebert_mass_external_pressure,
    bonnor_ebert_mass_pressure,
    bonnor_ebert_mass_radius,
    bonnor_ebert_ratio,
    effective_support_speed,
    external_pressure_from_density,
    implied_critical_boundary,
    mass_density_from_h2_number_density,
    mass_density_from_surface_density,
    nonthermal_velocity_dispersion,
    sound_speed,
)


DEFAULT_CONFIG = Path("tasks/bonnor_ebert/task.toml")


def load_task_config(path: Path) -> dict[str, object]:
    """Load task settings from TOML."""
    if not path.is_file():
        raise FileNotFoundError(f"Task configuration not found: {path}")
    with path.open("rb") as handle:
        return tomllib.load(handle)


def _mapping(
    config: Mapping[str, object],
    name: str,
) -> Mapping[str, object]:
    value = config[name]
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a TOML table")
    return value


def _project_path(project: ProjectPaths, value: object) -> Path:
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = project.root / path
    return path.resolve()


def resolve_task_paths(
    project: ProjectPaths,
    config: Mapping[str, object],
) -> dict[str, Path]:
    """Resolve the analysis catalogue and policy-safe base output root."""
    inputs = _mapping(config, "inputs")
    outputs = _mapping(config, "outputs")
    paths = {
        name: _project_path(project, value)
        for name, value in inputs.items()
    }
    paths["output_base"] = project.validate_output(str(outputs["root"]))
    for name in inputs:
        if not paths[name].is_file():
            raise FileNotFoundError(f"{name} not found: {paths[name]}")
    return paths


def select_profile(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    requested: str | None = None,
) -> tuple[str, Mapping[str, object], dict[str, Path]]:
    """Select a configured formula profile and isolate its output tree."""
    calculation = _mapping(config, "calculation")
    profiles = _mapping(config, "profiles")
    name = requested or str(calculation["active_profile"])
    if name not in profiles:
        raise KeyError(
            f"Unknown BE profile {name!r}; choose from "
            + ", ".join(str(value) for value in profiles)
        )
    profile = profiles[name]
    if not isinstance(profile, Mapping):
        raise TypeError(f"profiles.{name} must be a TOML table")
    selected_paths = dict(paths)
    selected_paths["output_root"] = (
        paths["output_base"] / name
    ).resolve()
    return name, profile, selected_paths


def catalogue_path(paths: Mapping[str, Path]) -> Path:
    return (
        paths["output_root"]
        / "tables"
        / "bonnor_ebert_stability.ecsv"
    )


def summary_path(paths: Mapping[str, Path]) -> Path:
    return (
        paths["output_root"]
        / "tables"
        / "bonnor_ebert_summary.json"
    )


def example_path(paths: Mapping[str, Path]) -> Path:
    return (
        paths["output_root"]
        / "tables"
        / "bonnor_ebert_examples.json"
    )


def _figure_paths(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
) -> list[Path]:
    figures = _mapping(config, "figures")
    formats = [str(value).lower() for value in figures["formats"]]
    unsupported = set(formats) - {"pdf", "png"}
    if unsupported:
        raise ValueError(
            "Unsupported figure formats: " + ", ".join(sorted(unsupported))
        )
    return [
        paths["output_root"]
        / "figures"
        / f"BonnorEbert_stability.{extension}"
        for extension in formats
    ]


def _prepare_outputs(
    outputs: Sequence[Path],
    *,
    overwrite: bool,
) -> None:
    for output in outputs:
        if output.exists() and not overwrite:
            raise FileExistsError(
                f"Output already exists: {output}; "
                "pass --overwrite to replace it"
            )
    for output in outputs:
        output.parent.mkdir(parents=True, exist_ok=True)


def _catalogue_settings(
    config: Mapping[str, object],
) -> Mapping[str, object]:
    return _mapping(config, "catalogue")


def _quantity_source(
    catalogue: Table,
    settings: Mapping[str, object],
    *,
    name: str,
) -> u.Quantity | None:
    """Resolve a profile quantity from a catalogue column or scalar."""
    source = str(settings["source"]).lower()
    if source == "none":
        return None
    unit_text = str(settings["unit"])
    if source == "column":
        column_name = str(settings["column"])
        if column_name not in catalogue.colnames:
            raise KeyError(
                f"{name} column {column_name!r} is missing from catalogue"
            )
        column = catalogue[column_name]
        if unit_text == "native":
            if column.unit is None:
                raise ValueError(
                    f"{name} requests native units but {column_name} "
                    "has no unit"
                )
            quantity = u.Quantity(column)
        else:
            quantity = np.asarray(column, dtype=float) * u.Unit(unit_text)
    elif source == "scalar":
        if unit_text == "native":
            raise ValueError(f"{name} scalar source requires an explicit unit")
        quantity = (
            np.full(len(catalogue), float(settings["value"]))
            * u.Unit(unit_text)
        )
    else:
        raise ValueError(
            f"{name}.source must be column, scalar, or none"
        )
    if np.asarray(quantity.value).shape != (len(catalogue),):
        raise ValueError(f"{name} must resolve to one value per catalogue row")
    return quantity


def _valid_quantity(
    quantity: u.Quantity,
    *,
    allow_zero: bool = False,
) -> np.ndarray:
    values = np.asarray(quantity.value, dtype=float)
    return np.isfinite(values) & (
        values >= 0 if allow_zero else values > 0
    )


def _source_description(settings: Mapping[str, object]) -> str:
    """Return a compact provenance label for a configured quantity."""
    source = str(settings["source"]).lower()
    if source == "column":
        return f"column:{settings['column']}"
    if source == "scalar":
        return (
            f"scalar:{settings['value']} "
            f"{settings['unit']}"
        )
    if source == "external_fraction":
        return f"external_nonthermal_fraction:{settings['fraction']}"
    return source


def _velocity_dispersion(
    catalogue: Table,
    settings: Mapping[str, object],
    temperature: u.Quantity,
    *,
    name: str,
) -> tuple[u.Quantity | None, str]:
    """Resolve a nonthermal or total one-dimensional dispersion."""
    if str(settings["source"]).lower() == "external_fraction":
        return None, "external_fraction"
    values = _quantity_source(catalogue, settings, name=name)
    if values is None:
        return None, "none"
    velocity = values.to(u.km / u.s)
    kind = str(settings["kind"]).lower()
    if kind == "nonthermal_sigma":
        return velocity, "nonthermal"
    if kind == "total_sigma":
        return velocity, "total"
    if kind not in {"observed_sigma", "observed_fwhm"}:
        raise ValueError(
            f"{name}.kind must be nonthermal_sigma, total_sigma, "
            "observed_sigma, or observed_fwhm"
        )
    observed_sigma = (
        velocity
        if kind == "observed_sigma"
        else velocity / np.sqrt(8.0 * np.log(2.0))
    )
    if not bool(settings["remove_tracer_thermal"]):
        return observed_sigma, "total"
    tracer_weight = float(settings["tracer_molecular_weight"])
    tracer_thermal = sound_speed(
        temperature,
        mean_particle_weight=tracer_weight,
    )
    difference = observed_sigma**2 - tracer_thermal**2
    result = np.full(len(catalogue), np.nan) * u.km / u.s
    valid = (
        np.isfinite(difference.to_value((u.km / u.s) ** 2))
        & (difference >= 0 * (u.km / u.s) ** 2)
    )
    result[valid] = np.sqrt(difference[valid])
    return result, "nonthermal"


def _base_catalogue_mask(
    catalogue: Table,
    settings: Mapping[str, object],
) -> np.ndarray:
    required = [
        str(settings[name])
        for name in (
            "id_column",
            "region_column",
            "mass_column",
            "mass_error_column",
        )
    ]
    missing = set(required) - set(catalogue.colnames)
    if missing:
        raise KeyError(
            "Catalogue is missing columns: " + ", ".join(sorted(missing))
        )
    mass = np.asarray(
        catalogue[str(settings["mass_column"])],
        dtype=float,
    )
    return np.isfinite(mass) & (mass > 0)


def calculate_catalogue(
    catalogue: Table,
    config: Mapping[str, object],
    profile_name: str | None = None,
) -> Table:
    """Calculate a configured radius- or pressure-based BE catalogue."""
    settings = _catalogue_settings(config)
    profiles = _mapping(config, "profiles")
    if profile_name is None:
        profile_name = str(
            _mapping(config, "calculation")["active_profile"]
        )
    profile_value = profiles[profile_name]
    if not isinstance(profile_value, Mapping):
        raise TypeError(f"profiles.{profile_name} must be a TOML table")
    profile = profile_value
    formula = str(profile["formula"]).lower()
    if formula not in {"radius", "external_pressure"}:
        raise ValueError(
            "Profile formula must be radius or external_pressure"
        )
    mean_particle_weight = float(profile["mean_particle_weight"])
    mean_mass_per_h2 = float(profile["mean_mass_per_h2"])

    temperature_settings = profile["support_temperature"]
    internal_settings = profile["internal_nonthermal"]
    if not isinstance(temperature_settings, Mapping) or not isinstance(
        internal_settings,
        Mapping,
    ):
        raise TypeError("Profile quantity definitions must be TOML tables")
    temperature = _quantity_source(
        catalogue,
        temperature_settings,
        name="support_temperature",
    )
    assert temperature is not None
    temperature = temperature.to(u.K)
    internal_nonthermal, internal_role = _velocity_dispersion(
        catalogue,
        internal_settings,
        temperature,
        name="internal_nonthermal",
    )
    if internal_role == "total":
        raise ValueError(
            "Internal support must be nonthermal or none; total_sigma "
            "would double-count the gas thermal term"
        )
    configured_support = str(profile["internal_support"]).lower()
    expected_support = (
        "thermal_plus_nonthermal"
        if internal_role == "external_fraction"
        else (
            "thermal"
            if internal_nonthermal is None
            else "thermal_plus_nonthermal"
        )
    )
    if configured_support != expected_support:
        raise ValueError(
            f"internal_support={configured_support!r} conflicts with "
            f"internal_nonthermal.source; expected {expected_support!r}"
        )

    mask = _base_catalogue_mask(catalogue, settings)
    mask &= _valid_quantity(temperature)
    if internal_nonthermal is not None:
        mask &= _valid_quantity(internal_nonthermal, allow_zero=True)

    extra_columns: list[Column] = []
    extra_names: list[str] = []
    formula_definition: str
    pressure_definition: str
    density_definition = "Not applicable"
    internal_turbulence_fraction: float | None = None
    if formula == "radius":
        if internal_role == "external_fraction":
            raise ValueError(
                "internal_nonthermal.source=external_fraction requires "
                "the external_pressure formula"
            )
        radius_settings = profile["radius"]
        if not isinstance(radius_settings, Mapping):
            raise TypeError("profiles.<name>.radius must be a TOML table")
        radius_base = _quantity_source(
            catalogue,
            radius_settings,
            name="radius",
        )
        assert radius_base is not None
        radius_factor = float(radius_settings["factor"])
        if not np.isfinite(radius_factor) or radius_factor <= 0:
            raise ValueError("radius.factor must be finite and positive")
        radius = radius_base.to(u.au) * radius_factor
        mask &= _valid_quantity(radius)
        if radius_settings.get("resolved_only", False):
            mask &= resolved_core_mask(catalogue)
        selected_temperature = temperature[mask]
        selected_internal = (
            None
            if internal_nonthermal is None
            else internal_nonthermal[mask]
        )
        selected_radius = radius[mask]
        critical_mass = bonnor_ebert_mass_radius(
            selected_temperature,
            selected_radius,
            nonthermal_dispersion=selected_internal,
            mean_particle_weight=mean_particle_weight,
        )
        implied = implied_critical_boundary(
            selected_temperature,
            selected_radius,
            nonthermal_dispersion=selected_internal,
            mean_particle_weight=mean_particle_weight,
            mean_mass_per_h2=mean_mass_per_h2,
        )
        pressure = implied["external_pressure"]
        extra_columns.extend(
            [
                Column(selected_radius.to_value(u.au), unit=u.au),
                Column(
                    implied["edge_number_density_h2"].to_value(u.cm**-3),
                    unit=u.cm**-3,
                ),
                Column(
                    pressure.to_value(u.dyne / u.cm**2),
                    unit=u.dyne / u.cm**2,
                ),
                Column(
                    (pressure / const.k_B).to_value(u.K / u.cm**3),
                    unit=u.K / u.cm**3,
                ),
            ]
        )
        extra_names.extend(
            [
                "Radius_BE",
                "Implied_edge_number_density_H2",
                "Implied_external_pressure",
                "Implied_external_pressure_over_k",
            ]
        )
        formula_definition = "2.4 R a_support^2 / G"
        pressure_definition = (
            "Implied by critical BE radius; not independently measured"
        )
    else:
        mode = str(profile["external_pressure_mode"]).lower()
        external_temperature_settings = profile["external_temperature"]
        if not isinstance(external_temperature_settings, Mapping):
            raise TypeError("external_temperature must be a TOML table")
        external_temperature = _quantity_source(
            catalogue,
            external_temperature_settings,
            name="external_temperature",
        )
        assert external_temperature is not None
        external_temperature = external_temperature.to(u.K)
        mask &= _valid_quantity(external_temperature)

        external_density_number: u.Quantity | None = None
        external_density_mass: u.Quantity | None = None
        external_surface_density: u.Quantity | None = None
        external_path_length: u.Quantity | None = None
        external_velocity: u.Quantity | None = None
        external_observed_sigma: u.Quantity | None = None
        external_velocity_role = "none"
        if mode == "density_velocity":
            density_settings = profile["external_density"]
            velocity_settings = profile["external_velocity"]
            if not isinstance(density_settings, Mapping) or not isinstance(
                velocity_settings,
                Mapping,
            ):
                raise TypeError(
                    "external density/velocity definitions must be TOML tables"
                )
            density = _quantity_source(
                catalogue,
                density_settings,
                name="external_density",
            )
            assert density is not None
            density_kind = str(density_settings["kind"]).lower()
            if density_kind == "number_density_h2":
                external_density_number = density.to(u.cm**-3)
                mask &= _valid_quantity(external_density_number)
                density_definition = (
                    "rho_external = mean_mass_per_H2 * m_H * n_H2"
                )
            elif density_kind == "mass_density":
                external_density_mass = density.to(u.g / u.cm**3)
                mask &= _valid_quantity(external_density_mass)
                density_definition = "Configured external mass density"
            elif density_kind == "surface_density_path_length":
                external_surface_density = density.to(u.g / u.cm**2)
                external_path_length = (
                    float(density_settings["path_length_pc"]) * u.pc
                )
                if (
                    not np.isfinite(external_path_length.to_value(u.pc))
                    or external_path_length <= 0 * u.pc
                ):
                    raise ValueError(
                        "external_density.path_length_pc must be finite "
                        "and positive"
                    )
                external_density_mass = mass_density_from_surface_density(
                    external_surface_density,
                    external_path_length,
                )
                mask &= _valid_quantity(external_density_mass)
                density_definition = (
                    "rho_external = Sigma_intercore / "
                    f"{external_path_length.to_value(u.pc):g} pc"
                )
            else:
                raise ValueError(
                    "external_density.kind must be number_density_h2, "
                    "mass_density, or surface_density_path_length"
                )
            velocity_input = _quantity_source(
                catalogue,
                velocity_settings,
                name="external_velocity",
            )
            velocity_kind = str(velocity_settings["kind"]).lower()
            if (
                velocity_input is not None
                and velocity_kind in {"observed_sigma", "observed_fwhm"}
            ):
                external_observed_sigma = velocity_input.to(u.km / u.s)
                if velocity_kind == "observed_fwhm":
                    external_observed_sigma = (
                        external_observed_sigma
                        / np.sqrt(8.0 * np.log(2.0))
                    )
            external_velocity, external_velocity_role = _velocity_dispersion(
                catalogue,
                velocity_settings,
                external_temperature,
                name="external_velocity",
            )
            if external_velocity is not None:
                mask &= _valid_quantity(
                    external_velocity,
                    allow_zero=external_velocity_role == "nonthermal",
                )
        elif mode == "direct":
            direct_settings = profile["direct_external_pressure"]
            if not isinstance(direct_settings, Mapping):
                raise TypeError(
                    "direct_external_pressure must be a TOML table"
                )
            direct_pressure = _quantity_source(
                catalogue,
                direct_settings,
                name="direct_external_pressure",
            )
            assert direct_pressure is not None
            representation = str(
                direct_settings["representation"]
            ).lower()
            if representation == "pressure_over_k":
                direct_pressure = (
                    direct_pressure.to(u.K / u.cm**3) * const.k_B
                ).to(u.dyne / u.cm**2)
            elif representation == "pressure":
                direct_pressure = direct_pressure.to(u.dyne / u.cm**2)
            else:
                raise ValueError(
                    "direct pressure representation must be pressure "
                    "or pressure_over_k"
                )
            mask &= _valid_quantity(direct_pressure)
        else:
            raise ValueError(
                "external_pressure_mode must be density_velocity or direct"
            )

        selected_temperature = temperature[mask]
        selected_external_temperature = external_temperature[mask]
        if mode == "direct":
            if internal_role == "external_fraction":
                raise ValueError(
                    "internal_nonthermal.source=external_fraction requires "
                    "external_pressure_mode=density_velocity"
                )
            selected_internal = (
                None
                if internal_nonthermal is None
                else internal_nonthermal[mask]
            )
            pressure = direct_pressure[mask]
            pressure_definition = "Direct configured external pressure"
            density_for_output = None
            velocity_for_output = None
            observed_velocity_for_output = None
            surface_density_for_output = None
            path_length_for_output = None
        else:
            number_density = (
                None
                if external_density_number is None
                else external_density_number[mask]
            )
            mass_density = (
                None
                if external_density_mass is None
                else external_density_mass[mask]
            )
            selected_external_velocity = (
                None
                if external_velocity is None
                else external_velocity[mask]
            )
            if internal_role == "external_fraction":
                if (
                    selected_external_velocity is None
                    or external_velocity_role != "nonthermal"
                ):
                    raise ValueError(
                        "Scaled internal turbulence requires a derived "
                        "external nonthermal dispersion"
                    )
                internal_turbulence_fraction = float(
                    internal_settings["fraction"]
                )
                if (
                    not np.isfinite(internal_turbulence_fraction)
                    or not 0 <= internal_turbulence_fraction <= 1
                ):
                    raise ValueError(
                        "internal_nonthermal.fraction must be between "
                        "0 and 1"
                    )
                selected_internal = (
                    internal_turbulence_fraction
                    * selected_external_velocity
                )
            else:
                selected_internal = (
                    None
                    if internal_nonthermal is None
                    else internal_nonthermal[mask]
                )
            pressure = external_pressure_from_density(
                temperature=selected_external_temperature,
                number_density_h2=number_density,
                mass_density=mass_density,
                nonthermal_dispersion=(
                    selected_external_velocity
                    if external_velocity_role == "nonthermal"
                    else None
                ),
                total_velocity_dispersion=(
                    selected_external_velocity
                    if external_velocity_role == "total"
                    else None
                ),
                mean_particle_weight=mean_particle_weight,
                mean_mass_per_h2=mean_mass_per_h2,
            )
            density_for_output = (
                mass_density_from_h2_number_density(
                    number_density,
                    mean_mass_per_h2=mean_mass_per_h2,
                )
                if mass_density is None
                else mass_density
            )
            velocity_for_output = selected_external_velocity
            observed_velocity_for_output = (
                None
                if external_observed_sigma is None
                else external_observed_sigma[mask]
            )
            surface_density_for_output = (
                None
                if external_surface_density is None
                else external_surface_density[mask]
            )
            path_length_for_output = external_path_length
            pressure_definition = (
                "P_external = rho_external * "
                "(c_s(environment)^2 + sigma_nt(environment)^2)"
                if external_velocity_role == "nonthermal"
                else (
                    "Measured/configured density times external effective "
                    "one-dimensional velocity dispersion squared"
                )
            )
        support = effective_support_speed(
            selected_temperature,
            nonthermal_dispersion=selected_internal,
            mean_particle_weight=mean_particle_weight,
        )
        critical_mass = bonnor_ebert_mass_external_pressure(
            support,
            pressure,
        )
        extra_columns.extend(
            [
                Column(
                    selected_external_temperature.to_value(u.K),
                    unit=u.K,
                ),
                Column(
                    pressure.to_value(u.dyne / u.cm**2),
                    unit=u.dyne / u.cm**2,
                ),
                Column(
                    (pressure / const.k_B).to_value(u.K / u.cm**3),
                    unit=u.K / u.cm**3,
                ),
            ]
        )
        extra_names.extend(
            [
                "External_temperature",
                "External_pressure",
                "External_pressure_over_k",
            ]
        )
        if density_for_output is not None:
            extra_columns.append(
                Column(
                    density_for_output.to_value(u.g / u.cm**3),
                    unit=u.g / u.cm**3,
                )
            )
            extra_names.append("External_mass_density")
            extra_columns.append(
                Column(
                    (
                        density_for_output
                        / (mean_mass_per_h2 * const.m_p)
                    ).to_value(u.cm**-3),
                    unit=u.cm**-3,
                )
            )
            extra_names.append("External_number_density_H2")
        if surface_density_for_output is not None:
            extra_columns.append(
                Column(
                    surface_density_for_output.to_value(u.g / u.cm**2),
                    unit=u.g / u.cm**2,
                )
            )
            extra_names.append("External_surface_density_intercore")
        if path_length_for_output is not None:
            extra_columns.append(
                Column(
                    np.full(
                        len(selected_temperature),
                        path_length_for_output.to_value(u.pc),
                    ),
                    unit=u.pc,
                )
            )
            extra_names.append("External_path_length")
        if observed_velocity_for_output is not None:
            extra_columns.append(
                Column(
                    observed_velocity_for_output.to_value(u.km / u.s),
                    unit=u.km / u.s,
                )
            )
            extra_names.append("External_observed_tracer_dispersion")
        if velocity_for_output is not None:
            extra_columns.append(
                Column(
                    velocity_for_output.to_value(u.km / u.s),
                    unit=u.km / u.s,
                )
            )
            extra_names.append("External_velocity_dispersion")
            if external_velocity_role == "nonthermal":
                extra_columns.append(
                    Column(
                        velocity_for_output.to_value(u.km / u.s),
                        unit=u.km / u.s,
                    )
                )
                extra_names.append("External_nonthermal_dispersion")
        if density_for_output is not None:
            extra_columns.append(
                Column(
                    np.sqrt(pressure / density_for_output).to_value(
                        u.km / u.s
                    ),
                    unit=u.km / u.s,
                )
            )
            extra_names.append("External_total_dispersion")
        if internal_turbulence_fraction is not None:
            extra_columns.append(
                Column(
                    np.full(
                        len(selected_temperature),
                        internal_turbulence_fraction,
                    ),
                )
            )
            extra_names.append("Internal_turbulence_fraction")
        formula_definition = (
            "1.18 a_support^4 / (G^(3/2) P_external^(1/2))"
        )

    selected = catalogue[mask]
    selected_temperature = temperature[mask]
    support = effective_support_speed(
        selected_temperature,
        nonthermal_dispersion=selected_internal,
        mean_particle_weight=mean_particle_weight,
    )
    mass = (
        np.asarray(
            selected[str(settings["mass_column"])],
            dtype=float,
        )
        * u.M_sun
    )
    ratio = bonnor_ebert_ratio(mass, critical_mass)
    mass_error = np.asarray(
        selected[str(settings["mass_error_column"])],
        dtype=float,
    )
    ratio_error = mass_error / critical_mass.to_value(u.M_sun)
    threshold = float(settings["instability_threshold"])

    output = selected[
        [
            str(settings["id_column"]),
            str(settings["region_column"]),
            str(settings["mass_column"]),
            str(settings["mass_error_column"]),
        ]
    ].copy()
    output.add_columns(
        [
            Column(selected_temperature.to_value(u.K), unit=u.K),
            Column(
                sound_speed(
                    selected_temperature,
                    mean_particle_weight=mean_particle_weight,
                ).to_value(u.km / u.s),
                unit=u.km / u.s,
            ),
            Column(
                (
                    np.zeros(len(selected))
                    if selected_internal is None
                    else selected_internal.to_value(u.km / u.s)
                ),
                unit=u.km / u.s,
            ),
            Column(
                support.to_value(u.km / u.s),
                unit=u.km / u.s,
            ),
            Column(critical_mass.to_value(u.M_sun), unit=u.M_sun),
            ratio,
            ratio_error,
            np.where(ratio >= threshold, "supercritical", "subcritical"),
            *extra_columns,
        ],
        names=[
            "Support_temperature",
            "Sound_speed",
            "Internal_nonthermal_dispersion",
            "Support_speed",
            "Mass_BE_critical",
            "Alpha_BE",
            "err_Alpha_BE_mass_only",
            "BE_state",
            *extra_names,
        ],
    )
    output.meta["profile"] = profile_name
    output.meta["formula"] = formula
    output.meta["critical_mass_definition"] = formula_definition
    output.meta["external_pressure_definition"] = pressure_definition
    output.meta["external_density_definition"] = density_definition
    output.meta["internal_support"] = str(profile["internal_support"])
    output.meta["support_temperature_source"] = _source_description(
        temperature_settings
    )
    output.meta["internal_nonthermal_source"] = _source_description(
        internal_settings
    )
    if internal_turbulence_fraction is not None:
        output.meta["internal_turbulence_fraction"] = (
            internal_turbulence_fraction
        )
    if formula == "external_pressure":
        output.meta["external_temperature_source"] = _source_description(
            external_temperature_settings
        )
        if mode == "density_velocity":
            output.meta["external_density_source"] = _source_description(
                density_settings
            )
            output.meta["external_velocity_source"] = _source_description(
                velocity_settings
            )
            output.meta["external_velocity_kind"] = str(
                velocity_settings["kind"]
            )
            output.meta["external_tracer_molecular_weight"] = float(
                velocity_settings["tracer_molecular_weight"]
            )
    output.meta["alpha_be_definition"] = "Mass_core / Mass_BE_critical"
    output.meta["instability_threshold"] = threshold
    output.meta["uncertainty_scope"] = (
        "Alpha_BE uncertainty propagates core-mass uncertainty only; "
        "support temperature, radius, density, velocity, and pressure "
        "uncertainties are not propagated"
    )
    output.meta["mean_particle_weight"] = mean_particle_weight
    output.meta["mean_mass_per_h2"] = mean_mass_per_h2
    output.meta["input_rows"] = len(catalogue)
    output.meta["valid_rows"] = len(output)
    return output


def summarize_catalogue(
    catalogue: Table,
    config: Mapping[str, object],
) -> dict[str, object]:
    """Summarize supercritical counts globally and by subregion."""
    settings = _catalogue_settings(config)
    threshold = float(settings["instability_threshold"])
    ratio = np.asarray(catalogue["Alpha_BE"], dtype=float)
    regions = np.asarray(
        catalogue[str(settings["region_column"])]
    ).astype(str)

    def one_summary(mask: np.ndarray) -> dict[str, object]:
        values = ratio[mask]
        unstable = values >= threshold
        return {
            "count": int(values.size),
            "supercritical_count": int(np.count_nonzero(unstable)),
            "supercritical_fraction": float(np.mean(unstable)),
            "median_alpha_be": float(np.median(values)),
            "p16_alpha_be": float(np.percentile(values, 16)),
            "p84_alpha_be": float(np.percentile(values, 84)),
        }

    summary: dict[str, object] = {
        "profile": catalogue.meta.get("profile"),
        "formula": catalogue.meta.get("formula"),
        "critical_mass_definition": catalogue.meta.get(
            "critical_mass_definition"
        ),
        "external_pressure_definition": catalogue.meta.get(
            "external_pressure_definition"
        ),
        "external_density_definition": catalogue.meta.get(
            "external_density_definition"
        ),
        "internal_support": catalogue.meta.get("internal_support"),
        "internal_turbulence_fraction": catalogue.meta.get(
            "internal_turbulence_fraction"
        ),
        "definition": "Alpha_BE = Mass_core / Mass_BE_critical",
        "instability_threshold": threshold,
        "all": one_summary(np.ones(len(catalogue), dtype=bool)),
        "subregions": {},
    }
    diagnostic_columns = (
        "External_surface_density_intercore",
        "External_number_density_H2",
        "External_observed_tracer_dispersion",
        "External_nonthermal_dispersion",
        "External_total_dispersion",
        "External_pressure_over_k",
        "Internal_nonthermal_dispersion",
        "Support_speed",
        "Mass_BE_critical",
    )
    summary["diagnostic_medians"] = {
        name: float(np.median(np.asarray(catalogue[name], dtype=float)))
        for name in diagnostic_columns
        if name in catalogue.colnames
    }
    summary["diagnostic_units"] = {
        name: str(catalogue[name].unit)
        for name in diagnostic_columns
        if name in catalogue.colnames
    }
    region_summary = summary["subregions"]
    assert isinstance(region_summary, dict)
    for code in np.unique(regions):
        region_summary[str(code)] = one_summary(regions == code)
    return summary


def inspect_inputs(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    profile_name: str,
    profile: Mapping[str, object],
) -> None:
    """Report the selected formula, sources, and output readiness."""
    catalogue = Table.read(paths["catalogue"], format="ascii.ecsv")
    output = calculate_catalogue(catalogue, config, profile_name)
    summary = summarize_catalogue(output, config)
    all_summary = summary["all"]
    assert isinstance(all_summary, Mapping)
    print("Bonnor-Ebert input:")
    print(
        f"  catalogue: {paths['catalogue']} "
        f"(rows={len(catalogue)}, valid={len(output)})"
    )
    print(f"Profile: {profile_name}")
    print(f"Formula: {profile['formula']}")
    print(f"Internal support: {profile['internal_support']}")
    print(
        f"Critical mass: {summary['critical_mass_definition']}"
    )
    print(
        "External pressure: "
        f"{summary['external_pressure_definition']}"
    )
    if str(profile["formula"]).lower() == "external_pressure":
        print(
            "External density: "
            f"{summary['external_density_definition']}"
        )
    if summary["internal_turbulence_fraction"] is not None:
        print(
            "Internal/environment nonthermal fraction: "
            f"{float(summary['internal_turbulence_fraction']):g}"
        )
    print(
        "Supercritical: "
        f"{all_summary['supercritical_count']}/{all_summary['count']} "
        f"({100 * float(all_summary['supercritical_fraction']):.1f}%) "
        f"with Alpha_BE >= {summary['instability_threshold']}"
    )
    print(f"Output root: {paths['output_root']}")


def build_catalogue(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    profile_name: str,
    *,
    overwrite: bool,
) -> Table:
    """Write the per-core stability table and JSON summary."""
    source = Table.read(paths["catalogue"], format="ascii.ecsv")
    output = calculate_catalogue(source, config, profile_name)
    summary = summarize_catalogue(output, config)
    products = [catalogue_path(paths), summary_path(paths)]
    _prepare_outputs(products, overwrite=overwrite)
    output.write(products[0], format="ascii.ecsv", overwrite=True)
    products[1].write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {products[0]} ({len(output)} cores)")
    print(f"Wrote {products[1]}")
    return output


def build_example(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool,
) -> dict[str, float]:
    """Evaluate the configured pressure- and radius-based examples."""
    example = _mapping(config, "example")
    temperature = float(example["temperature_k"]) * u.K
    fwhm = float(example["line_fwhm_kms"]) * u.km / u.s
    nonthermal = nonthermal_velocity_dispersion(
        fwhm,
        temperature=temperature,
        tracer_molecular_weight=float(
            example["tracer_molecular_weight"]
        ),
    )
    thermal_pressure_mass = bonnor_ebert_mass_pressure(
        temperature,
        number_density_h2=(
            float(example["number_density_h2_cm3"]) * u.cm**-3
        ),
        mean_particle_weight=float(example["mean_particle_weight"]),
        mean_mass_per_h2=float(example["mean_mass_per_h2"]),
    )
    turbulent_pressure_mass = bonnor_ebert_mass_pressure(
        temperature,
        number_density_h2=(
            float(example["number_density_h2_cm3"]) * u.cm**-3
        ),
        nonthermal_dispersion=nonthermal,
        mean_particle_weight=float(example["mean_particle_weight"]),
        mean_mass_per_h2=float(example["mean_mass_per_h2"]),
    )
    radius_mass = bonnor_ebert_mass_radius(
        temperature,
        float(example["radius_pc"]) * u.pc,
        mean_particle_weight=float(example["mean_particle_weight"]),
    )
    observed_mass = float(example["observed_mass_msun"]) * u.M_sun
    result = {
        "temperature_k": temperature.to_value(u.K),
        "number_density_h2_cm3": float(
            example["number_density_h2_cm3"]
        ),
        "line_fwhm_kms": fwhm.to_value(u.km / u.s),
        "nonthermal_dispersion_kms": nonthermal.to_value(u.km / u.s),
        "pressure_thermal_mass_msun": thermal_pressure_mass.to_value(
            u.M_sun
        ),
        "pressure_with_nonthermal_mass_msun": (
            turbulent_pressure_mass.to_value(u.M_sun)
        ),
        "radius_mass_msun": radius_mass.to_value(u.M_sun),
        "observed_to_radius_critical_ratio": float(
            bonnor_ebert_ratio(observed_mass, radius_mass)
        ),
    }
    output = example_path(paths)
    _prepare_outputs([output], overwrite=overwrite)
    output.write_text(
        json.dumps(result, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {output}")
    return result


def plot_catalogue(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool,
) -> list[Path]:
    """Plot critical versus measured mass and the stability-ratio distribution."""
    source = catalogue_path(paths)
    if not source.is_file():
        raise FileNotFoundError(
            f"Stability catalogue not found: {source}; run --stage catalogue"
        )
    catalogue = Table.read(source, format="ascii.ecsv")
    outputs = _figure_paths(paths, config)
    _prepare_outputs(outputs, overwrite=overwrite)
    figures = _mapping(config, "figures")
    settings = _catalogue_settings(config)
    regions = _mapping(config, "regions")

    figure, axes = plt.subplots(
        1,
        2,
        figsize=(
            float(figures["width_inches"]),
            float(figures["height_inches"]),
        ),
    )
    mass = np.asarray(catalogue[str(settings["mass_column"])], dtype=float)
    critical = np.asarray(catalogue["Mass_BE_critical"], dtype=float)
    ratio = np.asarray(catalogue["Alpha_BE"], dtype=float)
    region_values = np.asarray(
        catalogue[str(settings["region_column"])]
    ).astype(str)
    for _, region in regions.items():
        if not isinstance(region, Mapping):
            raise TypeError("Each regions entry must be a TOML table")
        code = str(region["code"])
        mask = region_values == code
        axes[0].scatter(
            critical[mask],
            mass[mask],
            s=float(figures["marker_size"]),
            facecolor=str(region["color"]),
            edgecolor="black",
            linewidth=float(figures["marker_edge_width"]),
            alpha=float(figures["marker_alpha"]),
            label=str(region["label"]),
            zorder=3,
        )
        axes[1].hist(
            ratio[mask],
            bins=np.logspace(
                np.log10(float(figures["ratio_min"])),
                np.log10(float(figures["ratio_max"])),
                int(figures["ratio_bins"]) + 1,
            ),
            histtype="step",
            linewidth=1.5,
            color=str(region["color"]),
            label=str(region["label"]),
        )

    bounds = np.array(
        [
            min(np.min(critical), np.min(mass)),
            max(np.max(critical), np.max(mass)),
        ]
    )
    axes[0].plot(
        bounds,
        bounds,
        color="black",
        linestyle="--",
        linewidth=1.2,
        label=r"$M_{\rm core}=M_{\rm BE}$",
    )
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlim(bounds)
    axes[0].set_ylim(bounds)
    axes[0].set_xlabel(r"$M_{\rm BE}$ ($M_\odot$)")
    axes[0].set_ylabel(r"$M_{\rm core}$ ($M_\odot$)")
    axes[0].legend(frameon=False, fontsize=9)
    axes[0].grid(alpha=0.18, which="both")

    threshold = float(settings["instability_threshold"])
    axes[1].axvline(
        threshold,
        color="black",
        linestyle="--",
        linewidth=1.2,
    )
    axes[1].set_xscale("log")
    axes[1].set_xlim(
        float(figures["ratio_min"]),
        float(figures["ratio_max"]),
    )
    axes[1].set_xlabel(r"$\alpha_{\rm BE}=M_{\rm core}/M_{\rm BE}$")
    axes[1].set_ylabel("Number of cores")
    axes[1].grid(alpha=0.18, which="both")
    profile_label = str(catalogue.meta.get("profile", "BE profile")).replace(
        "_",
        " ",
    )
    formula_label = str(
        catalogue.meta.get("formula", "unknown")
    ).replace("_", " ")
    fraction = catalogue.meta.get("internal_turbulence_fraction")
    fraction_label = (
        ""
        if fraction is None
        else rf", $f_{{\rm turb}}={float(fraction):g}$"
    )
    figure.suptitle(
        f"{profile_label}{fraction_label} ({formula_label} formula)",
        fontsize=11,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    dpi = int(figures["dpi"])
    for output in outputs:
        figure.savefig(output, dpi=dpi, bbox_inches="tight")
        print(f"Wrote {output}")
    plt.close(figure)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--stage",
        choices=("inspect", "catalogue", "example", "plot", "all"),
        default="inspect",
    )
    parser.add_argument(
        "--profile",
        help="Configured calculation profile; defaults to active_profile",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    project = ProjectPaths.discover()
    config_path = _project_path(project, args.config)
    config = load_task_config(config_path)
    base_paths = resolve_task_paths(project, config)
    profile_name, profile, paths = select_profile(
        base_paths,
        config,
        args.profile,
    )
    if args.stage == "inspect":
        inspect_inputs(paths, config, profile_name, profile)
    elif args.stage == "catalogue":
        build_catalogue(
            paths,
            config,
            profile_name,
            overwrite=args.overwrite,
        )
    elif args.stage == "example":
        build_example(paths, config, overwrite=args.overwrite)
    elif args.stage == "plot":
        plot_catalogue(paths, config, overwrite=args.overwrite)
    else:
        build_catalogue(
            paths,
            config,
            profile_name,
            overwrite=args.overwrite,
        )
        plot_catalogue(paths, config, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
