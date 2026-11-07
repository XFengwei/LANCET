"""Tests for Bonnor--Ebert stability calculations and task wiring."""

import copy
from pathlib import Path
import subprocess
import sys

from astropy import constants as const
from astropy import units as u
from astropy.table import Table
import numpy as np
import pytest

from lancet.stability import (
    bonnor_ebert_mass_external_pressure,
    bonnor_ebert_mass_pressure,
    bonnor_ebert_mass_radius,
    bonnor_ebert_ratio,
    effective_support_speed,
    external_pressure_from_density,
    implied_critical_boundary,
    mass_density_from_surface_density,
    nonthermal_velocity_dispersion,
    sound_speed,
)
from scripts.bonnor_ebert import calculate_catalogue, load_task_config


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_sound_speed_matches_isothermal_definition() -> None:
    result = sound_speed(20 * u.K)
    expected = np.sqrt(const.k_B * 20 * u.K / (2.33 * const.m_p))
    assert result.to_value(u.km / u.s) == pytest.approx(
        expected.to_value(u.km / u.s)
    )


def test_nonthermal_linewidth_removes_tracer_thermal_width() -> None:
    observed = nonthermal_velocity_dispersion(
        1.0 * u.km / u.s,
        temperature=30 * u.K,
        tracer_molecular_weight=30.0,
    )
    raw = 1.0 / np.sqrt(8.0 * np.log(2.0))
    assert 0 < observed.to_value(u.km / u.s) < raw


def test_radius_and_pressure_critical_masses_are_positive() -> None:
    radius_mass = bonnor_ebert_mass_radius(20 * u.K, 0.01 * u.pc)
    pressure_mass = bonnor_ebert_mass_pressure(
        20 * u.K,
        number_density_h2=1e5 * u.cm**-3,
    )
    assert radius_mass.to_value(u.M_sun) > 0
    assert pressure_mass.to_value(u.M_sun) > 0
    assert bonnor_ebert_ratio(
        2 * radius_mass,
        radius_mass,
    ) == pytest.approx(2.0)


def test_pressure_mass_requires_exactly_one_density() -> None:
    with pytest.raises(ValueError, match="exactly one"):
        bonnor_ebert_mass_pressure(20 * u.K)
    with pytest.raises(ValueError, match="exactly one"):
        bonnor_ebert_mass_pressure(
            20 * u.K,
            number_density_h2=1e5 * u.cm**-3,
            mass_density=1e-18 * u.g / u.cm**3,
        )


def test_measured_pressure_separates_internal_and_external_terms() -> None:
    temperature = 20 * u.K
    density = 2e4 * u.cm**-3
    nonthermal = 0.8 * u.km / u.s
    pressure = external_pressure_from_density(
        temperature=temperature,
        number_density_h2=density,
        nonthermal_dispersion=nonthermal,
    )
    support = effective_support_speed(temperature)
    mass = bonnor_ebert_mass_external_pressure(support, pressure)
    wrapper = bonnor_ebert_mass_pressure(
        temperature,
        number_density_h2=density,
        nonthermal_dispersion=nonthermal,
    )
    assert mass.to_value(u.M_sun) == pytest.approx(
        wrapper.to_value(u.M_sun)
    )


def test_surface_density_and_depth_give_volume_density() -> None:
    surface_density = 0.4 * u.g / u.cm**2
    depth = 0.4 * u.pc
    result = mass_density_from_surface_density(surface_density, depth)
    assert result.to_value(u.g / u.cm**3) == pytest.approx(
        (surface_density / depth).to_value(u.g / u.cm**3)
    )


def test_radius_form_records_its_implied_boundary_pressure() -> None:
    temperature = 20 * u.K
    radius = 0.01 * u.pc
    boundary = implied_critical_boundary(temperature, radius)
    radius_mass = bonnor_ebert_mass_radius(temperature, radius)
    pressure_mass = bonnor_ebert_mass_external_pressure(
        sound_speed(temperature),
        boundary["external_pressure"],
    )
    assert pressure_mass.to_value(u.M_sun) == pytest.approx(
        radius_mass.to_value(u.M_sun),
        rel=0.004,
    )
    assert boundary["edge_number_density_h2"].to_value(u.cm**-3) > 0


def test_bonnor_ebert_task_inspection_reports_radius_assumption() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/bonnor_ebert.py",
            "--stage",
            "inspect",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    assert "rows=231, valid=209" in completed.stdout
    assert "Profile: radius_thermal" in completed.stdout
    assert "Formula: radius" in completed.stdout
    assert "not independently measured" in completed.stdout


def test_bonnor_ebert_task_can_select_measured_external_profile() -> None:
    config = load_task_config(
        PROJECT_ROOT / "tasks" / "bonnor_ebert" / "task.toml"
    )
    fraction = config["profiles"]["measured_external"][
        "internal_nonthermal"
    ]["fraction"]
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/bonnor_ebert.py",
            "--stage",
            "inspect",
            "--profile",
            "measured_external",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    assert "Profile: measured_external" in completed.stdout
    assert "Formula: external_pressure" in completed.stdout
    assert "rho_external = Sigma_intercore / 0.4 pc" in completed.stdout
    assert (
        "Internal/environment nonthermal fraction: "
        f"{float(fraction):g}"
    ) in completed.stdout


def test_measured_external_rejects_unphysical_turbulence_fraction() -> None:
    config = load_task_config(
        PROJECT_ROOT / "tasks" / "bonnor_ebert" / "task.toml"
    )
    invalid = copy.deepcopy(config)
    invalid["profiles"]["measured_external"]["internal_nonthermal"][
        "fraction"
    ] = 1.1
    catalogue = Table.read(
        PROJECT_ROOT
        / "results"
        / "environment_measurements"
        / "tables"
        / "core_environment_catalogue.ecsv",
        format="ascii.ecsv",
    )
    with pytest.raises(ValueError, match="between 0 and 1"):
        calculate_catalogue(catalogue[:3], invalid, "measured_external")


def test_measured_external_scales_internal_environment_turbulence() -> None:
    config = load_task_config(
        PROJECT_ROOT / "tasks" / "bonnor_ebert" / "task.toml"
    )
    catalogue = Table.read(
        PROJECT_ROOT
        / "results"
        / "environment_measurements"
        / "tables"
        / "core_environment_catalogue.ecsv",
        format="ascii.ecsv",
    )
    output = calculate_catalogue(catalogue[:3], config, "measured_external")
    internal = np.asarray(
        output["Internal_nonthermal_dispersion"],
        dtype=float,
    )
    external = np.asarray(
        output["External_nonthermal_dispersion"],
        dtype=float,
    )
    fraction = float(
        config["profiles"]["measured_external"]["internal_nonthermal"][
            "fraction"
        ]
    )
    assert internal == pytest.approx(fraction * external)


def test_generated_profiles_record_pressure_provenance() -> None:
    root = PROJECT_ROOT / "results" / "bonnor_ebert"
    radius_path = (
        root
        / "radius_thermal"
        / "tables"
        / "bonnor_ebert_stability.ecsv"
    )
    pressure_path = (
        root
        / "measured_external"
        / "tables"
        / "bonnor_ebert_stability.ecsv"
    )
    if not radius_path.is_file() or not pressure_path.is_file():
        pytest.skip("Generated BE profile catalogues are unavailable")
    radius = Table.read(radius_path, format="ascii.ecsv")
    pressure = Table.read(pressure_path, format="ascii.ecsv")
    assert radius.meta["formula"] == "radius"
    assert "Implied_edge_number_density_H2" in radius.colnames
    assert pressure.meta["formula"] == "external_pressure"
    assert "External_pressure" in pressure.colnames
    assert "External_velocity_dispersion" in pressure.colnames
    assert "External_surface_density_intercore" in pressure.colnames
    assert "External_number_density_H2" in pressure.colnames
    assert "External_observed_tracer_dispersion" in pressure.colnames
    assert "External_nonthermal_dispersion" in pressure.colnames
    assert "External_total_dispersion" in pressure.colnames
    assert "Internal_turbulence_fraction" in pressure.colnames
