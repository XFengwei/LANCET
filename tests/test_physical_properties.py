"""Tests for core physical-property calculations."""

from pathlib import Path
import subprocess
import sys

from astropy import units as u
from astropy.table import Table
import numpy as np
import pytest

from lancet.physical_properties import (
    blackbody_cgs,
    core_optical_depth,
    propagate_uncertainty,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_core_opacity_recovers_heated_core_behind_cooler_environment() -> None:
    nu = 223.8 * u.GHz
    tcore, tenv = 50 * u.K, 20 * u.K
    tau_core, tau_half = 0.12, 0.015
    benv = blackbody_cgs(nu, tenv)
    background = benv * (1 - np.exp(-2 * tau_half))
    rear = benv * (1 - np.exp(-tau_half))
    peak = np.exp(-tau_half) * (1 - np.exp(-tau_core)) * (blackbody_cgs(nu, tcore) - rear)
    recovered = core_optical_depth(nu=nu, temp=tcore, Isr=peak, Ibg=background,
                                  environment_temp=tenv)
    np.testing.assert_allclose(recovered.value, tau_core, rtol=1e-12)


def test_explicit_equal_environment_preserves_old_opacity() -> None:
    kwargs = dict(nu=223.8*u.GHz, temp=20*u.K,
                  Isr=100*u.MJy/u.sr, Ibg=300*u.MJy/u.sr)
    np.testing.assert_allclose(core_optical_depth(**kwargs).value,
        core_optical_depth(**kwargs, environment_temp=kwargs["temp"]).value, rtol=1e-12)


def test_blackbody_increases_with_temperature() -> None:
    values = blackbody_cgs(
        223.8 * u.GHz,
        np.array([10.0, 20.0, 40.0]) * u.K,
    )

    assert np.all(np.diff(values.value) > 0)


def test_quantity_uncertainty_is_deterministic() -> None:
    parameters = {
        "value": {
            "dist": "normal",
            "mean": 10.0 * u.K,
            "std": 1.0 * u.K,
        }
    }
    first = propagate_uncertainty(
        lambda value: 2 * value,
        parameters,
        sample_count=1000,
        random_seed=42,
    )
    second = propagate_uncertainty(
        lambda value: 2 * value,
        parameters,
        sample_count=1000,
        random_seed=42,
    )

    assert first["mean"] == second["mean"]
    assert first["std"].unit == u.K


def test_physical_properties_inspection_reads_generated_catalogue() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/physical_properties.py",
            "--stage",
            "inspect",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    assert "rows=231" in completed.stdout
    assert "samples=10000, seed=42" in completed.stdout


def test_generated_physical_properties_schema_and_regression_values() -> None:
    generated = (
        PROJECT_ROOT
        / "results"
        / "physical_properties"
        / "tables"
        / "core_physical_properties.ecsv"
    )
    if not generated.is_file():
        pytest.skip("Generated physical-property catalogue is unavailable")
    actual = Table.read(generated, format="ascii.ecsv")

    assert actual.colnames == [
        "ID",
        "WCS_ACOOR",
        "WCS_DCOOR",
        "Tenv",
        "Subregion",
        "Resolvability",
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
    ]
    assert len(actual) == 231
    indices = [0, 115, 230]
    np.testing.assert_array_equal(actual["ID"][indices], [1, 116, 231])
    np.testing.assert_allclose(
        actual["Mass_core"][indices],
        [25.56380823008521, 1.4310267089414677, 0.7611526267089587],
        rtol=1e-12,
    )
    np.testing.assert_allclose(
        actual["nVolden_core"][indices],
        [3.3051022428170023, 0.01952763541595397, 0.013915990678481108],
        rtol=5e-9,  # Allow the CODATA proton-mass update in Astropy.
    )
    assert np.count_nonzero(actual["Resolvability"] == "N") == 22
    summary = Table.read(generated.with_name("core_property_summary.ecsv"))
    counts = {r["property"]: r["n"] for r in summary}
    assert counts["nVolden_core"] == counts["FWHM_linear"] == 209
    assert counts["Mass_core"] == counts["Surface_core"] == counts["Tau_core"] == 231
    np.testing.assert_array_equal(
        np.asarray(actual["Subregion"][indices], dtype=str),
        ["E", "I", "I"],
    )
