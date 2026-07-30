"""Tests for core-mass-function estimation and task wiring."""

from pathlib import Path
import subprocess
import sys

from astropy.table import Table
import numpy as np
import pytest

from lancet.core_mass_function import (
    diagnose_minimum_mass,
    differential_mass_function,
    powerlaw_alpha_mle,
    powerlaw_cdf,
    powerlaw_pdf,
    sample_alpha_posterior,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CATALOGUE = (
    PROJECT_ROOT
    / "results"
    / "physical_properties"
    / "tables"
    / "core_physical_properties.ecsv"
)


def test_powerlaw_mle_matches_closed_form() -> None:
    masses = np.array([1.0, 2.0, 4.0, 8.0])
    expected = 1.0 + masses.size / np.sum(np.log(masses))
    assert powerlaw_alpha_mle(masses, 1.0) == pytest.approx(expected)


def test_powerlaw_pdf_and_cdf_respect_lower_cutoff() -> None:
    masses = np.array([0.5, 1.0, 2.0, np.inf])
    pdf = powerlaw_pdf(masses, 2.0, 1.0)
    cdf = powerlaw_cdf(masses, 2.0, 1.0)
    np.testing.assert_allclose(pdf[:3], [0.0, 1.0, 0.25])
    np.testing.assert_allclose(cdf, [0.0, 0.0, 0.5, 1.0])


def test_minimum_mass_diagnosis_respects_minimum_count() -> None:
    masses = np.geomspace(1.0, 100.0, 40)
    diagnosis = diagnose_minimum_mass(masses, minimum_count=15)
    assert int(diagnosis["count"]) >= 15
    assert np.all(np.asarray(diagnosis["counts"]) >= 15)
    assert float(diagnosis["alpha"]) > 1


def test_posterior_sampler_is_reproducible() -> None:
    masses = np.geomspace(1.0, 30.0, 80)
    first = sample_alpha_posterior(
        masses,
        1.0,
        steps=1200,
        burn_in=200,
        seed=9,
    )
    second = sample_alpha_posterior(
        masses,
        1.0,
        steps=1200,
        burn_in=200,
        seed=9,
    )
    np.testing.assert_allclose(first["posterior"], second["posterior"])
    assert 0 < float(first["acceptance_fraction"]) < 1


def test_differential_mass_function_preserves_counts() -> None:
    result = differential_mass_function(
        [1.0, 1.5, 2.5, 3.5],
        [1.0, 2.0, 4.0],
    )
    np.testing.assert_array_equal(result["counts"], [2, 2])
    np.testing.assert_allclose(result["density"], [2.0, 1.0])


def test_current_catalogue_reproduces_opacity_consistent_global_fit() -> None:
    if not CATALOGUE.is_file():
        pytest.skip("Generated physical catalogue is unavailable")
    table = Table.read(CATALOGUE, format="ascii.ecsv")
    diagnosis = diagnose_minimum_mass(
        np.asarray(table["Mass_core"], dtype=float),
        minimum_count=30,
    )
    assert diagnosis["minimum_mass"] == pytest.approx(
        1.7506050248345124,
        rel=1e-10,
    )
    assert diagnosis["count"] == 92
    assert diagnosis["alpha"] == pytest.approx(
        2.4147247376365835,
        rel=1e-10,
    )


def test_core_mass_function_task_inspection_reports_current_catalogue() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/core_mass_function.py",
            "--stage",
            "inspect",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    assert "rows=231, valid positive masses=231" in completed.stdout
    assert "Mmin=1.75061 solMass" in completed.stdout
    assert "Ntail=92" in completed.stdout
    assert "Completeness floor: none" in completed.stdout
