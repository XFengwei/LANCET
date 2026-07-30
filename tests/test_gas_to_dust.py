"""Tests for the Giannetti et al. (2017) gas-to-dust relation."""

import json
from pathlib import Path
import subprocess
import sys

import pytest

from lancet.gas_to_dust import estimate_gas_to_dust_ratio


def test_g316_value_matches_the_reference_central_calculation() -> None:
    estimate = estimate_gas_to_dust_ratio(6.5, 0.2)

    assert estimate.log10_ratio == pytest.approx(2.0055)
    assert estimate.ratio == pytest.approx(10**2.0055)


def test_systematic_branches_match_the_published_paired_shifts() -> None:
    estimate = estimate_gas_to_dust_ratio(6.5, 0.2)

    assert estimate.systematic_log10_minus == pytest.approx(2.0055 - 1.848)
    assert estimate.systematic_log10_plus == pytest.approx(2.053 - 2.0055)
    assert estimate.ratio - estimate.systematic_ratio_minus == pytest.approx(
        10**1.848
    )
    assert estimate.ratio + estimate.systematic_ratio_plus == pytest.approx(
        10**2.053
    )


def test_statistical_propagation_includes_radius_uncertainty() -> None:
    estimate = estimate_gas_to_dust_ratio(6.5, 0.2)
    expected = ((6.5 * 0.007) ** 2 + (0.087 * 0.2) ** 2 + 0.03**2) ** 0.5

    assert estimate.statistical_log10_sigma == pytest.approx(expected)


def test_fit_range_is_enforced_unless_extrapolation_is_explicit() -> None:
    with pytest.raises(ValueError, match="2–20 kpc"):
        estimate_gas_to_dust_ratio(21.0)

    estimate = estimate_gas_to_dust_ratio(21.0, allow_extrapolation=True)
    assert estimate.radius_kpc == 21.0


def test_command_prints_machine_readable_result() -> None:
    project_root = Path(__file__).resolve().parents[1]
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/gas_to_dust_ratio.py",
            "--radius-kpc",
            "6.5",
            "--radius-error-kpc",
            "0.2",
            "--json",
        ],
        cwd=project_root,
        check=True,
        capture_output=True,
        text=True,
    )
    result = json.loads(completed.stdout)

    assert result["radius_kpc"] == 6.5
    assert result["ratio"] == pytest.approx(10**2.0055)
