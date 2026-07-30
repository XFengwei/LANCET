"""Tests for the background-emission sensitivity workflow."""

from pathlib import Path
import subprocess
import sys

from astropy import units as u
from astropy.table import Table
import numpy as np
import pytest

from lancet.physical_properties import (
    symmetric_environment_optical_depth,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_symmetric_environment_opacity_starts_at_zero_and_increases() -> None:
    background = np.array([0.0, 100.0, 1000.0]) * u.MJy / u.sr
    result = symmetric_environment_optical_depth(
        nu=223.8 * u.GHz,
        temp=20.0 * u.K,
        Ibg=background,
    )
    values = result.to_value(u.dimensionless_unscaled)
    assert values[0] == pytest.approx(0.0)
    assert np.all(np.diff(values) > 0)


def test_background_effects_inspection_is_ready() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/background_effects.py",
            "--stage",
            "inspect",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    assert "frequency=223.8 GHz, temperature=20.0 K" in completed.stdout
    assert "1, 10, 100" in completed.stdout


def test_generated_background_table_has_finite_corrections() -> None:
    path = (
        PROJECT_ROOT
        / "results"
        / "background_effects"
        / "tables"
        / "background_effects.ecsv"
    )
    if not path.is_file():
        pytest.skip("Generated background-effects table is unavailable")
    table = Table.read(path, format="ascii.ecsv")
    assert len(table) == 50
    assert table.meta["normalization"] == (
        "zero-background core optical depth"
    )
    for name in table.colnames:
        assert np.all(np.isfinite(np.asarray(table[name], dtype=float)))
