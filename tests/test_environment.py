"""Tests for core-environment aperture measurements."""

from pathlib import Path
import subprocess
import sys

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from photutils.aperture import SkyCircularAperture
import numpy as np
import pytest

from lancet.environment import (
    add_environment_measurements,
    measure_aperture_statistics,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _header() -> fits.Header:
    return fits.Header(
        {
            "NAXIS": 2,
            "NAXIS1": 11,
            "NAXIS2": 11,
            "CTYPE1": "RA---TAN",
            "CTYPE2": "DEC--TAN",
            "CRPIX1": 6.0,
            "CRPIX2": 6.0,
            "CRVAL1": 221.0,
            "CRVAL2": -59.0,
            "CDELT1": -1.0 / 3600.0,
            "CDELT2": 1.0 / 3600.0,
        }
    )


def test_aperture_statistics_exclude_nonpositive_pixels() -> None:
    data = np.ones((11, 11), dtype=float)
    data[5, 5] = 5.0
    data[5, 6] = 0.0
    data[6, 5] = np.nan
    aperture = SkyCircularAperture(
        SkyCoord(221.0 * u.deg, -59.0 * u.deg),
        r=2.5 * u.arcsec,
    )

    result = measure_aperture_statistics(
        data,
        aperture,
        WCS(_header()),
    )

    assert result["mean"].shape == (1,)
    assert result["mean"][0] > 1.0
    assert result["median"][0] == pytest.approx(1.0)
    assert 0.0 < result["validfrac"][0] < 1.0


def test_measurements_are_appended_with_units() -> None:
    catalogue = Table({"ID": [1, 2]})
    values = {
        name: np.array([1.0, 2.0])
        for name in ("mean", "median", "std", "mad", "validfrac")
    }

    output = add_environment_measurements(
        catalogue,
        {"Temp": values},
        units={"Temp": "K"},
    )

    assert output["Temp_env_mean"].unit == u.K
    assert output["Temp_env_validfrac"].unit is None
    np.testing.assert_array_equal(output["ID"], [1, 2])


def test_environment_inspection_reports_current_dependencies() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/environment_measurements.py",
            "--stage",
            "inspect",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    assert "rows=231" in completed.stdout
    assert "c18o.mom2.fits.gz" in completed.stdout
    assert "radius=15.0 arcsec" in completed.stdout


def test_generated_environment_catalogue_schema_and_regression_values() -> None:
    generated = (
        PROJECT_ROOT
        / "results"
        / "environment_measurements"
        / "tables"
        / "core_environment_catalogue.ecsv"
    )
    if not generated.is_file():
        pytest.skip("Generated environment catalogue is unavailable")
    actual = Table.read(generated, format="ascii.ecsv")

    expected_environment_columns = [
        f"{prefix}_env_{statistic}"
        for prefix in ("Temp", "Vsig", "Interco", "Surface")
        for statistic in ("mean", "median", "std", "mad", "validfrac")
    ]
    assert len(actual) == 231
    assert actual.colnames[-20:] == expected_environment_columns
    indices = [0, 115, 230]
    np.testing.assert_array_equal(actual["ID"][indices], [1, 116, 231])
    np.testing.assert_allclose(
        actual["Temp_env_mean"][indices],
        [28.121459688459122, 18.486155101231166, 19.997605419158937],
        rtol=1e-12,
    )
    np.testing.assert_allclose(
        actual["Vsig_env_mean"][indices],
        [1.5996262516294206, 1.3432199358940125, 1.2145213782787323],
        rtol=1e-12,
    )
    np.testing.assert_allclose(
        actual["Interco_env_mean"][indices],
        [0.6610533325924676, 0.4524789945258619, 0.4896424765607729],
        rtol=1e-12,
    )
    np.testing.assert_allclose(
        actual["Surface_env_mean"][indices],
        [0.529169282742909, 1.0387836070287795, 0.8509821206331253],
        rtol=1e-12,
    )
