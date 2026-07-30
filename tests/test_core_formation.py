"""Tests for the core-formation analysis."""

from pathlib import Path
import subprocess
import sys

from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import pytest

from lancet.core_formation import (
    add_catalogue_ellipses,
    counting_fractional_uncertainty,
    create_core_overlay_figure,
    dust_mass_opacity_scale,
    summarize_core_formation,
)


matplotlib.use("Agg")
PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _region_settings() -> dict[str, dict[str, object]]:
    return {
        "young": {
            "label": "Young",
            "code": "Y",
            "cloud_mass_msun": 100.0,
            "paper1_fragment_mass_msun": 10.0,
            "fragment_count": 3,
        },
        "evolved": {
            "label": "Evolved",
            "code": "E",
            "cloud_mass_msun": 200.0,
            "paper1_fragment_mass_msun": 40.0,
            "fragment_count": 4,
        },
    }


def _header() -> fits.Header:
    return fits.Header(
        {
            "NAXIS": 2,
            "NAXIS1": 41,
            "NAXIS2": 41,
            "CTYPE1": "RA---TAN",
            "CTYPE2": "DEC--TAN",
            "CRPIX1": 21.0,
            "CRPIX2": 21.0,
            "CRVAL1": 221.0,
            "CRVAL2": -59.0,
            "CDELT1": -1.0 / 3600.0,
            "CDELT2": 1.0 / 3600.0,
            "BUNIT": "MJy/sr",
        }
    )


def test_counting_uncertainty_handles_cloud_and_sources() -> None:
    assert counting_fractional_uncertainty(np.inf) == 0.0
    assert counting_fractional_uncertainty(5) == pytest.approx(0.25)
    with pytest.raises(ValueError, match="greater than one"):
        counting_fractional_uncertainty(1)


def test_dust_mass_opacity_scale_is_inverse() -> None:
    assert dust_mass_opacity_scale(0.9, 0.5) == pytest.approx(1.8)
    with pytest.raises(ValueError, match="positive"):
        dust_mass_opacity_scale(0.9, 0.0)


def test_summary_uses_generated_core_masses_and_counts() -> None:
    catalogue = Table(
        {
            "Subregion": ["Y", "Y", "E", "E", "E"],
            "Mass_core": [1.0, 2.0, 4.0, 5.0, 6.0] * u.M_sun,
            "err_Mass_core": [0.1, 0.2, 0.3, 0.4, 0.5] * u.M_sun,
        }
    )

    summary = summarize_core_formation(
        catalogue,
        _region_settings(),
        reference_opacity_cm2_per_g=0.9,
        adopted_opacity_cm2_per_g=0.5,
    )

    assert list(summary["region"]) == ["young", "evolved"]
    assert list(summary["core_count"]) == [2, 3]
    np.testing.assert_allclose(summary["core_mass"], [3.0, 15.0])
    np.testing.assert_allclose(
        summary["core_mass_error"],
        [np.hypot(0.1, 0.2), np.sqrt(0.3**2 + 0.4**2 + 0.5**2)],
    )
    np.testing.assert_allclose(
        summary["paper1_fragment_mass"],
        [10.0, 40.0],
    )
    np.testing.assert_allclose(
        summary["fragment_mass"],
        [18.0, 72.0],
    )
    np.testing.assert_allclose(
        summary["fragment_efficiency"],
        [0.18, 0.36],
    )
    np.testing.assert_allclose(
        summary["core_efficiency"],
        [0.03, 0.075],
    )
    np.testing.assert_allclose(
        summary["fragment_to_core_efficiency"],
        [3.0 / 18.0, 15.0 / 72.0],
    )
    np.testing.assert_allclose(
        summary["fragment_to_core_efficiency_core_mass_error"],
        [
            np.hypot(0.1, 0.2) / 18.0,
            np.sqrt(0.3**2 + 0.4**2 + 0.5**2) / 72.0,
        ],
    )


def test_catalogue_ellipses_respect_astronomical_angle() -> None:
    table = Table(
        {
            "ra": [221.0],
            "dec": [-59.0],
            "major": [2.0],
            "minor": [1.0],
            "pa": [30.0],
        }
    )
    figure = plt.figure()
    axis = figure.add_subplot(projection=WCS(_header()).celestial)

    artists = add_catalogue_ellipses(
        axis,
        WCS(_header()).celestial,
        table,
        ra_column="ra",
        dec_column="dec",
        major_column="major",
        minor_column="minor",
        position_angle_column="pa",
        position_angle_type="astronomical",
    )

    assert len(artists) == 1
    assert artists[0].angle == pytest.approx(120.0)
    assert artists[0].width == pytest.approx(2.0)
    plt.close(figure)


def test_overlay_figure_can_be_rendered(tmp_path: Path) -> None:
    data = np.linspace(0.0, 10.0, 41 * 41).reshape(41, 41)
    sources = Table(
        {
            "WCS_ACOOR": [221.0],
            "WCS_DCOOR": [-59.0],
            "AFWHM01": [2.0],
            "BFWHM01": [1.0],
            "THETA01": [20.0],
        }
    )
    fragments = Table(
        {
            "x_cen": [221.0],
            "y_cen": [-59.0],
            "major_sigma": [3.0],
            "minor_sigma": [2.0],
            "position_angle": [40.0],
        }
    )
    figure = create_core_overlay_figure(
        data,
        _header(),
        data,
        _header(),
        sources,
        fragments,
        center_ra_deg=221.0,
        center_dec_deg=-59.0,
        width_arcsec=30.0,
        height_arcsec=30.0,
        label="Test",
        image_vmin=0.0,
        image_vmax=10.0,
        contour_levels=(2.0, 5.0),
    )
    output = tmp_path / "overlay.png"

    figure.savefig(output, dpi=50)
    plt.close(figure)

    assert output.is_file()
    assert output.stat().st_size > 0


def test_core_formation_inspection_reports_optional_reference() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/core_formation.py",
            "--stage",
            "inspect",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    assert "physical_catalogue" in completed.stdout
    assert "rows=231" in completed.stdout
    assert "not used in calculations" in completed.stdout
    assert "Young: cores=16" in completed.stdout
