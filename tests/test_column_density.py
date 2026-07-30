"""Tests for dust column-density calculations."""

import numpy as np
import pytest

from lancet.column_density import (
    brightness_to_column_density,
    column_density_to_surface_density,
    dense_mass_msun,
    modified_blackbody_mjysr,
    planck_mjysr,
)


def test_planck_function_is_positive_and_increases_with_temperature() -> None:
    values = planck_mjysr(223.8, np.array([10.0, 20.0, 40.0]))

    assert np.all(values > 0)
    assert np.all(np.diff(values) > 0)


def test_brightness_conversion_inverts_modified_blackbody() -> None:
    expected = np.array([1.0e22, 1.0e23, 5.0e23])
    intensity = modified_blackbody_mjysr(
        expected,
        20.0,
        frequency_ghz=223.8,
        opacity_cm2_per_g=0.005,
    )

    actual = brightness_to_column_density(
        intensity,
        20.0,
        frequency_ghz=223.8,
        opacity_cm2_per_g=0.005,
    )

    np.testing.assert_allclose(actual, expected, rtol=2e-13)


def test_column_density_surface_density_conversion() -> None:
    surface_density = column_density_to_surface_density(1.0e23)

    assert surface_density == pytest.approx(0.4684, rel=2e-3)


def test_dense_mass_uses_pixel_area_and_threshold() -> None:
    surface_density = np.array([[0.5, 1.5], [2.0, np.nan]])
    selection = np.array([[True, True], [False, True]])

    mass = dense_mass_msun(
        surface_density,
        selection,
        pixel_scale_x_arcsec=1.0,
        pixel_scale_y_arcsec=1.0,
        distance_pc=1.0,
        threshold_g_cm2=1.0,
    )

    expected = 1.5 * (1.495978707e13**2) / 1.988409870698051e33
    assert mass == pytest.approx(expected, rel=1e-10)


def test_invalid_physical_parameters_are_rejected() -> None:
    with pytest.raises(ValueError, match="Frequency"):
        planck_mjysr(0.0, 20.0)
    with pytest.raises(ValueError, match="Opacity"):
        brightness_to_column_density(
            1.0,
            20.0,
            frequency_ghz=223.8,
            opacity_cm2_per_g=0.0,
        )
