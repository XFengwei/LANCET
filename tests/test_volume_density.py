"""Tests for resolution-aware volume-density safeguards."""

import numpy as np

from lancet.volume_density import (
    apply_resolution_floor,
    cube_storage_bytes,
    density_from_width,
    number_density_h2,
)


def test_resolution_floor_replaces_zero_and_unresolved_widths() -> None:
    widths = np.array([0.0, 1.0e16, 1.0e18])
    resolutions = np.array([0.65, 18.0, 0.65])

    effective, defined, floored = apply_resolution_floor(
        widths,
        resolutions,
        distance_pc=2800.0,
    )

    assert np.array_equal(defined, [False, True, True])
    assert np.array_equal(floored, [True, True, False])
    assert np.all(np.isfinite(effective))
    assert np.all(effective > 0)


def test_density_keeps_undefined_width_as_nan() -> None:
    density = density_from_width(np.array([1.0, 1.0]), np.array([0.0, 1.0e18]))

    assert np.isnan(density[0])
    assert np.isfinite(density[1])


def test_number_density_and_cube_estimate() -> None:
    number = number_density_h2(np.array([2.8 * 1.67262192369e-24]))

    np.testing.assert_allclose(number, [1.0], rtol=1.0e-8)
    assert cube_storage_bytes((1500, 1500)) == 13_500_000_000

