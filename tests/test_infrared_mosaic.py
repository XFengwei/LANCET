"""Tests for the VVV/infrared mosaicking helpers."""

import numpy as np
import pytest

from lancet.infrared_mosaic import (
    block_any,
    block_median,
    boost_saturation,
    feather_weight,
    fill_saturated_cores,
    solve_offset_surfaces,
    stretch_channel,
    upsample_surface,
)


def test_block_median_blanks_partly_covered_cells() -> None:
    image = np.arange(16, dtype=np.float32).reshape(4, 4)
    image[3, 3] = np.nan

    medians = block_median(image, 2)

    assert medians.shape == (2, 2)
    assert medians[0, 0] == pytest.approx(np.median([0, 1, 4, 5]))
    assert np.isnan(medians[1, 1])


def test_block_median_pads_incomplete_edge_cells() -> None:
    medians = block_median(np.ones((3, 3), np.float32), 2)

    assert medians.shape == (2, 2)
    assert medians[0, 0] == 1.0
    assert np.isnan(medians[1, 1])


def test_block_any_marks_cells_with_any_true_pixel() -> None:
    mask = np.zeros((4, 4), bool)
    mask[3, 0] = True

    cells = block_any(mask, 2)

    np.testing.assert_array_equal(cells, [[False, False], [True, False]])


def test_offset_surfaces_remove_a_constant_tile_offset() -> None:
    sky = np.full((4, 6), 10.0)
    first = sky.copy()
    second = sky + 3.0
    domain_first = np.zeros((4, 6), bool)
    domain_first[:, :4] = True
    domain_second = np.zeros((4, 6), bool)
    domain_second[:, 2:] = True

    surfaces = solve_offset_surfaces(
        [np.where(domain_first, first, np.nan), np.where(domain_second, second, np.nan)],
        [domain_first, domain_second],
        smoothness=4.0,
        decay_cells=1e6,
    )

    overlap = domain_first & domain_second
    difference = surfaces[0][overlap] - surfaces[1][overlap]
    np.testing.assert_allclose(difference, -3.0, atol=1e-6)


def test_upsample_surface_fills_missing_cells_and_matches_shape() -> None:
    surface = np.array([[1.0, np.nan], [1.0, 1.0]])

    full = upsample_surface(surface, 3, (5, 6))

    assert full.shape == (5, 6)
    np.testing.assert_allclose(full, 1.0)


def test_feather_weight_rises_from_edge_to_one() -> None:
    footprint = np.ones((21, 21), bool)

    weight = feather_weight(footprint, 5.0)

    assert weight[0, 10] < 0.2
    assert weight[10, 10] == pytest.approx(1.0)
    assert np.all((weight >= 0) & (weight <= 1))


def test_fill_saturated_cores_fills_enclosed_low_pixels() -> None:
    image = np.zeros((5, 5), np.float32)
    image[1:4, 1:4] = 10.0
    image[2, 2] = 1.0

    filled = fill_saturated_cores(image, 5.0)

    assert filled[2, 2] == 5.0
    assert filled[0, 0] == 0.0


def test_stretch_channel_maps_to_unit_interval() -> None:
    rng = np.random.default_rng(1)
    image = rng.normal(100.0, 1.0, (50, 50)).astype(np.float32)
    image[25, 25] = 1000.0
    image[0, 0] = np.nan

    scaled = stretch_channel(image, 99.5, 0.05, 3.0)

    assert scaled.dtype == np.float32
    assert np.all(np.isfinite(scaled))
    assert scaled.min() >= 0.0 and scaled.max() == pytest.approx(1.0)


def test_boost_saturation_keeps_grey_and_clips() -> None:
    rgb = np.array([[[0.5, 0.5, 0.5], [0.9, 0.5, 0.1]]])

    boosted = boost_saturation(rgb, 2.0)

    np.testing.assert_allclose(boosted[0, 0], 0.5)
    np.testing.assert_allclose(boosted[0, 1], [1.0, 0.5, 0.0])
