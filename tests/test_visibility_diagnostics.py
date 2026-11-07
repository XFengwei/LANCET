"""Scientific averaging and bin-accounting checks for continuum diagnostics."""

import numpy as np

from lancet.visibility_diagnostics import (
    average_visibility_rows,
    profile_columns,
    radial_visibility_sums,
)


def test_vector_average_precedes_amplitude_and_excludes_flagged_outlier():
    data = np.array([[[1 + 2j, 3 + 0j], [1 - 2j, 1000 + 0j]]])
    flags = np.array([[[False, False], [False, True]]])
    visibility, weight, frequency, count = average_visibility_rows(
        data, flags, np.ones((1, 1, 2)), [220e9, 230e9])
    np.testing.assert_allclose(visibility, [1, 3])
    np.testing.assert_allclose(weight, [2, 1])
    np.testing.assert_allclose(frequency, [225e9, 220e9])
    np.testing.assert_array_equal(count, [2, 1])


def test_invalid_weights_and_data_do_not_contaminate_good_channels():
    data = np.array([[[2 + 1j, complex(np.nan, 0)], [100, 4], [200, 6]]])
    weights = np.array([[[2, 5], [0, 3], [np.inf, -2]]])
    visibility, weight, frequency, count = average_visibility_rows(
        data, np.zeros(data.shape, bool), weights, [220e9, 225e9, 230e9])
    np.testing.assert_allclose(visibility, [2 + 1j, 4])
    np.testing.assert_allclose(weight, [2, 3])
    np.testing.assert_allclose(frequency, [220e9, 225e9])
    np.testing.assert_array_equal(count, [1, 1])


def test_radial_statistics_are_chunk_invariant_and_preserve_empty_bins():
    radius, visibility, weights = np.array([1., 2., 3.]), np.array([1., -3., 5.]), np.array([1., 3., 2.])
    edges = np.array([0., 2.5, 4., 5.])
    total = radial_visibility_sums(radius, visibility, weights, edges)
    chunks = sum(radial_visibility_sums(radius[s], visibility[s], weights[s], edges)
                 for s in [slice(0, 1), slice(1, 3)])
    np.testing.assert_allclose(total, chunks)
    values = profile_columns(total, edges)
    np.testing.assert_allclose(values['mean_amplitude_jy'][:2], [2.5, 5.])
    np.testing.assert_array_equal(values['averaged_row_count'], [2, 1, 0])
    assert np.isnan(values['mean_amplitude_jy'][2])
    assert values['sum_weight'].sum() == 6


def test_spectral_channel_split_preserves_vector_mean_and_total_weight():
    # Same integrated sensitivity represented as one channel or four channels.
    one = average_visibility_rows(np.array([[[2+3j]]]), np.zeros((1, 1, 1), bool),
                                  np.array([[[12.]]]), [225e9])
    four = average_visibility_rows(np.full((1, 4, 1), 2+3j), np.zeros((1, 4, 1), bool),
                                   np.full((1, 4, 1), 3.), [225e9]*4)
    for a, b in zip(one[:3], four[:3]):
        np.testing.assert_allclose(a, b)
