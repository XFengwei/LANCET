"""Scientific averaging and bin-accounting checks for continuum diagnostics."""

import numpy as np

from lancet.combination.visibility_diagnostics import average_visibility_rows


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


def test_spectral_channel_split_preserves_vector_mean_and_total_weight():
    # Same integrated sensitivity represented as one channel or four channels.
    one = average_visibility_rows(np.array([[[2+3j]]]), np.zeros((1, 1, 1), bool),
                                  np.array([[[12.]]]), [225e9])
    four = average_visibility_rows(np.full((1, 4, 1), 2+3j), np.zeros((1, 4, 1), bool),
                                   np.full((1, 4, 1), 3.), [225e9]*4)
    for a, b in zip(one[:3], four[:3]):
        np.testing.assert_allclose(a, b)
