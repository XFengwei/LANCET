import numpy as np
from scipy import fft

from lancet.combination.uvplane import gaussian_transfer, uv_distance_klambda
from lancet.combination.prep import mjysr_to_jybeam
from lancet.combination.validation import (
    beam_volume_outside,
    large_scale_ratio,
    smooth_to,
    synthetic_gate,
)


def _smooth_sky(shape, pix, seed=0):
    rng = np.random.default_rng(seed)
    u, v, _ = uv_distance_klambda(shape, pix)
    field = fft.irfft2(fft.rfft2(rng.random(shape)) * gaussian_transfer(u, v, 30.0, 30.0, 0.0), s=shape)
    return field - field.min() + 1.0


def test_large_scale_ratio_of_scaled_tp():
    shape, pix = (300, 300), 1.0
    sky = _smooth_sky(shape, pix)
    u, v, _ = uv_distance_klambda(shape, pix)
    high_beam = (2.0, 1.5, 20.0)
    product = mjysr_to_jybeam(1.03 * fft.irfft2(fft.rfft2(sky) * gaussian_transfer(u, v, *high_beam), s=shape),
                              2.0, 1.5)
    tp = fft.irfft2(fft.rfft2(sky) * gaussian_transfer(u, v, 18.1, 18.1, 0.0), s=shape)
    mask = np.zeros(shape, dtype=bool)
    mask[100:200, 100:200] = True
    result = large_scale_ratio(product, tp, pix, high_beam, 18.1, mask)
    assert abs(result["median"] - 1.03) < 0.005


def test_beam_volume_outside_jcomb_zero():
    shape = (128, 128)
    u, v, _ = uv_distance_klambda(shape, 0.5)
    high = gaussian_transfer(u, v, 2.0, 1.5, 20.0)
    assert beam_volume_outside(high, high, shape) < 1e-6
    feather = gaussian_transfer(u, v, 18.1, 18.1, 0.0)
    assert beam_volume_outside(feather, high, shape) > 0.1


def test_synthetic_gate_passes_on_small_grid():
    shape, pix = (400, 400), 0.5
    y, x = np.indices(shape)
    gain = ((x - 200) ** 2 + (y - 200) ** 2 <= 120**2).astype(np.float32)
    settings = {"k1": 5.0, "k2": 9.0, "order": 5, "tp_factor": 1.0, "minimum_tp_response": 0.1,
                "padding_arcsec": 60.0}
    result = synthetic_gate(gain, pix, (2.0, 1.5, 20.0), (18.1, 18.1, 0.0), settings,
                            aca_shortest_klambda=6.6, min_gain=0.5)
    assert result["flux_error"] < 0.02
    assert result["rms_fraction_of_peak"] < 0.01
    assert result["passed"] is True


def _gaussians(shape, pix, positions, fwhm_arcsec, amplitude=1.0):
    yy, xx = np.indices(shape)
    out = np.zeros(shape)
    for x, y in positions:
        out += amplitude * np.exp(-4 * np.log(2) * ((xx - x) ** 2 + (yy - y) ** 2) / (fwhm_arcsec / pix) ** 2)
    return out


def test_tp_at_high_resolution_recovers_high_beam_sky():
    from lancet.combination.validation import tp_at_high_resolution

    shape, pix = (300, 360), 1.0
    rng = np.random.default_rng(4)
    # 40" sources: negligible power at the 9 klambda cut, flat towards the edges
    sky = 1.0 + _gaussians(shape, pix, rng.uniform([150, 120], [210, 180], size=(6, 2)), 40.0)
    u, v, k = uv_distance_klambda(shape, pix)
    high_beam, low_beam = (2.0, 1.5, 20.0), (18.1, 18.1, 0.0)
    low = fft.irfft2(fft.rfft2(sky) * gaussian_transfer(u, v, *low_beam), s=shape)
    out = tp_at_high_resolution(low, pix, low_beam, high_beam, k_max=9.0, padding_arcsec=60.0)
    expected = fft.irfft2(fft.rfft2(sky) * gaussian_transfer(u, v, *high_beam) * (k <= 9.0), s=shape)
    inner = (slice(60, -60), slice(60, -60))
    assert out.shape == shape
    assert np.max(np.abs(out[inner] - expected[inner])) < 1e-3 * np.abs(expected).max()


def test_cross_spectrum_slopes_bracket_the_scale():
    from lancet.combination.validation import cross_spectrum_slopes

    shape, pix = (256, 256), 1.0
    b = _smooth_sky(shape, 4.0, seed=5)          # 7.5" structure: signal in every bin
    window = np.ones(shape)
    edges = [3.0, 5.0, 7.0, 9.0]
    s_ab, s_aa = cross_spectrum_slopes(1.25 * b, b, window, pix, edges)
    np.testing.assert_allclose(s_ab, 1.25, rtol=1e-6)
    np.testing.assert_allclose(s_aa, 1.25, rtol=1e-6)
    noisy = b + 0.05 * np.random.default_rng(6).standard_normal(shape)
    s_ab, s_aa = cross_spectrum_slopes(1.25 * b, noisy, window, pix, edges)
    assert np.all(s_ab < 1.25)                    # noise in the reference biases this one low
    np.testing.assert_allclose(s_aa, 1.25, rtol=0.05)   # unbiased on average; one noise draw


def test_windowed_amplitude_ratio_identity_and_scale():
    from lancet.combination.validation import windowed_amplitude_ratio

    b = _smooth_sky((128, 128), 1.0, seed=7)
    window = np.ones(b.shape)
    np.testing.assert_allclose(windowed_amplitude_ratio(b, b, window, 1.0, [3.0, 6.0, 9.0]), 1.0)
    np.testing.assert_allclose(windowed_amplitude_ratio(0.9 * b, b, window, 1.0, [3.0, 6.0, 9.0]), 0.9)


def test_mosaic_window_is_smooth_and_bounded():
    from lancet.combination.validation import mosaic_window

    gain = np.zeros((200, 200))
    gain[50:150, 50:150] = 0.8
    w = mosaic_window(gain, 0.5, 1.0, 18.1)
    assert w.min() >= 0.0 and w.max() <= 1.0
    assert w[100, 100] > 0.99 and w[5, 5] < 1e-3
    assert 0.3 < w[100, 50] < 0.7


def test_core_contrast_ignores_smooth_additions():
    from lancet.combination.validation import core_contrast_ratio

    shape, pix = (300, 300), 0.5
    positions = [(80, 90), (150, 160), (220, 200)]
    cores = _gaussians(shape, pix, positions, 1.5, amplitude=5e-3)
    smooth = 2e-3 * _smooth_sky(shape, 0.25, seed=8)   # 60" structure on the 0.5" grid
    same = core_contrast_ratio(cores + smooth, cores, positions, pix, annulus_arcsec=(3.0, 6.0))
    assert abs(same["median"] - 1.0) < 0.01 and same["n"] == 3
    scaled = core_contrast_ratio(1.2 * cores, cores, positions, pix, annulus_arcsec=(3.0, 6.0))
    assert abs(scaled["median"] - 1.2) < 1e-6


def test_smooth_to_does_not_wrap_around_the_grid():
    shape, pix = (120, 400), 1.0
    ramp = np.tile(np.linspace(0.0, 1.0, shape[1]), (shape[0], 1))
    out = smooth_to(ramp, pix, (2.0, 2.0, 0.0), 18.1)
    assert np.max(np.abs(out[:, 0] - ramp[:, 0])) < 0.05     # wrapping pulls this towards 0.5
    assert np.max(np.abs(out[:, -1] - ramp[:, -1])) < 0.05


def test_core_contrast_with_no_positions():
    from lancet.combination.validation import core_contrast_ratio

    r = core_contrast_ratio(np.ones((50, 50)), np.ones((50, 50)), [], 0.5)
    assert r["n"] == 0 and np.isnan(r["median"])


def test_tp_only_ratio_percentiles():
    from lancet.combination.validation import tp_only_ratio

    tp = _smooth_sky((80, 90), 1.0, seed=9)
    region = np.zeros(tp.shape, dtype=bool)
    region[10:70, 10:80] = True
    r = tp_only_ratio(mjysr_to_jybeam(1.02 * tp, 0.7, 0.6), tp, (0.7, 0.6, 30.0), region)
    assert abs(r["median"] - 1.02) < 1e-6 and r["n"] == region.sum()
    assert set(r) >= {"p01", "p16", "median", "p84", "p99"}


def test_tp_only_ratio_empty_region_returns_none_values():
    from lancet.combination.validation import tp_only_ratio

    tp = np.ones((20, 20))
    product = np.full((20, 20), np.nan)                 # product blank outside the mosaic
    r = tp_only_ratio(product, tp, (0.7, 0.6, 30.0), np.ones((20, 20), dtype=bool))
    assert r == {"p01": None, "p16": None, "median": None, "p84": None, "p99": None, "n": 0}


def test_tapered_window_is_one_inside_and_zero_at_the_edge():
    from lancet.combination.validation import tapered_window

    mask = np.zeros((60, 80), dtype=bool)
    mask[10:50, 10:70] = True
    w = tapered_window(mask, 1.0, 8.0)
    assert w[~mask].max() == 0.0 and w[30, 40] == 1.0
    assert 0.0 < w[10, 40] < 0.1 and np.all(np.diff(w[10:17, 40]) > 0)     # rises from the edge inwards
    full = tapered_window(np.ones((20, 20), dtype=bool), 1.0, 4.0)
    assert full[0, 10] < 0.2 and full[10, 10] == 1.0                       # the grid edge counts as an edge


def test_windowed_power_scales_with_amplitude_squared_and_ignores_the_mean():
    from lancet.combination.validation import windowed_power

    sky = _smooth_sky((128, 160), 2.0, seed=4)
    window = np.ones(sky.shape)
    edges = [2.0, 4.0, 8.0, 16.0]
    p = windowed_power(sky, window, 2.0, edges)
    assert np.allclose(windowed_power(1.4 * sky, window, 2.0, edges), 1.96 * p)
    assert np.allclose(windowed_power(sky + 50.0, window, 2.0, edges), p)
    half = np.zeros(sky.shape)
    half[:, :80] = 1.0
    tiled = np.hstack([sky[:, :80], sky[:, :80]])                          # same sky in both halves
    assert np.allclose(windowed_power(tiled, half, 2.0, edges), windowed_power(tiled, window, 2.0, edges), rtol=0.5)


def test_psd_match_factor_takes_the_largest_ratio_inside_the_band():
    import pytest

    from lancet.combination.validation import psd_match_factor

    centres = np.array([2.0, 4.0, 6.0, 8.0, 12.0])
    p_tp = np.ones(5)
    p_ref = np.array([9.0, 1.0, 1.96, 1.21, 25.0])       # amplitude ratios 3, 1, 1.4, 1.1, 5
    r = psd_match_factor(p_ref, p_tp, centres, [4.0, 10.0])
    assert r["factor"] == pytest.approx(1.4) and r["at_klambda"] == 6.0
    assert r["bins_klambda"] == [4.0, 6.0, 8.0]
    with pytest.raises(ValueError, match="No uv bin"):
        psd_match_factor(p_ref, p_tp, centres, [20.0, 30.0])



def test_zero_spacing_check_reports_flux_and_pixel_ratios():
    import pytest

    from lancet.combination.validation import zero_spacing_check

    tp = np.full((20, 30), 10.0)
    tp[:, 15:] = 30.0
    combined = tp.copy()
    combined[:, 15:] *= 1.1                      # 10% high in the bright half
    combined[0, 0] = np.nan
    regions = {"left": np.arange(30)[None, :] < 15 + np.zeros((20, 1), bool), "all": np.ones(tp.shape, bool),
               "empty": np.zeros(tp.shape, bool)}
    out = zero_spacing_check(combined, tp, regions)
    assert out["left"]["flux_ratio"] == pytest.approx(1.0) and out["left"]["n"] == 299   # the NaN pixel is left out
    assert out["all"]["pixel_ratio"]["median"] == pytest.approx(1.1, abs=0.06)
    assert out["all"]["flux_ratio"] == pytest.approx((299 * 10 + 300 * 33) / (299 * 10 + 300 * 30))
    assert out["empty"] == {"n": 0}
