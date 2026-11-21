import numpy as np
import pytest
from scipy import fft

from lancet.combination.uvplane import (
    butterworth_taper,
    combine_plane,
    feather_weights,
    gaussian_transfer,
    jcomb_weights,
    uv_distance_klambda,
)

ARCSEC = np.pi / 180 / 3600
LOW_BEAM = (18.1, 18.1, 0.0)


def _convolve(sky, pixel_arcsec, beam):
    u, v, _ = uv_distance_klambda(sky.shape, pixel_arcsec)
    return fft.irfft2(fft.rfft2(sky) * gaussian_transfer(u, v, *beam), s=sky.shape)


def _sky(shape, seed=1):
    rng = np.random.default_rng(seed)
    return rng.random(shape) + 0.5


def test_uv_grid_matches_fftfreq():
    u, v, k = uv_distance_klambda((64, 81), 0.5)
    d = 0.5 * np.pi / 648000.0
    assert u.shape == (64, 41)
    np.testing.assert_allclose(u[0, :], np.fft.rfftfreq(81, d) / 1e3, rtol=1e-12)
    np.testing.assert_allclose(v[:, 0], np.fft.fftfreq(64, d) / 1e3, rtol=1e-12)
    np.testing.assert_allclose(k, np.hypot(u, v))


def _image_gaussian(shape, pixel_arcsec, bmaj, bmin, bpa_deg):
    """Major axis rotated from +y toward -x (east, CDELT1 < 0) by BPA."""
    ny, nx = shape
    y = (np.arange(ny) + ny // 2) % ny - ny // 2
    x = (np.arange(nx) + nx // 2) % nx - nx // 2
    xx, yy = np.meshgrid(x, y)
    theta = np.deg2rad(bpa_deg)
    major = -xx * np.sin(theta) + yy * np.cos(theta)
    minor = xx * np.cos(theta) + yy * np.sin(theta)
    to_sigma = 1 / np.sqrt(8 * np.log(2)) / pixel_arcsec
    return np.exp(-0.5 * ((major / (bmaj * to_sigma)) ** 2 + (minor / (bmin * to_sigma)) ** 2))


def test_gaussian_transfer_matches_numerical_fft():
    shape, pix = (128, 129), 0.5
    numerical = np.abs(fft.rfft2(_image_gaussian(shape, pix, 6.0, 3.0, 30.0)))
    numerical /= numerical[0, 0]
    u, v, _ = uv_distance_klambda(shape, pix)
    assert np.max(np.abs(numerical - gaussian_transfer(u, v, 6.0, 3.0, 30.0))) < 1e-3
    assert np.max(np.abs(numerical - gaussian_transfer(u, v, 6.0, 3.0, -30.0))) > 0.05


def test_butterworth_taper_limits():
    k = np.linspace(0, 12, 1201)
    w = butterworth_taper(k, 5.0, 9.0, 5)
    assert np.all(w[k <= 5] == 1.0)
    assert np.all(w[k >= 9] == 0.0)
    assert np.all(np.diff(w) <= 0)
    assert np.max(np.abs(np.diff(w))) < 0.02
    with pytest.raises(ValueError):
        butterworth_taper(k, 9.0, 5.0, 5)
    with pytest.raises(ValueError):
        butterworth_taper(k, 5.0, 9.0, 0)


def test_guard_refuses_weak_tp_response():
    u, v, _ = uv_distance_klambda((200, 200), 0.5)
    high = (0.71, 0.61, 47.4)
    with pytest.raises(ValueError):
        jcomb_weights(u, v, LOW_BEAM, high, 5.0, 12.0, 5)
    w, r = jcomb_weights(u, v, LOW_BEAM, high, 5.0, 12.0, 5, minimum_tp_response=0.0)
    assert w.shape == r.shape == u.shape


def test_fully_sampled_identity_jcomb():
    shape, pix = (255, 301), 0.5
    sky = _sky(shape)
    high_beam = (2.0, 1.5, 20.0)
    low = _convolve(sky, pix, LOW_BEAM)
    high = _convolve(sky, pix, high_beam)
    u, v, _ = uv_distance_klambda(shape, pix)
    w, r = jcomb_weights(u, v, LOW_BEAM, high_beam, 5.0, 9.0, 5)
    out = combine_plane(low, high, w, r)
    assert out.shape == shape
    assert out.dtype == np.float32
    assert np.max(np.abs(out - high)) < 1e-4 * high.max()


def test_feather_keeps_tp_total_flux():
    shape, pix = (255, 301), 0.5
    sky = _sky(shape)
    low = _convolve(sky, pix, LOW_BEAM)
    high = _convolve(sky, pix, (2.0, 1.5, 20.0))
    high = high - high.mean()  # remove the k = 0 term, as an interferometer does
    u, v, _ = uv_distance_klambda(shape, pix)
    w, r = feather_weights(u, v, LOW_BEAM)
    out = combine_plane(low, high, w, r)
    assert abs(out.sum() - low.sum()) < 1e-4 * abs(low.sum())


def test_tp_factor_scales_low_k_only():
    shape, pix = (255, 301), 0.5
    sky = _sky(shape)
    high_beam = (2.0, 1.5, 20.0)
    low = _convolve(sky, pix, LOW_BEAM)
    high = _convolve(sky, pix, high_beam)
    u, v, _ = uv_distance_klambda(shape, pix)
    w, r = jcomb_weights(u, v, LOW_BEAM, high_beam, 5.0, 9.0, 5)
    diff = fft.rfft2(combine_plane(low, high, w, r, 2.0).astype(np.float64)
                     - combine_plane(low, high, w, r, 1.0).astype(np.float64))
    expected = fft.rfft2(low) * r * w
    scale = np.abs(fft.rfft2(low)).max()
    assert np.max(np.abs(diff[w == 0])) < 1e-5 * scale
    assert np.max(np.abs(diff - expected)) < 1e-5 * scale


def test_feather_fully_sampled_transfer():
    """Feather adds the TP's FT unchanged and down-weights the interferometer by 1 - G_low."""
    shape, pix = (255, 301), 0.5
    sky = _sky(shape)
    high_beam = (2.0, 1.5, 20.0)
    low = _convolve(sky, pix, LOW_BEAM)
    high = _convolve(sky, pix, high_beam)
    u, v, _ = uv_distance_klambda(shape, pix)
    w, r = feather_weights(u, v, LOW_BEAM)
    out = combine_plane(low, high, w, r)
    g_low, g_high = gaussian_transfer(u, v, *LOW_BEAM), gaussian_transfer(u, v, *high_beam)
    expected = fft.irfft2(fft.rfft2(sky) * (g_low + (1 - g_low) * g_high), s=shape)
    assert np.max(np.abs(out - expected)) < 1e-4 * expected.max()
