import numpy as np
from astropy.io import fits

from lancet.aca_forward_test import binned_rms, pair_weights, scale_for_unit_ratio, sky_model


def test_scale_for_unit_ratio_interpolates_inside_range():
    assert np.isclose(scale_for_unit_ratio([0, 1, 2, 3], [0.4, 0.8, 1.2, 1.6]), 1.5)
    assert np.isnan(scale_for_unit_ratio([0, 1], [1.2, 1.5]))
    assert np.isnan(scale_for_unit_ratio([0, 1, 2], [0.5, 1.4, 1.1]))


def test_bootstrap_weights_repeat_pairs():
    assert np.array_equal(pair_weights(3, np.array([0, 0, 2])), [2, 0, 1])
    power = np.array([1.0, 4.0, 9.0])
    use = np.ones(3, bool)
    assert np.isclose(binned_rms(power, use, np.arange(3), pair_weights(3)), np.sqrt(14/3))
    assert np.isclose(binned_rms(power, use, np.arange(3), np.array([2.0, 0, 1])), np.sqrt(11/3))


def test_jy_per_beam_is_converted_to_mjy_per_sr(tmp_path):
    header = fits.Header({"BUNIT": "Jy/beam", "BMAJ": 1/3600, "BMIN": 1/3600, "CDELT1": -0.1/3600,
                          "CDELT2": 0.1/3600, "CRPIX1": 5.5, "CRPIX2": 5.5, "CTYPE1": "RA---SIN",
                          "CTYPE2": "DEC--SIN", "CRVAL1": 0.0, "CRVAL2": 0.0})
    path = tmp_path/"image.fits"
    fits.PrimaryHDU(np.ones((10, 10), np.float32), header).writeto(path)
    data, out = sky_model(path, 0.5)
    area = np.pi/(4*np.log(2))*(np.pi/648000)**2
    assert data.shape == (2, 2) and out["BUNIT"] == "MJy/sr"
    assert np.allclose(data, 1e-6/area, rtol=1e-6) and np.isclose(out["CDELT2"]*3600, 0.5)


def test_gaussian_transfer_and_limit_are_consistent():
    from lancet.amplitude_alignment import gaussian_transfer, response_limit
    limit = response_limit(18.1, 0.1)
    assert np.isclose(gaussian_transfer(limit, 18.1), 0.1)
    half = response_limit(18.1, 0.5)
    assert np.isclose(half, 5.0287, rtol=1e-3)  # matches the TP summary half-response radius
