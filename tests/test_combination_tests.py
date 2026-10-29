"""Controls for empirical noise, time grouping and common-sky Fourier tests."""

import numpy as np
from astropy.io import fits

from lancet.combination_tests import (
    adjacent_time_pairs, empirical_noise, radial_test_profile,
    sample_image_visibilities, matching_cells, common_cell_ratios,
)


def _rows(n):
    return {"x": np.ones(n, complex), "y": np.ones(n, complex),
            "wx": np.ones(n), "wy": np.ones(n), "frequency": np.full(n, 225e9),
            "uvw": np.tile([20., 0., 0.], (n, 1)), "time": np.arange(n, dtype=float),
            **{k: np.zeros(n, int) for k in ["field", "ddid", "scan", "observation", "antenna1"]},
            "antenna2": np.ones(n, int)}


def test_time_pairs_never_cross_scan_or_observation_and_do_not_overlap():
    r = _rows(9)
    r["scan"] = np.array([0, 0, 0, 1, 1, 1, 1, 1, 1])
    r["observation"] = np.array([0, 0, 0, 0, 0, 0, 1, 1, 1])
    i, j = adjacent_time_pairs(r, 3.)
    np.testing.assert_array_equal(i, [0, 3, 6])
    np.testing.assert_array_equal(j, [1, 4, 7])
    assert len(np.unique(np.r_[i, j])) == 2 * len(i)


def test_empirical_noise_and_cross_power_recover_known_signal():
    r = _rows(20000)
    rng = np.random.default_rng(123)
    sigma = 0.2
    for key in ["x", "y"]:
        r[key] = 0.3 + sigma * (rng.normal(size=20000) + 1j * rng.normal(size=20000))
    r["wx"][:] = r["wy"][:] = 1 / sigma**2
    kappa, estimates, _ = empirical_noise(r, 3.)
    assert abs(estimates['0']['kappa_time'] - 1) < 0.05
    assert abs(estimates['0']['kappa_hand_difference'] - 1) < 0.05
    p = radial_test_profile(r, [0., 100.], kappa)
    np.testing.assert_allclose(p['debiased_power_jy2'], [0.3**2], rtol=0.03)
    np.testing.assert_allclose(p['cross_power_jy2'], [0.3**2], rtol=0.03)
    assert p['raw_power_jy2'][0] > p['debiased_power_jy2'][0]


def test_forward_gaussian_sky_and_pb_match_analytic_fourier_transform():
    n, pixel_arcsec = 512, 0.25
    x = (np.arange(n) - n//2) * pixel_arcsec
    xx, yy = np.meshgrid(x, x)
    source_fwhm, beam_fwhm, peak = 4., 28., 100.
    image = peak * np.exp(-4*np.log(2)*(xx**2 + yy**2)/source_fwhm**2)
    header = fits.Header({'NAXIS': 2, 'NAXIS1': n, 'NAXIS2': n,
        'CTYPE1': 'RA---SIN', 'CTYPE2': 'DEC--SIN', 'CRPIX1': n//2+1,
        'CRPIX2': n//2+1, 'CRVAL1': 221.25, 'CRVAL2': -59.8,
        'CDELT1': -pixel_arcsec/3600, 'CDELT2': pixel_arcsec/3600})
    u, v = np.array([0., 8000., 20000., 30000.]), np.array([0., 5000., 6000., -10000.])
    actual = sample_image_visibilities(image, header, 221.25, -59.8, u, v,
        fwhm_arcsec=beam_fwhm, pixel_arcsec=0.5, patch_arcsec=128.)
    theta = 1/np.sqrt(source_fwhm**-2 + beam_fwhm**-2) * np.pi/648000
    expected = peak * 1e6 * np.pi*theta**2/(4*np.log(2)) * np.exp(
        -np.pi**2*theta**2*(u*u+v*v)/(4*np.log(2)))
    np.testing.assert_allclose(actual.real, expected, rtol=0.003)
    np.testing.assert_allclose(actual.imag, 0., atol=1e-9)


def test_noise_subtraction_keeps_negative_power_in_low_signal_bin():
    r = _rows(4)
    r['x'][:] = r['y'][:] = 0
    p = radial_test_profile(r, [0., 100.], np.ones(4))
    assert p['debiased_power_jy2'][0] < 0


def test_forward_transform_preserves_off_axis_east_north_phase():
    n, pixel = 512, 0.25
    grid = (np.arange(n)-n//2)*pixel
    east, north = np.meshgrid(-grid, grid)  # native image RA increases to the left
    sx, sy, source, beam, peak = 2., -3., 2., 28., 100.
    image = peak*np.exp(-4*np.log(2)*((east-sx)**2+(north-sy)**2)/source**2)
    header = fits.Header({'NAXIS': 2, 'NAXIS1': n, 'NAXIS2': n,
        'CTYPE1': 'RA---SIN', 'CTYPE2': 'DEC--SIN', 'CRPIX1': n//2+1,
        'CRPIX2': n//2+1, 'CRVAL1': 221.25, 'CRVAL2': -59.8,
        'CDELT1': -pixel/3600, 'CDELT2': pixel/3600})
    u, v = np.array([0., 10000., 20000.]), np.array([0., 6000., -12000.])
    actual = sample_image_visibilities(image, header, 221.25, -59.8, u, v,
        fwhm_arcsec=beam, pixel_arcsec=0.5, patch_arcsec=128.)
    conversion = np.pi/648000
    theta = 1/np.sqrt(source**-2+beam**-2)*conversion
    ratio = beam**2/(beam**2+source**2)
    flux = peak*1e6*np.pi*theta**2/(4*np.log(2))*np.exp(-4*np.log(2)*(sx*sx+sy*sy)/(source*source+beam*beam))
    expected = flux*np.exp(-np.pi**2*theta**2*(u*u+v*v)/(4*np.log(2)))*np.exp(-2j*np.pi*ratio*(u*sx+v*sy)*conversion)
    np.testing.assert_allclose(actual, expected, rtol=0.003)


def test_common_cells_remove_unequal_pointing_population_weights():
    a, b = _rows(104), _rows(104)
    a['field'] = np.r_[np.zeros(100, int), np.ones(4, int)]
    b['field'] = np.r_[np.zeros(4, int), np.ones(100, int)]
    for r in (a, b):
        r['x'] = r['y'] = np.where(r['field'] == 0, 1., 3.).astype(complex)
    pairs = [{'pair': i, 'aca_field': i, 'alma_field': i} for i in range(2)]
    edges = np.array([0., 100.])
    ca, ga = matching_cells(a, pairs, 'aca', edges, 1)
    cb, gb = matching_cells(b, pairs, 'alma', edges, 1)
    records, counts = common_cell_ratios(a, b, ca, ga, cb, gb,
        np.zeros(104), np.zeros(104), pairs, edges,
        {'orientation_sectors': 1, 'random_seed': 1, 'bootstrap_repeats': 20})
    assert records[0]['power_ratio'] == 1.
    np.testing.assert_allclose(records[0]['power_aca_rms_jy'], np.sqrt(5.))
    np.testing.assert_allclose(records[0]['power_alma_rms_jy'], np.sqrt(5.))
    assert counts['retained_aca_rows'] == counts['retained_alma_rows'] == 104
