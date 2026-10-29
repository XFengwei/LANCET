"""Regression checks for complementary coverage and TP aperture normalization."""

import numpy as np
from astropy.io import fits

from lancet.full_sampling import combine_sampling_cells, separate_array_profiles
from lancet.combination_figure import tp_amplitude_reference
from lancet.tp_beam_test import smooth_to_circular_beam


def test_unique_short_and_long_modes_are_retained_outside_common_overlap():
    a_count, b_count = np.array([8, 8, 0, 0, 0, 0]), np.array([0, 8, 8, 0, 0, 0])
    a = {'count': a_count, 'radius': np.array([5., 15., np.nan, np.nan, np.nan, np.nan]),
         'power': np.full(6, 4.), 'cross_power': np.full(6, 4.)}
    b = {'count': b_count, 'radius': np.array([np.nan, 15., 100., np.nan, np.nan, np.nan]),
         'power': np.ones(6), 'cross_power': np.ones(6)}
    for model in ['allscale', 'twelve_m']:
        a[model+'_native'], a[model+'_target'] = np.full(6, 4.), np.ones(6)
    rows = combine_sampling_cells(a, b, 1, 3, 1)
    assert np.isfinite(rows[0]['aca_rms_mjy']) and np.isnan(rows[0]['alma_rms_mjy'])
    assert rows[1]['matched_overlap']
    assert not rows[0]['matched_overlap'] and not rows[2]['matched_overlap']
    assert np.isfinite(rows[2]['alma_rms_mjy']) and np.isnan(rows[2]['aca_rms_mjy'])
    np.testing.assert_allclose(rows[0]['allscale_aca_12m_pb_mjy'], 1000.)
    np.testing.assert_allclose(rows[1]['allscale_aca_12m_pb_mjy'], rows[1]['alma_rms_mjy'])
    assert np.isnan(rows[2]['allscale_aca_12m_pb_mjy'])


def test_uniform_tp_sky_gets_the_12m_pb_area_without_a_fitted_scale(tmp_path):
    n = 101
    header = fits.Header({'NAXIS': 2, 'NAXIS1': n, 'NAXIS2': n,
        'CTYPE1': 'RA---SIN', 'CTYPE2': 'DEC--SIN', 'CRPIX1': 51., 'CRPIX2': 51.,
        'CRVAL1': 221.25, 'CRVAL2': -59.8, 'CDELT1': -3/3600, 'CDELT2': 3/3600,
        'BMAJ': 18.1/3600, 'BMIN': 18.1/3600, 'BUNIT': 'Jy/beam'})
    path = tmp_path/'tp.fits'
    fits.writeto(path, np.ones((n, n)), header)
    statistics = tmp_path/'cells.npz'
    np.savez(statistics, aca_count=np.array([4, 0, 4, 0]))
    summary = {'settings': {'orientation_sectors': 1, 'pb_12m_fwhm_arcsec': 26.,
        'pb_reference_frequency_ghz': 225., 'model_pixel_arcsec': .5,
        'model_patch_arcsec': 256., 'model_padding_factor': 2},
        'uv_edges_klambda': [4.8, 6., 10.],
        'target_pointings': [{'pair': 0, 'ra_deg': 221.25, 'dec_deg': -59.8}]}
    radius, amplitude, meta = tp_amplitude_reference(path, statistics, summary, .1)
    curves = []
    for frequency in [216.74, 230.95]:
        beam = 26*225/frequency
        theta = beam*np.pi/648000
        curves.append((beam/18.1)**2*np.exp(-np.pi**2*theta**2*(radius*1000)**2/(4*np.log(2))))
    expected_mjy = 1000*np.sqrt(np.mean(np.array(curves)**2, axis=0))
    np.testing.assert_allclose(amplitude, expected_mjy, rtol=.004)
    assert meta['reference_aca_cells'] == 2


def test_separate_arrays_use_common_cells_and_preserve_unique_modes_and_limits():
    def cells(count, power, radius):
        return {'count': np.array(count), 'cross_power': np.array(power),
                'cross_variance': np.ones(4), 'radius': np.array(radius)}
    aca = cells([8, 8, 0, 0], [4., 100., 0., 0.], [5., 6., 0., 0.])
    aca['allscale_native'], aca['allscale_target'] = np.full(4, 4.), np.ones(4)
    c2 = cells([8, 0, 0, 0], [1., 0., 0., 0.], [5.1, 0., 0., 0.])
    c5 = cells([8, 8, 8, 0], [1., 100., -1., 0.], [5.2, 6.2, 100., 0.])
    rows = separate_array_profiles({'aca_7m': aca, 'configuration_c2': c2,
                                    'configuration_c5': c5}, 2, 2)
    common = [r for r in rows if r['bin'] == 0]
    assert len(common) == 3 and all(r['matched_arrays'] and r['cells'] == 1 for r in common)
    np.testing.assert_allclose([r['rms_mjy'] for r in common], 1000.)
    unique = [r for r in rows if r['bin'] == 1]
    assert len(unique) == 1 and unique[0]['array'] == 'configuration_c5'
    assert not unique[0]['matched_arrays'] and not unique[0]['detected']
    assert np.isnan(unique[0]['rms_mjy']) and np.isfinite(unique[0]['upper_rms_mjy'])


def test_beam_matching_preserves_flux_and_produces_target_covariance():
    n, pixel, target = 151, .5, 10.
    header = fits.Header({'NAXIS': 2, 'NAXIS1': n, 'NAXIS2': n,
        'CTYPE1': 'RA---SIN', 'CTYPE2': 'DEC--SIN', 'CRPIX1': 76., 'CRPIX2': 76.,
        'CRVAL1': 221.25, 'CRVAL2': -59.8, 'CDELT1': -pixel/3600, 'CDELT2': pixel/3600,
        'BMAJ': 3/3600, 'BMIN': 2/3600, 'BPA': 30., 'BUNIT': 'MJy/sr'})
    x, y = np.meshgrid(np.arange(n)-n//2, np.arange(n)-n//2)
    east, north = -x*pixel, y*pixel
    angle = np.deg2rad(30.)
    major = east*np.sin(angle)+north*np.cos(angle)
    minor = east*np.cos(angle)-north*np.sin(angle)
    image = np.exp(-4*np.log(2)*((major/3)**2+(minor/2)**2))
    smoothed, info = smooth_to_circular_beam(image, header, target)
    np.testing.assert_allclose(smoothed.sum(), image.sum(), rtol=1e-8)
    moment = np.array([[np.sum(smoothed*x*x), np.sum(smoothed*x*y)],
                       [np.sum(smoothed*x*y), np.sum(smoothed*y*y)]])/smoothed.sum()
    expected = np.eye(2)*(target/pixel)**2/(8*np.log(2))
    np.testing.assert_allclose(moment, expected, atol=.04)
    np.testing.assert_allclose(info['kernel_sum'], 1.)
