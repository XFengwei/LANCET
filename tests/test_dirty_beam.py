import numpy as np

from lancet.dirty_beam import cumulative_volume, gaussian_beam, radial_profile


def test_gaussian_volume_matches_analytic_beam_area():
    beam = gaussian_beam((401, 401), 0.05, 0.7, 0.6, 47.0, (200, 200))
    volume = cumulative_volume(beam, (200, 200), 0.05, np.array([9.0]))[0]
    assert np.isclose(volume, np.pi*0.7*0.6/(4*np.log(2)), rtol=1e-3)


def test_position_angle_is_east_of_north():
    beam = gaussian_beam((201, 201), 0.05, 1.0, 0.2, 90.0, (100, 100))
    # PA 90 deg puts the major axis east-west: along image rows, not columns.
    assert beam[100, 100+8] > 0.5 > beam[100+8, 100]


def test_radial_profile_peaks_at_centre():
    beam = gaussian_beam((101, 101), 0.1, 0.7, 0.7, 0.0, (50, 50))
    mean, minimum, radius = radial_profile(beam, (50, 50), 0.1, np.arange(0, 3.0, 0.1))
    assert mean[0] == 1.0 and radius[0] == 0.0
    assert np.all(np.diff(mean) < 0) and np.all(minimum <= mean)
