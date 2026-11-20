import numpy as np
import pytest

from lancet.pointing_coverage import (
    gain_against_pb,
    pb_at_gain_edge,
    primary_beam_fwhm_arcsec,
    unique_positions,
)


def test_primary_beam_fwhm_at_224ghz():
    assert primary_beam_fwhm_arcsec(12.0, 224.4) == pytest.approx(26.0, abs=0.2)
    assert primary_beam_fwhm_arcsec(7.0, 224.4) == pytest.approx(44.5, abs=0.3)


def test_unique_positions_drops_repeats_only():
    dec0 = -60.0
    step = 10.0 / 3600 / np.cos(np.radians(dec0))  # 10 arcsec in RA
    ra = np.array([100.0, 100.0 + 0.1 / 3600, 100.0 + step])
    dec = np.array([dec0, dec0, dec0])
    r, d = unique_positions(ra, dec, tolerance_arcsec=1.0)
    assert r.size == 2
    assert r[1] == pytest.approx(100.0 + step)


def test_gain_against_pb_capped_relation():
    pb = np.linspace(0.05, 1.0, 2000)
    gain = np.minimum(1.0, pb / 0.5)
    rows = gain_against_pb(gain, pb, [0.1, 0.2, 0.6, 1.01])
    assert rows[0]["gain_median"] == pytest.approx(0.3, abs=0.01)
    assert rows[2]["gain_p5"] == 1.0
    edge = pb_at_gain_edge(gain, pb)
    assert edge["pb_median"] == pytest.approx(0.495, abs=0.01)


def test_gain_against_pb_ignores_nan_and_empty_bins():
    rows = gain_against_pb(np.array([np.nan, 1.0]), np.array([0.5, 0.9]), [0.0, 0.1, 1.0])
    assert rows[0]["n"] == 0 and "gain_median" not in rows[0]
    assert rows[1]["n"] == 1
