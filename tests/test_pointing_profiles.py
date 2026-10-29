"""Pointing isolation and matched TP population checks."""

import numpy as np

from lancet.pointing_profiles import per_pointing_profiles


def test_pointing_profiles_do_not_mix_fields_and_tp_uses_aca_population():
    def cells(count, power):
        return {"count": np.array(count), "cross_power": np.array(power),
                "cross_variance": np.full(4, .01), "radius": np.full(4, 5.5)}
    aca = cells([8, 8, 8, 8], [1., 1., 9., 9.])
    aca["allscale_native"] = aca["allscale_target"] = np.ones(4)
    c2 = cells([0, 0, 8, 0], [0., 0., 9., 0.])
    tp = cells([8, 8, 8, 8], [4., 4., 16., 16.])
    pairs = [{"pair": i, "aca_field": 10+i, "alma_field": 20+i,
              "ra_deg": 221.+i, "dec_deg": -59.8} for i in range(2)]
    rows = per_pointing_profiles({"aca_7m": aca, "configuration_c2": c2, "tp": tp},
                                 pairs, np.array([4.8, 6.]), 1)
    lookup = {(r["pair"], r["array"]): r for r in rows}
    assert (0, "configuration_c2") not in lookup
    for p, aca_rms, tp_rms, ncells in [(0, 1000., 2000., 2), (1, 3000., 4000., 1)]:
        a, t = lookup[p, "aca_7m"], lookup[p, "tp"]
        np.testing.assert_allclose(a["rms_mjy"], aca_rms)
        np.testing.assert_allclose(t["rms_mjy"], tp_rms)
        assert a["cells"] == t["cells"] == ncells
        assert a["rows"] == t["rows"]
        assert a["aca_field"] == 10+p and t["alma_field"] == 20+p
