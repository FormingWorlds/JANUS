from janus.utils.height import gravity
from janus.utils.planets import Earth
import numpy as np
from janus.utils.cp_funcs import cpv

def test_gravity_decreases_with_radius():
    m = 5.972e24
    g_surface = gravity(m, 6.371e6)
    g_high = gravity(m, 2 * 6.371e6)
    assert g_surface > g_high

def test_earth_constants_match_database():
    assert Earth.g == 9.798
    assert Earth.albedo == 0.306
    assert Earth.Tsbar == 288.0

def test_cpv_constant_positive():
    for vol in ['H2O', 'CO2', 'N2', 'H2', 'CH4', 'O2']:
        cp = cpv(vol, 300.0, cp_mode='constant')
        assert np.isfinite(cp)
        assert cp > 0.0