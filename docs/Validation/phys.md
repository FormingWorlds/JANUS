# Validation: `src/janus/utils/phys.py`

This page tracks the `reference_pinned` tests that anchor the thermodynamic
constants and saturation-vapour-pressure relations in `janus.utils.phys`
against published data and analytical limits.

| Test id | Reference | Anchor type | Scope |
|---|---|---|---|
| `tests/utils/test_phys.py::test_satvps_returns_reference_pressure_at_reference_temperature` | Clausius-Clapeyron reference identity plus off-reference exponent | Analytical limit | Pins the simplified Clausius-Clapeyron law `satvps` to its reference pressure `e0` at `T = T0`, where the exponent vanishes, and to 53.0 Pa at 30 K below `T0`, where the exponent is far from zero and pins the `L / Rv` coefficient. |
| `tests/utils/test_phys.py::test_satvpw_matches_one_atm_at_steam_point` | Smithsonian Meteorological Tables[^cite-smithsonian] | Published benchmark | Pins the liquid-water saturation curve `satvpw` to 1.0151e5 Pa at the 373.16 K steam point, 0.19% above one standard atmosphere, and brackets the value to rule out a Pa/hPa unit slip. |
| `tests/utils/test_phys.py::test_planck_radiance_positive_and_increases_with_temperature` | Planck closed-form radiance | Analytical limit | Pins the Planck function `B(nu, T)` to its closed form 6.193e-13 W m^-2 sr^-1 Hz^-1 at 5e13 Hz and 300 K, fixing the `nu**3` prefactor and the `h nu / k T` exponent, and brackets the 300-to-1500 K Wien-tail ratio. |

## Re-derivation note

`satvps(T, T0, e0, M, L)` evaluates the integrated Clausius-Clapeyron relation
`e(T) = e0 * exp(-(L / Rv) * (1/T - 1/T0))` with `Rv = Rstar / M`. At the
reference temperature `T = T0` the exponent is identically zero, so `e(T0) = e0`
to machine precision independent of `M` and `L`; this is the degenerate identity
point. Thirty kelvin below `T0` the exponent is `-2.44` and `e = 53.0 Pa` for the
water-like coefficient, so the test pins that value: a doubled `L / Rv` collapses
it to about 4.6 Pa and a halved one lifts it to about 180 Pa, both outside the
`10 Pa < e < 100 Pa` band. A dropped minus sign in the exponent inverts the
temperature response, caught by `e(T0 - 30) < e0 < e(T0 + 30)`.

`satvpw(T)` is the Smithsonian polynomial fit over liquid water. It computes in
dyn/cm^2 and multiplies by 0.1 to return Pa. At 373.16 K it returns 101514.5 Pa;
a forgotten unit conversion would leave the value near 1.015e6 Pa (ten times too
large), which the scale guard `1.00e5 < p < 1.02e5` rejects.

`B(nu, T) = (2 h nu**3 / c**2) / (exp(h nu / k T) - 1)` is the Planck function of
frequency. At 5e13 Hz and 300 K the closed form is 6.193e-13 W m^-2 sr^-1 Hz^-1;
the pin fixes both the `nu**3` prefactor and the `h nu / k T` exponent. Warming
from 300 K to 1500 K at this frequency multiplies the radiance by about 750, far
from the factor of 5 a radiance linear in temperature would give, so the ratio
guard separates the Wien-tail exponential from a linear-in-`T` slip. The overflow
clamp caps `h nu / (k T)` at 500 so an extreme frequency returns a finite,
vanishingly small radiance rather than overflowing.

## Anchor type

Analytical limits (`satvps(T0) = e0` with an off-reference coefficient pin, and
the Planck closed form) plus a published benchmark (the Smithsonian liquid-water
fit at the steam point). Positivity, Clausius-Clapeyron monotonicity, and
`dB/dT > 0` are asserted as the physics invariants.

## Cross-references

- `src/janus/utils/phys.py`: `satvps`, `satvpw`, `satvpi`, and the Planck
  function `B(nu, T)`; the physical constants (`R_gas`, `molar_mass`) that the
  height integrator and adiabat setup import.
- `tests/utils/test_height.py` reuses `phys.R_gas` and `phys.molar_mass` in the
  isothermal scale-height cross-check.

## References

[^cite-smithsonian]: R. J. List, *Smithsonian Meteorological Tables*, 6th revised edition, Smithsonian Institution Press, 1951. The liquid-water (`satvpw`) and ice (`satvpi`) saturation formulae in `phys.py` follow the Smithsonian fits.
