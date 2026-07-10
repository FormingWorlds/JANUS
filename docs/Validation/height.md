# Validation: `src/janus/utils/height.py`

This page tracks the `reference_pinned` tests that anchor the surface-gravity
law and the hydrostatic altitude integrator in `janus.utils.height` against
Earth reference values and the analytic isothermal scale-height limit.

| Test id | Reference | Anchor type | Scope |
|---|---|---|---|
| `tests/utils/test_height.py::test_surface_gravity_matches_earth_inverse_square` | Earth mass and mean radius | Analytical limit | Pins the Newtonian law `gravity(M, r) = G M / r^2` to ~9.82 m/s^2 at Earth's mass and mean radius, and verifies the inverse-square and mass-linear scalings. |
| `tests/utils/test_height.py::test_integrate_heights_recovers_isothermal_scale_height` | Isothermal hydrostatic column | Analytical limit | Pins the integrated column height to the analytic scale-height integral `H ln(p_surface / p)` with `H = R_gas T / (M g)` to a few percent on a fine grid. |

## Re-derivation note

`gravity(M, r) = G M / r^2` is the Newtonian surface field. At Earth's mass
5.972e24 kg and mean radius 6.371e6 m it returns ~9.82 m/s^2, slightly above the
standard 9.81 because the mean radius is smaller than the equatorial radius.
Doubling the radius quarters `g` (inverse-square), which discriminates an
inverse-linear `1/r` bug that would only halve it; doubling the mass doubles `g`.

`integrate_heights` sums the hydrostatic step `dz = -R_gas T / (M g) * dln(p)`
up the column. For an isothermal ideal-gas column this integrates in closed form
to `z(p) = H ln(p_surface / p)` with scale height `H = R_gas T / (M g)`. The test
pins the total span against that analytic integral to `rel=0.05`; a missing `1/p`
factor or a wrong molar mass would throw the scale height far outside 5%. The
routine reverses the pressure array internally, so the test also asserts the
input array is restored on return.

## Anchor type

Analytical limit (Earth inverse-square gravity; isothermal scale-height
integral). Monotonicity (altitude anti-monotone with pressure) and the
blow-up error contract (a diverging integration is flagged and returns a bounded
dummy grid) are asserted as the physics invariants.

## Cross-references

- `src/janus/utils/height.py`: `gravity`, `integrate_heights`.
- `src/janus/utils/phys.py`: `R_gas`, `molar_mass` are imported here for the
  scale-height evaluation; see [phys.md](phys.md).
