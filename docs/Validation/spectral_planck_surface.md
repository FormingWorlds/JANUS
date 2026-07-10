# Validation: `src/janus/modules/spectral_planck_surface.py`

This page tracks the `reference_pinned` test that anchors the band-integrated
surface thermal emission in `janus.modules.spectral_planck_surface` against the
grey Stefan-Boltzmann limit.

| Test id | Reference | Anchor type | Scope |
|---|---|---|---|
| `tests/modules/test_spectral_planck_surface.py::test_band_integrated_surface_flux_matches_stefan_boltzmann` | Stefan-Boltzmann law | Analytical limit | Pins the sum of the per-band surface fluxes to the grey black-body value `(1 - a) sigma Ts^4` and verifies the quartic temperature scaling. |

## Re-derivation note

The surface routine fills each spectral band with the Planck radiance at the
surface temperature and integrates over the band grid, weighting by
`(1 - surface_albedo)`. Summed over a grid that spans the thermal range, the
band-integrated flux converges to the Stefan-Boltzmann result
`F = (1 - a) sigma Ts^4` with `sigma = 5.670374419e-8 W m^-2 K^-4`. The test pins
the integral to that closed form and, by evaluating two surface temperatures,
confirms the ratio tracks `(Ts2 / Ts1)^4`.

The quartic scaling is the discriminating check: a linearised or cubic emission
law would match at a single temperature but fail the ratio. The `(1 - a)`
prefactor is exercised by the companion albedo test, which drives the flux to
zero at unit albedo and asserts the linear `(1 - albedo)` response between the
end points.

## Anchor type

Analytical limit (grey Stefan-Boltzmann black-body flux). Positivity of the
band fluxes and the `(1 - albedo)` linear response are the physics invariants.

## Cross-references

- `src/janus/modules/spectral_planck_surface.py`: the surface-flux assembly.
- `src/janus/utils/phys.py`: the Planck function `B(nu, T)` the band integral
  samples; see [phys.md](phys.md).
