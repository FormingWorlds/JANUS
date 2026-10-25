# Validation: `src/janus/modules/spectral_planck_surface.py`

## What this page checks

This page checks the thermal emission of the planet surface, summed over the spectral bands of JANUS. `janus.modules.spectral_planck_surface` fills each band with the Planck radiance at the surface temperature. Validation means that a test runs the code on a case with a known right answer and checks the result within a stated tolerance. If this part of the code is wrong, the surface emission that `plot_flux_balance` draws is wrong. The case below says why the check exists, what it tests, how we know the right answer and why it matters, and ends with a box "Test details".

## Band-integrated surface flux and the Stefan-Boltzmann law

**Why this check exists.** The surface flux is a sum over bands of a Planck radiance with an albedo factor. A wrong band weight or a wrong temperature power gives a wrong flux without an error message. This test runs with albedo 0; the companion albedo test checks the albedo factor.

**What it tests.** The routine fills each band with the Planck radiance at the surface temperature $T_s$, integrates over the band grid and weights by $(1 - a)$, with $a$ the surface albedo. Summed over a grid that spans the thermal range, the band fluxes must give the Stefan-Boltzmann result $F = (1 - a) \sigma T_s^4$, with $\sigma$ = 5.670374419e-8 W m$^{-2}$ K$^{-4}$. The test evaluates 2 surface temperatures and checks that the ratio of the fluxes follows $(T_{s2} / T_{s1})^4$.

**How we know the right answer.** The integral of the Planck function over all frequencies is the Stefan-Boltzmann law. The quartic scaling tells the right law apart: a linearised or cubic emission law could match at one temperature but fails the ratio. A companion albedo test drives the flux to zero at unit albedo and checks the linear $(1 - a)$ response between the end points. The physics invariants are positive band fluxes and the linear albedo response.

**Why it matters.** `plot_flux_balance` uses `surf_Planck_nu` to draw the surface emission. The band integral samples the Planck function $B(\nu, T)$ of `janus.utils.phys` (see [phys.md](phys.md)).

??? info "Test details"
    - Reference: Stefan-Boltzmann law, grey black-body flux (analytical limit).
    - Tolerance: relative 2e-3 on the band-integrated flux; relative 1e-2 on the flux ratio of 16 for a doubled temperature.
    - Test id: `tests/modules/test_spectral_planck_surface.py::test_band_integrated_surface_flux_matches_stefan_boltzmann`
    - Command: `pytest tests/modules/test_spectral_planck_surface.py::test_band_integrated_surface_flux_matches_stefan_boltzmann`

## Related code

- `src/janus/modules/spectral_planck_surface.py`: the surface-flux assembly.
- `src/janus/utils/phys.py`: the Planck function `B(nu, T)` that the band integral samples; see [phys.md](phys.md).
