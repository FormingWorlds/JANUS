# Validation: `src/janus/modules/spectral_planck_surface.py`

## What this page checks

This page checks the thermal emission of the planet surface, summed over spectral bands. If it is wrong, a direct user of `surf_Planck_nu` gets a wrong surface flux. `surf_Planck_nu` in `janus.modules.spectral_planck_surface` fills each band with the Planck emission at the surface temperature. Validation means that a test runs the code on a case with a known right answer and checks the result within a stated tolerance. The case below says why the check exists, what it tests, how we know the right answer and why it matters. A box "Test details" ends the case.

## Band-integrated surface flux and the Stefan-Boltzmann law

**Why this check exists.** The surface flux is a sum over bands of a Planck emission with an albedo factor. A wrong band weight or a wrong temperature power gives a wrong flux without an error message.

**What it tests.** The routine fills each band with the Planck emission at the surface temperature $T_s$. It weights the emission by $(1 - a)$, with $a$ the surface albedo, and multiplies it by $\pi$ and by the band width divided by 1000 (a unit conversion). Summed over a grid that spans the thermal range, the band fluxes must give the Stefan-Boltzmann result $F = (1 - a) \sigma T_s^4$, with $\sigma$ = 5.670374419e-8 W m$^{-2}$ K$^{-4}$. The test uses 60000 bands from 1 to 30000 cm$^{-1}$ (wavenumber), $T_s$ = 1500 K and albedo 0. Halving $T_s$ must lower the flux by $2^4 = 16$.

**How we know the right answer.** The integral of the Planck function over all frequencies is the Stefan-Boltzmann law. The quartic scaling tells the right law apart: a linearised or cubic emission law could match at one temperature but fails the ratio.

**Why it matters.** `surf_Planck_nu` is a standalone helper: no module under `src/janus` calls it, and `plot_flux_balance` only imports it.

![Band-summed surface flux against surface temperature with the Stefan-Boltzmann law](../assets/validation/janus_spectral_stefan_boltzmann.svg#only-light){ width="80%" }
![Band-summed surface flux against surface temperature with the Stefan-Boltzmann law](../assets/validation/janus_spectral_stefan_boltzmann_dark.svg#only-dark){ width="80%" }

**How to read the plot.** The line is $(1-a)\,\sigma T^4$, the flux of a grey surface with albedo $a$, and gives context. The marks are the band sums of `surf_Planck_nu` at 1500 K and 750 K. The test asserts the sum at 1500 K against the Stefan-Boltzmann value (square); the tolerance, a relative 2e-3, is too small to see. It asserts the sum at 750 K only through the ratio of the 2 sums, which must be 16 to within 1 %. A pass means that the mark at 1500 K sits in the square and is 16 times the mark at 750 K.

**Known limit.** The test runs with albedo 0, so it does not check the albedo factor; the companion test does. The test grid is a uniform fine grid, not the band grid of a SOCRATES spectral file.

??? info "Test details"
    - Reference: Stefan-Boltzmann law, grey black-body flux (analytical limit).
    - Tolerance: relative 2e-3 on the band-integrated flux; relative 1e-2 on the flux ratio of 16 between $T_s$ and $T_s / 2$.
    - Test id: `tests/modules/test_spectral_planck_surface.py::test_band_integrated_surface_flux_matches_stefan_boltzmann`
    - Command: `pytest tests/modules/test_spectral_planck_surface.py::test_band_integrated_surface_flux_matches_stefan_boltzmann`
    - Companion test: `test_surface_flux_scales_with_albedo_and_vanishes_at_unit_albedo` drives the flux to zero at unit albedo and checks the linear $(1 - a)$ response between the end points.

## Related code

- `src/janus/modules/spectral_planck_surface.py`: the surface-flux assembly, with its own Planck formula in wavenumber.
