# Validation: `src/janus/utils/phys.py`

## What this page checks

This page checks the thermodynamic relations of `janus.utils.phys`: the saturation vapour pressure of a gas and the Planck function. The saturation vapour pressure is the partial pressure at which a gas is in equilibrium with its condensate. Validation means that a test runs the code on a case with a known right answer and checks the result within a stated tolerance. If `satvps` is wrong, the condensation curves of the moist adiabats in `GeneralAdiabat` are wrong. Each case below says why the check exists, what it tests, how we know the right answer and why it matters. A box "Test details" ends each case.

## Clausius-Clapeyron saturation pressure

**Why this check exists.** `satvps` evaluates an exponential law. A wrong coefficient or a dropped minus sign in the exponent changes the saturation pressure by large factors without an error message.

**What it tests.** `satvps(T, T0, e0, M, L)` evaluates the integrated Clausius-Clapeyron relation $e(T) = e_0 \exp(-(L / R_v)(1/T - 1/T_0))$, with $R_v = R^* / M$. Here $R^*$ is the universal gas constant, $e_0$ the saturation pressure at the reference temperature $T_0$, $L$ the latent heat (J/kg) and $M$ the molar mass (g/mol). At $T = T_0$ the exponent is exactly zero, so $e(T_0) = e_0$ for any $M$ and $L$. 30 K below $T_0$ the exponent is -2.44, and the test expects $e$ = 53.0 Pa for the water-like coefficient. The test also checks $e(T_0 - 30) < e_0 < e(T_0 + 30)$.

**How we know the right answer.** Both values follow from the closed form. The reference point alone does not test $L / R_v$, because the exponent is zero there, so the second point does. A doubled $L / R_v$ gives about 4.6 Pa and a halved one about 180 Pa, both outside the range check of 10 to 100 Pa. A dropped minus sign inverts the temperature response, which the ordering check catches.

**Why it matters.** `GeneralAdiabat.p_sat` builds the saturation curve of each condensing volatile (water, methane, CO$_2$, CO, N$_2$, O$_2$, O$_3$, He, NH$_3$) with `satvps_function`, which calls `satvps`. H$_2$ is forced dry. For water with `water_lookup=True` (off by default), steam tables replace `satvps` above the triple point.

![Saturation pressure from satvps against temperature, with the values at T0 and 30 K below and above it](../assets/validation/janus_phys_satvps.svg#only-light){ width="80%" }
![Saturation pressure from satvps against temperature, with the values at T0 and 30 K below and above it](../assets/validation/janus_phys_satvps_dark.svg#only-dark){ width="80%" }

**How to read the plot.** The line is `satvps` with the inputs of the test and gives context. The inputs are $T_0$ = 273.16 K, $e_0$ = 611 Pa, a molar mass of 18 g/mol and a latent heat of 2.5e6 J/kg. The grey band is the guard band of 10 to 100 Pa. The test asserts 3 marks. At $T_0$ the code value (circle) equals $e_0$ (square) to a relative 1e-12. At 30 K below $T_0$ the code value (circle) is 53.02 Pa (square, with its error bar of 1 %). At 30 K above $T_0$ the code value (circle) is higher than $e_0$. A pass means that the circles at $T_0$ and 30 K below sit in their squares. The circle at 30 K below must also lie in the grey band. The circle at 30 K above must lie above $e_0$.

??? info "Test details"
    - Reference: Clausius-Clapeyron reference identity plus an off-reference exponent (analytical limit).
    - Tolerance: relative 1e-12 at $T_0$; relative 1e-2 on 53.02 Pa at $T_0 - 30$ K; range check 10 to 100 Pa.
    - Test id: `tests/utils/test_phys.py::test_satvps_returns_reference_pressure_at_reference_temperature`
    - Command: `pytest tests/utils/test_phys.py::test_satvps_returns_reference_pressure_at_reference_temperature`

## Saturation pressure of liquid water at the steam point

**Why this check exists.** `satvpw` computes in dyn/cm$^2$ and multiplies by 0.1 to return Pa. A missing conversion gives a pressure 10 times too large.

**What it tests.** `satvpw(T)` is the Smithsonian empirical formula for liquid water. At the steam point, 373.16 K, it returns 101514.5 Pa, 0.19 % above one standard atmosphere. A range check, $1.00 \times 10^5 < p < 1.02 \times 10^5$ Pa, rejects a value near 1.015e6 Pa from a forgotten conversion. The test also checks that the curve rises from 363.16 K to 383.16 K.

**How we know the right answer.** The Smithsonian Meteorological Tables[^cite-smithsonian] give the formula, and at the steam point water boils at about one standard atmosphere. The exact value of 101514.5 Pa is a regression value of the formula; the independent anchor is the range check around one atmosphere.

**Why it matters.** `satvpg` blends the liquid-water and ice formulae (`satvpw`, `satvpi`) between -20 and 0 °C. No other JANUS module calls these 3 functions, so this check protects them for direct use.

![Saturation pressure of liquid water from satvpw against temperature near the steam point](../assets/validation/janus_phys_satvpw.svg#only-light){ width="80%" }
![Saturation pressure of liquid water from satvpw against temperature near the steam point](../assets/validation/janus_phys_satvpw_dark.svg#only-dark){ width="80%" }

**How to read the plot.** The line is `satvpw` from 355 to 390 K and gives context. The test asserts the code value at 373.16 K (circle) against 101514.5 Pa (square, with its error bar of relative 1e-4, too small to see). The 2 triangles are the values 10 K below and above, which must bracket the circle. A pass means that the circle sits in the square and lies between the triangles.

??? info "Test details"
    - Reference: Smithsonian Meteorological Tables[^cite-smithsonian] (published benchmark).
    - Tolerance: relative 1e-4 on 101514.5 Pa; range check 1.00e5 to 1.02e5 Pa.
    - Test id: `tests/utils/test_phys.py::test_satvpw_matches_one_atm_at_steam_point`
    - Command: `pytest tests/utils/test_phys.py::test_satvpw_matches_one_atm_at_steam_point`

## Planck radiance

**Why this check exists.** The Planck function has a $\nu^3$ prefactor and an exponential in $h\nu / kT$. A wrong power or a wrong constant in either part changes the radiance without an error message.

**What it tests.** $B(\nu, T) = (2 h \nu^3 / c^2) / (\exp(h \nu / k T) - 1)$ is the Planck function of frequency $\nu$, with $h$ the Planck constant, $c$ the speed of light and $k$ the Boltzmann constant. At 5e13 Hz and 300 K the test expects 6.193e-13 W m$^{-2}$ sr$^{-1}$ Hz$^{-1}$. Warming from 300 K to 1500 K at this frequency must raise the radiance by more than 100 times. An overflow clamp caps $h \nu / (k T)$ at 500, so an extreme frequency must return a finite, very small radiance.

**How we know the right answer.** The closed form fixes both the $\nu^3$ prefactor and the exponent. Warming from 300 K to 1500 K multiplies the radiance by about 750, far from the factor of 5 that a radiance linear in temperature would give. So the ratio check tells the exponential Wien tail (the high-frequency side of the Planck curve) apart from a linear error.

**Why it matters.** No JANUS source calls `B`, so this check protects it for direct use.

![Planck radiance at 5e13 Hz against temperature, with the values at 300 and 1500 K](../assets/validation/janus_phys_planck.svg#only-light){ width="80%" }
![Planck radiance at 5e13 Hz against temperature, with the values at 300 and 1500 K](../assets/validation/janus_phys_planck_dark.svg#only-dark){ width="80%" }

**How to read the plot.** Both axes are logarithmic. The line is `B` at 5e13 Hz from 200 to 2000 K and gives context. The test asserts the code value at 300 K (circle) against 6.193e-13 W m$^{-2}$ sr$^{-1}$ Hz$^{-1}$ (square). Its error bar of relative 1e-3 is too small to see. The value at 1500 K (circle) must be more than 100 times higher; it is 753 times higher. A pass means that the circle at 300 K sits in the square. The circle at 1500 K must be more than 100 times higher. The test also checks a value at 1e16 Hz, which is not on this axis.

??? info "Test details"
    - Reference: Planck closed-form radiance (analytical limit).
    - Tolerance: relative 1e-3 on 6.193e-13 W m$^{-2}$ sr$^{-1}$ Hz$^{-1}$; range check 5e-13 to 7e-13; ratio above 100 for 300 to 1500 K.
    - Test id: `tests/utils/test_phys.py::test_planck_radiance_positive_and_increases_with_temperature`
    - Command: `pytest tests/utils/test_phys.py::test_planck_radiance_positive_and_increases_with_temperature`
    - Companion test: `test_planck_temperature_derivative_positive_and_matches_finite_difference` checks $dB/dT > 0$.

## Related code

- `src/janus/utils/phys.py`: `satvps`, `satvpw`, `satvpi`, and the Planck function `B(nu, T)`; the physical constants `R_gas` and `molar_mass` that the height integrator imports.
- `tests/utils/test_height.py` uses `phys.R_gas` and `phys.molar_mass` in the isothermal scale-height check.

## References

[^cite-smithsonian]: R. J. List, *Smithsonian Meteorological Tables*, 6th revised edition, Smithsonian Institution Press, 1951. The liquid-water (`satvpw`) and ice (`satvpi`) saturation formulae in `phys.py` follow the Smithsonian fits.
