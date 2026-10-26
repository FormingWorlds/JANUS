# Validation: `src/janus/utils/cp_funcs.py`

## What this page checks

This page checks the heat capacity of a gas species in JANUS. If it is wrong, the adiabats and the stratosphere of a run that uses this branch are wrong too. The molar heat capacity $c_p$ is the heat needed to warm 1 mol of a gas by 1 K at constant pressure. `cpv` in `janus.utils.cp_funcs` gives it, either as a flat constant or as a function of temperature. Validation means that a test runs the code on a case with a known right answer and checks the result within a stated tolerance. The case below says why the check exists, what it tests, how we know the right answer and why it matters. A box "Test details" ends the case.

## Heat capacity of water vapour at 500 K

**Why this check exists.** `cpv(vol, tmp, cp_mode)` has 2 branches: a flat constant (`cp_mode = "constant"`) and the NIST Shomate polynomial (`cp_mode = "T-dependent"`). A wrong coefficient of the low-temperature water fit, a per-gram value where a per-mole value belongs, or the constant branch in place of the Shomate branch gives a wrong heat capacity without an error message.

**What it tests.** The Shomate polynomial is a standard fit of a heat capacity against temperature: $c_p = A + B t + C t^2 + D t^3 + E / t^2$ with $t = T / 1000$. The test uses the tabulated coefficients of water vapour. At 500 K the code returns 35.218 J/mol/K, and the test expects 35.22 J/mol/K. It also checks that the value at 1500 K is higher than at 500 K, because vibrational modes activate as the gas warms. A range check, $25 < c_p < 55$ J/mol/K, rejects a per-gram value (which would be near 2 J/g/K) or a factor-of-ten error.

**How we know the right answer.** The NIST-JANAF thermochemical tables[^cite-janaf] give 35.22 J/mol/K for gaseous water at 500 K. This value tells the Shomate branch apart from the flat constant branch, which does not rise with temperature.

**Why it matters.** The heat capacity enters the adiabats and the stratosphere of JANUS through `cpv`. The ratio of the gas constant to $c_p$ sets the slope of the dry adiabat (see [dry_adiabat_setup.md](dry_adiabat_setup.md)).

![Heat capacity of water vapour against temperature from both branches of cpv, with the values at 500 and 1500 K](../assets/validation/janus_cp_h2o_500k.svg#only-light){ width="80%" }
![Heat capacity of water vapour against temperature from both branches of cpv, with the values at 500 and 1500 K](../assets/validation/janus_cp_h2o_500k_dark.svg#only-dark){ width="80%" }

**How to read the plot.** The lines are `cpv` for water vapour from 300 to 2000 K and give context: the solid line is the `T-dependent` branch and the dashed line the constant branch. Below 500 K the solid line is flat, because the code never uses a temperature below 500 K, the lower limit of the fit. The test asserts the code value at 500 K (circle) against 35.22 J/mol/K (square, with its error bar of relative 2e-3). It also asserts that the value at 1500 K (triangle) is higher. A pass means that the circle lies inside the error bar and the triangle lies above the circle.

**Known limit.** No JANUS call passes `cp_mode`, so the model uses the constant branch by default. This test checks only the `T-dependent` branch, and only the low-temperature coefficients of water at 500 K.

??? info "Test details"
    - Reference: NIST-JANAF Thermochemical Tables[^cite-janaf], Shomate coefficients of water vapour (published benchmark).
    - Tolerance: relative 2e-3 at 500 K; range check 25 to 55 J/mol/K.
    - Test id: `tests/utils/test_cp_funcs.py::test_water_heat_capacity_matches_nist_shomate_at_500K`
    - Command: `pytest tests/utils/test_cp_funcs.py::test_water_heat_capacity_matches_nist_shomate_at_500K`
    - Companion test: `test_unknown_species_returns_zero_and_known_species_positive`. It checks that an unknown species returns 0.0 in both branches. It also checks that the constant and Shomate values of CO$_2$ at 500 K are positive and agree within 35 %.

## Related code

- `src/janus/utils/cp_funcs.py`: `cpv`.

## References

[^cite-janaf]: M. W. Chase, *NIST-JANAF Thermochemical Tables*, 4th edition, Journal of Physical and Chemical Reference Data, Monograph 9, 1998. Machine-readable Shomate coefficients for individual species are available through the [NIST Chemistry WebBook](https://webbook.nist.gov).
