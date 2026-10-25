# Validation: `src/janus/utils/cp_funcs.py`

## What this page checks

This page checks the heat capacity of each gas species that JANUS uses. The molar heat capacity $c_p$ is the heat needed to warm 1 mol of a gas by 1 K at constant pressure. `cpv` in `janus.utils.cp_funcs` gives it, either as a flat constant or as a function of temperature. Validation means that a test runs the code on a case with a known right answer and checks the result within a stated tolerance. If this part of the code is wrong, the adiabats and the temperature profile of the atmosphere are wrong too. The case below says why the check exists, what it tests, how we know the right answer and why it matters, and ends with a box "Test details".

## Heat capacity of water vapour at 500 K

**Why this check exists.** `cpv(vol, tmp, cp_mode)` has 2 branches: a flat constant heat capacity (`cp_mode = "constant"`) and the NIST Shomate polynomial (`cp_mode = "T-dependent"`). A wrong coefficient, a per-gram value where a per-mole value belongs, or the wrong branch gives a heat capacity that is wrong without an error message.

**What it tests.** The Shomate branch evaluates $c_p = A + B t + C t^2 + D t^3 + E / t^2$ with $t = T / 1000$ and the tabulated coefficients of water vapour. At 500 K the code returns 35.218 J/mol/K, and the test expects 35.22 J/mol/K. The test also checks that the heat capacity at 1500 K is higher than at 500 K, because vibrational modes activate as the gas warms. A scale guard, $25 < c_p < 55$ J/mol/K, rejects a per-gram value (which would be near 3.5 J/g/K) or a factor-of-ten error.

**How we know the right answer.** The NIST-JANAF thermochemical tables[^cite-janaf] give 35.22 J/mol/K for gaseous water at 500 K. This value tells the Shomate branch apart from the flat constant branch, which does not rise with temperature. Two physics invariants are also asserted: the heat capacity of a known species is positive, and the constant and Shomate values of CO$_2$ roughly agree. An unknown species with no table entry returns 0.0 in both branches; a separate test checks this path.

**Why it matters.** The ratio of the gas constant to $c_p$ sets the slope of the dry adiabat (see [dry_adiabat_setup.md](dry_adiabat_setup.md)). A caller must name the `cp_mode` it needs: if it relies on the default, a change of the default changes its results without a message.

??? info "Test details"
    - Reference: NIST-JANAF Thermochemical Tables[^cite-janaf], Shomate coefficients of water vapour (published benchmark).
    - Tolerance: relative 2e-3 at 500 K; scale guard 25 to 55 J/mol/K.
    - Test id: `tests/utils/test_cp_funcs.py::test_water_heat_capacity_matches_nist_shomate_at_500K`
    - Command: `pytest tests/utils/test_cp_funcs.py::test_water_heat_capacity_matches_nist_shomate_at_500K`

## Related code

- `src/janus/utils/cp_funcs.py`: `cpv`.

## References

[^cite-janaf]: M. W. Chase, *NIST-JANAF Thermochemical Tables*, 4th edition, Journal of Physical and Chemical Reference Data, Monograph 9, 1998. Machine-readable Shomate coefficients for individual species are available through the [NIST Chemistry WebBook](https://webbook.nist.gov).
