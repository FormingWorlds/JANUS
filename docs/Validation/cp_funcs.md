# Validation: `src/janus/utils/cp_funcs.py`

This page tracks the `reference_pinned` test that anchors the per-species molar
heat capacity `cpv` in `janus.utils.cp_funcs` against the NIST-JANAF
thermochemical tables.

| Test id | Reference | Anchor type | Scope |
|---|---|---|---|
| `tests/utils/test_cp_funcs.py::test_water_heat_capacity_matches_nist_shomate_at_500K` | NIST-JANAF Thermochemical Tables[^cite-janaf] | Published benchmark | Pins the temperature-dependent Shomate branch for gaseous water to 35.22 J/mol/K at 500 K and verifies the heat capacity rises with temperature as vibrational modes activate. |

## Re-derivation note

`cpv(vol, tmp, cp_mode)` dispatches between a flat constant heat capacity
(`cp_mode = "constant"`) and the NIST Shomate polynomial (`cp_mode =
"T-dependent"`). The Shomate branch evaluates
`Cp = A + B t + C t^2 + D t^3 + E / t^2` with `t = T / 1000` and the tabulated
water-vapour coefficients. At 500 K the tabulated molar heat capacity of gaseous
water is 35.22 J/mol/K; the implementation returns 35.218, pinned to `rel=2e-3`.

The value 35.22 J/mol/K distinguishes the Shomate branch from the flat constant
branch, which does not rise with temperature. The scale guard `25 < cp < 55`
J/mol/K rejects a per-gram versus per-mole unit slip (which would land near 3.5
J/g/K numerically) or a factor-of-ten error. An unknown species with no table
entry returns 0.0 in both branches, the guard path exercised separately.

## Anchor type

Published benchmark (NIST-JANAF Shomate water-vapour coefficients). Positivity of
a known species' heat capacity and rough constant-versus-Shomate agreement for
CO2 are asserted as the physics invariants.

## Cross-references

- `src/janus/utils/cp_funcs.py`: `cpv`.
- The `cp_mode` string selects the branch; a caller that relies on the default
  without naming the mode it needs is a review red flag, since a default flip
  silently redirects every implicit caller.

## References

[^cite-janaf]: M. W. Chase, *NIST-JANAF Thermochemical Tables*, 4th edition, Journal of Physical and Chemical Reference Data, Monograph 9, 1998. Machine-readable Shomate coefficients for individual species are available through the [NIST Chemistry WebBook](https://webbook.nist.gov).
