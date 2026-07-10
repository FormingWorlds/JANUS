# Validation: `src/janus/modules/dry_adiabat_setup.py`

This page tracks the `reference_pinned` test that anchors the dry-adiabat profile
builder in `janus.modules.dry_adiabat_setup` against the Poisson potential-
temperature identity.

| Test id | Reference | Anchor type | Scope |
|---|---|---|---|
| `tests/modules/test_dry_adiabat_setup.py::test_dry_adiabat_conserves_potential_temperature` | Poisson adiabat identity | Analytical limit | Pins the dry-adiabatic profile to constant potential temperature `theta = T (p0 / p)^(R/cp)` along the column, equal to the surface temperature at the reference pressure. |

## Re-derivation note

A dry adiabat conserves potential temperature
`theta = T (p_ref / p)^(R / cp)`. With the reference pressure taken at the
surface, `theta` equals the surface temperature `ts` at every level of an ideal
dry-adiabatic profile. The test constructs the profile and asserts `theta` is
constant along the column and equal to `ts`, which pins the constant-potential-
temperature form and the reference-pressure convention, and brackets the exponent
`R / cp` to the physical `(0, 0.4)` band for a diatomic-dominated column.

The constant-`theta` identity is the discriminating check: a profile built with a
non-Poisson lapse (for example an isothermal column, a "forgot the adiabat" bug)
gives a `theta` that drifts with height instead of staying flat. The test pins
the counterfactual explicitly: an isothermal profile would place `theta` more
than 50% above `ts` at the low-pressure top. The companion test asserts the
temperature rises monotonically from the top of the column to the surface, ruling
out a sign error in the lapse direction.

## Anchor type

Analytical limit (Poisson potential-temperature conservation). Monotonicity of
temperature with pressure along the adiabat is the accompanying physics
invariant.

## Cross-references

- `src/janus/modules/dry_adiabat_setup.py`: the profile builder.
- `src/janus/utils/cp_funcs.py`: the heat capacity `cpv` that sets the `R / cp`
  exponent; see [cp_funcs.md](cp_funcs.md).
