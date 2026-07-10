# Validation: `src/janus/modules/moist_adjustment_H2O.py`

This page tracks the `reference_pinned` test that anchors the water moist-
adjustment relaxation in `janus.modules.moist_adjustment_H2O` against the water
triple point and the constant-latent-heat Clausius-Clapeyron dew-point relation.

| Test id | Reference | Anchor type | Scope |
|---|---|---|---|
| `tests/modules/test_moist_adjustment_H2O.py::test_dew_point_target_matches_clausius_clapeyron_reference` | Water triple point (273.15 K at 611.657 Pa) and the constant-L Clausius-Clapeyron inversion | Analytical limit | Pins the dew-point relaxation target to the published water triple-point temperature, verifies the inversion recovers a hand-computed Clausius-Clapeyron forward pressure, and confirms the module drives a triple-point-pressure level to 273.15 K. |

## Re-derivation note

The moist adjustment relaxes a supersaturated water layer toward its local dew
point over a relaxation time `tau`: the temperature tendency is
`(Tdew(pp) - T) / tau`, where `pp` is the local water partial pressure and the
target `Tdew` is `GeneralAdiabat.Tdew("H2O", pp)`. The dew-point curve is the
closed-form inverse of the constant-latent-heat Clausius-Clapeyron relation,
`Tdew = Tref / (1 - (Tref R / L) ln(p / pref))`, with the water triple point
(`Tref = 273.15 K`, `pref = 611.657 Pa`) as its low-pressure reference and the
boiling point (`Tref = 373.15 K`, `pref = 1e5 Pa`) above the triple point.

Three checks anchor the target. First, at the triple-point pressure the dew point
must equal 273.15 K, the published thermodynamic fixed point; this pins the
reference. That single point is degenerate in the latent heat (the logarithm
vanishes at `pref`), so second, an off-reference round-trip pins the latent heat:
a hand-computed forward pressure `p(T) = pref exp((L / R)(1/Tref - 1/T))` fed back
through `Tdew` must recover `T` to nine significant figures, and doubling `L`
shifts the recovered temperature by about 49 K. Third, the published anchor is
tied to the module output: a cold column with an interior level at the triple-
point pressure relaxes that level's target to 273.15 K through `moist_adj`.

The relaxation wiring is verified separately in
`test_supersaturated_column_relaxes_toward_local_dew_point`: the tendency on the
adjusted interior levels equals `(Tdew - T) / tau`, latent heating warms the cold
column, and halving `tau` doubles the tendency, pinning the `1 / tau` dependence.
Two degenerate limits (no water, or an already-subsaturated column) return zero
adjustment and are asserted in
`test_moist_adjustment_zero_without_water_or_when_subsaturated`.

## Anchor type

Analytical limit: the water triple point as a published fixed point of the
dew-point curve, plus the constant-latent-heat Clausius-Clapeyron inversion
verified by a forward round-trip. The sign of the tendency (supersaturated layers
warm toward the dew point) and the `1 / tau` scaling are the physics invariants.

## Cross-references

- `src/janus/modules/moist_adjustment_H2O.py`: the relaxation routine.
- `src/janus/utils/GeneralAdiabat.py`: `Tdew`, the dew-point evaluation the
  adjustment relaxes toward.
- `src/janus/utils/phys.py`: the saturation-vapour-pressure relations and the
  latent heat and gas constant the Clausius-Clapeyron round-trip uses; see
  [phys.md](phys.md).
