# Validation: `src/janus/modules/moist_adjustment_H2O.py`

## What this page checks

This page checks how JANUS relaxes a layer that holds more water vapour than it can hold at its temperature. Such a supersaturated layer warms toward its dew point, the temperature at which its water vapour is just saturated. `janus.modules.moist_adjustment_H2O` computes this temperature change. Validation means that a test runs the code on a case with a known right answer and checks the result within a stated tolerance. If this part of the code is wrong, `compute_dry_adiabat`, which calls `moist_adj`, moves supersaturated layers to a wrong temperature. The case below says why the check exists, what it tests, how we know the right answer and why it matters, and ends with a box "Test details".

## Dew-point target and the Clausius-Clapeyron relation

**Why this check exists.** The relaxation target comes from an inverted Clausius-Clapeyron relation with a reference point and a latent heat. A wrong reference point, or a latent heat in `Tdew` that differs from the one in `phys`, moves the target temperature without an error message.

**What it tests.** The moist adjustment relaxes a supersaturated water layer toward its local dew point over a relaxation time $\tau$. The temperature tendency is $(T_\mathrm{dew}(p_p) - T) / \tau$, with $p_p$ the local partial pressure of water and the target $T_\mathrm{dew}$ = `GeneralAdiabat.Tdew("H2O", pp)`. The dew-point curve is the closed-form inverse of the Clausius-Clapeyron relation with a constant latent heat $L$:

$$T_\mathrm{dew} = \frac{T_\mathrm{ref}}{1 - (T_\mathrm{ref} R / L) \ln(p / p_\mathrm{ref})}.$$

The reference is the water triple point ($T_\mathrm{ref}$ = 273.15 K, $p_\mathrm{ref}$ = 611.657 Pa) at low pressure, and the boiling point ($T_\mathrm{ref}$ = 373.15 K, $p_\mathrm{ref}$ = 1e5 Pa) above the triple point. The test makes 3 checks:

- At the triple-point pressure the dew point must be 273.15 K.
- A forward pressure computed by hand, $p(T) = p_\mathrm{ref} \exp((L / R)(1/T_\mathrm{ref} - 1/T))$, fed back through $T_\mathrm{dew}$ must give $T$ again to 9 significant figures.
- A cold column with an interior level at the triple-point pressure must relax that level's target to 273.15 K through `moist_adj`.

**How we know the right answer.** The triple point of water is a published thermodynamic fixed point. At that point the logarithm is zero, so it does not test the latent heat; the round trip does. Doubling $L$ shifts the recovered temperature by about 49 K. The third check ties the published fixed point to the output of the module. The physics invariants are the sign of the tendency (supersaturated layers warm toward the dew point) and its $1 / \tau$ scaling.

**Why it matters.** `compute_dry_adiabat` calls `moist_adj`, and `RadConvEqm` with `cp_dry=True` calls `compute_dry_adiabat` (see [solve_pt.md](solve_pt.md)). Other tests check the relaxation itself: in `test_supersaturated_column_relaxes_toward_local_dew_point` the tendency on the adjusted interior levels equals $(T_\mathrm{dew} - T) / \tau$, latent heating warms the cold column, and halving $\tau$ doubles the tendency. In `test_moist_adjustment_zero_without_water_or_when_subsaturated`, a column with no water or an already subsaturated column has no adjustment.

**Known limit.** The round trip takes the latent heat of water from `phys.water`, the same source that `Tdew` uses, so a wrong value there does not fail this test.

??? info "Test details"
    - Reference: water triple point (273.15 K at 611.657 Pa) and the constant-$L$ Clausius-Clapeyron inversion (analytical limit).
    - Tolerance: absolute 0.05 K at the triple point; relative 1e-9 for the round trip; absolute 0.1 K for the `moist_adj` target; a doubled $L$ must miss by more than 1 K.
    - Test id: `tests/modules/test_moist_adjustment_H2O.py::test_dew_point_target_matches_clausius_clapeyron_reference`
    - Command: `pytest tests/modules/test_moist_adjustment_H2O.py::test_dew_point_target_matches_clausius_clapeyron_reference`

## Related code

- `src/janus/modules/moist_adjustment_H2O.py`: the relaxation routine.
- `src/janus/utils/GeneralAdiabat.py`: `Tdew`, the dew-point evaluation the adjustment relaxes toward.
- `src/janus/utils/phys.py`: the saturation-vapour-pressure relations and the latent heat and gas constant the Clausius-Clapeyron round trip uses; see [phys.md](phys.md).
