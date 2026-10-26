# Validation: `src/janus/modules/moist_adjustment_H2O.py`

## What this page checks

This page checks how JANUS relaxes a layer that holds more water vapour than it can hold at its temperature. Such a supersaturated layer warms toward its dew point, the temperature at which its water vapour is just saturated. If this code is wrong, those layers move to a wrong temperature. `moist_adj` in `janus.modules.moist_adjustment_H2O` computes this temperature change. Validation means that a test runs the code on a case with a known right answer and checks the result within a stated tolerance. The case below says why the check exists, what it tests, how we know the right answer and why it matters. A box "Test details" ends the case.

## Dew-point target and the Clausius-Clapeyron relation

**Why this check exists.** The relaxation target comes from an inverted Clausius-Clapeyron relation with a reference point and a latent heat. A wrong reference point, or a latent heat in `Tdew` that differs from the one in `phys`, moves the target temperature without an error message.

**What it tests.** The moist adjustment relaxes a supersaturated water layer toward its local dew point over a relaxation time $\tau$. The temperature tendency is $(T_\mathrm{dew}(p_w) - T) / \tau$, with $p_w$ the local partial pressure of water and $T_\mathrm{dew}$ = `GeneralAdiabat.Tdew("H2O", p_w)`. The dew-point curve is the closed-form inverse of the Clausius-Clapeyron relation with a constant molar latent heat $L$ (J/mol) and the gas constant $R$:

$$T_\mathrm{dew} = \frac{T_\mathrm{ref}}{1 - (T_\mathrm{ref} R / L) \ln(p_w / p_\mathrm{ref})}.$$

At or below 611.657 Pa the code takes $T_\mathrm{ref}$ = 273.15 K, $p_\mathrm{ref}$ = 611.657 Pa and the latent heat of sublimation. Above it, it takes the boiling point, $T_\mathrm{ref}$ = 373.15 K and $p_\mathrm{ref}$ = 1e5 Pa, and the latent heat of vaporization. The test makes 3 checks:

- At the triple-point pressure the dew point must be 273.15 K.
- A forward pressure computed by hand, $p(T) = p_\mathrm{ref} \exp((L / R)(1/T_\mathrm{ref} - 1/T))$, fed back through $T_\mathrm{dew}$ must give $T$ = 340 K again to 9 significant figures.
- A cold column with an interior level at the triple-point pressure must relax that level's target to 273.15 K through `moist_adj`.

**How we know the right answer.** The water triple point, 273.16 K at 611.657 Pa, is a published fixed point. The code uses 273.15 K, inside the 0.05 K tolerance. At that point the logarithm is zero, so it does not test the latent heat; the round trip does. A doubled $L$ shifts the recovered temperature by about 28 K, and the test asserts a shift of more than 1 K. The third check ties the fixed point to the output of the module.

**Why it matters.** `compute_dry_adiabat` calls `moist_adj` when `pure_steam_adj=True`, and `RadConvEqm` calls `compute_dry_adiabat` when `cp_dry=True` (see [solve_pt.md](solve_pt.md)).

**Known limit.** The round trip takes the latent heat of water from `phys.water`, the same source that `Tdew` uses, so a wrong value there does not fail this test. The round trip runs above the triple-point pressure only, so the sublimation branch below it has no round-trip check.

??? info "Test details"
    - Reference: water triple point (273.16 K at 611.657 Pa; the code uses 273.15 K) and the constant-$L$ Clausius-Clapeyron inversion (analytical limit).
    - Tolerance: absolute 0.05 K at the triple point; relative 1e-9 for the round trip; absolute 0.1 K for the `moist_adj` target; a doubled $L$ must miss by more than 1 K.
    - Test id: `tests/modules/test_moist_adjustment_H2O.py::test_dew_point_target_matches_clausius_clapeyron_reference`
    - Command: `pytest tests/modules/test_moist_adjustment_H2O.py::test_dew_point_target_matches_clausius_clapeyron_reference`
    - Companion test: `test_supersaturated_column_relaxes_toward_local_dew_point`. The tendency on the adjusted interior levels equals $(T_\mathrm{dew} - T) / \tau$, latent heating warms the cold column, and doubling $\tau$ halves the tendency.
    - Companion test: `test_moist_adjustment_zero_without_water_or_when_subsaturated`. There is no adjustment without water or in a subsaturated column.

## Related code

- `src/janus/modules/moist_adjustment_H2O.py`: the relaxation routine.
- `src/janus/utils/GeneralAdiabat.py`: `Tdew`, the dew-point evaluation the adjustment relaxes toward.
- `src/janus/utils/phys.py`: the saturation-vapour-pressure relations and the latent heat and gas constant of the Clausius-Clapeyron round trip; see [phys.md](phys.md).
