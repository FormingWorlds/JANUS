# Validation: `src/janus/modules/dry_adiabat_setup.py`

## What this page checks

This page checks how JANUS builds a dry adiabat, the temperature profile of a well-mixed dry atmosphere. If it is wrong, the dry-adiabat solve starts from a wrong profile. `dry_adiabat_atm` in `janus.modules.dry_adiabat_setup` builds this profile from the surface temperature and pressure. Validation means that a test runs the code on a case with a known right answer and checks the result within a stated tolerance. The case below says why the check exists, what it tests, how we know the right answer and why it matters. A box "Test details" ends the case.

## Constant potential temperature along a dry adiabat

**Why this check exists.** A profile that is not an adiabat can still look plausible. An isothermal column, for example from a missing adiabat step, gives a wrong temperature at every level without an error message.

**What it tests.** The potential temperature is $\theta = T (p_\mathrm{ref} / p)^{R / c_p}$. It is the temperature that a gas parcel at pressure $p$ and temperature $T$ would have at the reference pressure $p_\mathrm{ref}$ after a dry-adiabatic move. $R / c_p$ is the gas constant over the heat capacity of the gas mixture. With the reference pressure at the surface, $\theta$ must equal the surface temperature $T_s$ at every level of a dry adiabat. The test builds the profile and checks that $\theta$ is constant along the column and equal to $T_s$. It also checks that $R / c_p$ lies between 0 and 0.4 (2/7 for a diatomic gas, 2/5 for a monatomic gas).

**How we know the right answer.** A dry adiabat conserves potential temperature (the Poisson relation), which follows from the first law of thermodynamics for an ideal gas. The test also states the wrong answer: an isothermal profile would put $\theta$ more than 50 % above $T_s$ at the low-pressure top.

**Why it matters.** `compute_dry_adiabat` in `janus.modules.dry_adiabat_timestep` builds its start column with `dry_adiabat_atm`, and `RadConvEqm` calls `compute_dry_adiabat` when `cp_dry=True`. The exponent $R / c_p$ comes from the heat capacity `cpv` (see [cp_funcs.md](cp_funcs.md)).

**Known limit.** The test computes $\theta$ with the `Rcp` that the code itself returns, and checks only that it lies between 0 and 0.4. A wrong heat capacity, and so a wrong exponent inside that range, does not fail this test.

??? info "Test details"
    - Reference: Poisson adiabat identity, conservation of potential temperature (analytical limit).
    - Tolerance: relative 1e-10 for $\theta$ against $T_s$ along the column; relative 1e-12 for the surface temperature.
    - Test id: `tests/modules/test_dry_adiabat_setup.py::test_dry_adiabat_conserves_potential_temperature`
    - Command: `pytest tests/modules/test_dry_adiabat_setup.py::test_dry_adiabat_conserves_potential_temperature`
    - Companion test: `test_dry_adiabat_temperature_rises_monotonically_toward_surface` checks that the temperature rises from the top of the column to the surface, which rules out a lapse of the wrong sign.

## Related code

- `src/janus/modules/dry_adiabat_setup.py`: the profile builder.
- `src/janus/utils/cp_funcs.py`: the heat capacity `cpv` that sets the $R / c_p$ exponent; see [cp_funcs.md](cp_funcs.md).
