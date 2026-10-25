# Validation: `src/janus/modules/dry_adiabat_setup.py`

## What this page checks

This page checks how JANUS builds a dry adiabat, the temperature profile of a well-mixed dry atmosphere. `janus.modules.dry_adiabat_setup` builds this profile from the surface temperature and pressure. Validation means that a test runs the code on a case with a known right answer and checks the result within a stated tolerance. If this part of the code is wrong, the dry-adiabat solve (`compute_dry_adiabat`), which builds its column with `dry_adiabat_atm`, starts from a wrong temperature profile. The case below says why the check exists, what it tests, how we know the right answer and why it matters, and ends with a box "Test details".

## Constant potential temperature along a dry adiabat

**Why this check exists.** A profile built with the wrong lapse rate looks plausible but is not an adiabat. An isothermal column, for example from a missing adiabat step, or a lapse with the wrong sign gives a wrong temperature at every level.

**What it tests.** The potential temperature $\theta = T (p_\mathrm{ref} / p)^{R / c_p}$ is the temperature that a gas parcel at pressure $p$ and temperature $T$ would have if it were moved dry-adiabatically to the reference pressure $p_\mathrm{ref}$. Here $R / c_p$ is the gas constant over the heat capacity. With the reference pressure at the surface, $\theta$ must equal the surface temperature $T_s$ at every level of a dry adiabat. The test builds the profile and checks that $\theta$ is constant along the column and equal to $T_s$. It also checks that the exponent $R / c_p$ lies between 0 and 0.4, the physical range for a column made mostly of diatomic gas.

**How we know the right answer.** A dry adiabat conserves potential temperature (the Poisson relation); this follows from the first law of thermodynamics for an ideal gas. The test also states the wrong answer explicitly: an isothermal profile would put $\theta$ more than 50 % above $T_s$ at the low-pressure top. A companion test checks that the temperature rises monotonically from the top of the column to the surface, which rules out a sign error in the lapse direction.

**Why it matters.** `RadConvEqm` with `cp_dry=True` calls `compute_dry_adiabat`, which builds its column with this profile (see [solve_pt.md](solve_pt.md)). Its exponent $R / c_p$ comes from the heat capacity `cpv` (see [cp_funcs.md](cp_funcs.md)).

??? info "Test details"
    - Reference: Poisson adiabat identity, conservation of potential temperature (analytical limit).
    - Tolerance: relative 1e-10 for $\theta$ against $T_s$ along the column; relative 1e-12 for the surface temperature.
    - Test id: `tests/modules/test_dry_adiabat_setup.py::test_dry_adiabat_conserves_potential_temperature`
    - Command: `pytest tests/modules/test_dry_adiabat_setup.py::test_dry_adiabat_conserves_potential_temperature`

## Related code

- `src/janus/modules/dry_adiabat_setup.py`: the profile builder.
- `src/janus/utils/cp_funcs.py`: the heat capacity `cpv` that sets the $R / c_p$ exponent; see [cp_funcs.md](cp_funcs.md).
