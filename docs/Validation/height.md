# Validation: `src/janus/utils/height.py`

## What this page checks

This page checks how JANUS computes the gravity at the surface and the altitude of each level of the atmosphere. `janus.utils.height` holds the surface-gravity law and the hydrostatic altitude integrator. Validation means that a test runs the code on a case with a known right answer and checks the result within a stated tolerance. If this part of the code is wrong, the heights and the local gravity that JANUS writes to its netCDF output are wrong. Each case below says why the check exists, what it tests, how we know the right answer and why it matters, and ends with a box "Test details".

## Surface gravity of the Earth

**Why this check exists.** A wrong power of the radius, for example $1/r$ in place of $1/r^2$, or a wrong constant gives a wrong gravity without an error message.

**What it tests.** `gravity(M, r)` returns the Newtonian surface field $g = G M / r^2$, with $G$ the gravitational constant, $M$ the planet mass and $r$ the radius. At the mass of the Earth, 5.972e24 kg, and its mean radius, 6.371e6 m, the test expects about 9.82 m/s$^2$. Doubling the radius must quarter $g$, and doubling the mass must double it.

**How we know the right answer.** The value follows from Newton's law with the Earth's mass and mean radius. It is slightly above the standard 9.81 m/s$^2$ because the mean radius is smaller than the equatorial radius. An inverse-linear $1/r$ error would only halve $g$ when the radius doubles, so the test tells the 2 laws apart.

**Why it matters.** The surface gravity sets the scale height and with it the altitude of every level (next case).

??? info "Test details"
    - Reference: Earth mass and mean radius, Newtonian inverse-square law (analytical limit).
    - Tolerance: absolute 0.02 m/s$^2$ at the Earth; relative 1e-12 for the radius and mass scalings.
    - Test id: `tests/utils/test_height.py::test_surface_gravity_matches_earth_inverse_square`
    - Command: `pytest tests/utils/test_height.py::test_surface_gravity_matches_earth_inverse_square`

## Height of an isothermal column

**Why this check exists.** The altitude integrator sums a hydrostatic step over the pressure levels. A missing $1/p$ factor in the step moves every altitude.

**What it tests.** `integrate_heights` sums the hydrostatic step $dz = -(R_\mathrm{gas} T / (M g))\, d\ln p$ up the column. Here $R_\mathrm{gas}$ is the gas constant, $T$ the temperature, $M$ the molar mass and $g$ the gravity. For an isothermal ideal-gas column this integrates in closed form to $z(p) = H \ln(p_\mathrm{surface} / p)$, with the scale height $H = R_\mathrm{gas} T / (M g)$. The test compares the total span of the column with this integral on a fine grid. The routine reverses the pressure array internally, so the test also checks that the input array is the same on return.

**How we know the right answer.** The closed-form integral of the hydrostatic equation for an isothermal ideal gas. A missing $1/p$ factor would put the span far outside the 5 % tolerance. Further physics invariants are asserted:

- the altitude falls as the pressure rises, and the local gravity falls with height;
- the stored local gravity agrees with the returned heights, $g(z) = g_s R^2 / (R + z)^2$ on the same levels, also through the netCDF writer;
- a diverging integration is flagged and returns a bounded dummy grid, with the gravity evaluated on it.

**Why it matters.** `write_ncdf` calls `integrate_heights`, so these heights are the altitudes of the JANUS output. The routine uses `R_gas` and `molar_mass` from `janus.utils.phys` (see [phys.md](phys.md)).

**Known limit.** The test computes the analytic span with the same `R_gas`, `molar_mass` and `gravity` that the code uses, so a wrong value of one of these does not fail this test. The gravity case above checks `gravity` on its own.

??? info "Test details"
    - Reference: isothermal hydrostatic column, analytic scale-height integral (analytical limit).
    - Tolerance: relative 0.05 on the total span; the input pressure array identical on return.
    - Test id: `tests/utils/test_height.py::test_integrate_heights_recovers_isothermal_scale_height`
    - Command: `pytest tests/utils/test_height.py::test_integrate_heights_recovers_isothermal_scale_height`

## Related code

- `src/janus/utils/height.py`: `gravity`, `integrate_heights`.
- `src/janus/utils/phys.py`: `R_gas`, `molar_mass`, imported for the scale-height evaluation; see [phys.md](phys.md).
