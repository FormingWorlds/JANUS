# Validation: `src/janus/utils/height.py`

## What this page checks

This page checks how JANUS computes the gravity at the surface and the altitude of each level of the atmosphere. `janus.utils.height` holds the surface-gravity law and the hydrostatic altitude integrator. Validation means that a test runs the code on a case with a known right answer and checks the result within a stated tolerance. If this part of the code is wrong, the heights and the local gravity that JANUS writes to its netCDF output are wrong. Each case below says why the check exists, what it tests, how we know the right answer and why it matters. A box "Test details" ends each case.

## Surface gravity of the Earth

**Why this check exists.** A wrong power of the radius, for example $1/r$ in place of $1/r^2$, or a wrong constant gives a wrong gravity without an error message.

**What it tests.** `gravity(M, r)` returns the Newtonian surface field $g = G M / r^2$, with $G$ the gravitational constant, $M$ the planet mass and $r$ the radius. At the mass of the Earth, 5.972e24 kg, and its mean radius, 6.371e6 m, the test expects about 9.82 m/s$^2$. Doubling the radius must quarter $g$. Doubling the mass must double it.

**How we know the right answer.** The value follows from Newton's law with the Earth's mass and mean radius. It is slightly above 9.81 m/s$^2$. An inverse-linear $1/r$ error would only halve $g$ when the radius doubles. The test tells the 2 laws apart.

**Why it matters.** The surface gravity sets the scale height. That sets the altitude of every level (next case).

![Gravity against distance from the centre of the Earth, with the values at 1 and 2 Earth radii](../assets/validation/janus_height_gravity.svg#only-light){ width="80%" }
![Gravity against distance from the centre of the Earth, with the values at 1 and 2 Earth radii](../assets/validation/janus_height_gravity_dark.svg#only-dark){ width="80%" }

**How to read the plot.** The line is `gravity` for the mass of the Earth from 1 to 3 Earth radii and gives context. The test asserts the code value at 1 Earth radius (left circle) against 9.82 m/s$^2$ (square, with its error bar of 0.02). At 2 Earth radii it asserts the code value (right circle) against $g/4$ (square) to a relative 1e-12. It also asserts that this value lies below $g/3$ (dashed line), which rules out an inverse-linear law. A pass means that each circle sits in its square and the right circle lies below the dashed line. That twice the mass gives twice $g$ is asserted too but not drawn.

??? info "Test details"
    - Reference: Earth mass and mean radius, Newtonian inverse-square law (analytical limit).
    - Tolerance: absolute 0.02 m/s$^2$ at the Earth; relative 1e-12 for the radius and mass scalings.
    - Test id: `tests/utils/test_height.py::test_surface_gravity_matches_earth_inverse_square`
    - Command: `pytest tests/utils/test_height.py::test_surface_gravity_matches_earth_inverse_square`

## Height of an isothermal column

**Why this check exists.** The altitude integrator sums a hydrostatic step over the pressure levels. A missing $1/p$ factor in the step moves every altitude.

**What it tests.** `integrate_heights` sums the hydrostatic step $dz = -(R_\mathrm{gas} T / (M g))\, d\ln p$ up the column. Here $R_\mathrm{gas}$ is the gas constant, $T$ the temperature, $M$ the molar mass and $g$ the gravity. For an isothermal ideal-gas column this integrates to $z(p) = H \ln(p_\mathrm{surface} / p)$, with the scale height $H = R_\mathrm{gas} T / (M g)$. The test compares the total span of a 200-level N$_2$ column at 288 K, from 1e3 to 1e5 Pa, with this integral; $M$ is in kg/mol. It also checks that the altitude falls as the pressure rises. The routine reverses the pressure array internally, so the test also checks that the input array is the same on return.

**How we know the right answer.** The span follows from the closed-form integral of the hydrostatic equation for an isothermal ideal gas. A missing $1/p$ factor would put the span far outside the 5 % tolerance. The code span is 0.53 % shorter than the integral. Its discrete step alone gives 1.15 % less, and the gravity, which falls with height in the code, adds 0.62 % back.

**Why it matters.** `write_ncdf` calls `integrate_heights`, so these heights are the altitudes of the JANUS output. The routine uses `R_gas` and `molar_mass` from `janus.utils.phys` (see [phys.md](phys.md)).

![Pressure against height of an isothermal nitrogen column, from the code and from the scale height](../assets/validation/janus_height_isothermal_column.svg#only-light){ width="80%" }
![Pressure against height of an isothermal nitrogen column, from the code and from the scale height](../assets/validation/janus_height_isothermal_column_dark.svg#only-dark){ width="80%" }

**How to read the plot.** The y-axis is the pressure on a log scale, with the surface at the bottom. The line is $H \ln(p_s/p)$ with the scale height $H$ of the test and gives context. The marks are the heights that `integrate_heights` gives, measured from the lowest level, for every 5th level and for the lowest and highest levels. The test asserts the full height of the column (diamond, 39.87 km) against the height from the scale height (square, 40.09 km). The error bar of the square is the tolerance of 5 %. A pass means that the diamond lies inside the error bar.

**Known limit.** The test computes the analytic span with the same `R_gas`, `molar_mass` and `gravity` that the code uses, so a wrong value of one of these does not fail this test. The gravity case above checks `gravity` on its own.

??? info "Test details"
    - Reference: isothermal hydrostatic column, analytic scale-height integral (analytical limit).
    - Tolerance: relative 0.05 on the total span; the input pressure array identical on return.
    - Test id: `tests/utils/test_height.py::test_integrate_heights_recovers_isothermal_scale_height`
    - Command: `pytest tests/utils/test_height.py::test_integrate_heights_recovers_isothermal_scale_height`
    - Companion test: `test_integrate_heights_flags_blowup_and_returns_bounded_dummy`. A diverging integration is flagged and returns a bounded dummy grid, with the gravity evaluated on it.
    - Companion test: `test_integrate_heights_writes_local_gravity_on_returned_heights`. The local gravity $g(z) = g_s R^2 / (R + z)^2$ is stored on the returned heights and falls with height. Here $g_s$ is the surface gravity, $R$ the planet radius and $z$ the height.
    - Companion test: `tests/utils/test_atmosphere_column.py::test_write_ncdf_gravity_matches_written_heights`. The same relation holds through the netCDF writer.

## Related code

- `src/janus/utils/height.py`: `gravity`, `integrate_heights`.
- `src/janus/utils/phys.py`: `R_gas`, `molar_mass`, imported for the scale-height evaluation; see [phys.md](phys.md).
