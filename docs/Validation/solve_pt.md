# Validation: `src/janus/modules/solve_pt.py`

## What this page checks

This page checks 2 of the atmosphere solves of JANUS, `RadConvEqm` and `MCPA_CBL` in `janus.modules.solve_pt`. If they are wrong, the outgoing fluxes, the surface temperature and the tropopause temperature that they return are wrong. Each one builds a full atmosphere column and returns its fluxes. Validation means that a test runs the code on a case with a known right answer and checks the result within a stated tolerance. Here the right answer is a stored result of the same code (a regression pin). The tests catch any change of the result, not an error that was there from the start. They run the real SOCRATES radiative-transfer code over a full column solve, so they are integration tests that run every night. Each case below says why the check exists, what it tests, how we know the right answer and why it matters; a box "Test details" ends each case.

## Outgoing radiation along the runaway greenhouse

**Why this check exists.** `RadConvEqm` combines the column setup, the moist adiabat and the SOCRATES radiation. A change in any of them moves the outgoing radiation, and the test catches it.

**What it tests.** `RadConvEqm` builds a radiative-convective equilibrium column and returns the top-of-atmosphere fluxes for a given surface temperature. For a pure-steam atmosphere, the outgoing longwave radiation (OLR) flattens toward the Simpson-Nakajima limit as the surface warms. This limit bounds the OLR while the condensing steam atmosphere is optically thick; at high enough surface temperature the OLR rises again. The test runs the column with `cp_dry=False` at 3 surface temperatures, 200, 1705.3 and 2800 K, and compares the OLR with stored values. 1705.3 K lies near the end of the plateau, and 2800 K beyond it. Before that comparison it checks that the returned flux is finite and positive.

**How we know the right answer.** The stored values are the JANUS + SOCRATES result of the runaway-greenhouse curve of Nakajima et al. (1992) and Goldblatt et al. (2013) (see References). They lock the current solve. A solve that returns `nan` or a negative flux fails the finiteness check first.

**Why it matters.** `RadConvEqm` returns the atmosphere state. When a physics change moves these values on purpose, the stored values change in the same commit, with the reason in the commit message.

![Outgoing longwave radiation against surface temperature, with the 3 pinned values and the code value](../assets/validation/janus_runaway_olr_200k.svg#only-light){ width="80%" }
![Outgoing longwave radiation against surface temperature, with the 3 pinned values and the code value](../assets/validation/janus_runaway_olr_200k_dark.svg#only-dark){ width="80%" }

![Outgoing longwave radiation against surface temperature, with the 3 pinned values and the code value](../assets/validation/janus_runaway_olr_1705k.svg#only-light){ width="80%" }
![Outgoing longwave radiation against surface temperature, with the 3 pinned values and the code value](../assets/validation/janus_runaway_olr_1705k_dark.svg#only-dark){ width="80%" }

![Outgoing longwave radiation against surface temperature, with the 3 pinned values and the code value](../assets/validation/janus_runaway_olr_2800k.svg#only-light){ width="80%" }
![Outgoing longwave radiation against surface temperature, with the 3 pinned values and the code value](../assets/validation/janus_runaway_olr_2800k_dark.svg#only-dark){ width="80%" }

**How to read the plot.** There is one plot for each of the 3 surface temperatures of the test. The grey squares are the 3 pinned values: 90.73, 278.99 and 6583.79 W/m$^2$. In each plot the large square is the pin of that temperature and the circle is the code value, which the test asserts to a relative 1e-5, too small to see. No line is drawn, because each point needs a full radiative-convective solve. A pass means that the circle sits in the large square.

**Known limit.** Only 3 surface temperatures are checked; the test file holds further reference rows on the plateau, but they are commented out. With `cp_dry=False`, the test does not reach the dry-adiabat or the moist-adjustment code.

??? info "Test details"
    - Reference: stored OLR along the pure-steam runaway-greenhouse curve (regression pin of the JANUS + SOCRATES pipeline): 90.7314, 278.990 and 6583.79 W/m$^2$ at 200, 1705.3 and 2800 K.
    - Tolerance: relative 1e-5, absolute 0.
    - Test id: `tests/test_runaway_greenhouse.py::test_runaway_greenhouse`
    - Command: `pytest tests/test_runaway_greenhouse.py::test_runaway_greenhouse` (needs SOCRATES)

## Fluxes against orbital distance

**Why this check exists.** `MCPA_CBL` finds the surface temperature at which the net atmospheric flux equals the heat conducted through a surface skin (Fourier's law). A change in the solve or in the stellar input moves its fluxes and temperatures, and the test catches it.

**What it tests.** The test runs `MCPA_CBL` at orbital distances of 0.3 and 1.4 au (astronomical units). At each distance it compares 5 values with stored values: the downward shortwave flux, the upward longwave flux and the net flux at the top of the atmosphere, the surface temperature and the tropopause temperature. Before that comparison it checks that the surface temperature and the tropopause temperature are positive and that every value is finite.

**How we know the right answer.** Like the case above, this is a regression pin of the JANUS + SOCRATES pipeline. The stored downward shortwave flux follows the inverse square of the distance: 2342.97 / 107.585 = 21.78 = (1.4 / 0.3)$^2$.

**Why it matters.** `MCPA_CBL` returns the atmosphere state with a conductive surface layer.

![Ratio of the code value to the pinned value for the 5 quantities of the instellation test](../assets/validation/janus_instellation_0p3au.svg#only-light){ width="80%" }
![Ratio of the code value to the pinned value for the 5 quantities of the instellation test](../assets/validation/janus_instellation_0p3au_dark.svg#only-dark){ width="80%" }

![Ratio of the code value to the pinned value for the 5 quantities of the instellation test](../assets/validation/janus_instellation_1p4au.svg#only-light){ width="80%" }
![Ratio of the code value to the pinned value for the 5 quantities of the instellation test](../assets/validation/janus_instellation_1p4au_dark.svg#only-dark){ width="80%" }

**How to read the plot.** There is one plot for each of the 2 orbital distances, 0.3 and 1.4 AU. Each mark is one of the 5 quantities divided by its pinned value: downward shortwave flux, upward longwave flux and net flux at the top of the atmosphere, surface temperature and tropopause temperature. The test asserts each to a relative 1e-5, too small to see on this axis. A pass means that every mark lies on the dotted line at 1.

**Known limit.** Only the downward shortwave flux has an analytical scaling; the other 4 values are a regression pin only. `MCPA`, the third solve that `janus.modules` exports, has only a mocked unit test (`tests/modules/test_solve_pt.py`), no stored SOCRATES result.

??? info "Test details"
    - Reference: stored top-of-atmosphere fluxes, surface temperature and tropopause temperature at 0.3 and 1.4 au (regression pin of the JANUS + SOCRATES pipeline).
    - Tolerance: relative 1e-5, absolute 0.
    - Test id: `tests/test_instellation.py::test_instellation`
    - Command: `pytest tests/test_instellation.py::test_instellation` (needs SOCRATES)

## Related code

- `src/janus/modules/solve_pt.py`: `RadConvEqm`, `MCPA_CBL`.
- `src/janus/modules/dry_adiabat_setup.py`, `src/janus/modules/moist_adjustment_H2O.py`: column steps that `RadConvEqm` reaches only with `cp_dry=True` (and, for `moist_adj`, `pure_steam_adj=True`); see [dry_adiabat_setup.md](dry_adiabat_setup.md) and [moist_adjustment_H2O.md](moist_adjustment_H2O.md).
- `tests/helpers/`: the shared column-setup fixtures of these tests. They import the MORS package, so each test module skips when MORS is not installed.

## References

The runaway-greenhouse outgoing-radiation limit follows Nakajima et al. (1992), *Journal of the Atmospheric Sciences* 49, 2256, and the pure-steam formulation of Goldblatt et al. (2013), *Nature Geoscience* 6, 661. The stored values are the JANUS + SOCRATES result of that curve, not a direct transcription of either paper.
