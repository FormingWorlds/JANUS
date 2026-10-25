# Validation: `src/janus/modules/solve_pt.py`

## What this page checks

This page checks the 2 top-level atmosphere solves of JANUS, `RadConvEqm` and `MCPA_CBL` in `janus.modules.solve_pt`. Each one builds a full atmosphere column and returns its fluxes. Validation means that a test runs the code on a case with a known right answer and checks the result within a stated tolerance. Here the right answer is a stored result of the same code, so these tests catch any change of the result, not an error that was there from the start. They run the real SOCRATES radiative-transfer code over a full column solve, are in the integration tier and run every night. Each case below says why the check exists, what it tests, how we know the right answer and why it matters, and ends with a box "Test details".

## Outgoing radiation along the runaway greenhouse

**Why this check exists.** `RadConvEqm` combines the column setup, the convection and the SOCRATES radiation. A change in any of them can move the outgoing radiation without an error message.

**What it tests.** `RadConvEqm` builds a radiative-convective equilibrium column and returns the top-of-atmosphere fluxes for a given surface temperature. For a pure-steam atmosphere the outgoing longwave radiation (OLR) flattens toward the Simpson-Nakajima limit as the surface warms. The test runs surface temperatures on both sides of that plateau and compares the OLR with stored values. Before that comparison it checks that the returned fluxes are finite and positive.

**How we know the right answer.** The stored values are the JANUS + SOCRATES result of the runaway-greenhouse curve of Nakajima et al. (1992) and Goldblatt et al. (2013) (see References). They are a regression pin: they lock the current solve, and they bracket the analytical Simpson-Nakajima limit. A change of the plateau value or of its onset temperature fails the test. A solve that returns `nan` or a negative flux fails the finiteness check first.

**Why it matters.** `RadConvEqm` is one of the 2 entry points that return the atmosphere state. When a physics change moves these values on purpose, the stored values change in the same commit, with the reason in the commit message.

??? info "Test details"
    - Reference: stored OLR along the pure-steam runaway-greenhouse curve (regression pin of the JANUS + SOCRATES pipeline).
    - Tolerance: relative 1e-5, absolute 0.
    - Test id: `tests/test_runaway_greenhouse.py::test_runaway_greenhouse`
    - Command: `pytest tests/test_runaway_greenhouse.py::test_runaway_greenhouse` (needs SOCRATES)

## Fluxes against orbital distance

**Why this check exists.** `MCPA_CBL` solves the column with a moist convective boundary layer. A change in the solve or in the instellation input moves its fluxes without an error message.

**What it tests.** `MCPA_CBL` returns the surface energy balance for a given instellation, the stellar flux at the planet. The test steps the orbital distance and compares the returned fluxes with stored values. Before that comparison it checks that the surface temperature and the tropopause temperature are positive and that every returned value is finite.

**How we know the right answer.** The stellar flux falls with the inverse square of the distance, and the stored values hold both the magnitude and this scaling. Like the case above, this is a regression pin of the JANUS + SOCRATES pipeline.

**Why it matters.** `MCPA_CBL` is the second entry point that returns the atmosphere state.

??? info "Test details"
    - Reference: stored surface fluxes against orbital distance (regression pin of the JANUS + SOCRATES pipeline).
    - Tolerance: relative 1e-5, absolute 0.
    - Test id: `tests/test_instellation.py::test_instellation`
    - Command: `pytest tests/test_instellation.py::test_instellation` (needs SOCRATES)

## Related code

- `src/janus/modules/solve_pt.py`: `RadConvEqm`, `MCPA_CBL`.
- `src/janus/modules/dry_adiabat_setup.py`, `src/janus/modules/moist_adjustment_H2O.py`: column steps; see [dry_adiabat_setup.md](dry_adiabat_setup.md) and [moist_adjustment_H2O.md](moist_adjustment_H2O.md).
- `tests/helpers/`: the shared column-setup fixtures of these tests, which import `mors` (hence the `importorskip` at the top of each module).

## References

The runaway-greenhouse outgoing-radiation limit follows Nakajima et al. (1992), *Journal of the Atmospheric Sciences* 49, 2256, and the pure-steam formulation of Goldblatt et al. (2013), *Nature Geoscience* 6, 661. The stored values are the JANUS + SOCRATES result of that curve, not a direct transcription of either paper.
