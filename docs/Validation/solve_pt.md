# Validation: `src/janus/modules/solve_pt.py`

This page tracks the `reference_pinned` pipeline tests that anchor the two
top-level atmosphere-structure entry points, `RadConvEqm` and `MCPA_CBL`, in
`janus.modules.solve_pt` against pinned outgoing-radiation references. These are
integration-tier tests: they run the real SOCRATES binary over a full column
solve and are exercised nightly.

| Test id | Reference | Anchor type | Scope |
|---|---|---|---|
| `tests/test_runaway_greenhouse.py::test_runaway_greenhouse` | Pinned OLR along the pure-steam runaway-greenhouse curve | Cross-implementation regression pin | Runs `RadConvEqm` at surface temperatures spanning the runaway-greenhouse branch and pins the outgoing longwave radiation to its reference values, bracketing the Simpson-Nakajima limit. |
| `tests/test_instellation.py::test_instellation` | Pinned surface flux versus orbital separation | Cross-implementation regression pin | Runs `MCPA_CBL` across orbital separations and pins the returned fluxes to their reference values, verifying the inverse-square instellation response. |

## Re-derivation note

`RadConvEqm` builds a radiative-convective-equilibrium column and returns the
top-of-atmosphere fluxes for a given surface temperature. Along the pure-steam
runaway-greenhouse branch the outgoing longwave radiation flattens toward the
Simpson-Nakajima limit as the surface warms; the test pins the OLR at surface
temperatures on either side of that plateau, so a change that shifts the plateau
value or its onset temperature is caught. Before the pinned comparison the test
asserts the returned flux array is finite and positive, so a solver that returns
`nan` or a non-physical negative flux fails on the invariant rather than slipping
through the numeric pin.

`MCPA_CBL` solves the moist convective boundary-layer column and returns the
surface energy balance for a given instellation. The instellation test steps the
orbital separation and pins the returned fluxes; the physical expectation is the
inverse-square falloff of stellar flux with distance, and the pinned values
encode both the magnitude and that scaling. The pre-pin sanity block asserts the
surface temperature and tropopause temperature are positive and every returned
value is finite.

Both pins are regression anchors: they lock the current coupled SOCRATES + JANUS
solve against silent drift. When a physics change legitimately moves these
values, the reference numbers are updated in the same change and the motivation
recorded in the commit message.

## Anchor type

Cross-implementation regression pin (the full JANUS + SOCRATES pipeline against
stored outgoing-radiation references), bracketing the analytical
Simpson-Nakajima runaway-greenhouse limit and the inverse-square instellation
law. Finiteness and positivity of the returned fluxes are the accompanying
physics invariants.

## Cross-references

- `src/janus/modules/solve_pt.py`: `RadConvEqm`, `MCPA_CBL`.
- `src/janus/modules/dry_adiabat_setup.py`, `src/janus/modules/moist_adjustment_H2O.py`:
  the column-construction steps the pipeline composes; see
  [dry_adiabat_setup.md](dry_adiabat_setup.md) and
  [moist_adjustment_H2O.md](moist_adjustment_H2O.md).
- `tests/helpers/`: the shared column-setup fixtures these pipeline tests import,
  which pull in `mors` transitively (hence the module-top `importorskip`).

## References

The runaway-greenhouse outgoing-radiation limit follows Nakajima et al. (1992),
*Journal of the Atmospheric Sciences* 49, 2256, and the pure-steam formulation of
Goldblatt et al. (2013), *Nature Geoscience* 6, 661. The pinned values are the
JANUS + SOCRATES realisation of that curve, not a direct transcription of either
paper.
