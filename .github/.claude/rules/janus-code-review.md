---
description: JANUS-specific code review criteria for the generator-evaluator pattern. Applies domain expertise (radiative transfer, saturation thermodynamics, spectral unit conventions, PROTEUS coupling) to all code review in this repo.
---

# JANUS Code Review Criteria

When reviewing JANUS code (either your own or via code-reviewer agents), apply these domain-specific checks in addition to standard code quality review.

> **Discovery note.** JANUS keeps its Claude-Code rule files under `.github/.claude/rules/` (not the conventional repo-root `.claude/`) so they can be tracked in git and shared across collaborators. Claude does NOT auto-discover them at this path; the repo-root `CLAUDE.md` (symlinked to `.github/copilot-instructions.md`) names this file and `janus-tests.md` explicitly. **Before opening any review pass, read both this file and `janus-tests.md`.**

## Physics plausibility

- Temperature must be positive everywhere (Kelvin). Flag any code path where `T` could reach zero or go negative (an unclamped moist-adjustment step, a divergent Poisson exponent, an inversion routine that overshoots).
- Pressure must be positive everywhere (Pa). Flag any column-construction path that lets `P <= 0` reach the integrator.
- Saturation vapour pressures (`satvps`, `satvpw`, `satvpi`) must be non-negative and finite. Flag any Clausius-Clapeyron evaluation that could return `nan` / `inf` (a zero or negative temperature reaching the exponential, an overflowing argument that is not capped).
- Planck emission `B(nu, T)` must be non-negative. The source caps the exponent argument at 500 to prevent overflow; flag any new Planck evaluation that drops the cap.
- The band-integrated surface flux must equal the grey Stefan-Boltzmann limit `(1 - albedo) sigma ts**4` as the spectral grid widens. Flag any change to `spectral_planck_surface.py` that breaks this limit.
- The dry-adiabat Poisson exponent `Rcp = R / cp` must be strictly positive. Flag any path that could make it zero, negative, or `nan` (a zero heat capacity in the denominator).
- Integrated height must increase monotonically upward through the column. Flag any hydrostatic-integration change that could produce a non-monotone or blowing-up height array; the inverse-square `gravity(m, r) = G m / r**2` must guard against `r -> 0`.

## Unit convention boundaries

JANUS mixes SI bulk quantities with per-band spectral quantities:

- **Bulk atmosphere state**: `T` in K, `P` in Pa, mass in kg, height in m, gravity in m/s^2.
- **Spectral grid**: wavenumbers in cm^-1, band widths and per-band fluxes normalized per band. `spectral_planck_surface.py` multiplies the Planck function by `band_widths / 1000.0`; that `/1000.0` is the band-width normalization, not a stray constant.
- **Saturation thermodynamics**: the Smithsonian formulae (`satvpw`, `satvpi`) compute in dyn/cm^2 and multiply by 0.1 to return Pa. `satvps` returns Pa directly from `e0` in Pa.

When reviewing code that crosses these boundaries (a new spectral routine, a new saturation formula, a new PROTEUS-side caller), verify the unit is correct at each conversion site. The two recurring traps are (1) **wavenumber vs wavelength** (cm^-1 vs microns) at the spectral grid boundary and (2) **per-band vs per-wavenumber** normalization: a dropped or doubled `/1000.0` band-width factor silently rescales the surface flux while keeping its sign, so it survives a sign-only test.

## Dry-adiabat and saturation argument safety

Several JANUS routines take multiple same-typed float arguments where a swap is plausible and silent:

- `satvps(T, T0, e0, MolecularWeight, LatentHeat)`: `MolecularWeight` and `LatentHeat` are both positive floats. The exponent coefficient `LatentHeat * MolecularWeight / Rstar` is symmetric under the swap, so a swapped call returns an identical value at every temperature and no test can catch it; the call site is the only line of defence. Flag any call site whose argument order is not obviously correct. A value pin at `T != T0` is still worth having, but it guards a different failure: it pins the `L M / Rstar` coefficient and so catches a factor slip in `L` or `M`, not a clean swap.
- `cpv(vol, tmp, cp_mode=...)`: the `cp_mode` string dispatches between a constant heat capacity and the T-dependent Shomate polynomial. Flag any new caller that relies on the default `cp_mode` without stating which mode it needs, and any change to the default value (it propagates to every implicit caller; see the `janus-tests.md` Section 16 trap).

## Config and atmosphere-state mutability

The atmosphere-column object carries the evolving column state and is mutated in place by the adjustment and solver steps; that is by design. What must NOT be mutated at runtime after initialization is user-supplied configuration (planet constants, orbital distance, spectral-file selection). Flag any code that overwrites a configuration field mid-run instead of threading a local variable. When a routine both reads and writes the column arrays, confirm the array lengths and ordering are preserved on return (`p`, `tmp`, `tmpl`, height arrays stay consistent).

## Cross-module constant duplication

Physical constants (`R` / `Rstar`, `G`, `sigma`, `k`, `h`, `c`, molecular weights, latent heats) are defined in `src/janus/utils/phys.py`. When reviewing code that uses a physical constant, check that it is imported from `janus.utils.phys` and not re-derived. A new constant introduced as a literal in a function body (e.g. a hard-coded `5.67e-8` for the Stefan-Boltzmann constant, or a re-typed gas constant) is a red flag: it will drift out of sync with `phys.py`.

## PROTEUS coupling

JANUS is called by PROTEUS through the `atmos_clim` JANUS wrapper during the atmosphere step. PROTEUS-side tests pin the keys of the atmosphere-state result that JANUS returns (surface and top-of-atmosphere fluxes, the temperature profile, the tropopause temperature `trppT`). Flag any change that:

- Renames, drops, or changes the units of a returned state key (`SW_flux_down`, `LW_flux_up`, `net_flux`, `ts`, `trppT`, the temperature / pressure profile arrays) that PROTEUS reads back.
- Silently substitutes a clamped or fallback value for one of those keys without surfacing the substitution.
- Changes the sign convention of a flux (up-positive vs down-positive) at the wrapper boundary.

A JANUS PR that alters the returned-dict contract but does not anticipate the PROTEUS-side fallout is a red flag during review. Announce such a change and confirm the PROTEUS-side pinned values are updated in the same cycle.

## Test marker discipline

Every test file must begin with a module-level `pytestmark = [pytest.mark.<tier>, pytest.mark.timeout(<budget>)]` (unit/30 s, smoke/60 s, integration/300 s, slow/3600 s). Per-function markers are additive but do not replace the module-level marker; CI runs `pytest -m "(unit or smoke) and not skip"` and any file missing the tier marker ships untested.

## Test quality (cross-reference)

Test-content rules (anti-happy-path, discriminating-value guards, physics-invariant tiering, `physics_invariant` / `reference_pinned` certification markers, adversarial-review trigger, mocking discipline, `importorskip` + module-constant-monkeypatch traps, `cp_mode` default flip, `satvps` argument-swap, hypothesis seed stability, pipeline solver intermediate-state assertions) live in [`janus-tests.md`](janus-tests.md). When reviewing tests, apply both files: this one for marker discipline and review-pass gate, the deep-dive for the content contract.

## Sister rules (cross-link)

- [`.github/copilot-instructions.md`](../../copilot-instructions.md) "Testing Standards" -- high-level rules visible to all readers. Repo-root `CLAUDE.md` is a symlink to this file.
- [`janus-tests.md`](janus-tests.md) -- test quality deep-dive; the canonical source for anti-happy-path patterns and the validation certification markers.
