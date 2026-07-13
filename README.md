[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Docs](https://img.shields.io/github/actions/workflow/status/FormingWorlds/JANUS/docs.yaml?branch=main&label=Docs)](https://proteus-framework.org/JANUS/)
[![codecov](https://img.shields.io/codecov/c/github/FormingWorlds/JANUS?label=coverage&logo=codecov)](https://app.codecov.io/gh/FormingWorlds/JANUS)
[![Unit Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/JANUS/tests.yaml?branch=main&label=Unit%20Tests)](https://github.com/FormingWorlds/JANUS/actions/workflows/tests.yaml)
[![Integration Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/JANUS/nightly.yml?branch=main&label=Integration%20Tests)](https://github.com/FormingWorlds/JANUS/actions/workflows/nightly.yml)

[![total tests](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/FormingWorlds/JANUS/badges/tests-total.json)](https://proteus-framework.org/validation)
[![unit tests](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/FormingWorlds/JANUS/badges/tests-unit.json)](https://proteus-framework.org/validation)
[![integration tests](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/FormingWorlds/JANUS/badges/tests-integration.json)](https://proteus-framework.org/validation)

# JANUS

**JANUS** is the 1D convective-radiative atmosphere module of the [PROTEUS](https://proteus-framework.org/PROTEUS) coupled atmosphere-interior evolution framework. It builds the vertical temperature-pressure structure of a rocky-planet or magma-ocean atmosphere and computes its radiative fluxes.

Given a surface temperature and pressure, a set of volatile partial pressures, the instellation, and a spectral file, JANUS constructs the profile along a generalised multi-species moist pseudoadiabat capped by a prescribed stratosphere, then evaluates the shortwave, longwave, and net fluxes with [SOCRATES](https://github.com/FormingWorlds/SOCRATES). It returns the temperature-pressure profile, the tropopause temperature, and the top-of-atmosphere and surface fluxes. Within PROTEUS it is called each step through the `atmos_clim` interface; it is also usable standalone for atmosphere-structure studies.

Named for [Janus](https://en.wikipedia.org/wiki/Janus), the two-faced Roman god of transitions and doorways. Pronounced *jan-us*: *Jan* as in "january", *us* as in the collective pronoun.

## Atmosphere model

- **Generalised moist pseudoadiabat.** The troposphere follows the multi-species pseudoadiabat of Graham et al. (2021), which lets several volatiles condense together, capped by a prescribed stratosphere above a diagnosed tropopause.
- **Radiative transfer.** Shortwave, longwave, and net fluxes are computed with the SOCRATES correlated-k code over a user-selected spectral file and stellar spectrum.
- **Two entry points.** `RadConvEqm` sets the atmosphere to a temperature profile along the general adiabat, while `MCPA_CBL` builds the multiple-condensible pseudoadiabat and steps the surface temperature to conserve energy across a conductive boundary layer.

## Documentation

Full documentation is at **[proteus-framework.org/JANUS](https://proteus-framework.org/JANUS/)**, including:

- [Getting started](https://proteus-framework.org/JANUS/getting_started.html): installation and the quickest path to a first run.
- [Tutorial](https://proteus-framework.org/JANUS/Tutorials/first_run.html): a first atmosphere-structure calculation.
- [How-to guides](https://proteus-framework.org/JANUS/How-to/installation.html): install, run the tests, build the documentation.
- [Explanations](https://proteus-framework.org/JANUS/Explanations/model.html): model overview and the testing suite.
- [Validation anchors](https://proteus-framework.org/JANUS/Validation/phys.html): the per-source reference-pinned test inventory.
- [Publications](https://proteus-framework.org/JANUS/Reference/publications.html): the papers that developed and applied JANUS.

## Installation

> **Note:** The standard way of installing JANUS is within the PROTEUS framework, as described in the [PROTEUS installation guide](https://proteus-framework.org/PROTEUS/How-to/installation.html#10-install-submodules-as-editable). The steps below are for a standalone installation.

JANUS has one compiled dependency, [SOCRATES](https://github.com/FormingWorlds/SOCRATES), which a helper script clones and builds. **Prerequisites:** `git`, Python 3.10+, a Fortran/C build toolchain (`gfortran`, `gcc`, `make`), and NetCDF tools with NetCDF-Fortran development headers.

```console
git clone https://github.com/FormingWorlds/JANUS.git
cd JANUS

# 1. Build SOCRATES and point RAD_DIR at it (omit the path to use ./socrates)
bash tools/get_socrates.sh /path/to/socrates
export RAD_DIR=/path/to/socrates

# 2. Install JANUS (add ,docs for a local documentation build)
pip install -e .[develop]

# 3. Point FWL_DATA at a data directory and download spectral + stellar data
export FWL_DATA=/path/to/fwl_data
janus download spectral
janus download stellar
```

`janus env` prints the resolved locations of `RAD_DIR` and `FWL_DATA`. The `docs` extra pulls in [Zensical](https://zensical.org/) so you can build this documentation locally with `zensical serve`.

## Quick start

The two entry points are importable directly from `janus.modules`:

```python
from janus.modules import RadConvEqm, MCPA_CBL
```

The `examples/` folder holds runnable scripts that set up an atmosphere, prepare a spectral file, and solve for the profile and fluxes:

```console
python examples/demo_runaway_greenhouse.py   # runaway-greenhouse OLR curve
python examples/demo_instellation.py         # profile vs instellation
```

See the [first-run tutorial](https://proteus-framework.org/JANUS/Tutorials/first_run.html) for the full walkthrough.

## Citation

If you use JANUS in published work, please cite the methods papers below. The full reference list is on the [Publications page](https://proteus-framework.org/JANUS/Reference/publications.html).

- Graham, R.J., Lichtenberg, T., Boukrouche, R., & Pierrehumbert, R.T. (2021). *A multispecies pseudoadiabat for simulating condensable-rich exoplanet atmospheres.* **PSJ** 2, 207. [\[DOI\]](https://doi.org/10.3847/PSJ/ac214c)
- Lichtenberg, T., Bower, D.J., Hammond, M., Boukrouche, R., Sanan, P., Tsai, S.-M., & Pierrehumbert, R.T. (2021). *Vertically resolved magma ocean-protoatmosphere evolution.* **JGR Planets** 126, e2020JE006711. [\[DOI\]](https://doi.org/10.1029/2020JE006711)
- Boukrouche, R., Lichtenberg, T., & Pierrehumbert, R.T. (2021). *Beyond runaway: initiation of the post-runaway greenhouse state on rocky exoplanets.* **ApJ** 919, 130. [\[DOI\]](https://doi.org/10.3847/1538-4357/ac1345)

## License

[Apache License 2.0](LICENSE.md). JANUS is part of the [PROTEUS framework](https://proteus-framework.org/).
