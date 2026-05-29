[![Documentation](https://github.com/FormingWorlds/JANUS/actions/workflows/docs.yaml/badge.svg)](https://proteus-framework.org/JANUS/)
[![Tests](https://github.com/FormingWorlds/JANUS/actions/workflows/tests.yaml/badge.svg?branch=main)](https://github.com/FormingWorlds/JANUS/actions/workflows/tests.yaml)
![Coverage](https://gist.githubusercontent.com/stefsmeets/99391a66bb9229771504c3a4db611d05/raw/covbadge.svg)

## JANUS (1D convective atmosphere model)

Generates a temperature profile using the generalised moist pseudoadiabat and a prescribed stratosphere. Calculates radiative fluxes using SOCRATES.

Pronounced *jan-us*. *Jan* as in "january", and *us* as in the collective pronoun.

### Documentation

The documentation for JANUS can be found [here](https://proteus-framework.org/JANUS/).

### Repository structure

* `README.md`           - This file
* `src/janus/data/`     - Janus data files
* `src/janus/modules/`  - Utility python scripts
* `src/janus/utils/`    - Utility python scripts
* `examples/`           - Typical use scripts
* `tools/`              - Useful tools

### Installation

> **Note:** The standard way of installing JANUS is within the PROTEUS Framework, as described in the [PROTEUS installation guide](https://proteus-framework.org/PROTEUS/How-to/installation.html#10-install-submodules-as-editable). The steps below are for a standalone developer installation only.

**Prerequisites:** `git`, Python 3.10+ (3.11 recommended), a Fortran/C build toolchain (`gfortran`, `gcc`, `make`), and NetCDF tools with NetCDF-Fortran development headers. See [detailed installation instructions](https://proteus-framework.org/JANUS/How-to/installation.html) for more information.

#### 0. Create a Conda environment (optional)
```console
conda create -n janus python=3.11 -y
conda activate janus
```

#### 1. Install SOCRATES

A helper script clones and builds [SOCRATES](https://github.com/FormingWorlds/SOCRATES):
```console
bash tools/get_socrates.sh /path/to/socrates
```
If the path argument is omitted, SOCRATES is cloned into a `socrates/` subdirectory of the current working directory. Once built, add the following to your `~/.bashrc` or `~/.zshrc`:
```console
export RAD_DIR=/path/to/socrates
```

#### 2. Install JANUS
```console
git clone git@github.com:FormingWorlds/JANUS.git
cd JANUS
pip install -e .
```

#### 3. Download JANUS data

Set `FWL_DATA` to define where spectral and stellar data are stored:
```console
export FWL_DATA=/path/to/fwl_data
```
Download the default datasets:
```console
janus download spectral
janus download stellar
```
To download a specific spectral dataset with a given number of bands:
```console
janus download spectral -n Frostflow -b 4096
```

#### 4. Verify installation
```console
janus env
```
This prints the resolved locations of `RAD_DIR` and `FWL_DATA`.

### Run instructions

In the `examples/` folder you can find Python scripts showing typical use cases and workflows for atmosphere modelling with JANUS.