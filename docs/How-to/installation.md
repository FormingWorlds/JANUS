# Installation

This page describes a manual developer installation of JANUS and SOCRATES.

!!! note
    The standard way of installing JANUS is within the PROTEUS Framework, as described in the [PROTEUS installation guide](https://proteus-framework.org/PROTEUS/How-to/installation.html#10-install-submodules-as-editable). Only use the guide below for a standalone installation of JANUS. 

!!! info "Prerequisites"
    - `git`
    - Python 3.10+ (recommended: 3.11)
    - a Fortran/C build toolchain (`gfortran`, `gcc`, `make`)
    - NetCDF tools and NetCDF-Fortran development headers/libraries
    - optional but recommended: Conda

    Typical installation commands:

    === "Ubuntu / Debian"

        ```bash
        sudo apt install libnetcdff-dev gfortran
        ```

    === "Fedora / RedHat"

        ```bash
        sudo dnf install gcc gcc-gfortran gcc-c++ netcdf netcdf-fortran netcdf-fortran-devel \
            lapack lapack-devel lapack-static
        ```

    === "macOS (Homebrew)"

        ```bash
        brew install git gcc netcdf netcdf-fortran
        ```

    === "macOS (MacPorts)"

        ```bash
        sudo port install netcdf-fortran +gcc8 wget gcc13
        ```

## 0. Create a Conda environment [optional]

Using a dedicated Conda environment helps avoid dependency conflicts.

```console
conda create -n janus python=3.11 -y
conda activate janus
```

## 1. Install SOCRATES

JANUS uses the [SOCRATES](https://github.com/nichollsh/SOCRATES) radiative transfer
code as an external dependency. A helper script is provided to clone and build it:

```console
bash tools/get_socrates.sh /path/to/socrates
```

This clones the repository, runs `./configure`, and compiles the code. Once done,
add the following to your `~/.bashrc` (or `~/.zshrc`):

```console
export RAD_DIR=/path/to/socrates
```

!!! tip "Custom install path"
    The script accepts an optional path argument. If omitted, SOCRATES is cloned
    into a `socrates/` subdirectory of the current working directory.

!!! tip "Running tests with a custom SOCRATES path"
    You can override `RAD_DIR` for a single command without modifying your shell profile:

    ```console
    RAD_DIR=/path/to/socrates pytest
    ```

## 2. Install JANUS

```console
git clone git@github.com:FormingWorlds/JANUS.git
cd JANUS
pip install -e .
```

## 3. Download JANUS data

Data is downloaded from the [OSF repository](https://osf.io/vehxg/).

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

## 4. Verify installation

Inspect active environment paths with:

```console
janus env
```

This prints the resolved locations of `RAD_DIR` and `FWL_DATA`.

---

## Environment variables

### `RAD_DIR`

Points to the compiled SOCRATES directory. JANUS will not function without this.
Set it in your shell profile as shown in step 1, or export it inline:

```console
export RAD_DIR=/path/to/socrates
```

### `FWL_DATA`

Controls where JANUS stores downloaded spectral and stellar data. By default the
location follows the [XDG base directory specification](https://specifications.freedesktop.org/basedir-spec/latest/).
Override it with:

```console
export FWL_DATA=/path/to/fwl_data
```