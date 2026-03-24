# Installation

This page describes a manual developer installation of JANUS and SOCRATES. JANUS contains a small CLI tool to help get set up with JANUS.

!!! info "Prerequisites"
	- `git`
	- Python 3.10+ (recommended: 3.11)
	- a Fortran/C build toolchain (`gfortran`, `gcc`, `make`)
	- NetCDF tools and NetCDF-Fortran development headers/libraries
	- optional but recommended: Conda 

	Typical installation commands:

	=== "Ubuntu / Debian"

		```console
		sudo apt install libnetcdff-dev gfortran
		```

	=== "Fedora / RedHat"
		```console
		sudo dnf install gcc gcc-gfortran gcc-c++ netcdf netcdf-fortran netcdf-fortran-devel \
    	lapack lapack-devel lapack-static sundials-mpich openmpi openmpi-devel f2c f2c-libs
		```

	=== "macOS (Homebrew)"

		```console
		brew install git gcc netcdf netcdf-fortran
		```

	=== "macOS (MacPorts)"

		```console
		sudo port install netcdf-fortran +gcc8 wget gcc13 openmpi
		```

## 0. Create a Conda environment [optional]

Using a dedicated Conda environment helps avoid dependency conflicts.

```console
conda create -n janus python=3.11 -y
conda activate janus
```

If your system does not already provide NetCDF-Fortran development tools,
install them before building SOCRATES as shown above.

## 1. Install SOCRATES

```console
git clone git@github.com:nichollsh/SOCRATES.git
cd SOCRATES
./configure
./build-code
source set_rad_env
cd ..
```

	If you install and compile [SOCRATES](https://github.com/nichollsh/SOCRATES) yourself,
	you can override the path using the `SOCRATES` environment variable, e.g.

	```console
	RAD_DIR=/home/user/path/to/SOCRATES pytest
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

To download a specific spectral dataset (optionally with a band count):

```console
janus download spectral -n Frostflow -b 4096
```

## 4. Verify installation

You can inspect active JANUS environment paths with:

```console
janus env
```

### `SOCRATES`



### `FWL_DATA`

Set this variable to modify where janus stores its stellar and spectral data. By default this is based on the [XDG specification](https://specifications.freedesktop.org/basedir-spec/latest/).
You can override the path using the `FWL_DATA` environment variable, e.g.

```console
FWL_DATA=/home/user/path/to/fwl_data pytest
```
