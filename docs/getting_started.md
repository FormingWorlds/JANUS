# Install & Getting Started

JANUS contains a small CLI tool to help get set up with JANUS.

## Install SOCRATES

Download and install SOCRATES:

```console
janus download socrates
```

Make sure you have the netcdf fortran libraries installed:

```
sudo apt install libnetcdff-dev netcdf-bin
```

## Download data

Download spectral and stellar data:

```console
janus download spectral
janus download stellar
```

## Environment variables

To see all environment variables and locations, use:

```console
janus env
```

### `SOCRATES`

By default, SOCRATES is installed to the default location based on the [XDG specification](https://specifications.freedesktop.org/basedir-spec/latest/).

If you install and compile [SOCRATES](https://github.com/nichollsh/SOCRATES) yourself,
you can override the path using the `SOCRATES` environment variable, e.g.

```console
RAD_DIR=/home/user/path/to/SOCRATES pytest
```

### `FWL_DATA`

Set this variable to modify where janus stores its stellar and spectral data. By default this is based on the [XDG specification](https://specifications.freedesktop.org/basedir-spec/latest/).
You can override the path using the `FWL_DATA` environment variable, e.g.

```console
FWL_DATA=/home/user/path/to/fwl_data pytest
```
