# Getting started

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

### `SOCRATES`

By default, SOCRATES is installed to the default location based on the [XDG specification](https://specifications.freedesktop.org/basedir-spec/latest/).

If you install and compile [SOCRATES](https://github.com/nichollsh/SOCRATES) yourself,
you can override the path using the `SOCRATES` environment variable, e.g.

```console
SOCRATES=/home/user/path/to/SOCRATES pytest
```
