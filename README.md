[![Documentation Status](https://readthedocs.org/projects/fwl-janus/badge/?version=latest)](https://fwl-janus.readthedocs.io/en/latest/?badge=latest)
![Coverage](https://gist.githubusercontent.com/stefsmeets/99391a66bb9229771504c3a4db611d05/raw/covbadge.svg)

## JANUS (1D convective atmosphere model)

Generates a temperature profile using the generalised moist pseudoadiabat and a prescribed stratosphere. Calculates radiative fluxes using SOCRATES.

Pronounced *jan-us*. *Jan* as in "january", and *us* as in the collective pronoun.

### Documentation
https://fwl-proteus.readthedocs.io

## Contributors

| Name  | Email address |
| -     | -             |
Tim Lichtenberg         | tim.lichtenberg[at]rug.nl |
Harrison Nicholls       | harrison.nicholls[at]physics.ox.ac.uk |
Laurent Soucasse        | l.soucasse[at]esciencecenter.nl |
Stef Smeets             | s.smeets[at]esciencecenter.nl |
Mark Hammond            | mark.hammond[at]physics.ox.ac.uk |
RJ Graham               | arejaygraham[at]uchicago.edu |
Raymond Pierrehumbert   | raymond.pierrehumbert[at]physics.ox.ac.uk |
Ryan Boukrouche         | ryan.boukrouche[at]astro.su.se |
Hamish Innes            | hamish.innes[at]fu-berlin.de |


### Repository structure
* `README.md`           - This file
* `src/janus/data/`     - Janus data files
* `src/janus/modules/`  - Utility python scripts
* `src/janus/utils/`    - Utility python scripts
* `examples/`           - Typical use scripts
* `tools/`              - Useful tools

### Developer installation instructions
1. Download and install Socrates
```console
git clone git@github.com:nichollsh/SOCRATES.git
cd SOCRATES
./configure
./build-code
source set_rad_env
cd ..
```
2. Download and install Janus
```console
git clone git@github.com:FormingWorlds/JANUS.git
cd JANUS
pip install -e .
```
3. Download data from the [OSF repository](https://osf.io/vehxg/)
    * Set the environment variable FWL_DATA to define where the spectral data files will be stored
        * `export FWL_DATA=...`
    * Run the following commands to download all basic data
        * `janus download spectral`
        * `janus download stellar`
    * Alternatively, you can specify which spectral data you want to download, and optionally the number of bands
        * `janus download spectral -n Frostflow -b 4096`

### Run instructions
In the examples folder you can find python scripts showing typical usecases/workflows of atmosphere modelling with Janus.
