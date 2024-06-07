## JANUS (temperature structure generator)

Generates a temperature profile using the generalised moist pseudoadiabat and a prescribed stratosphere. Calculates radiative fluxes using SOCRATES.   

Pronounced *jan-us*. *Jan* as in "january", and *us* as in the collective pronoun.

### Documentation
https://proteus-code.readthedocs.io

### Contributors (abbreviations & email addresses):
* TL - Tim Lichtenberg (tim.lichtenberg@rug.nl)
* MH - Mark Hammond (mark.hammond@physics.ox.ac.uk)
* RB – Ryan Boukrouche (ryan.boukrouche@astro.su.se)
* RJG – RJ Graham (arejaygraham@uchicago.edu)
* HN - Harrison Nicholls (harrison.nicholls@physics.ox.ac.uk)
* HII - Hamish Innes (hamish.innes@physics.ox.ac.uk)
* LS - Laurent Soucasse (l.soucasse@esciencecenter.nl)

### Repository structure
* `README.md`           - This file
* `src/janus/data/`     - Janus data files
* `src/janus/modules/`  - Utility python scripts
* `src/janus/utils/`    - Utility python scripts
* `examples/`           - Typical use scripts
* `tools/`              - Useful tools

### Developer installation instructions
1. Download and install Socrates
    * `git clone git@github.com:FormingWorlds/SOCRATES.git`
    * `cd SOCRATES`
    * `./configure`
    * `./build-code`
    * `source set_rad_env`
    * `cd ..`
2. Download and install Janus
    * `git clone git@github.com:FormingWorlds/JANUS.git`
    * `cd JANUS`
    * `pip install -e .`

### Run instructions
In the examples folder you can find python scripts showing typical usecases/workflows of atmosphere modelling with Janus.
