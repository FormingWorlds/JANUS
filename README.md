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
* `README.md`                          - This file
* `JANUS.env`                          - Sets environment flags to run the code
* `src/janus/SocRadConv.py`            - Main JANUS Python script
* `src/janus/data/luminosity_tracks/`  - Stellar evolution data
* `src/janus/data/spectral_files/`     - Spectral files for SOCRATES
* `src/janus/modules/`                 - Utility python scripts
* `src/janus/plotting_tools/`          - Plotting scripts
* `src/janus/utils/`                   - Utility python scripts
* `tests/demo_runaway_greenhouse.py`   - Demonstrate pure-steam runaway greenhouse OLR curve
* `tests/demo_instellation.py`         - Calculate fluxes (and temperatures) for different instellations
* `tools/`                             - Useful tools

### Developer installation instructions
1. Download and install Socrates
    * `git clone git@github.com:FormingWorlds/SOCRATES.git ./my_socrates`
    * `cd my_socrates`
    * `./configure`
    * `./build-code`
    * `source set_rad_env`
2. Download and install Janus
    * `git clone git@github.com:FormingWorlds/JANUS.git ./my_janus`
    * `cd my_janus`
    * `pip install -e .`

### Run instructions
Only attempt to run JANUS after you have followed all of the instructions in INSTALL.md    
If using a fresh shell, it is necessary to perform the following steps:     
1. `source JANUS.env`    
2. `conda activate janus`    
Then you can run the code by running: `python SocRadConv.py`      
