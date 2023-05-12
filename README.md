## AEOLUS (radiative-convective solver)

Runs until the OLR changes by less than a threshold value in W/m<sup>2</sup>, or stops after a fixed number of iterations.

### Documentation
https://proteus-code.readthedocs.io

### Contributors (abbreviations & email addresses):
* TL - Tim Lichtenberg (tim.lichtenberg@rug.nl)
* MH - Mark Hammond (mark.hammond@physics.ox.ac.uk)
* RB – Ryan Boukrouche (ryan.boukrouche@astro.su.se)
* RJG – RJ Graham (arejaygraham@uchicago.edu)
* HN - Harrison Nicholls (harrison.nicholls@physics.ox.ac.uk)

### Repository structure

* `SocRadConv.py`               - Main AEOLUS Python script
* `README.md`                   - This file
* `INSTALL.md`                  - Installation instructions
* `AEOLUS.env`                  - Sets environment flags to run the code
* `MoistAdiabat_RayCode.py`     - Legacy file
* `utils/`                      - Utility python scripts
* `modules/`                    - Utility python scripts
* `output/`                     - Output folder
* `luminosity_tracks/`          - Stellar evolution data
* `plotting_tools/`             -  Plotting scripts
* `rad_trans/`                  - Subfolder with adiative transfer code/s


* `spectral_files`         - Spectral files for SOCRATES

### Installation instructions
See `INSTALL.md` for steps. If you encounter issues check `INSTALL.md` and `TROUBLESHOOTING.md` in [PROTEUS](https://github.com/FormingWorlds/PROTEUS) repository. `INSTALL.md` will soon be superseded by https://proteus-code.readthedocs.io.

### Run instructions
Only attempt to run AEOLUS after you have followed all of the instructions in INSTALL.md    
If using a fresh shell, it is necessary to perform the following steps:     
1. `source AEOLUS.env`
2. `conda activate aeolus`
Then you can run the code by running: `python SocRadConv.py`      
