## AEOLUS (Radiative-convective solver)

Runs until the OLR changes by less than a threshold value in W/m<sup>2</sup>, or stops after a fixed number of iterations.

### Contributors (abbreviations & email addresses):
* TL - Tim Lichtenberg (tim.lichtenberg@rug.nl)
* MH - Mark Hammond (mark.hammond@physics.ox.ac.uk)
* RTP - Ray Pierrehumbert (raymond.pierrehumbert@physics.ox.ac.uk)
* RB – Ryan Boukrouche (ryan.boukrouche@astro.su.se)
* HN - Harrison Nicholls (harrison.nicholls@physics.ox.ac.uk)

### Repository structure

* GGRadConv.py - Grey gas radiative convective model
* GreyHeat.py - Calculates heating rates in grey model
* SocRadConv.py - Radiative convective model w/ SOCRATES
* SocRadModel.py - Calculates heating rates in for SOCRATES version
* atmosphere_column.py - Class for atmospheric column data
* nctools.py - Some useful functions for netCDF
* phys.py - Constants
* planets.py - Planet-specific constants
* README.md - This file
* INSTALL.md - Installation instructions
* surfaceT.txt - Surface temperature form interior model, coupler-file

### Installation instructions
See INSTALL.md for steps.

### Run instructions
Only attempt to run AEOLUS after you have followed all of the instructions in INSTALL.md    
If using a fresh shell, it is necessary to perform the following steps:     
1. `source AEOLUS.env`
2. `conda activate aeolus`
Then you can run the code by running: `python SocRadConv.py`      
