
## Access & installation instructions

1. Install dependencies


    * MacOS: Ensure installation of command line developer tools
        * Open Xcode application and install
        * Type `xcode-select --install`

    * Install FORTRAN NetCDF library via the most appropriate method for you
        * `brew install netcdf`    
        * `brew install netcdf-fortran`     
        OR    
        * `sudo port install netcdf-fortran +gcc8`    
        OR     
        * `sudo apt install libnetcdff-dev`
    
    * Setup a Python environment:
        * Option A: Using the `Anaconda` package manager (careful, probably breaks on ARM machines/newer Macs)
            * Install `conda`:
                * Download the appropriate Miniconda installer from their website
            https://docs.conda.io/en/latest/miniconda.html#id36
                * Create a conda environment for PROTEUS:
                * `conda create -n aeolus python=3.10.9`    
                * `conda activate aeolus`
            * `conda install netcdf4 matplotlib numpy pandas scipy sympy natsort`
            * `conda install -c conda-forge f90nml`
        * Option B: using the `brew` package manager (*recommended*)
            * Delete all traces of Anaconda package manager from your system and switch to a different Python environment, for example brew/pip
            * Follow the steps at https://docs.anaconda.com/free/anaconda/install/uninstall/
            * Delete all Anaconda-related entries from your .bash_profile (Intel) or .zshrc (ARM)
            * Install Python via `brew`: 
                * `brew install python`
                * Update to the latest stable version: `brew upgrade python`
                * Install `tkinter`: `brew install python-tk@3.11`
                * Refresh your shell / `source ~/.zsrhrc` (ARM) / `source ~/.bash_profile` (Intel)
                * Install all necessary packages: `pip3 install matplotlib pandas netcdf4 matplotlib numpy pandas scipy sympy natsort`
            * Make the new Python version the system default (check what `brew` tells you during/after the `brew install python` step):
                * ARM: `export PATH="/opt/homebrew/opt/python/libexec/bin:$PATH"`
                * Intel: `export PATH="/usr/local/opt/python/libexec/bin:$PATH"`

    * *Optionally* register your public SSH key with Github:
        * Generate key with `ssh-keygen -t rsa`
        * Copy to clipboard using `cat ~/.ssh/id_rsa.pub`
        * Add "New SSH key" in the settings at https://github.com/settings/keys 

2. Ensure that you have access to all of the following codes
    * Radiative-convective scheme: **AEOLUS** 
        * URL: https://github.com/FormingWorlds/AEOLUS/
        * Contact: TL, MH, RB

    * Radiation transport: **SOCRATES** 
        * Main development URL: https://code.metoffice.gov.uk/trac/socrates
        * Contact: james.manners@metoffice.gov.uk
        * Obtain the SOCRATES source code from: https://simplex.giss.nasa.gov/gcm/ROCKE-3D/ (Latest released version of SOCRATES code)
        * Latest tested version: *socrates_2211.tar.xz*

3. Setup codes and modules in the following order (ignore their individual README files)

    1. Download AEOLUS (*AEOLUS, SPIDER, VULCAN*)
        * `git clone  git@github.com:FormingWorlds/AEOLUS.git`

    2. Download and extract the SOCRATES archive to the correct location
        * `cd AEOLUS/rad_trans/socrates_code/`
        * `curl -L -o ../socrates_2211.tar.xz https://www.dropbox.com/sh/ixefmrbg7c94jlj/AAChibnZU9PRi8pXdVxVbdj3a/socrates_2211.tar.xz?dl=1`
        * `tar --strip-components 1 -xvf PATH_TO_ARCHIVE -C ./`
        * `cp -f ../build_code_modified build_code`

    4. Overwrite the `Mk_cmd` file with the right setup for your machine
        * `cp -rf ../Mk_cmd_SYSTEM make/Mk_cmd`    
        * Options are:
            * `cp -rf ../Mk_cmd_MAC_INTEL make/Mk_cmd`
            * `cp -rf ../Mk_cmd_MAC_APPLESILICON make/Mk_cmd`
            * `cp -rf ../Mk_cmd_AOPP_CLUSTER make/Mk_cmd`
            
        The command `nf-config` might be helpful if none of these options work for you.

    5. Setup SOCRATES 
        * `./build_code`
        * `type ksh >/dev/null 2>&1 || sed -i 's/ksh/bash/g' sbin/* `
        * `cd ../../../`

    6. Setup environment variables
        * `source AEOLUS.env`
 
**Done!**
