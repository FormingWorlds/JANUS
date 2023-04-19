
## Access & installation instructions

1. Install dependencies

    * Setup a conda environment

        * Install conda
            Download the appropriate Miniconda installer from their website
            https://docs.conda.io/en/latest/miniconda.html#id36

        * Create a conda environment for PROTEUS
            `conda create -n aeolus python=3.10.9`    
            `conda activate aeolus`
            
    * Install FORTRAN NetCDF library via the most appropriate method for you
        * `brew install netcdf`    
        * `brew install netcdf-fortran`     
        OR    
        * `sudo port install netcdf-fortran +gcc8`    
        OR     
        * `sudo apt install libnetcdff-dev`
    
    * Install Python libraries:
        * `conda install netcdf4 matplotlib numpy pandas scipy seaborn natsort`
        * `conda install -c conda-forge f90nml`

    * *Optionally* register your public SSH key with Github:
        * Generate key with `ssh-keygen -t rsa`
        * Copy to clipboard using `cat ~/.ssh/id_rsa.pub`
        * Add "New SSH key" in the settings at https://github.com/settings/keys 

2. Ensure that you have access to all of the following codes
    * Radiative-convective scheme: **AEOLUS** 
        * URL: https://github.com/FormingWorlds/AEOLUS/
        * Contact: TL, MH, RB

    * Radiation transport: **SOCRATES** 
        * URL: https://code.metoffice.gov.uk/trac/socrates
        * Contact: james.manners@metoffice.gov.uk
        * Latest tested version: *socrates_2211.tar.xz*

3. Setup codes and modules in the following order (ignore their individual README files)

    1. Download AEOLUS (*AEOLUS, SPIDER, VULCAN*)
        * `git clone  git@github.com:FormingWorlds/AEOLUS.git`

    2. Enter into AEOLUS folder
        * `cd AEOLUS`

    3. Extract SOCRATES archive to the correct location
        * `cd rad_trans/socrates_code/`
        * `tar --strip-components 1 -xvf PATH_TO_ARCHIVE -C ./`
        * `cp -f ../build_code_modified build_code`

    4. Overwrite the `Mk_cmd` file with the right setup for your machine
        * `cp -rf ../Mk_cmd_SYSTEM make/Mk_cmd`    
        Options are: *Mk_cmd_MAC_INTEL*, *Mk_cmd_MAC_APPLESILICON*, *Mk_cmd_AOPP_CLUSTER*.    
        The command `nf-config` might be helpful if none of these options work for you.

    5. Setup SOCRATES 
        * `./build_code`
        * `type ksh >/dev/null 2>&1 || sed -i 's/ksh/bash/g' sbin/* `
        * `cd ../../../`

    6. Setup environment variables
        * `source AEOLUS.env`
 
**Done!**
