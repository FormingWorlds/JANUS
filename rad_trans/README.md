- Socrates source goes here. Path should be /rad_trans/socrates_code/
- Unzip, e.g., socrates_2002/ in ./socrates_code/
- Overwrite ./socrates_code/make/Mk_cmd with the appropriate Mk_cmd from this folder

Machine/cluster-specific instructions:

AOPP Cluster:
- Replace all /bin/ksh references with /bin/bash if you use bash shell (enjoy!)
  To find all such files, do (in socrates_code/) 
    $ grep -r /bin/ksh *
- Compile with the following modules
    $ module load intel-compilers netcdf/netcdf-c-4.7.3 netcdf/netcdf-fortran-4.5.2 openmpi/4.0.1-intel
- If you want to run Ccorr_k (spectral file generation) with openMP parallel support (-np option), make sure to compile with:
    $ OMPARG = -qopenmp

