#!/bin/bash
#SBATCH --job-name=hitr_03_16
#SBATCH --output=hitr_03_16.out
#SBATCH --error=hitr_03_16.err
#SBATCH --time=336:00:00
#SBATCH --ntasks=16
#SBATCH --mem=64000
#SBATCH -p priority-rp
#SBATCH --cpus-per-task=1
date
export OMP_NUM_THREADS=16
module load intel-compilers netcdf/netcdf-c-4.7.3 netcdf/netcdf-fortran-4.5.2 openmpi/4.0.1-intel
source /home/lichtenberg/codes/socrates/socrates_2002/set_rad_env
python3 sp_gen_script.py
date
