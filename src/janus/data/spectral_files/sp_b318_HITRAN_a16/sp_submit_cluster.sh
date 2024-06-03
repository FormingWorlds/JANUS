# Run prefix
prefix=${PWD##*/}

# Submit filename
SUBMIT="exec_submit.sh"

# Delete old submit file
rm -rf $SUBMIT

# Number of CPUs/talks
nprocs=16

echo "#!/bin/bash" >> $SUBMIT

# SLURM commands (man sbatch / http://www.arc.ox.ac.uk/content/slurm-job-scheduler)
echo "#SBATCH --job-name=$prefix" >> $SUBMIT
echo "#SBATCH --output=$prefix.out" >> $SUBMIT
echo "#SBATCH --error=$prefix.err" >> $SUBMIT
echo "#SBATCH --time=336:00:00" >> $SUBMIT
echo "#SBATCH --ntasks=$nprocs" >> $SUBMIT
echo "#SBATCH --mem=64000" >> $SUBMIT # memory per node
#echo "#SBATCH --mem-per-cpu=8000" >> $SUBMIT # memory per cpu/task
echo "#SBATCH -p priority-rp" >> $SUBMIT # Pierrehumbert priority queue
#echo "#SBATCH -p shared" >> $SUBMIT # Physics/AOPP shared queue
echo "#SBATCH --cpus-per-task=1" >> $SUBMIT

# Print date and time before job runs
echo "date" >> $SUBMIT

# Multithreading, should be equal to number of tasks
# Sets the number of threads in openMP (== number of tasks/CPUs)
echo "export OMP_NUM_THREADS=$nprocs" >> $SUBMIT

# Load environment
echo "module load intel-compilers netcdf/netcdf-c-4.7.3 netcdf/netcdf-fortran-4.5.2 openmpi/4.0.1-intel" >> $SUBMIT
echo "source /home/lichtenberg/codes/socrates/socrates_2002/set_rad_env" >> $SUBMIT

# Generate initial conditions and run the code
echo "python3 sp_gen_script.py" >> $SUBMIT

#Print date and time after job finished
echo "date" >> $SUBMIT

# Submit created script
sbatch ./$SUBMIT

