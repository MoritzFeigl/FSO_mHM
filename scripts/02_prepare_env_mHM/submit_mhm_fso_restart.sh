#!/bin/bash

# ----------------------------------------------------------------------
# qsub arguments
# ----------------------------------------------------------------------

#$ -N mhm_spinup
#$ -S /bin/bash
#$ -wd /work/$USER
#$ -o /work/$USER/stdout/$JOB_NAME-$JOB_ID-$TASK_ID.log
#$ -e /work/$USER/stderr/$JOB_NAME-$JOB_ID-$TASK_ID.log
#$ -j y
#$ -l h_rt=00:20:00
#$ -l h_vmem=8G
#$ -binding linear:1
#$ -cwd

# call it via 'qsub -t 1-$(wc -l /path/to/input.files) array-job.sh /path/to/input.files'
echo qsub arguments

# ----------------------------------------------------------------------
# load required modules
# ----------------------------------------------------------------------
echo load required modules


# load python environment
sys_name=$(uname -n)
if [[ ${sys_name:0:8} != 'frontend' ]]; then
  module purge
  module load foss/2019b
  module load netCDF-Fortran
  module load CMake
  module use /global/apps/modulefiles
  module load python_env_mpr
  module load git-subrepo
  module load pFUnit/4.0.0_foss_2019b
  export NETCDF_DIR="$EBROOTNETCDF"
  export NETCDF_FORTRAN_DIR="$EBROOTNETCDFMINFORTRAN"
  export FC=mpifort
fi

# ----------------------------------------------------------------------
# main script
# ----------------------------------------------------------------------
echo main script


# get all basins
# find <path>/sub* -type d -maxdepth 1 -printf "%f\n" | tee input.files

# get current basin by indexing
# get list of input files from first argument
job_list=$1

# get the specific input file for this array job task
# this is the n-th (where n is current task ID) line of the file
folder=$(awk "NR==$SGE_TASK_ID" "$job_list")

# cd
cd /work/ottor/FSO_mHM_major_basins/config/"${folder}" || exit
# execute
./mhm
