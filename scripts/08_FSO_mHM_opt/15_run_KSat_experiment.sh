#!/bin/sh
#SBATCH -J ksat_exp
#SBATCH -N 1
#SBATCH --qos=vsc3plus_0256
#SBATCH --partition=vsc3plus_0256
#SBATCH -L intel@vsc
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=15_run_KSat_experiment.out
module purge
module load intel/18 intel-mkl/2018 pcre2/10.35 R/4.0.2
module load intel/18 intel-mpi/2018 hdf5/1.8.12-MPI pnetcdf/1.10.0 netcdf_C/4.4.1.1 netcdf_Fortran/4.4.4
cd /home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt
Rscript 15_run_KSat_experiment.R
