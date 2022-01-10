#!/bin/sh
#SBATCH -J val
#SBATCH -N 1
#SBATCH --qos=vsc3plus_0064
#SBATCH --partition=vsc3plus_0064
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=FAIL
#SBATCH --output=val.out
module purge
module load intel/18 intel-mpi/2018 hdf5/1.8.12-MPI pnetcdf/1.10.0 netcdf_C/4.4.1.1 netcdf_Fortran/4.4.4 cmake/3.9.6
