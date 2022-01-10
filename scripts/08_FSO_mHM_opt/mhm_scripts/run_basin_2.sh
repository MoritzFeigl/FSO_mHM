#!/bin/bash
PROJECTPATH=/home/lv71468/mfeigl/FSO_mHM/mHM_setup
I=6335600
cd $PROJECTPATH/FSO_mHM_major_basins/config/sub_$I
module load intel/18 intel-mpi/2018 hdf5/1.8.12-MPI pnetcdf/1.10.0 netcdf_C/4.4.1.1 netcdf_Fortran/4.4.4 cmake/3.9.6
./mhm > output.txt
