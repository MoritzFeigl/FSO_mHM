#!/bin/bash
PROJECTPATH=/home/lv71468/mfeigl/FSO_mHM/mHM_setup
I=6335600
cd $PROJECTPATH/FSO_mHM_major_basins/config/sub_$I
module load gcc/7.3.0 intel-mpi/2018 pnetcdf/1.10.0 hdf5/1.8.18-MPI cmake/3.9.6 netcdf_C/4.4.1.1 netcdf_Fortran/4.4.4
./mhm > output.txt
