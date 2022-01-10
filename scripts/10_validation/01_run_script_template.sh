#!/bin/bash
PROJECTPATH=/home/lv71468/mfeigl/FSO_mHM/mHM_setup/FSO_mHM_all_basins
I=6335304
cd $PROJECTPATH/config/sub_$I
echo -e "running"
module load intel/18 intel-mpi/2018 hdf5/1.8.12-MPI pnetcdf/1.10.0 netcdf_C/4.4.1.1 netcdf_Fortran/4.4.4 cmake/3.9.6
./mhm > output.txt
