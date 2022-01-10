#!/bin/sh
# Preparational steps for FSO mHM optimization
PROJECTPATH=/home/lv71468/mfeigl/FSO_mHM/mHM_setup
# Load modules
module load intel/18 intel-mpi/2018 hdf5/1.8.12-MPI pnetcdf/1.10.0 netcdf_C/4.4.1.1 netcdf_Fortran/4.4.4 cmake/3.9.6
# Download correct mHM version
cd $PROJECTPATH
#rm -r mhm
echo -e "Cloning mhm from gitlab"
git clone -b FSO_project https://git.ufz.de/ottor/mhm.git
# compile mhm
echo -e "Compiling mhm"
cd $PROJECTPATH/mhm
mkdir build
cd build
cmake -DCMAKE_WITH_OpenMP:STRING=ON -DCMAKE_BUILD_TYPE=Release -DNETCDF_LIBRARY=/opt/sw/x86_64/glibc-2.17/ivybridge-ep/netcdf_C/4.4.1.1/intel/18/intel-mpi/2018/hdf5/1.8.12-MPI/pnetcdf/1.10.0/lib/libnetcdf.so -DNETCDF_F90_LIBRARY=/opt/sw/x86_64/glibc-2.17/ivybridge-ep/netcdf_Fortran/4.4.4/intel/18/intel-mpi/2018/hdf5/1.8.12-MPI/pnetcdf/1.10.0/netcdf_C/4.4.1.1/lib/libnetcdff.so ..
make -j 8
# save original source
cp -R $PROJECTPATH/mhm $PROJECTPATH/mhm_original
# run test basins
cp $PROJECTPATH/mhm/build/mhm-0.5.9 $PROJECTPATH/mhm/mhm-0.5.9
cp $PROJECTPATH/mhm/build/mhm $PROJECTPATH/mhm/mhm
cd $PROJECTPATH/mhm
./mhm
#cp -R $PROJECTPATH/mhm_original/mhm $PROJECTPATH
