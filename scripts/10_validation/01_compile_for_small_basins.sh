#!/bin/sh
# FSO mHM optimization mHM runs
# define working directory
PROJECTPATH=/home/lv71468/mfeigl/FSO_mHM/Runs/small_basins_validation

# get tf update
# get python env
cd $PROJECTPATH
source fso_mhm_py3/bin/activate
cd $PROJECTPATH/mhm/deps/mpr
python -m src_python.pre_proc.update_tfs_in_fortran_source -c /home/lv71468/mfeigl/FSO_mHM/mHM_setup/FSO_mHM_all_basins/config/sub_6335304/mpr.nml -p /home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_master_namelists/mpr_global_parameter_mhm.nml --clean
deactivate

# Load modules
module load intel/18 intel-mpi/2018 hdf5/1.8.12-MPI pnetcdf/1.10.0 netcdf_C/4.4.1.1 netcdf_Fortran/4.4.4 cmake/3.9.6
# compile mhm
now=$(date +"%T")
echo -e "Compiling time start: $now"
cd $PROJECTPATH/mhm
rm -r build
mkdir build
cd build
cmake -DCMAKE_WITH_OpenMP:STRING=ON -DCMAKE_BUILD_TYPE=Release -DNETCDF_LIBRARY=/opt/sw/x86_64/glibc-2.17/ivybridge-ep/netcdf_C/4.4.1.1/intel/18/intel-mpi/2018/hdf5/1.8.12-MPI/pnetcdf/1.10.0/lib/libnetcdf.so -DNETCDF_F90_LIBRARY=/opt/sw/x86_64/glibc-2.17/ivybridge-ep/netcdf_Fortran/4.4.4/intel/18/intel-mpi/2018/hdf5/1.8.12-MPI/pnetcdf/1.10.0/netcdf_C/4.4.1.1/lib/libnetcdff.so ..
make -j 8
now=$(date +"%T")
echo -e "Compiling time end: $now"
