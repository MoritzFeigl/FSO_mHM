#!/bin/sh
# FSO mHM optimization mHM runs

PROJECTPATH=/home/lv71468/mfeigl/FSO_mHM/mHM_setup

# get python env
cd $PROJECTPATH
virtualenv --python python3 fso_mhm_py3
source fso_mhm_py3/bin/activate
pip install f90nml

# update tfs
cd $PROJECTPATH/mhm/deps/mpr
python -m src_python.pre_proc.update_tfs_in_fortran_source -c $PROJECTPATH/FSO_mHM_major_basins/config/sub_6335304/mpr.nml -p /home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_master_namelists/mpr_global_parameter_mhm.nml --clean
deactivate

# compile mhm
module load gcc/7.3.0 intel-mpi/2018 pnetcdf/1.10.0 hdf5/1.8.18-MPI cmake/3.9.6 netcdf_C/4.4.1.1 netcdf_Fortran/4.4.4
now=$(date +"%T")
echo -e "Compiling time start: $now"
cd $PROJECTPATH/mhm
rm -r build
mkdir build
cd build
cmake -DCMAKE_WITH_OpenMP:STRING=ON -DCMAKE_BUILD_TYPE=Release -DNETCDF_LIBRARY=/opt/sw/x86_64/glibc-2.17/ivybridge-ep/netcdf_C/4.4.1.1/gcc/7.3.0/intel-mpi/2018/hdf5/1.8.18-MPI/pnetcdf/1.10.0/lib/libnetcdf.so -DNETCDF_F90_LIBRARY=/opt/sw/x86_64/glibc-2.17/ivybridge-ep/netcdf_Fortran/4.4.4/gcc/7.3.0/intel-mpi/2018/hdf5/1.8.18-MPI/netcdf_C/4.4.1.1/lib/libnetcdff.so ..
make -j 8
now=$(date +"%T")
echo -e "Compiling time end: $now"
# copy mhm into basin folders
for I in 6335304 6335600 6337200 6338100 6340300 6340600 6342800
do
  ln -sf $PROJECTPATH/mhm/build/mhm-0.5.9 $PROJECTPATH/FSO_mHM_major_basins/config/sub_$I/mhm-0.5.9
  ln -sf $PROJECTPATH/mhm/build/mhm $PROJECTPATH/FSO_mHM_major_basins/config/sub_$I/mhm
done
