#!/bin/sh
# FSO mHM optimization mHM runs
# define working directory
PROJECTPATH=/home/lv71468/mfeigl/FSO_mHM/mHM_setup

# get tf update
# get python env
cd $PROJECTPATH
source fso_mhm_py3/bin/activate
cd $PROJECTPATH/mhm/deps/mpr
python -m src_python.pre_proc.update_tfs_in_fortran_source -c $PROJECTPATH/FSO_mHM_major_basins/config/sub_6335304/mpr.nml -p /home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_master_namelists/mpr_global_parameter_mhm.nml --clean   > /dev/null 2>&1
deactivate

# Load modules
module load gcc/7.3.0 intel-mpi/2018 pnetcdf/1.10.0 hdf5/1.8.18-MPI cmake/3.9.6 netcdf_C/4.4.1.1 netcdf_Fortran/4.4.4  > /dev/null 2>&1
# compile mhm
now=$(date +"%T")
echo -e "Compiling time start: $now"
cd $PROJECTPATH/mhm/build
make -j 2  > /dev/null 2>&1
now=$(date +"%T")
echo -e "Compiling time end: $now"

# copy mhm into basin folders
for I in 6335304 6335600 6337200 6338100 6340300 6340600 6342800
do
  ln -sf $PROJECTPATH/mhm/build/mhm-0.5.9 $PROJECTPATH/FSO_mHM_major_basins/config/sub_$I/mhm-0.5.9
  ln -sf $PROJECTPATH/mhm/build/mhm $PROJECTPATH/FSO_mHM_major_basins/config/sub_$I/mhm
done
