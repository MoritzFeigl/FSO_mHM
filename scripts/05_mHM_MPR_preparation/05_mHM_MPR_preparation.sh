#!/bin/sh

# Script for compiling mHM and MPR and generating spatial predictors for FSO
# Refer to workflow.txt for more details.

# Colors
GREEN='\033[1;32m'
NC='\033[0m'
RED='\033[0;31m'
# Paths
PROJECTPATH=/home/lv71468/mfeigl/FSO_mHM/mHM_setup
MPR_PATH=/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/04_prepare_predictors_for_FSO
MHM_PATH=/home/lv71468/mfeigl/FSO_mHM/mHM_setup
echo -e $GREEN"***** Starting compilation tests and Spatial predictor generation *****"$NC
echo -e $RED"This code will be run in:"$NC
echo -e $PROJECTPATH
echo -e $RED"Major Basins data is assumed to be in:"$NC
echo -e "$MHM_PATH/FSO_mHM_major_basins"
cd $PROJECTPATH

# 1. Download mHM and MPR from gitlab
echo -e $GREEN"1. Download mHM and MPR from gitlab"$NC
if [ ! -d "mhm" ]
then
  git clone -b 22-couple-mhm-with-mpr-1-0 https://git.ufz.de/ottor/mhm.git
fi

if [ ! -d "MPR" ]
then
  git clone -b develop https://git.ufz.de/chs/MPR.git
fi

# 2. Check if TFs are registered
echo -e $GREEN"2. Check if TFs are registered"$NC
module purge
module load anaconda3/5.3.0
conda init bash
source ~/.bashrc
module load anaconda3/5.3.0
ENV_NAME="local"
if conda info --envs | grep ${ENV_NAME}; then
  echo -e "activating conda environment ${ENV_NAME}"
else
  echo -e "New conda environment ${ENV_NAME} gets initialized"
  conda create -n ${ENV_NAME} python=3.7
  return 1
fi
conda activate ${ENV_NAME}
pip install f90nml > /dev/null 2>&1
# source python skript for checking and updating TFS
cd $PROJECTPATH/MPR
python -m src_python.pre_proc.update_tfs_in_fortran_source -c $MPR_PATH/mpr_predictor_FSO_only.nml -p $MPR_PATH/mpr_global_parameter_mhm.nml
cd $PROJECTPATH/mhm/deps/mpr
python -m src_python.pre_proc.update_tfs_in_fortran_source -c $MPR_PATH/mpr_predictor_FSO_only.nml -p $MPR_PATH/mpr_global_parameter_mhm.nml
conda deactivate

# 3. Compile mHM and MPR
echo -e $GREEN"3. Compile mHM and MPR"$NC
module load intel/18 intel-mpi/2018 hdf5/1.8.12-MPI pnetcdf/1.10.0 netcdf_C/4.4.1.1 netcdf_Fortran/4.4.4 cmake/3.9.6

# Compile mHM
cd $PROJECTPATH/mhm # works for mHM
if [ -d "build" ]
then
    rm -r build
fi
mkdir build
cd build
echo -e "cmake mHM..."
cmake -DNETCDF_LIBRARY=/opt/sw/x86_64/glibc-2.17/ivybridge-ep/netcdf_C/4.4.1.1/intel/18/intel-mpi/2018/hdf5/1.8.12-MPI/pnetcdf/1.10.0/lib/libnetcdf.so -DNETCDF_F90_LIBRARY=/opt/sw/x86_64/glibc-2.17/ivybridge-ep/netcdf_Fortran/4.4.4/intel/18/intel-mpi/2018/hdf5/1.8.12-MPI/pnetcdf/1.10.0/netcdf_C/4.4.1.1/lib/libnetcdff.so .. > /dev/null 2>&1
make -j 8 > /dev/null 2>&1
echo -e "done!"
# Compile MPR
cd $PROJECTPATH/MPR
if [ -d "build" ]
then
    rm -r build
fi
mkdir build
cd build
echo -e "cmake MPR..."
cmake -DNETCDF_LIBRARY=/opt/sw/x86_64/glibc-2.17/ivybridge-ep/netcdf_C/4.4.1.1/intel/18/intel-mpi/2018/hdf5/1.8.12-MPI/pnetcdf/1.10.0/lib/libnetcdf.so -DNETCDF_F90_LIBRARY=/opt/sw/x86_64/glibc-2.17/ivybridge-ep/netcdf_Fortran/4.4.4/intel/18/intel-mpi/2018/hdf5/1.8.12-MPI/pnetcdf/1.10.0/netcdf_C/4.4.1.1/lib/libnetcdff.so .. > /dev/null 2>&1
make -j 8 > /dev/null 2>&1
echo -e "done!"

# 4. Create only FSO spatial predictors with MPR
echo -e $GREEN"4. Create FSO spatial predictors with MPR (this takes a while)"$NC

# create folder for spatial predictors
if [ ! -d "$PROJECTPATH/FSO_mHM_major_basins/FSO_spatial_predictors" ]
then
    mkdir $PROJECTPATH/FSO_mHM_major_basins/FSO_spatial_predictors
fi
# copy MPR into folder
cp $PROJECTPATH/MPR/build/MPR-0.6.7 $PROJECTPATH/FSO_mHM_major_basins/FSO_spatial_predictors/MPR-0.6.7
cp $PROJECTPATH/MPR/build/MPR $PROJECTPATH/FSO_mHM_major_basins/FSO_spatial_predictors/MPR

# Create spatial predictors with MPR
cd $PROJECTPATH
for I in 6335304 6335600 6337200 6338100 6340300 6340600 6342800
do
  echo "Compute spatial predictors for catchment $I: FSO_spatial_predictors/sub_${I}_mHM_parameters.nc"
  cp $MPR_PATH/mpr_predictor_FSO_only.nml $MPR_PATH/mpr_predictor_FSO_only_tmp.nml
  sed -i -e "s%out_filename = '<path>/sub_<id>_mHM_parameters.nc'%out_filename = '$PROJECTPATH/FSO_mHM_major_basins/FSO_spatial_predictors/sub_<id>_mHM_parameters.nc'%g" $MPR_PATH/mpr_predictor_FSO_only_tmp.nml
  sed -i -e "s/<id>/$I/g" $MPR_PATH/mpr_predictor_FSO_only_tmp.nml
  sed -i -e "s%<path>%$MHM_PATH%g" $MPR_PATH/mpr_predictor_FSO_only_tmp.nml
  # Run MPR
  ./FSO_mHM_major_basins/FSO_spatial_predictors/MPR -c $MPR_PATH/mpr_predictor_FSO_only_tmp.nml -p $MPR_PATH/mpr_global_parameter_mhm.nml  > /dev/null 2>&1
  rm $MPR_PATH/mpr_predictor_FSO_only_tmp.nml
done

echo -e $GREEN"Done!"$NC
