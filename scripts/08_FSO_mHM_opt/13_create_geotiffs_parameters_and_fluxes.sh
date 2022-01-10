#!/bin/sh
# Create geoTIFFs from netcdf

PROJECTPATH=/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt

cd $PROJECTPATH

module purge
module load anaconda3/5.3.0
conda init bash
source ~/.bashrc
module load anaconda3/5.3.0
ENV_NAME="para_map"
conda create -n ${ENV_NAME} python=3.7
conda activate ${ENV_NAME}

conda install typing --yes
conda install cartopy --yes
conda install matplotlib --yes
conda install xarray --yes
conda install rasterio --yes
conda install -c conda-forge rioxarray --yes
 conda install -c anaconda netcdf4
#conda install glob --yes

python3.7 -m 13_create_geotiffs_parameters_and_fluxes.py
