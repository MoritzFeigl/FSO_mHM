#!/usr/bin/env python
# encoding: utf-8

"""
File Name   : modify_soil_depth
Project Name: 2020_FSO_mHM
Description : insert your description here, if applicable
Author      : ottor
Created     : 08.11.20 22:43
"""

# IMPORTS
import xarray as xr
import pathlib as pl

# GLOBAL VARIABLES
SUBBASINS = ['9305022', '9305031', '9316014', '9316037']
FILES = ['bd', 'clay', 'sand']
# FUNCTIONS

# CLASSES

# SCRIPT
if __name__ == '__main__':
    for basin in SUBBASINS:
        for file in FILES:
            path = pl.Path(f'static/sub_{basin}/mpr/{file}.nc')
            with xr.open_dataset(path) as ds:
                ds.load()
            ds['horizon'].data[-1] = 2.0
            ds['horizon_bnds'].data[-1, 1] = 2.0
            ds.to_netcdf(path, encoding={var: {
                'zlib': True,
                'complevel': 5,
            } for var in ds.data_vars.keys()})

