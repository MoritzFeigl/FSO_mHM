#!/usr/bin/env python3
# encoding: utf-8

"""
File Name   : plot_parameter_field
Project Name: 2020_FSO_mHM
Description : generate geotiff files for mHM parameter and fluxes
Author      : ottor, Moritz Feigl
Created     : 08.02.21 10:35
"""

# IMPORTS
from typing import Union, Iterable, Optional
# all external modules installed with conda
import cartopy.crs as ccrs  # used version 0.17.0
import matplotlib.pyplot as plt  # 3.2.1
import numpy as np  # 1.18.1
import xarray as xr  # 0.15.0
from rasterio.warp import transform  # 1.1.3
import rioxarray
import os
import re
from itertools import repeat
#os.environ['PROJ_LIB'] = 'C:\\Users\\morit\\anaconda3\\envs\\pytorch\\Library\\share\\proj'

# GLOBAL VARIABLES
RENAME_DIMS = {
    'lon_out': 'x',
    'lat_out': 'y',
}


# FUNCTIONS

def _all_dims_contained(dims: Iterable[str], obj: Union[xr.DataArray, xr.Dataset]):
    """check if all dimensions of a list are containing in an xarray object"""
    all([dim in obj.dims for dim in dims])


def read(filename, select_variable: Optional[Union[str, Iterable[str]]] = None, select_slice: Optional[dict] = None):
    """open a netcdf file from disk, optionally selecting variables or slices of dimensions"""
    with xr.open_dataset(filename) as ds:
        # select one or more variables
        if select_variable is not None:
            ds = ds[select_variable]
        # select a part of the coordinates
        if select_slice is not None and _all_dims_contained(select_slice.keys(), ds):
            ds = ds.sel(**select_slice)
        # rename the coordinates contained in RENAME_DIMS
        if _all_dims_contained(RENAME_DIMS.keys(), ds):
            ds = ds.rename(RENAME_DIMS)
        # load only the part needed and close the file handle
        ds.load()
    return ds


def select(ds: Union[xr.DataArray, xr.Dataset], isel: Optional[dict] = None, sel: Optional[dict] = None):
    """again, perform a selection for coordinates based on labels (sel) or indices (isel)"""
    selected = ds
    if isel:
        selected = selected.isel(**isel)
    if sel:
        selected = selected.sel(**sel)
    return selected


def transform_coord(ds: Union[xr.DataArray, xr.Dataset], source_epsg: int, target_epsg: int):
    """transform the 2D coordinates using its EPSG codes, assumes 1D input coordinates named x,y and creates lon, lat"""
    # Compute the lon/lat coordinates with rasterio.warp.transform
    ny, nx = len(ds['y']), len(ds['x'])
    x, y = np.meshgrid(ds['x'], ds['y'])
    # Rasterio works with 1D arrays
    lon, lat = transform({'init': f'EPSG:{source_epsg}'}, {'init': f'EPSG:{target_epsg}'},
                         x.flatten(), y.flatten())
    lon = np.asarray(lon).reshape((ny, nx))
    lat = np.asarray(lat).reshape((ny, nx))
    ds.coords['lon'] = (('y', 'x'), lon)
    ds.coords['lat'] = (('y', 'x'), lat)
    ds['lat'].attrs = {'standard_name': 'latitude', 'long_name': 'Latitude',
                       'units': 'degrees_north', '_CoordinateAxisType': 'Lat', 'axis': 'y'}
    ds['lon'].attrs = {'standard_name': 'longitude', 'long_name': 'Longitude',
                       'units': 'degrees_east', '_CoordinateAxisType': 'Lon', 'axis': 'x'}
    return ds


def print_simple(filename: str, ds: Union[xr.DataArray, xr.Dataset], x_coord: str = 'lon', y_coord: str = 'lat'):
    """print using the native xarray plotting, provide a 2D xarray object as ds"""
    ds.plot(x=x_coord, y=y_coord)
    plt.savefig(filename)
    plt.close()


def export_to_geotiff(filename: str, ds: Union[xr.DataArray, xr.Dataset]):
    ds.rio.to_raster(filename)


def fso_mhm_parameter_to_geotiffs(run: str):
    # for each basin: load and add to overall xarray
    files = [run + '/' + dir for dir in os.listdir(run)]
    variables = ['KSat_notill', 'FieldCap_notill', 'SoilMoistureSaturationDeficit', 'VarKSat_horizontal_relative',
                 'VarKSat_vertical_relative', 'L1_FieldCap', 'L1_Alpha', 'L1_Kperco', 'L1_SatSoilMoisture']
    horizon_notill = ['KSat_notill', 'FieldCap_notill']
    horizon_all_latlon = ['SoilMoistureSaturationDeficit', 'VarKSat_horizontal_relative', 'VarKSat_vertical_relative']
    horizon_all_xy = ['L1_Kperco', 'L1_Alpha']
    horizon_out = ['L1_FieldCap', 'L1_SatSoilMoisture']
    for var_select in variables:
        for id, file in enumerate(files):
            print(f'{var_select} : {file}: ')
            # load the dataset, optionally select variables or slices of dimensions
            my_ds = read(file, select_variable=var_select)
            if var_select in horizon_notill:
                my_selected = select(my_ds, isel=dict(horizon_notill=0, land_cover_period=1))
                my_selected = my_selected.rename({'lat': 'y', 'lon': 'x'})
            elif var_select in horizon_all_latlon:
                my_selected = select(my_ds, isel=dict(horizon_all=0, land_cover_period=1))
                my_selected = my_selected.rename({'lat': 'y', 'lon': 'x'})
            elif var_select in horizon_all_xy:
                my_selected = select(my_ds, isel=dict(horizon_all=0, land_cover_period_out=1))
                my_selected = my_selected.rename({'lat_out': 'y', 'lon_out': 'x'})
            elif var_select in horizon_out:
                my_selected = select(my_ds, isel=dict(horizon_out=2, land_cover_period_out=1))
                my_selected = my_selected.rename({'lat_out': 'y', 'lon_out': 'x'})
            # transform and save as geotiff
            my_transformed = transform_coord(my_selected, 31468, 4326)
            if id == 0:
                full_ds = my_transformed
            else:
                full_ds = full_ds.combine_first(my_transformed)
        export_to_geotiff(filename=f'{run}_{var_select}.tiff', ds=full_ds)



def fso_mhm_time_parameter_to_geotiffs(run: str):
    # for each basin: load and add to overall xarray
    files = [run + '/' + dir for dir in os.listdir(run)]
    variables = ['QD', 'QIf', 'QIs', 'QB']
    for var_select in variables:
        for id, file in enumerate(files):
            print(f'{var_select} : {file}: ')
            # get total runoff 'Q' for scaling
            Q = read(file, select_variable='Q')
            # load the dataset, optionally select variables or slices of dimensions
            my_ds = read(file, select_variable=var_select)
            scaled_ds = my_ds / Q
            final_ds = scaled_ds.mean(dim="time")
            final_ds = final_ds.rename({'northing': 'y', 'easting': 'x'})
            # transform and save as geotiff
            my_transformed = transform_coord(final_ds, 31468, 4326)
            if id == 0:
                full_ds = my_transformed
            else:
                full_ds = full_ds.combine_first(my_transformed)
        export_to_geotiff(filename=f'{run}_{var_select}.tiff', ds=full_ds)


def fso_mhm_sw_fraction_to_geotiffs(run: str):
    # for each basin: load and add to overall xarray
    files = [run + '/' + dir for dir in os.listdir(run)]
    nr_files = len(files)
    restart_files = ['../restarts/' + re.sub("fulx_states", "restart", file) for file in files]
    for depth in [0, 1, 2]:
        for id in range(nr_files):
            print(f'depth {depth} : {files[id]}: ')
            # get SWC_Lxx from flux_states.nc
            SWC_Lxx = read(files[id], select_variable='SWC_L0' + str(depth+1))
            # get L1_soilMoistFC
            L1_soilMoistFC = read(restart_files[id], select_variable='L1_soilMoistFC')
            L1_soilMoistFC = select(L1_soilMoistFC, isel=dict(L1_SoilHorizons=depth, L1_LandCoverPeriods=1))
            for t in range(0, 1827):
                diff = SWC_Lxx.values[t, :, :] - L1_soilMoistFC.values
                SWC_Lxx.values[t, :, :] = (diff < 0)*1
                SWC_Lxx.values[t, np.isnan(diff)] = np.nan
            SWC_Lxx_agg = SWC_Lxx.sum(dim="time")
            final_ds = SWC_Lxx_agg / SWC_Lxx.values.shape[0]
            final_ds.values[np.isnan(SWC_Lxx.values[t, :, :])] = np.nan
            final_ds = final_ds.rename({'northing': 'y', 'easting': 'x'})
            # transform and save as geotiff
            my_transformed = transform_coord(final_ds, 31468, 4326)
            if id == 0:
                full_ds = my_transformed
            else:
                full_ds = full_ds.combine_first(my_transformed)
        export_to_geotiff(filename=f'{run}_SWC_fraction_0{depth+1}.tiff', ds=full_ds)




# CLASSES

# SCRIPT
if __name__ == '__main__':
    os.chdir("/home/lv71468/mfeigl/FSO_mHM/Results/parameter_maps")
    runs = next(os.walk('.'))[1]
    # variables to be exported as geotiff
    list(map(fso_mhm_parameter_to_geotiffs, runs))
    # fluxes
    os.chdir("/home/lv71468/mfeigl/FSO_mHM/Results/flux_states")
    runs = next(os.walk('.'))[1]
    list(map(fso_mhm_time_parameter_to_geotiffs, runs))
    list(map(fso_mhm_sw_fraction_to_geotiffs, runs))
