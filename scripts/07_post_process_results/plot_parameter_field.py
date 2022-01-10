#!/usr/bin/env python
# encoding: utf-8

"""
File Name   : plot_parameter_field
Project Name: 2020_FSO_mHM
Description : insert your description here, if applicable
Author      : ottor
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

# GLOBAL VARIABLES
RENAME_DIMS = {
    'lon_out': 'x',
    'lat_out': 'y',
}


# FUNCTIONS
def _all_dims_contained(dims: Iterable[str], obj: Union[xr.DataArray, xr.Dataset]):
    """check if all dimensions of a list are containing in an xarray object"""
    return all([dim in obj.dims for dim in dims])


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


def print_map(filename: str, ds: Union[xr.DataArray, xr.Dataset], x_coord: str = 'lon', y_coord: str = 'lat'):
    """print using the cartopy library, provide a 2D xarray object as ds"""
    # set the map projection
    ax = plt.axes(projection=ccrs.PlateCarree())
    # add coastlines to map, add other items if you like
    ax.coastlines()

    ds.plot.pcolormesh(ax=ax, x=x_coord, y=y_coord,
                       # does args are not required if the transform is done elsewhere...
                       # transform=ccrs.PlateCarree(),
                       # subplot_kws=dict(projection=ccrs.Orthographic(), facecolor="gray"),
                       )
    # play around with the spatial extent
    #ax.set_global()
    #ax.set_xlim([ds[x_coord].min() - 20, ds[x_coord].max() + 20])
    #ax.set_ylim([ds[y_coord].min() - 20, ds[y_coord].max() + 20])

    plt.savefig(filename)
    plt.close()


def export_to_geotiff(filename: str, ds: Union[xr.DataArray, xr.Dataset]):
    ds.rio.to_raster(filename)

# CLASSES

# SCRIPT
if __name__ == '__main__':
    # load the dataset, optionally select variables or slices of dimensions
    my_filename = '~/Downloads/FSO_SCE_NSE_run_2/sub_6340600_parameter.nc'
    var_select = 'L1_FieldCap'
    my_ds = read(my_filename, select_variable=var_select)
    my_selected = select(my_ds, isel=dict(horizon_out=0, land_cover_period_out=0))
    my_transformed = transform_coord(my_selected, 31468, 4326)
    print_simple(filename='simple.png', ds=my_transformed)
    print_map(filename='map.png', ds=my_transformed)
    export_to_geotiff(filename=f'{var_select}.tiff', ds=my_transformed)
