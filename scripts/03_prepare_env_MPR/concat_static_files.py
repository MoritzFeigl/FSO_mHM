#!/usr/bin/env python
# encoding: utf-8

"""
File Name   : concat_static_files
Project Name: 2020_FSO_mHM
Description : combines all static files into one big file
Author      : ottor
Created     : 15.07.20 15:20
"""

# IMPORTS
import xarray as xr
from pathlib import Path
import numpy as np

# GLOBAL VARIABLES
# output paths (simulation environment)
ROOT_TARGET_FOLDER = '/work/ottor/FSO_mHM_major_basins'


# FUNCTIONS
def get_coord_values(dfs, coord_name, step=100.0):
    """
    get a 1D array containing all coordinate values present in each of dfs

    Parameters
    ----------
    dfs: all datasets to use for global coordinate name creation
    coord_name: coordinate name to get values from
    step: step size of coordinate values

    Returns
    -------
    1D np.array
    """
    min_value = min([float(df[coord_name].min()) for df in dfs])
    max_value = max([float(df[coord_name].max()) for df in dfs])
    count = int((max_value - min_value) / step) + 1
    return np.linspace(min_value, max_value, count)


# CLASSES


# SCRIPT
if __name__ == '__main__':
    # get all subbasins
    paths = sorted(Path(ROOT_TARGET_FOLDER, 'static').glob('sub_*'))
    # get all static files
    files = [path.parts[-2:] for path in sorted(paths[0].rglob('*.nc'))]
    for file in files:
        # skip if target file already exists
        output_path = Path(ROOT_TARGET_FOLDER, 'static_concat', *file)
        if output_path.exists():
            continue
        if not output_path.parent.exists():
            output_path.parent.mkdir(parents=True, exist_ok=True)
        # load all the datasets
        dfs = [xr.open_dataset(Path(path, *file)) for path in paths]
        # special case for latlon.nc
        if 'xc_l0' in dfs[0].dims:
            # determine the min and max boundary of the xc*, yc* coords
            coords = {
                'xc_l0': get_coord_values(dfs, 'xc_l0'),
                'yc_l0': get_coord_values(dfs, 'yc_l0'),
                'xc': get_coord_values(dfs, 'xc'),
                'yc': get_coord_values(dfs, 'yc'),
            }
        else:
            # determine the min and max boundary of the lat, lon coords
            coords = {
                'lat': get_coord_values(dfs, 'lat'),
                'lon': get_coord_values(dfs, 'lon'),
            }
        # add additional coords to the dictionary
        coords.update({k: v for k, v in dfs[0].coords.items() if k not in coords.keys()})

        # get all existing data_vars and init to new size and set to np.nan
        data_vars = {k: (v.dims, np.full([len(coords.get(dim, dfs[0][dim])) for dim in v.dims], np.nan)) for k, v in dfs[0].data_vars.items()}
        # build the complete DataSet
        df = xr.Dataset(data_vars=data_vars, coords=coords)
        # iteratively enter the non-missing values
        for i_df, _df in enumerate(dfs):
            print('working on', file[-1], 'for domain', paths[i_df].name)
            df = df.combine_first(_df)
        # create output folder
        # add the attrs
        for item in (*df.coords, *df.data_vars):
            if dfs[0][item].attrs:
                df[item].attrs = dfs[0][item].attrs
        if dfs[0].attrs:
            df.attrs = dfs[0].attrs

        print('writing to', output_path)
        # write to file
        df.to_netcdf(
            output_path,
            encoding={var: {
                'zlib': True,
                'complevel': 5,
            } for var in df.data_vars.keys()}
        )

