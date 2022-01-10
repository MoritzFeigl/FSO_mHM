#!/usr/bin/env python
# encoding: utf-8

"""
File Name   : retrieve_time_series_subbasins_zink
Project Name: 2020_FSO_mHM
Description : insert your description here, if applicable
Author      : ottor
Created     : 12.11.20 10:42
"""

# IMPORTS
import pandas as pd
import pathlib as pl


# GLOBAL VARIABLES
ROOT_PATH = '/data/stohyd/mHM_project/germany/sub_basins/'
PATTERN = 'sub_*'
PATH_ADDON = 'zink/default'
FILE = 'daily_discharge.out'
OUTPUT_FILE = '/home/ottor/temp/zink_qobs_qsim_all_basins'


# FUNCTIONS
def read_discharge_out(filename):
    return pd.read_fwf(filename,
                       parse_dates={'date': [0, 1, 2]},
                       usecols=[1, 2, 3, 4, 5],
                       index_col='date')
# CLASSES

# SCRIPT
if __name__ == '__main__':
    folders = sorted(pl.Path(ROOT_PATH).glob(PATTERN))
    dfs = {}
    failed = []
    for i, folder in enumerate(folders):
        basin_id = folder.name.split('_')[1]
        print(i, basin_id)
        path = pl.Path(folder, PATH_ADDON, FILE)
        if not path.exists():
            print(f'{basin_id} does not exist')
            failed.append(basin_id)
            continue
        dfs[basin_id] = read_discharge_out(path)
    df = pd.concat(dfs.values(), axis=1)
    print(failed)
    try:
        df.to_pickle(OUTPUT_FILE + '.pickle')
    except:
        df.to_csv(OUTPUT_FILE + '.csv')

