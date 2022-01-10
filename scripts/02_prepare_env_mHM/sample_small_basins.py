#!/usr/bin/env python
# encoding: utf-8

"""
File Name   : sample_small_basins
Project Name: 2020_FSO_mHM
Description : insert your description here, if applicable
Author      : ottor
Created     : 10.11.20 13:06
"""

# IMPORTS
import pandas as pd
from random import choice

# GLOBAL VARIABLES
BASIN_PROPERTY_FILE = 'LUT_german_basins_with_properties_and_water_balance.txt'
AREA_CUTOFF = 7500


# FUNCTIONS
def choose(arg):
    if len(arg) == 0:
        return None
    else:
        return choice(arg)


# CLASSES

# SCRIPT
if __name__ == '__main__':
    ds = pd.read_csv(BASIN_PROPERTY_FILE, sep=';')
    updated_ds = ds[ds['catArea'] < AREA_CUTOFF].\
        assign(pet_div_p=ds['PET[mm_a-1]'] / ds['Precip[mm_a-1]'],
               et_div_p=ds['AET[mm_a-1]'] / ds['Precip[mm_a-1]'])
    updated_ds.plot.scatter(y='et_div_p', x='pet_div_p')
    budyko_bins = pd.cut(updated_ds['pet_div_p'], bins=21)
    selection = [choose(updated_ds[budyko_bins == cat]['Stat_ID'].values) for cat in budyko_bins.cat.categories]
    print([_ for _ in selection if _ is not None])
