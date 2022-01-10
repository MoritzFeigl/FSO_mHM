#!/usr/bin/env python
# encoding: utf-8

"""
File Name   : merge_domain
Project Name: 2020_FSO_mHM
Description : insert your description here, if applicable
Author      : ottor
Created     : 26.08.20 11:09
"""

# IMPORTS
import pathlib as pl
import xarray as xr
import numpy as np
from utils.xarray_combine_patched import combine_by_coords
from typing import Optional, Tuple, Union, Iterable, Generator
from textwrap import dedent
import argparse
# GLOBAL VARIABLES

# FUNCTIONS
def parse_args(args: Optional[str] = None):
    """parse the arguments passed from the console and display help with '-h' switch"""
    parser = argparse.ArgumentParser(
        description=dedent('''\
        Merge spatial subdomains for netcdf files

        author: Robert Schweppe
        created: Aug 2020'''))
    parser.add_argument('--target_file', dest='target_file', type=str,
                        help="path of target file to be split merged")
    parser.add_argument('--source_dir', dest='source_dir', default=None, help="folder with source files")
    parser.add_argument('--pattern', dest='pattern', default=None, help=dedent('''\
        pattern of filename, use "*" for variable-length wildcards and "?" for character wildcard
        '''))
    parser.add_argument('--used_mask', dest='used_mask', action='store_true',
                        help="whether to use to slow merge algorithm for custom masks")


    return parser.parse_args(args)


def merge_subdomains(source_dir, filename_pattern, target_file, used_mask=False):
    paths = sorted(pl.Path(source_dir).expanduser().glob(filename_pattern))
    if used_mask:
        # this works for all mask
        all_ds = [xr.open_dataset(path) for path in paths]
        ds = all_ds[0]
        for i, ds_to_merge in enumerate(all_ds[1:]):
            if i == 49:
                print(i)
            ds = xr.merge([ds, ds_to_merge], fill_value=np.nan)
    else:
        try:
            # this works for all non-mask-based subdomains with SubdomainCreater.split(use_mask=False)
            ds = xr.open_mfdataset(paths, combine='by_coords', parallel=True)
        except ValueError:
            # this works for all non-mask-based subdomains with SubdomainCreater.split(use_mask=True)
            all_ds = [xr.open_dataset(path) for path in paths]
            # this is a patched version from xarray with a more sophisticated tiling approach
            ds = combine_by_coords(all_ds, fill_value=np.nan)

    ds.to_netcdf(target_file)

# CLASSES

# SCRIPT
if __name__ == '__main__':
    # example_args = '--target_file ~/temp/out/merged.nc --pattern subdomain_*.nc --source_dir ~/temp/out --used_mask'
    args = parse_args(
        # example_args.split()
    )
    merge_subdomains(args.source_dir, args.pattern, args.target_file, used_mask=args.used_mask)
