#!/usr/bin/env python
# encoding: utf-8

"""
File Name   : check_restart_runs
Project Name: 2020_FSO_mHM
Description : insert your description here, if applicable
Author      : ottor
Created     : 02.10.20 08:28
"""

# IMPORTS
import pathlib as pl
import pandas as pd
import argparse
from typing import Optional
from textwrap import dedent

# GLOBAL VARIABLES

# FUNCTIONS
def parse_console_output(path, i):
    d = {'failed': True, 'KGE': None, 'NSE': None, 'duration': None, 'basin': str(i), 'path': path.name}
    start = None
    end = None
    with open(path) as file_object:
        for line in file_object.readlines():
            if line.strip().startswith('mHM: Finished!'):
                d['failed'] = False
            if line.strip().startswith('Start at'):
                start = pd.Timestamp(line.split('at')[1].strip().strip('.'))
            if line.strip().startswith('Precipitation directory'):
                d['basin'] = pl.Path(line.split()[-1]).parts[-1]
            if line.strip().startswith('Finished at'):
                end = pd.Timestamp(line.split('at')[1].strip().strip('.'))
            if line.strip().startswith('KGE'):
                d['KGE'] = float(line.split()[-1])
            if line.strip().startswith('NSE'):
                d['NSE'] = float(line.split()[-1])
    if start is not None and end is not None:
        d['duration'] = (end - start).seconds
    return d

def parse_args(args: Optional[str] = None):
    """parse the arguments passed from the console and display help with '-h' switch"""
    parser = argparse.ArgumentParser(
        description=dedent('''\
        Analyze mhm console output
        scans
        author: Robert Schweppe
        created: Oct 2020'''))
    parser.add_argument('--pattern', dest='pattern', type=str, default='',
                        help="search in current dir for files with '<pattern>*.log', write file to 'all-{pattern}.log'")

    return parser.parse_args(args)


# CLASSES

# SCRIPT
if __name__ == '__main__':
    # pattern: mhm_spinup-7078401
    args = parse_args()
    paths = sorted(pl.Path('').glob(f'{args.pattern}*.log'))
    parsed = pd.DataFrame([parse_console_output(path, i) for i, path in enumerate(paths)]).set_index('basin').sort_index()
    parsed.to_string(f'all_{args.pattern}.info')
    print(sum(parsed['failed'])/len(parsed))
