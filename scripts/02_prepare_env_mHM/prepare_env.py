#!/usr/bin/env python
# encoding: utf-8

"""
File Name   : prepare_env
Project Name: 2020_FSO_mHM
Description : insert your description here, if applicable
Author      : ottor
Created     : 26.06.20 08:51
"""

# IMPORTS
from pathlib import Path
import xarray as xr
import shutil
from datetime import date
import f90nml
import re
import subprocess
from os import sep
import numpy as np
from pyproj import Transformer

# GLOBAL VARIABLES
# parts of the script to execute - available are ['forcing', 'static', 'config']
PROCESSES = ['forcing', 'static', 'config']
# whether to overwrite existing files
OVERWRITE = True
# a list of the ids for each basin
# basin_ids="Mulde Ems Neckar Saale Main Weser Danube"
BASIN_SELECTION = ['6340600', '6338100', '6335600', '6340300', '6335304', '6337200', '6342800']
# all the remaining others
BASIN_SELECTION = ['0579085', '6321100', '6335115', '6335125', '6335290', '6335300',
                   #'6335304', #
                   '6335350',
                   #'6335600', #
                   '6335601', '6335602', '6335800',
                   #'6337200', #
                   '6337500', '6337501', '6337502', '6337503', '6337504',
                   '6337505', '6337506', '6337507', '6337508', '6337509', '6337510', '6337511', '6337512', '6337513',
                   #'6338100', #
                   '6338110', '6338120', '6338130', '6338800', '6340050', '6340200',
                   #'6340300', #
                   '6340302', '6340501',
                   #'6340600', #
                   '6340610', '6340700', '6340800', '6342200', '6342520',
                   #'6342800', #
                   '9304001',
                   '9304003', '9304014', '9304015', '9304016', '9304027', '9304028', '9304029', '9304030', '9304031',
                   '9304032', '9304033', '9304034', '9304037', '9304038', '9304039', '9304042', '9304045', '9304047',
                   '9304050', '9304051', '9304052', '9304054', '9304056', '9304057', '9304059', '9304060', '9304061',
                   '9304062', '9304063', '9304064', '9304065', '9304066', '9304067', '9304069', '9304070', '9304071',
                   '9304072', '9304074', '9304076', '9304078', '9304079', '9304084', '9304085', '9304086', '9304087',
                   '9304088', '9304089', '9304090', '9304091', '9304098', '9304102', '9304105', '9304114', '9304119',
                   '9304122', '9304130', '9304131', '9304134', '9304135', '9304136', '9304137', '9304145', '9305004',
                   '9305006', '9305009', '9305022', '9305024', '9305031', '9305036', '9305038', '9305043', '9305044',
                   '9305045', '9305046', '9305047', '9305051', '9305052', '9305055', '9305056', '9305059', '9305062',
                   '9305063', '9305064', '9305065', '9305066', '9305068', '9305073', '9305075', '9305077', '9305078',
                   '9305081', '9305086', '9305087', '9305088', '9305090', '9305091', '9305092', '9305094', '9305095',
                   '9305096', '9305105', '9305106', '9305107', '9305111', '9305114', '9305116', '9305118', '9305128',
                   '9305129', '9305130', '9305131', '9305132', '9305133', '9305134', '9305135', '9305136', '9305137',
                   '9305138', '9305139', '9305140', '9305142', '9305144', '9305145', '9305151', '9315005', '9315006',
                   '9315010', '9315014', '9315015', '9316002', '9316003', '9316006', '9316008', '9316009', '9316014',
                   '9316016', '9316018', '9316023', '9316024', '9316025', '9316036', '9316037', '9316051', '9316052',
                   '9316074', '9316086', '9316089', '9316091', '9316094', '9316095', '9316096', '9316097', '9316099',
                   '9316108', '9316109', '9316111', '9316113', '9316129', '9316131', '9316132', '9316135', '9316136',
                   '9316138', '9316139', '9316142', '9316146', '9316150', '9316152', '9316154', '9316155', '9316158',
                   '9316160', '9316162', '9316163', '9316166', '9316170', '9316172', '9316175', '9316176', '9316177',
                   '9316180', '9316181', '9316184', '9316185', '9316186', '9316188', '9316189', '9316190', '9316191',
                   '9316194', '9316195', '9316197', '9316199', '9316200', '9316204', '9316206', '9316208', '9316210',
                   '9316211', '9316212', '9316213', '9316214', '9316215', '9316216', '9316218', '9316220', '9316221',
                   '9316222', '9316224', '9316227', '9316228', '9316229', '9316231', '9316235', '9316236', '9316238',
                   '9316239', '9316240', '9316241', '9316242', '9316243', '9316246', '9316249', '9316252', '9316253',
                   '9316260', '9316261', '9316264', '9316265', '9316266', '9316268', '9316270', '9316271', '9316273',
                   '9316277', '9316283', '9316284', '9316288', '9316289', '9316290', '9316291', '9316293', '9316297',
                   '9316298', '9316302', '9316303', '9316304', '9316310', '9316311', '9316315', '9316316', '9316318',
                   '9316319', '9316320', '9316321', '9316325', '9316327', '9316330', '9316337', '9316342', '9316343',
                   '9316344', '9316345', '9316348', '9316350', '9316352', '9316354', '9316355', '9316356', '9316359',
                   '9316361', '9316367', '9316378', '9316380', '9316383', '9316389', '9316392', '9316393', '9316394',
                   '9316405', '9316407', '9316408', '9316410', '9316412', '9316414', '9316415', '9316418', '9316419',
                   '9316420', '9316421', '9316426', '9316430', '9316433', '9316434', '9316436', '9316437', '9316438',
                   '9316439', '9316440', '9316441', '9316444', '9316445', '9316447', '9316452', '9316453', '9316457',
                   '9316458', '9316460', '9316463', '9316464', '9316466', '9316470', '9316471', '9316473', '9316474',
                   '9316475', '9316476', '9316479', '9316480', '9316483', '9316485', '9316487', '9316489', '9316494',
                   '9316495', '9316503', '9316505', '9316509', '9316515', '9316520', '9316522', '9316523', '9316524',
                   '9316526', '9316529', '9316534', '9324001', '9324002', ]
USE_SYMLINKS = False
# output paths (simulation environment)
ROOT_TARGET_FOLDER = '/work/ottor/FSO_mHM_major_basins'
# ROOT_TARGET_FOLDER = './output'
STATIC_TARGET_FOLDER = 'static'
FORCING_TARGET_FOLDER = 'forcing'
CONFIG_TARGET_FOLDER = 'config'
OUTPUT_TARGET_FOLDER = 'output'
TARGET_MPR_NML_NAME = 'mpr.nml'
# input paths (environment from Zink et al., 2017)
ROOT_SOURCE_FOLDER = '/data/stohyd/data/processed/Germany/basins/'
FORCING_SOURCE_FOLDER = 'meteo'
FORCING_SOURCE_SEL = ['pre', 'tavg', 'pet']
# local mHM repository
ROOT_REPO_FOLDER = '/home/ottor/lib/mhm_mpr/'
# ROOT_REPO_FOLDER = '/Users/ottor/nc/Home/local_libs/fortran/mhm_mpr/'
CONFIG_SOURCE_SEL = [
    'mhm.nml', 'mhm_outputs.nml', 'mhm_parameter.nml', 'mrm_outputs.nml',
    # 'mrm.nml',  'mrm_parameter.nml',
    'mpr_fso.nml'
]
EXECUTABLE_PATH = 'build_foss2019b_debug/mhm'
# EXECUTABLE_PATH = 'build_gnu/mhm'

# all predictor variables, whose paths need to be adapted for MPR
MPR_READ_VARS = ['bd', 'sand', 'clay', 'land_cover', 'slope', 'karstic', 'classunit', 'lai_class', 'aspect', 'dem']
# all default mHM settings for namelist
# commented lines denote keywords for upcoming mHM version
MHM_NML_REPLACE_DICT = {
    ('project_description', 'project_details'): 'FSO mHM coupling',
    ('project_description',
     'setup_description'): 'forward run setup for major German basins based on Zink et al. (2017)',
    ('project_description', 'simulation_type'): 'historical simulation',
    ('project_description', 'contact'): 'robert.schweppe@ufz.de',
    ('project_description', 'history'): 'first test run',
    ('mainconfig', 'iflag_cordinate_sys'): 0,
    #('mainconfig', 'ndomains'): 1,
    ('mainconfig', 'nbasins'): 1,
    ('mainconfig', 'resolution_hydrology'): [4000],
    #('mainconfig', 'l0domain'): [1],
    ('mainconfig', 'l0basin'): [1],
    ('mainconfig', 'write_restart'): True,
    #('mainconfig', 'read_opt_domain_data'): [0],
    ('mainconfig_mhm_mrm', 'dir_restartin'): [''],
    ('mainconfig_mhm_mrm', 'resolution_routing'): [4000],
    ('mainconfig_mhm_mrm', 'timestep'): 1,
    ('mainconfig_mhm_mrm', 'read_restart'): False,
    ('mainconfig_mhm_mrm', 'optimize'): False,
    ('mainconfig_mhm_mrm', 'optimize_restart'): False,
    ('mainconfig_mhm_mrm', 'opti_method'): 1,
    ('mainconfig_mhm_mrm', 'opti_function'): 10,
    ('mainconfig_mrm', 'alma_convention'): True,
    ('mainconfig_mrm', 'varnametotalrunoff'): 'total_runoff',
    ('mainconfig_mrm', 'filenametotalrunoff'): 'total_runoff',
    ('mainconfig_mrm', 'gw_coupling'): False,
    ('directories_general', 'dirconfigout'): OUTPUT_TARGET_FOLDER + sep,
    ('directories_general', 'dir_restartout'): [OUTPUT_TARGET_FOLDER + sep],
    ('directories_general', 'dir_out'): [OUTPUT_TARGET_FOLDER + sep],
    ('directories_mhm', 'path_mpr_nml'): [TARGET_MPR_NML_NAME],
    ('directories_mhm', 'dir_precipitation'): [str(Path(ROOT_TARGET_FOLDER, FORCING_TARGET_FOLDER)) + sep],
    ('directories_mhm', 'dir_temperature'): [str(Path(ROOT_TARGET_FOLDER, FORCING_TARGET_FOLDER)) + sep],
    ('directories_mhm', 'dir_referenceet'): [str(Path(ROOT_TARGET_FOLDER, FORCING_TARGET_FOLDER)) + sep],
    ('directories_mhm', 'dir_mintemperature'): [str(Path(ROOT_TARGET_FOLDER, FORCING_TARGET_FOLDER)) + sep],
    ('directories_mhm', 'dir_maxtemperature'): [str(Path(ROOT_TARGET_FOLDER, FORCING_TARGET_FOLDER)) + sep],
    ('directories_mhm', 'dir_netradiation'): [str(Path(ROOT_TARGET_FOLDER, FORCING_TARGET_FOLDER)) + sep],
    ('directories_mhm', 'dir_absvappressure'): [str(Path(ROOT_TARGET_FOLDER, FORCING_TARGET_FOLDER)) + sep],
    ('directories_mhm', 'dir_windspeed'): [str(Path(ROOT_TARGET_FOLDER, FORCING_TARGET_FOLDER)) + sep],
    ('directories_mhm', 'time_step_model_inputs'): [0],
    ('directories_mhm', 'timestep_lai_input'): [0],
    ('directories_mrm', 'dir_gauges'): [str(Path(ROOT_TARGET_FOLDER, STATIC_TARGET_FOLDER, 'routing')) + sep],
    ('directories_mrm', 'dir_total_runoff'): [str(Path(ROOT_TARGET_FOLDER, STATIC_TARGET_FOLDER, 'routing')) + sep],
    ('directories_mrm', 'dir_bankfull_runoff'): [str(Path(ROOT_TARGET_FOLDER, STATIC_TARGET_FOLDER, 'routing')) + sep],
    ('optional_data', 'dir_soil_moisture'): [str(Path(ROOT_TARGET_FOLDER, STATIC_TARGET_FOLDER, 'optional_data')) + sep],
    ('optional_data', 'nsoilhorizons_sm_input'): 1,
    ('optional_data', 'timestep_sm_input'): -2,
    ('optional_data', 'dir_neutrons'): [str(Path(ROOT_TARGET_FOLDER, STATIC_TARGET_FOLDER, 'optional_data')) + sep],
    ('optional_data', 'dir_evapotranspiration'): [str(Path(ROOT_TARGET_FOLDER, STATIC_TARGET_FOLDER, 'optional_data')) + sep],
    ('optional_data', 'timestep_et_input'): -2,
    #('optional_data', 'dir_tws'): [str(Path(ROOT_TARGET_FOLDER, STATIC_TARGET_FOLDER, 'optional_data')) + sep],
    ('optional_data', 'file_tws'): [str(Path(ROOT_TARGET_FOLDER, STATIC_TARGET_FOLDER, 'optional_data', 'tws_basin_1.txt'))],
    #('optional_data', 'timestep_tws_input'): -2,
    #('processselection', 'processcase'): [1, 1, 1, 1, 0, 1, 1, 3, 1, 0],
    ('processselection', 'processcase'): [1, 1, 1, 1, 0, 1, 1, 2, 1, 1],
    ('time_periods', 'warming_days'): [0],
    ('time_periods', 'eval_per'): [
        {
            'ystart': 1990,
            'mstart': 1,
            'dstart': 1,
            'yend': 1991,
            'mend': 12,
            'dend': 31,
        }
    ],
    ('lcover', 'nlandcoverperiods'): 3,
    ('evaluation_gauges', 'ngaugestotal'): 1,
    #('evaluation_gauges', 'nogauges_domain'): [1],
    ('evaluation_gauges', 'nogauges_basin'): [1],
    ('evaluation_gauges', 'gauge_id'): [[1234]],
    ('evaluation_gauges', 'gauge_filename'): [['1234.txt']],
    ('inflow_gauges', 'ninflowgaugestotal'): 0,
    #('inflow_gauges', 'noinflowgauges_domain'): [0],
    ('inflow_gauges', 'noinflowgauges_basin'): [0],
    ('inflow_gauges', 'inflowgauge_id'): [[-9]],
    ('inflow_gauges', 'inflowgauge_filename'): [['']],
    ('inflow_gauges', 'inflowgauge_headwater'): [[False]],
    ('panevapo', 'evap_coeff'): [1.3, 1.2, 0.72, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.5],
    ('nightdayratio', 'read_meteo_weights'): False,
    ('nightdayratio', 'fnight_prec'): [0.46, 0.5, 0.52, 0.51, 0.48, 0.5, 0.49, 0.48, 0.52, 0.56, 0.5, 0.47],
    ('nightdayratio', 'fnight_pet'): [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
    ('nightdayratio', 'fnight_temp'): [-0.76, -1.3, -1.88, -2.38, -2.72, -2.75, -2.74, -3.04, -2.44, -1.6, -0.94,
                                       -0.53],
    ('optimization', 'niterations'): 7,
    ('optimization', 'seed'): 1235876,
    ('optimization', 'dds_r'): 0.2,
    ('optimization', 'sa_temp'): -9.0,
    ('optimization', 'sce_ngs'): 2,
    ('optimization', 'sce_npg'): -9,
    ('optimization', 'sce_nps'): -9,
    ('optimization', 'mcmc_opti'): False,
    ('optimization', 'mcmc_error_params'): [0.01, 0.6],
}
COMPRESSION = {'zlib': True, 'complevel': 5}


# FUNCTIONS
def _read_namelist(path):
    parser = f90nml.Parser()
    parser.global_start_index = 1
    return parser.read(path)


def _configure_mhm_nml(in_path, out_path, output_path=None, forcing_path=None, routing_path=None, gauge_id=None,
                       start_date=None, end_date=None, spinup=None):
    local_mhm_nml_replace_dict = {}
    if output_path is not None:
        local_mhm_nml_replace_dict[('directories_general', 'dirconfigout')] = str(output_path) + sep
        local_mhm_nml_replace_dict[('directories_general', 'dir_restartout')] = [str(output_path) + sep]
        local_mhm_nml_replace_dict[('directories_general', 'dir_out')] = [str(output_path) + sep]
    if forcing_path is not None:
        local_mhm_nml_replace_dict[('directories_mhm', 'dir_precipitation')] = [str(forcing_path) + sep]
        local_mhm_nml_replace_dict[('directories_mhm', 'dir_temperature')] = [str(forcing_path) + sep]
        local_mhm_nml_replace_dict[('directories_mhm', 'dir_referenceet')] = [str(forcing_path) + sep]
        local_mhm_nml_replace_dict[('directories_mhm', 'dir_mintemperature')] = [str(forcing_path) + sep]
        local_mhm_nml_replace_dict[('directories_mhm', 'dir_maxtemperature')] = [str(forcing_path) + sep]
        local_mhm_nml_replace_dict[('directories_mhm', 'dir_netradiation')] = [str(forcing_path) + sep]
        local_mhm_nml_replace_dict[('directories_mhm', 'dir_absvappressure')] = [str(forcing_path) + sep]
        local_mhm_nml_replace_dict[('directories_mhm', 'dir_windspeed')] = [str(forcing_path) + sep]
    if routing_path is not None:
        local_mhm_nml_replace_dict[('directories_mrm', 'dir_gauges')] = [str(routing_path) + sep]
        local_mhm_nml_replace_dict[('directories_mrm', 'dir_total_runoff')] = [str(routing_path) + sep]
        local_mhm_nml_replace_dict[('directories_mrm', 'dir_bankfull_runoff')] = [str(routing_path) + sep]
    if gauge_id is not None:
        local_mhm_nml_replace_dict[('evaluation_gauges', 'gauge_id')] = [[gauge_id]]
        local_mhm_nml_replace_dict[('evaluation_gauges', 'gauge_filename')] = [['{}.txt'.format(gauge_id)]]
    if start_date is not None and end_date is not None:
        local_mhm_nml_replace_dict[('time_periods', 'warming_days')] = [spinup]
        local_mhm_nml_replace_dict[('time_periods', 'eval_per')] = [
            {
                'ystart': start_date.year,
                'mstart': start_date.month,
                'dstart': start_date.day,
                'yend': end_date.year,
                'mend': end_date.month,
                'dend': end_date.day,
            }
        ]

    nml = _read_namelist(in_path)
    # set global settings
    for key, value in MHM_NML_REPLACE_DICT.items():
        nml[key[0]][key[1]] = value
    # set "local" settings
    for key, value in local_mhm_nml_replace_dict.items():
        nml[key[0]][key[1]] = value
    print('Creating {} file: {}'.format('config', out_path))
    nml.write(out_path, force=True)


def _get_item_index(nml, item_type, name, item):
    names = nml[item_type][name]
    # it fails if item not in names
    return names.index(item)


def _configure_mpr_nml(in_path, out_path, coords=None, data_arrays=None, out_filename=''):
    if data_arrays is None:
        data_arrays = {}
    if coords is None:
        coords = {}
    nml = _read_namelist(in_path)
    for item, item_type, item_label in zip(
            [coords, data_arrays],
            ['coordinates', 'data_arrays'],
            ['coord_name', 'name']):
        for item_name, item_props in item.items():
            item_index = _get_item_index(nml, item_type=item_type, name=item_label, item=item_name)
            for prop, value in item_props.items():
                nml[item_type][prop][item_index] = value
    if out_filename:
        nml['main']['out_filename'] = out_filename
    print('Creating {} file: {}'.format('config', out_path))
    nml.write(out_path, force=True)


def read_header_file(file, cellsize, nodata=-9999):
    with open(file) as f_object:
        header = dict(line.split() for line in f_object.readlines())
    cols = int(int(header['ncols']) / (cellsize / 100))
    rows = int(int(header['nrows']) / (cellsize / 100))
    cellsize = float(cellsize)
    xllcorner = float(header['xllcorner'])
    yllcorner = float(header['yllcorner'])
    return cellsize, cols, nodata, rows, xllcorner, yllcorner


def prepare_forcings(basins):
    for basin in basins:
        is_projected = False
        forcing_source_root = Path(ROOT_SOURCE_FOLDER, basin, FORCING_SOURCE_FOLDER)
        forcing_target_root = Path(ROOT_TARGET_FOLDER, FORCING_TARGET_FOLDER, basin)

        if not forcing_target_root.exists():
            forcing_target_root.mkdir(parents=True, exist_ok=True)

        for var in FORCING_SOURCE_SEL:
            file_name = '{}.nc'.format(var)
            forcing_path = Path(forcing_source_root, file_name)
            if not forcing_path.exists():
                forcing_path = Path(forcing_source_root, var, file_name)
            target_file = Path(forcing_target_root, file_name)
            if target_file.exists() and not OVERWRITE:
                continue
            # sort y dimension
            print('Creating {} file: {}'.format('forcing', target_file))
            ds = xr.open_dataset(forcing_path)
            if all(coord in ds.coords for coord in ['xc', 'yc']):
                ds = ds.rename({'xc': 'easting', 'yc': 'northing'})
            is_projected = all(coord in ds.coords for coord in ['easting', 'northing'])
            if is_projected:
                x_coord = 'easting'
                y_coord = 'northing'
            else:
                x_coord = 'x'
                y_coord = 'y'
            ds = ds.sortby(y_coord, ascending=False)
            if is_projected:
                cellsize, cols, nodata, rows, xllcorner, yllcorner = get_header_information(ds, x_coord, y_coord, var=var)
            ds.to_netcdf(target_file, encoding={var: dict(_FillValue=-9999.0, **COMPRESSION) for var in ds.data_vars.keys()})
        # link the header.txt file containing grid information on the meteo files
        header_file = Path(forcing_source_root, 'header.txt')
        target_file = Path(forcing_target_root, 'header.txt')
        fallback_file = Path(ROOT_SOURCE_FOLDER, basin, 'morph', 'latlon.nc')
        fallback_file2 = Path(ROOT_SOURCE_FOLDER, basin, 'latlon', 'header_100m.txt')
        if not OVERWRITE and target_file.exists():
            continue
        if header_file.exists():
            _copy_or_link_file(header_file,
                               Path(forcing_target_root, 'header.txt'))
        elif is_projected:
            write_header_file(cellsize, cols, nodata, rows, target_file, xllcorner, yllcorner)
        elif fallback_file.exists():
            print(f'#### using fallback option 1 for {basin}')
            ds = xr.open_dataset(fallback_file)
            cellsize, cols, nodata, rows, xllcorner, yllcorner = get_header_information(ds, 'xc', 'yc')
            write_header_file(cellsize, cols, nodata, rows, target_file, xllcorner, yllcorner)
        elif fallback_file2.exists():
            # for some basins, not even the latlon file exists
            print(f'#### using fallback option 2 for {basin}')
            cellsize, cols, nodata, rows, xllcorner, yllcorner = read_header_file(fallback_file2, cellsize=4000)
            write_header_file(cellsize, cols, nodata, rows, target_file, xllcorner, yllcorner)
        else:
            print(f'#### no header file exists nor information to write properly: {basin}')


def write_header_file(cellsize, cols, nodata, rows, target_file, xllcorner, yllcorner):
    with open(target_file, 'w') as f_obj:
        f_obj.write('{:<14}{}\n'.format('cols', cols))
        f_obj.write('{:<14}{}\n'.format('rows', rows))
        f_obj.write('{:<14}{}\n'.format('xllcorner', xllcorner))
        f_obj.write('{:<14}{}\n'.format('yllcorner', yllcorner))
        f_obj.write('{:<14}{}\n'.format('cellsize', cellsize))
        f_obj.write('{:<14}{}\n'.format('NODATA_value', nodata))


def get_header_information(ds, x_coord, y_coord, var=None, nodata=-9999.0):
    cols = len(ds[x_coord])
    rows = len(ds[y_coord])
    cellsize = abs(float(ds[x_coord][1] - ds[x_coord][0]))
    xllcorner = float(ds[x_coord][0] - cellsize * 0.5)
    yllcorner = float(ds[y_coord][0] - cellsize * 0.5)
    if var is not None:
        nodata = ds[var].encoding['_FillValue']
    return cellsize, cols, nodata, rows, xllcorner, yllcorner


def _copy_or_link_file(source, target):
    if USE_SYMLINKS:
        if target.is_symlink():
            target.unlink()
        target.symlink_to(source)
    else:
        if target.exists():
            target.unlink()
        shutil.copy(source, target)


def prepare_config(basins):
    exec_path = Path(ROOT_REPO_FOLDER, EXECUTABLE_PATH)
    for basin in basins:
        config_target_root = Path(ROOT_TARGET_FOLDER, CONFIG_TARGET_FOLDER, basin)
        output_target_root = Path(ROOT_TARGET_FOLDER, OUTPUT_TARGET_FOLDER, basin)
        for path in [config_target_root, output_target_root]:
            make_dir(path)
        _copy_or_link_file(exec_path, Path(config_target_root, exec_path.name))

        for file_name in CONFIG_SOURCE_SEL:
            config_path = Path(ROOT_REPO_FOLDER, file_name)
            if file_name == 'mhm.nml':
                _configure_mhm_nml(config_path, Path(config_target_root, file_name),
                                   output_path=output_target_root,
                                   forcing_path=Path(ROOT_TARGET_FOLDER, FORCING_TARGET_FOLDER, basin),
                                   routing_path=Path(ROOT_TARGET_FOLDER, STATIC_TARGET_FOLDER, basin, 'routing'),
                                   # this is a custom hack, make more flexible if needed
                                   gauge_id=int(re.search(r'(\d)+', basin).group()),
                                   start_date=date(2000, 1, 1),
                                   end_date=date(2004, 12, 31),
                                   spinup=1825,
                                   )
            elif file_name == 'mpr_fso.nml':
                data_arrays_props = {var: {
                    'from_file': str(Path(ROOT_TARGET_FOLDER, STATIC_TARGET_FOLDER, basin, 'mpr', '{}.nc'.format(var)))}
                         for var in MPR_READ_VARS}
                data_arrays_props['lat_l0'] = {
                    'from_file': str(Path(ROOT_TARGET_FOLDER, STATIC_TARGET_FOLDER, basin, 'mpr', 'latlon.nc'))}
                _configure_mpr_nml(config_path, Path(config_target_root, TARGET_MPR_NML_NAME),
                                   coords={
                                       'land_cover_period_out': {
                                           'coord_from_values': [2000, 2006, 2011],
                                           'coord_from_values_bound': 1950,
                                           'coord_stagger': 'end',
                                       },
                                       'lat_out': {
                                           'coord_from_range_step': 4000
                                       },
                                       'lon_out': {
                                           'coord_from_range_step': 4000
                                       },
                                       'horizon_out': {
                                           'coord_from_values': [0.05, 0.25, 2.0],
                                           'coord_from_values_bound': 0.0,
                                           'coord_stagger': 'end',
                                       },
                                       'horizon_till': {
                                           'coord_from_values': [0.05, 0.25],
                                           'coord_from_values_bound': 0.0,
                                           'coord_stagger': 'end',
                                       },
                                       'horizon_notill': {
                                           'coord_from_values': [2.0],
                                           'coord_from_values_bound': 0.25,
                                           'coord_stagger': 'end',
                                       },
                                   },
                                   data_arrays=data_arrays_props,
                                   out_filename=str(Path(output_target_root, 'mHM_parameters.nc'))
                                   )
            else:
                # create new file
                print('Creating {} file: {}'.format('config', Path(config_target_root, file_name)))
                _copy_or_link_file(config_path, Path(config_target_root, file_name))


def make_dir(path):
    if not path.exists():
        path.mkdir(parents=True, exist_ok=True)


def prepare_statics(basins):
    exec_path = Path(ROOT_REPO_FOLDER, 'pre-proc', 'asc2nc.py')
    for basin in basins:
        static_source_root = Path(ROOT_SOURCE_FOLDER, basin)
        static_target_root = Path(ROOT_TARGET_FOLDER, STATIC_TARGET_FOLDER, basin)
        make_dir(static_target_root)

        # execute
        command = ['python', exec_path, '-i', str(static_source_root), '-o', str(static_target_root)]
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        if out:
            print('output', out.decode("utf-8"))
        if err:
            print('error', err.decode("utf-8"))

        _manipulate_land_cover_period(static_target_root)
        _manipulate_gauge_timeseries(basin, static_source_root, static_target_root)
        _manipulate_latlon_file(static_target_root)


def _manipulate_latlon_file(static_target_root):
    # hack to fix the latlon bug
    # the lat and lat_l0 values are to high ~70, so they are recalculated based on yc, yc_l0
    lc_file = Path(static_target_root, 'latlon.nc')
    if lc_file.exists():
        lc_file.unlink()
    tmp_lc_file = Path(static_target_root, 'latlon_v5.5.nc')
    lc_file = Path(static_target_root, 'mpr', 'latlon.nc')
    ds = xr.open_dataset(tmp_lc_file)
    average_xc = float(ds['xc'].mean())
    # get the correct epsg code for the basin
    in_epsg = None
    if average_xc > 3000000 and average_xc < 3900000:
        in_epsg = 31467
    elif average_xc > 3900000 and average_xc < 4700000:
        in_epsg = 31468
    elif average_xc > 4700000 and average_xc < 6000000:
        in_epsg = 31469
    if in_epsg is not None:
        out_epsg = 4326
        transformer = Transformer.from_crs(in_epsg, out_epsg, always_xy=True)
        for level in ['', '_l0']:
            xc = np.tile(ds['xc{}'.format(level)], (len(ds['yc{}'.format(level)]), 1))
            yc = np.tile(ds['yc{}'.format(level)], (len(ds['xc{}'.format(level)]), 1)).T
            lons, lats = transformer.transform(xc, yc)
            ds['lon{}'.format(level)].data = lons
            ds['lat{}'.format(level)].data = lons
            for var in ['lon{}'.format(level), 'lat{}'.format(level)]:
                if 'missing_value' in ds[var].encoding:
                    ds[var].encoding.pop('missing_value')
        ds.to_netcdf(lc_file, encoding={var: COMPRESSION for var in ds.data_vars.keys()})
        tmp_lc_file.unlink()
    else:
        ds.close()
        shutil.move(tmp_lc_file, lc_file)


def _manipulate_land_cover_period(static_target_root):
    # hack to manipulate the land_cover_period coordinate
    lc_file = Path(static_target_root, 'mpr', 'land_cover.nc')
    tmp_lc_file = Path(static_target_root, 'mpr', 'land_cover.nc.bak')
    shutil.move(lc_file, tmp_lc_file)
    ds = xr.open_dataset(tmp_lc_file)
    attrs = ds['land_cover_period'].attrs
    ds['land_cover_period'] = [1950, 2000, 2006]
    ds['land_cover_period'].attrs = attrs
    attrs = ds['land_cover_period_bnds'].attrs
    ds['land_cover_period_bnds'].data = [[1950, 2000], [2000, 2006], [2006, 2011]]
    ds['land_cover_period_bnds'].attrs = attrs
    ds.to_netcdf(lc_file, encoding={var: COMPRESSION for var in ds.data_vars.keys()})
    tmp_lc_file.unlink()


def _manipulate_gauge_timeseries(basin, static_source_root, static_target_root):
    # hack to copy the respective gauge file
    filename = '{}.txt'.format(int(re.search(r'(\d)+', basin).group()))
    _copy_or_link_file(Path(static_source_root, 'gauge', filename),
                       Path(static_target_root, 'routing', filename))


# CLASSES

# SCRIPT


if __name__ == '__main__':
    # get all the paths for each basin existing on data archive
    basin_names = [path.name for path in sorted(Path(ROOT_SOURCE_FOLDER).glob('sub*')) if
                   any([sel in path.name for sel in BASIN_SELECTION])]
    if 'forcing' in PROCESSES:
        prepare_forcings(basin_names)
    if 'static' in PROCESSES:
        prepare_statics(basin_names)
    if 'config' in PROCESSES:
        prepare_config(basin_names)
