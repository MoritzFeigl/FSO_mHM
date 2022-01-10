#!/usr/bin/env python
# encoding: utf-8

"""
File Name   : convert_parameter_values
Project Name: 2020_FSO_mHM
Description : insert your description here, if applicable
Author      : ottor
Created     : 07.12.20 10:12
"""

# IMPORTS
import pandas as pd
import f90nml

# GLOBAL VARIABLES
N_PARAM = 68
NEW_NAMES_FILE = 'scripts/04_prepare_predictors_for_FSO/mpr_global_parameter_mhm.nml'
NAMES_SORTED = [
    'canopyInterceptionFactor',
    'snowThresholdTemperature',
    'degreeDayFactor_forest',
    'degreeDayFactor_impervious',
    'degreeDayFactor_pervious',
    'increaseDegreeDayFactorByPrecip',
    'maxDegreeDayFactor_forest',
    'maxDegreeDayFactor_impervious',
    'maxDegreeDayFactor_pervious',
    'orgMatterContent_forest',
    'orgMatterContent_impervious',
    'orgMatterContent_pervious',
    'PTF_lower66_5_constant',
    'PTF_lower66_5_clay',
    'PTF_lower66_5_Db',
    'PTF_higher66_5_constant',
    'PTF_higher66_5_clay',
    'PTF_higher66_5_Db',
    'infiltrationShapeFactor',
    'imperviousStorageCapacity',
    None,  # a scalar factor applied to pwp, default to 1, gone in mHM v.5
    None,  # a scalar factor applied to lower threshold, where aET is reduced linearly till PWP, default to 1, gone in mHM v.5
    'rootFractionCoefficient_forest',
    'rootFractionCoefficient_impervious',
    'rootFractionCoefficient_pervious',
    None,  # this is now interflowStorageCapacityFactor but the TF changed
    'interflowRecession_slope',  # 27 and 29 are same in mHM v.5 TF changed
    None,  # this parameter is gone in mHM v.5
    'interflowRecession_slope',  # 27 and 29 are same in mHM v.5 TF changed
    'slowInterflowRecession_Ks',  # TF changed
    'PTF_Ks_constant',  # TF changed
    'PTF_Ks_sand',  # TF changed
    'PTF_Ks_clay',  # TF changed
    None,  # this is now exponentSlowInterflow but the TF changed
    'rechargeCoefficient',
    None,  # TF changed
    'rechargeFactor_karstic', # TF changed but this part is retained
    'minCorrectionFactorPET',
    'maxCorrectionFactorPET',
    'aspectThresholdPET',
    None,  # routing
    None,  # routing
    None,  # routing
    None,  # routing
    None,  # routing
    'base_flow_geo_unit_01',
    'base_flow_geo_unit_02',
    'base_flow_geo_unit_03',
    'base_flow_geo_unit_04',
    'base_flow_geo_unit_05',
    'base_flow_geo_unit_06',
    'base_flow_geo_unit_07',
    'base_flow_geo_unit_08',
    'base_flow_geo_unit_09',
    'base_flow_geo_unit_10',
    'base_flow_geo_unit_11',
    'base_flow_geo_unit_12',
    'base_flow_geo_unit_13',
    'base_flow_geo_unit_14',
    'base_flow_geo_unit_15',
    'base_flow_geo_unit_16',
    'base_flow_geo_unit_17',
    'base_flow_geo_unit_18',
    'base_flow_geo_unit_19',
    'base_flow_geo_unit_20',
    'base_flow_geo_unit_21',
    'base_flow_geo_unit_22',
    'base_flow_geo_unit_23',
]
"""
# unused parameters
'exponentSlowInterflow',
'rootFractionCoefficient_sand',
'rootFractionCoefficient_clay',
'jarvis_sm_threshold_c1',
'PTF_Ks_curveSlope',  # TF changed
'interflowStorageCapacityFactor',
'fastInterflowRecession_forest',
'HargreavesSamaniCoeff',
'canopyheight_forest',
'canopyheight_impervious',
'canopyheight_pervious',
'displacementheight_coeff',
'roughnesslength_momentum_coeff',
'roughnesslength_heat_coeff',
'stomatal_resistance',
'PriestleyTaylorCoeff',
'PriestleyTaylorLAIcorr',
'PET_a_forest',
'PET_a_impervious',
'PET_a_pervious',
'PET_b',
'PET_c'
"""


# FUNCTIONS
def read_orig_param(filename):
    col_names = [f'gamma({id:03})' for id in range(1, N_PARAM + 1)]
    new_col_names = [f'gamma{id:03}' for id in range(1, N_PARAM + 1)]
    df = pd.read_fwf(filename, usecols=col_names).rename(columns=dict(zip(col_names, new_col_names)))
    return df


# CLASSES


# SCRIPTÎ
if __name__ == '__main__':
    sel_index = 0
    filename = './zink_original_parameters/filtered_param_sm_idchange.sets'
    target_filename = 'zink_parameters_new_format/mpr_global_parameter_predictor_FSO_only_mhmv43.nml'
    # filename = './scripts/02_prepare_env_mHM/zink_original_parameters/filtered_param_sm_idchange.sets'
    df = read_orig_param(filename)
    df_sel = df.iloc[0, :]
    nml = f90nml.Namelist({'Parameters':
                               f90nml.Namelist({
                                   'parameter_names': list(df_sel.index),
                                   'parameter_values': list(df_sel),
                               })})
    nml.write(nml_path=target_filename, force=True)
