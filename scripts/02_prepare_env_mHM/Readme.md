# Parameter sets

The original mHM source code as used in the Zink et al. paper is contained on folder [mHM_v4.3_Basin_inflow](mHM_v4.3_Basin_inflow).
It uses different transfer functions and parameters for some processes. An attempt to reverse-engineer these was done by establishing an [MPR namelist](zink_parameters_new_format/mpr_predictor_FSO_only_mhmv43.nml).
Note that another TF is also used for Ks in mHM v.43 and FieldCap is not directly calculated (only L1_LP might serve as an equivalent)! 

Original parameter sets as in Zink et al. as stored by Matthias Zink on Eve are in folder [zink_original_parameters](zink_original_parameters).
The processed one (e.g. parameter set with id 1 selected by [convert_parameter_values.py](convert_parameter_values.py)) to be used by the new MPR/mHM are in folder [mpr_predictor_FSO_only_mhmv43.nml](zink_parameters_new_format/mpr_predictor_FSO_only_mhmv43.nml).

# Data processing

In order to create the simulation environment the script [prepare_env.py](prepare_env.py) is run on Eve. The stuff to be edited needs to be hard-coded in the global parameters.

After that the script [prepare_env.py](modify_soil_depth.py) needs to be run, which alters the soil data in basins where only soils with soil depths up to 1.2m exist. They are set to 2.0 to be in line with the previous version and mHM layer depth.

After that the addon_danube folder needs to be copied to static/sub_6342800/routing containing inflow data for the Danube basin. The file [mhm](mhm.nml) needs to copied accordingly.

The script [submit_mhm_fso_restart.sh](submit_mhm_fso_restart.sh) sends all jobs to the eve cluster to execute one forward run. After that [check_restart_runs.py](check_restart_runs.py) scans the stdout and collects some information on the status.

 