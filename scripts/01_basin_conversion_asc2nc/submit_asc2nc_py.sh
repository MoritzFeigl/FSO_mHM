#!/bin/bash

# ----------------------------------------------------------------------
# qsub arguments
# ----------------------------------------------------------------------
echo qsub arguments

#$ -N asc2nc.py
#$ -S /bin/bash
#$ -o stdout/
#$ -e stderr/
#$ -l h_rt=00:10:00
#$ -l h_vmem=16G
#$ -t 1-7
#$ -j yes
#$ -binding linear:1
#$ -cwd 

# ----------------------------------------------------------------------
# load required modules
# ----------------------------------------------------------------------
echo load required modules


# load python environment
sys_name=$(uname -n)
if [[ ${sys_name:0:8} != 'frontend' ]]; then
  # unload all modules
  module purge

  # remove all module directories
  unset MODULEPATH

  # use the EasyBuild module directory
  module use /software/easybuild-broadwell/modules/all/Core
  module use /software/modulefiles
  module use ~/modules
  module use /global/apps/modulefiles
  module load python_env_geo
fi

# ----------------------------------------------------------------------
# main script
# ----------------------------------------------------------------------
echo main script


# get all basins
#basin_ids="Mulde Ems Neckar Saale Main Weser Danube"
basin_ids="6340600 6338100 6335600 6340300 6335304 6337200 6342800"

# get current basin by indexing
task_id=$SGE_TASK_ID
if [[ "${task_id}z" == "z" ]]; then
    task_id=1
fi
basin_id=$(echo "${basin_ids}" | cut -d " " -f ${task_id})

cd "$HOME"/lib/mhm_mpr/pre-proc || exit
in_path=/data/stohyd/data/processed/Germany/basins/sub_${basin_id}
out_path=/work/ottor/FSO_mHM_major_basins/static/sub_${basin_id}

# execute
python asc2nc.py -i "${in_path}" -o "${out_path}"
