#!/bin/sh
#SBATCH -J run_mhm
#SBATCH -N 1
#SBATCH --qos=normal_0128
#SBATCH --partition=mem_0128
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=run_mhm.out
# FSO mHM optimization mHM runs

# Run mhm in parallel
cd /home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_scripts
module load intel/18 intel-mpi/2018 hdf5/1.8.12-MPI pnetcdf/1.10.0 netcdf_C/4.4.1.1 netcdf_Fortran/4.4.4 cmake/3.9.6
#parallel -j0 bash :::: <(ls run_basin_{1..7}.sh)
bash run_basin_1.sh &
bash run_basin_2.sh &
bash run_basin_3.sh &
bash run_basin_4.sh &
bash run_basin_5.sh &
bash run_basin_6.sh &
bash run_basin_7.sh
#bash run_basin_${SLURM_ARRAY_TASK_ID}.sh
