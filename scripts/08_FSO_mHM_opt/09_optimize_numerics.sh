#!/bin/sh
#SBATCH -J numeric_opt
#SBATCH -N 1
#SBATCH --qos=normal_0128
#SBATCH --partition=mem_0128
#SBATCH -L intel@vsc
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=optimize_numerics.out
module purge
module load intel/18 intel-mkl/2018 pcre2/10.35 R/4.0.2
cd /home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt
Rscript 09_optimize_numerics.R
