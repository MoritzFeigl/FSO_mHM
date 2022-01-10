#!/bin/sh
#SBATCH -J samp_val
#SBATCH -N 1
#SBATCH --qos=vsc3plus_0256
#SBATCH --partition=vsc3plus_0256
#SBATCH -L intel@vsc
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=02_sampled_basin_validation.out
module purge
module load intel/18 intel-mkl/2018 pcre2/10.35 R/4.0.2
cd /home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/10_validation
Rscript 01_sampled_basins_validation.R
