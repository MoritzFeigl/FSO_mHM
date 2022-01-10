#!/bin/sh
#SBATCH -J VAE_data
#SBATCH -N 1
#SBATCH --qos=normal_0064
#SBATCH --partition=mem_0064
#SBATCH --mem=64000M
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=05_VAE_data.out
module purge
module load intel/18 intel-mkl/2018 pcre2/10.35 R/4.0.2 gcc/7.3
cd /home/lv71468/mfeigl/WT_Project/Modelling/final_models
Rscript 05_VAE_data.R
