#!/bin/sh
#SBATCH -J prepare_tfs
#SBATCH -N 1
#SBATCH --qos=normal_0064
#SBATCH --partition=mem_0064
#SBATCH --mem=64000M
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=03_prepare_tfs.out
module purge
module load intel/18 intel-mkl/2018 pcre2/10.35 R/4.0.2 gcc/7.3
Rscript 03_prepare_tfs.R
