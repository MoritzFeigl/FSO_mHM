#!/bin/sh
#SBATCH -J tf_distribution
#SBATCH -N 1
#SBATCH --qos=normal_0064
#SBATCH --partition=mem_0064
#SBATCH --mem=64000M
#SBATCH --mail-user=moritz.feigl@boku.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=04_estimate_distribution.out
module purge
module load intel/18 intel-mkl/2018 pcre2/10.35 R/4.0.2 gcc/7.3
Rscript 04_estimate_distributions.R
