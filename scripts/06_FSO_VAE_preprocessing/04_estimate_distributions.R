# Estimate TF distributions for FSO
# FSO-mHM project
# Moritz Feigl, Sep 2020


# wd depending on system
if(grepl("Ubuntu", Sys.info()[["version"]])){
  setwd("/home/cfgrammar/Dropbox/Diss/FSO_mHM/FSO_Data")
}
if(!grepl("Ubuntu", Sys.info()[["version"]])){
  setwd("/home/lv71468/mfeigl/FSO_mHM/FSO_Data")
  .libPaths("/home/lv71468/mfeigl/R/x86_64-pc-linux-gnu-library/4.0")
  chooseCRANmirror(ind = 6)
}
# Load FSO package from github
if("FSO" %in% list.files() == FALSE) {
  system("git clone https://github.com/MoritzFeigl/FSO.git")
  install.packages("FSO",
                   repos = NULL,
                   type = "source")
}
library(FSO)

setwd("prepared_TFs/prepared_functions")


# Create VSC skripts and folder for distribution estimation
lines <- readLines("../../../2020_fso_mhm/scripts/06_FSO_VAE_preprocessing/TEMPLATE_distribution_estimation.R")
if(!dir.exists("../scripts/distribution_estimation_scripts")) dir.create("../scripts/distribution_estimation_scripts")

folders <- list.files()
batch <- 0
for(folder in folders){
  files <- list.files(folder)
  for(file in files){
    batch <- batch + 1
    lines[18] <- paste0('functions_for_dist <- "', folder, '/', file, '"')
    writeLines(lines,
               con = paste0(
                 "../scripts/distribution_estimation_scripts/distribution_estimation_batch",
                 batch, ".R"))
  }
}

# Write bash file
setwd("../scripts/distribution_estimation_scripts")
fileConn <- file("distribution_estimation.sh")
writeLines(
  c("#!/bin/sh",
    "#SBATCH -J tf_distribution",
    "#SBATCH -N 1",
    "#SBATCH --qos=normal_0064",
    "#SBATCH --partition=mem_0064",
    "#SBATCH --mem=64000M",
    paste0("#SBATCH --array=1-", batch, "%", min(c(60, batch))),
    "#SBATCH --mail-user=moritz.feigl@boku.ac.at",
    "#SBATCH --mail-type=BEGIN,END,FAIL",
    "module purge",
    "module load intel/18 intel-mkl/2018 pcre2/10.35 R/4.0.2 gcc",
    "Rscript distribution_estimation_batch${SLURM_ARRAY_TASK_ID}.R"),
  fileConn)
close(fileConn)

system("sbatch distribution_estimation.sh")

