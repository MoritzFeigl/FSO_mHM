# Prepare transfer functions for FSO
# FSO-mHM project
# Moritz Feigl, Aug 2020


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
if("FSO" %in% list.files("..") == FALSE) {
  system("git clone https://github.com/MoritzFeigl/FSO.git")
  install.packages("FSO",
                   repos = NULL,
                   type = "source")
}
library(FSO)

if(!dir.exists("prepared_TFs")) dir.create("prepared_TFs")
setwd("prepared_TFs")

# Split sampled grammar into multiple batches
create_fst_batch_files <- function(files, number_of_batches_per_file){
  # Script to split available fst files in smaller batches for distributed computing
  # split all files
  count <- 0
  file_names <- gsub(".fst", "",
                     tail(unlist(strsplit(files[1], split = "/")), 1), fixed = TRUE)
  for(k in seq_along(files)){
    if(!dir.exists(paste0(file_names[k], "_batches"))){
      dir.create(paste0(file_names[k], "_batches"))
    }
    functions_para <- fst::read_fst(files[k])
    batch_cuts <- cut(1:nrow(functions_para),
                      breaks = number_of_batches_per_file)
    functions_batches <- split.data.frame(functions_para, batch_cuts)
    for(i in 1:number_of_batches_per_file){
      count <- count + 1
      fst::write.fst(
        functions_batches[[i]],
        path = paste0(file_names[k], "_batches/", file_names[k], "_batch",
                      count, ".fst"))
    }
  }
}

create_fst_batch_files("../Grammar_samples/sampled_grammar.fst", 60)



# Create grammar sampling scripts and corresponding bash files
lines <- readLines("../../2020_fso_mhm/scripts/06_FSO_VAE_preprocessing/TEMPLATE_prepare_tfs.R")
n_batches <- length(list.files("sampled_grammar_batches"))
seed_start <- sample(1000:10000, 1)
if(!dir.exists("scripts")) dir.create("scripts")
for(batch in 1:n_batches){
  # KSat variable inputs and simplify
  lines[21] <- paste0('variable_input(functions = "', "sampled_grammar_batches/",
                      "sampled_grammar_batch",
                      batch, '.fst",')
  lines[22] <- paste0("variable_list = KSat_variable_list,")
  lines[25] <- paste0("n_iter = 10, no_cores = 32, seed = ", seed_start+batch, ",")
  lines[26] <- paste0('         file_name = "prepared_functions/KSat_functions_batch_', batch,  '")')
  lines[27] <- paste0('setwd("prepared_functions/KSat_functions_batch_', batch, '")')
  lines[28] <- paste0('simplify_functions(files = paste0("KSat_functions_batch_', batch, '_batch", 1:10, ".fst"),')
  writeLines(lines, con = paste0("scripts/KSat_variable_input_batch_", batch, ".R"))
  # FieldCap variable inputs and simplify
  lines[21] <- paste0('variable_input(functions = "', "sampled_grammar_batches/",
                      "sampled_grammar_batch",
                      batch, '.fst",')
  lines[22] <- paste0("variable_list = FieldCap_variable_list,")
  lines[25] <- paste0("n_iter = 10, no_cores = 32, seed = ", seed_start+batch, ",")
  lines[26] <- paste0('         file_name = "prepared_functions/FieldCap_functions_batch_', batch,  '")')
  lines[27] <- paste0('setwd("prepared_functions/FieldCap_functions_batch_', batch, '")')
  lines[28] <- paste0('simplify_functions(files = paste0("FieldCap_functions_batch_', batch, '_batch", 1:10, ".fst"),')
  writeLines(lines, con = paste0("scripts/FieldCap_variable_input_batch_", batch, ".R"))
}

# create folder for prepared TFs
if(!dir.exists("prepared_functions")) dir.create("prepared_functions")

# Write bash file
setwd("scripts")
fileConn <- file("prepare_tfs_KSat.sh")
writeLines(
  c("#!/bin/sh",
    "#SBATCH -J Prepare_KSat_tfs",
    "#SBATCH -N 1",
    "#SBATCH --qos=normal_0064",
    "#SBATCH --partition=mem_0064",
    "#SBATCH --mem=64000M",
    paste0("#SBATCH --array=1-", n_batches, "%", min(c(20, n_batches))),
    "#SBATCH --mail-user=moritz.feigl@boku.ac.at",
    "#SBATCH --mail-type=BEGIN,END,FAIL",
    "module purge",
    "module load intel/18 intel-mkl/2018 pcre2/10.35 R/4.0.2 gcc",
    "Rscript KSat_variable_input_batch_${SLURM_ARRAY_TASK_ID}.R"),
  fileConn)
close(fileConn)

fileConn <- file("prepare_tfs_FieldCap.sh")
writeLines(
  c("#!/bin/sh",
    "#SBATCH -J Prepare_FieldCap_tfs",
    "#SBATCH -N 1",
    "#SBATCH --qos=normal_0064",
    "#SBATCH --partition=mem_0064",
    "#SBATCH --mem=64000M",
    paste0("#SBATCH --array=1-", n_batches, "%", min(c(20, n_batches))),
    "#SBATCH --mail-user=moritz.feigl@boku.ac.at",
    "#SBATCH --mail-type=BEGIN,END,FAIL",
    "module purge",
    "module load intel/18 intel-mkl/2018 pcre2/10.35 R/4.0.2 gcc",
    "Rscript FieldCap_variable_input_batch_${SLURM_ARRAY_TASK_ID}.R"),
  fileConn)
close(fileConn)

# Run grammar sampling batches
system("sbatch prepare_tfs_KSat.sh")
system("sbatch prepare_tfs_FieldCap.sh")

# remove all non simplified fst files
setwd("../prepared_functions")
folders <- list.files()
for(folder in folders){
  files <- list.files(folder)
  remove_files <- files[!grepl("_simplif", files)]
  file.remove(paste0(folder, "/", remove_files))
}


