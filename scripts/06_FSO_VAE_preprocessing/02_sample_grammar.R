# Sample Grammar for FSO
# FSO-mHM project
# Moritz Feigl, Aug 2020


# wd depending on system
if(grepl("Ubuntu", Sys.info()[["version"]])){
  setwd("/home/cfgrammar/Dropbox/Diss/FSO_mHM/")
}
if(!grepl("Ubuntu", Sys.info()[["version"]])){
  setwd("/home/lv71468/mfeigl/FSO_mHM")
  .libPaths("/home/lv71468/mfeigl/R/x86_64-pc-linux-gnu-library/4.0")
  chooseCRANmirror(ind = 6)
}
# Load FSO package from github
if("FSO" %in% list.files() == FALSE) {
  system("git clone https://github.com/MoritzFeigl/FSO.git")
  install.packages("FSO", 
                   repos = NULL, 
                   type = "source", dependencies = TRUE)
}
# create Data and grammar directories
if(!dir.exists("FSO_Data")) dir.create("FSO_Data")
setwd("FSO_Data")
if(!dir.exists("Grammar_samples")) dir.create("Grammar_samples")
setwd("Grammar_samples")
library(FSO)

# Create grammar sampling scripts and corresponding bash files
lines <- readLines("../../2020_fso_mhm/scripts/06_FSO_VAE_preprocessing/TEMPLATE_sample_grammar.R")
n_batches <- 10
n_grammar_samples <- 36000000
seed_start <- sample(1000:10000, 1)
if(!dir.exists("scripts")) dir.create("scripts")
for(batch in 1:n_batches){
  lines[24] <- paste0("funs <- grammar_sampler(n = ", 
                      ceiling(n_grammar_samples/n_batches), 
                      ", grammar = grammar, max_depth = 5,")
  lines[25] <- paste0("no_cores = 32, seed = ", seed_start+batch, ", save = TRUE,")
  lines[26] <- paste0('         file_name = "sampled_grammar_batch_', batch,  '")')
  writeLines(lines,
             con = paste0(
               "scripts/sample_grammar_batch", batch, ".R"
             )
  )
}

# Write bash file
setwd("scripts")
fileConn <- file("sample_grammar_batches.sh")
writeLines(
  c("#!/bin/sh",
    "#SBATCH -J Grammar_sampling",
    "#SBATCH -N 1",
    "#SBATCH --qos=normal_0064",
    "#SBATCH --partition=mem_0064",
    "#SBATCH --mem=64000M",
    paste0("#SBATCH --array=1-", n_batches, "%", min(c(10, n_batches))),
    "#SBATCH --mail-user=moritz.feigl@boku.ac.at",
    "#SBATCH --mail-type=BEGIN,END,FAIL",
    "module purge",
    "module load intel/18 intel-mkl/2018 pcre2/10.35 R/4.0.2 gcc",
    "Rscript sample_grammar_batch${SLURM_ARRAY_TASK_ID}.R"),
  fileConn)
close(fileConn)

# Run grammar sampling batches
system("sbatch sample_grammar_batches.sh")

setwd("..")
for(batch in 1:n_batches){
  sampled_batch <- fst::read_fst(paste0("sampled_grammar_batch_", batch, ".fst"))
  if(exists("grammar_samples")){
    grammar_samples <- rbind(grammar_samples, sampled_batch)
  } else grammar_samples <- sampled_batch
}

grammar_samples <- data.frame(functions = unique(grammar_samples$functions))
fst::write_fst(grammar_samples, "sampled_grammar.fst")
file.remove(paste0("sampled_grammar_batch_", 1:n_batches, ".fst"))








