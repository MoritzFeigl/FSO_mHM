# Prepare TFs for VAE 
# FSO-mHM project
# Moritz Feigl, Sep 2020


setwd("/home/lv71468/mfeigl/FSO_mHM/FSO_Data/prepared_TFs/prepared_functions")
.libPaths("/home/lv71468/mfeigl/R/x86_64-pc-linux-gnu-library/4.0")

library(tensorflow)
library(tokenizers)
library(magrittr)
# 1. estimate bounds for distribution values ---------------------------------------------
folders <- list.files()
files <- list.files(folders, full.names = TRUE)
files <- files[grepl("distribution", files)]
KSat_files <- files[grepl("KSat", files)]
FieldCap_files <- files[grepl("FieldCap", files)]

for(var in c("KSat", "FieldCap")){
  if(!file.exists(paste0("../../VAE_data/", var, "/distribution_cutoff_points.csv"))){
    count <- 0
    file_set <- get(paste0(var, "_files"))
    for(file in file_set){
      count <- count + 1
      cat(paste0(count,"/", length(file_set)), "\r")
      file_data <- fst::read_fst(file)
      try(rm(quant_data))
      if(!exists("quant_data")) {
        quant_data <- file_data[, c("10%", "90%")]
      } else {
        quant_data <- rbind(quant_data, file_data[, c("10%", "90%")])
      }
    }
    cat("\n")
    # estimate cutoff points
    cutoff <- c(median(quant_data$`10%`, na.rm = TRUE) - mad(quant_data$`10%`, na.rm = TRUE)*3, 
                median(quant_data$`90%`, na.rm = TRUE) + mad(quant_data$`90%`, na.rm = TRUE)*3)
    try(dir.create("../../VAE_data"))
    try(dir.create(paste0("../../VAE_data/", var)))
    write.csv(data.frame("lower_bound" = cutoff[1], "upper_bound" = cutoff[2]), 
              paste0("../../VAE_data/", var, "/distribution_cutoff_points.csv"), 
              row.names = FALSE)
  }
}

# 2. Estimate means and sds for distribution scaling -------------------------------------

# sample distribution parameters with 200 samples
for(var in c("KSat", "FieldCap")){
  if(!file.exists(paste0("../../VAE_data/", var, "/distribution_scale_parameters.csv"))){
  cat("*** Estimating distribution parameters of", var, "***\n")
    count <- 0
    file_set <- get(paste0(var, "_files"))
    file_samples <- file_set[sample(length(file_set), 200)]
    # concatenate samples
    for(file in file_samples){
      count <- count + 1
      cat(paste0(count,"/", length(file_samples)), "\r")
      file_data <- fst::read_fst(file)
      try(rm(quant_data))
      if(!exists("quant_data")) {
        quant_data <- file_data[, -1]
      } else {
        quant_data <- rbind(quant_data, file_data[, -1])
      }
    }
    cat("\n")
    # remove TFs projecting outside the cutoff range
    cutoff <- read.csv(paste0("../../VAE_data/", var, "/distribution_cutoff_points.csv"))
    quant_data <- quant_data[quant_data$`10%` > cutoff$lower_bound & 
                               quant_data$`90%`< cutoff$upper_bound, ]
    # estimate distribution parameters
    means <- sapply(quant_data, mean)
    sds <- sapply(quant_data, sd)
    write.csv(data.frame("mean" = means, "sd" = sds), 
              paste0("../../VAE_data/", var, "/distribution_scale_parameters.csv"), 
              row.names = FALSE)
  }
}
# 3. Define dictionary -------------------------------------------------------------------
dict_from_list <- function(...){
  dict <- c(...)
  dictionary <- as.list(1:length(dict))
  names(dictionary) <- dict
  return(dictionary)
}
f <- c('exp', 'log10', 'log', 'sin', 'cos', 'tan', 'abs', 'acos',
       'asin', 'atan', 'cosh', 'sinh', 'sqrt')
var_KSat <- c("bd", "sand", "clay", "slope", "aspect", "dem")
var_FieldCap <- c("bd", "sand", "clay", "slope", "aspect", "dem", "thetas", "ksat", "vgenu_n")
operators <- c('+', '-', '*', '/', "^", "(", ")")
numbers <- seq(0.05, 3, 0.05) %>% format(nsmall=2)

# 4. Assess max seq length of data -------------------------------------------------------
dictionary_KSat <- dict_from_list(var_KSat, f, operators, numbers)
dictionary_FieldCap <- dict_from_list(var_FieldCap, f, operators, numbers)
dictionaries <- list(KSat = dictionary_KSat,
                     FieldCap = dictionary_FieldCap)

prepare_text_for_tokenize <- function(text, operators = "\\^\\+\\-\\/\\*\\(\\)"){
  # Function to prepare TFs for tokenizing
  # input: strings of transfer functions
  # output: strings with ";" seperator
  text <- stringr::str_replace_all(text, paste0('([', operators, '])'), " \\1 ")
  text <- tokenizers::tokenize_words(text, strip_punct = FALSE) %>%
    sapply(FUN = function(x) paste0(c("", paste0(x, collapse = ";"), ""), collapse = ";"))
  return(text)
}

# KSat and FieldCap are based on the same sampled grammar -> padding_length should be equal
var <- "KSat"
if(!file.exists(paste0("../../VAE_data/", var, "/padding_length.csv"))){
cat("*** Estimate Padding length for varible", var, "***\n")
  files <- get(paste0(var, "_files"))
  sample_files <- files[sample(length(files), 80)]
  count <- 0
  try(rm(sample_tfs), silent = TRUE)
  for(file in sample_files){
    count <- count + 1
    cat(paste0("Concatenated ", var, " function files: ", count,"/", length(sample_files)), "\r")
    file_data <- fst::read_fst(file)
    if(!exists("sample_tfs")) {
      sample_tfs <- file_data[, "function"]
    } else {
      sample_tfs <- c(sample_tfs, file_data[, "function"])
    }
  }
  cat("\n")
  # keras tokenizer with pre-defined dictionary
  tt <- keras::text_tokenizer(filters = "\\;")
  tt$word_index <- dictionaries[[var]]
  cat("Preparing", var, "functions for tokenizing:\n")
  prep_funs <- pbmcapply::pbmclapply(sample_tfs, prepare_text_for_tokenize, mc.cores = 32)
  cat("Tokenize functions...\n")
  funs_tokenized <- keras::texts_to_sequences(tt, prep_funs)
  max_length <- max(lengths(funs_tokenized))
  padding_length <- max_length + 5
  cat("Files max seqence length is", max_length, "and will be padded to", padding_length, "\n")
  cat("Saving padding length information to", 
      paste0("VAE_data/", var, "/padding_length.csv"), "\n")
  write.csv(
    data.frame("variable" = var, "max_seq_length" = max_length, "padding_length" = padding_length),
    paste0("../../VAE_data/", var, "/padding_length.csv"),
    row.names = FALSE)
  # save FieldCap padding_length
  var <- "FieldCap"
  write.csv(
    data.frame("variable" = var, "max_seq_length" = NA, "padding_length" = padding_length),
    paste0("../../VAE_data/", var, "/padding_length.csv"),
    row.names = FALSE)
}
# 5. Pad sequences of all files ----------------------------------------------------------
# Create VSC skripts and folder for VAE data generation
lines <- readLines("../../../2020_fso_mhm/scripts/06_FSO_VAE_preprocessing/TEMPLATE_VAE_data.R")
try(dir.create("../../VAE_data/KSat/batches"), silent = TRUE)
try(dir.create("../../VAE_data/FieldCap/batches"), silent = TRUE)
try(dir.create("../../VAE_data/scripts"), silent = TRUE)
# KSat
batch <- 0
for(file in KSat_files){
  batch <- batch + 1
  lines[41] <- 'var <- "KSat"'
  lines[42] <- paste0('file <- "', file, '"') 
  lines[43] <- paste0('batch <- ', batch)
  writeLines(lines, paste0("../../VAE_data/scripts/KSat_batch", batch, ".R"))
}
# FieldCap
batch <- 0
for(file in FieldCap_files){
  batch <- batch + 1
  lines[41] <- 'var <- "FieldCap"'
  lines[42] <- paste0('file <- "', file, '"') 
  lines[43] <- paste0('batch <- ', batch)
  writeLines(lines, paste0("../../VAE_data/scripts/FieldCap_batch", batch, ".R"))
}


# 6. Write bash files and run preprocessing ----------------------------------------------
setwd("../../VAE_data/scripts")
# KSat
batches <- sum(grepl("KSat", list.files()))
fileConn <- file("VAE_data_KSat.sh")
writeLines(
  c("#!/bin/sh",
    "#SBATCH -J VAE_data_KSat",
    "#SBATCH -N 1",
    "#SBATCH --qos=normal_0064",
    "#SBATCH --partition=mem_0064",
    "#SBATCH --mem=64000M",
    paste0("#SBATCH --array=1-", batches, "%", min(c(30, batches))),
    "#SBATCH --mail-user=moritz.feigl@boku.ac.at",
    "#SBATCH --mail-type=BEGIN,END,FAIL",
    "#SBATCH --output=KSat_batch%a.out",
    "module purge",
    "module load intel/18 intel-mkl/2018 pcre2/10.35 R/4.0.2 gcc",
    "Rscript KSat_batch${SLURM_ARRAY_TASK_ID}.R"),
  fileConn)
close(fileConn)
# FieldCap
batches <- sum(grepl("FieldCap", list.files()))
fileConn <- file("VAE_data_FieldCap.sh")
writeLines(
  c("#!/bin/sh",
    "#SBATCH -J VAE_data_FieldCap",
    "#SBATCH -N 1",
    "#SBATCH --qos=normal_0064",
    "#SBATCH --partition=mem_0064",
    "#SBATCH --mem=64000M",
    paste0("#SBATCH --array=1-", batches, "%", min(c(30, batches))),
    "#SBATCH --mail-user=moritz.feigl@boku.ac.at",
    "#SBATCH --mail-type=BEGIN,END,FAIL",
    "#SBATCH -o FieldCap_batch${SLURM_ARRAY_TASK_ID}.out",
    "#SBATCH --output=FieldCap_batch%a.out",
    "module purge",
    "module load intel/18 intel-mkl/2018 pcre2/10.35 R/4.0.2 gcc",
    "Rscript FieldCap_batch${SLURM_ARRAY_TASK_ID}.R"),
  fileConn)
close(fileConn)

system("sbatch VAE_data_FieldCap.sh")
system("sbatch VAE_data_KSat.sh")






