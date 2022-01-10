setwd("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt")
script <- readLines("10_optimize_numerics_for_each_run.R")
name_line <- grep('run_name <- ', script, fixed = TRUE)
sbatch_script <- readLines("10_optimize_numerics_for_each_run.sh")

run_names <- c(paste0("FSO_SCE_NSE_run_", 1:4), "FSO_DDS_NSE_run_5")

for(run_name in run_names){
  run_folder <- paste0("/home/lv71468/mfeigl/FSO_mHM/Runs/", run_name)
  script[name_line] <- paste0('run_name <- "', run_name, '"')
  queue_name <- tail(unlist(strsplit(run_name, "_")), 1)
  sbatch_script[2] <- paste0("#SBATCH -J num_opt_", queue_name)
  sbatch_script[9] <- paste0("#SBATCH --output=num_opt_", queue_name, ".out")
  sbatch_script[12] <- paste0("cd ", run_folder)
  
  writeLines(script, paste0(run_folder, "/10_optimize_numerics_for_each_run.R"))
  writeLines(sbatch_script, paste0(run_folder, "/10_optimize_numerics_for_each_run.sh"))
  setwd(run_folder)
  system("sbatch 10_optimize_numerics_for_each_run.sh")
}
