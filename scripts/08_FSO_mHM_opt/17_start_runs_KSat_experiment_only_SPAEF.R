
setwd("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt")
dir.create("/home/lv71468/mfeigl/FSO_mHM/Runs", showWarnings = FALSE)
# Create run scripts for KSat experiment with only SPAEF runs
runs <- 1:5

for(run in runs){
  run_name <- paste0("KSat_experiment_spaef_", run)
  run_folder <- paste0("/home/lv71468/mfeigl/FSO_mHM/Runs/", run_name)
  # create FSO opt script
  dir.create(run_folder, showWarnings = FALSE)

    file.copy("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/16_run_KSat_experiment_only_SPAEF.R",
      paste0(run_folder, "/16_run_KSat_experiment_only_SPAEF.R"), overwrite = TRUE)
  script <- readLines(paste0(run_folder, "/16_run_KSat_experiment_only_SPAEF.R"))
  script <- gsub('run_name <- "KSat_experiment_spaef"',
                 paste0('run_name <- "', run_name, '"'),
                 script, fixed = TRUE)

  writeLines(script, paste0(run_folder, "/16_run_KSat_experiment_only_SPAEF.R"))
  # get sbatch script
  file.copy("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/16_run_KSat_experiment_only_SPAEF.sh",
            paste0(run_folder, "/16_run_KSat_experiment_only_SPAEF.sh"), overwrite = TRUE)
  start_fso_mhm <- readLines(paste0(run_folder, "/16_run_KSat_experiment_only_SPAEF.sh"))
  start_fso_mhm <- gsub("#SBATCH -J ksat_exp_spaef",
                        paste0("#SBATCH -J ksat_exp_spaef_", run), start_fso_mhm)
  start_fso_mhm <- gsub("cd /home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt",
                        paste0("cd ", run_folder), start_fso_mhm)
  writeLines(start_fso_mhm, paste0(run_folder, "/16_run_KSat_experiment_only_SPAEF.sh"))
  
  
  # start optimization
  setwd(run_folder)
  system("chmod 775 16_run_KSat_experiment_only_SPAEF.sh")
  system("sbatch 16_run_KSat_experiment_only_SPAEF.sh")
}
