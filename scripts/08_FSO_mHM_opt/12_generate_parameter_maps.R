# mHM FSO parameter maps
# FSO-mHM project
# Moritz Feigl, Mar 2021

# define run name to create environment automatically
runs <- c(paste0("FSO_SCE_NSE_run_", 1:4), "FSO_DDS_NSE_run_5")
for(run_name in runs){
  cat("*** start Run ", run_name, " ***\n")
  gcc <- ""
  dir.create("/home/lv71468/mfeigl/FSO_mHM/Runs", showWarnings = FALSE)
  run_folder <- paste0("/home/lv71468/mfeigl/FSO_mHM/Runs/", run_name)
  result_folder <- paste0("/home/lv71468/mfeigl/FSO_mHM/Results/", run_name)
  setwd(run_folder)
  
  # utils
  source("04_FSO_mhm_utils.R")
  # Load VAE generators
  .libPaths("/home/lv71468/mfeigl/R/x86_64-pc-linux-gnu-library/4.0")
  library(magrittr)
  
  # Change dates to training dates
  path <- paste0(run_folder, "/FSO_mHM_major_basins/config")
  config_paths <- list.files(path, pattern = "sub_", full.names = TRUE)
  for(cpath in config_paths){
    mhm_nml <- readLines(paste0(cpath, "/mhm.nml"))
    start_date <- as.POSIXct("2000-01-01") # period start 
    # end is always 31.12.1999
    mhm_nml[grep("eval_per(1)%yend = ", mhm_nml, fixed = TRUE)] <- 
      "    eval_per(1)%yend = 2004"
    # start depends on catchment data
    mhm_nml[grep("eval_per(1)%ystart = ", mhm_nml, fixed = TRUE)] <- 
      paste0("    eval_per(1)%ystart = ", format(start_date, "%Y"))
    mhm_nml[grep("eval_per(1)%dstart = ", mhm_nml, fixed = TRUE)] <- 
      paste0("    eval_per(1)%dstart = ", as.integer(format(start_date, "%d"))) 
    mhm_nml[grep("eval_per(1)%mstart = ", mhm_nml, fixed = TRUE)] <- 
      paste0("    eval_per(1)%mstart = ", as.integer(format(start_date, "%m")))
    writeLines(mhm_nml, paste0(cpath, "/mhm.nml"))
  }
  
  # prepare model parameters
  parameters <- suppressWarnings(
    read_nml("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_master_namelists/mhm_parameter.nml"))
  parameter_df <- do.call(rbind, parameters$mhm_parameters)
  ind_of_num_paras <- which(parameter_df[, 4] == 1)
  num_para_bounds <- parameter_df[ind_of_num_paras, c(1, 2)]
  standard_parameters <- parameter_df[ind_of_num_paras, 3]
  
  # Run best model in training time and generate parameter maps ----------------------------
  if("numerics_tracker.csv" %in% list.files(result_folder)){
    results <- read.csv(paste0(result_folder, "/numerics_tracker.csv"))
  } else {
    results <- read.csv(paste0(result_folder, "/result_tracker.csv"))
  }
  training_results <- results[which.max(results$wmulti_loss), ]
  start_col_num_paras <- which(colnames(training_results) == "canopyInterceptionFactor")
  opt_numeric_parameters <- training_results[1, start_col_num_paras:ncol(training_results)]
  # numerics of TFs
  KSat_tf <- training_results$KSat
  FieldCap_tf <- training_results$FieldCap
  
  
  mpr_nml <- readLines("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_master_namelists/mpr_fso.nml")
  original_mpr_nml <- mpr_nml
  to_file_start <- grep("to_file(", mpr_nml, fixed = TRUE)
  mpr_nml[to_file_start] # KSat: 28, 29, FieldCap: 37, 38
  mpr_nml[to_file_start + 4] <- "                     .true., .true., .false., .false., .false., .false.,"
  mpr_nml[to_file_start + 5] <- "                     .false., .false., .false., .true., .true., .false.,"
  mpr_nml[to_file_start + 7] <- "                     .false., .false., .false., .false., .true., .true.,"
  mpr_nml[to_file_start + 8] <- "                     .true., .false., .false., .false., , .true., .false."
  mpr_nml[to_file_start + 9] <- "                     .false., .false., , , .true., .false., .false., , , , , ,"
  mpr_nml[to_file_start + 11] <- "                     , , , .true., .true., , .false., , , , .false., .false.,"
  
  writeLines(mpr_nml, 
             "/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_master_namelists/mpr_fso.nml")
  
  # change mhm_output.nml of all basins
  path <- paste0(run_folder, "/FSO_mHM_major_basins/config")
  basin_folders <- list.files(path, pattern = "sub_", full.names = TRUE)
  for(basin in basin_folders){
    outputs <- readLines(paste0(basin, "/mhm_outputs.nml"))
    output_option <- grep("timeStep_model_outputs = ", outputs, fixed = TRUE)
    outputs[output_option] <- "timeStep_model_outputs = -1"
    writeLines(outputs, paste0(basin, "/mhm_outputs.nml"))
  }
  
  
  
  # Training time
  # run model and compute losses
  mhm_results <- run_mhm_from_tf(KSat_tf = KSat_tf, FieldCap_tf = FieldCap_tf,
                                 numeric_parameters = opt_numeric_parameters,
                                 run_folder = run_folder)
  writeLines(original_mpr_nml, 
             "/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_master_namelists/mpr_fso.nml")
  
  # reset output to monthly variables
  for(basin in basin_folders){
    outputs <- readLines(paste0(basin, "/mhm_outputs.nml"))
    output_option <- grep("timeStep_model_outputs = ", outputs, fixed = TRUE)
    outputs[output_option] <- "timeStep_model_outputs = -2"
    writeLines(outputs, paste0(basin, "/mhm_outputs.nml"))
  }
  
  # save parameter fields
  dir.create("/home/lv71468/mfeigl/FSO_mHM/Results/parameter_maps", showWarnings = FALSE)
  dir.create(paste0("/home/lv71468/mfeigl/FSO_mHM/Results/parameter_maps/", run_name))
  basin_folders <- list.files(paste0(run_folder, "/FSO_mHM_major_basins/output"))
  for(basin in basin_folders){
    file.copy(paste0(run_folder, "/FSO_mHM_major_basins/output/", basin, "/mHM_parameters.nc"),
              paste0("/home/lv71468/mfeigl/FSO_mHM/Results/parameter_maps/", run_name, "/", basin, "_parameter.nc"),
              overwrite = TRUE)
  }
  # save fluxstates
  dir.create("/home/lv71468/mfeigl/FSO_mHM/Results/flux_states", showWarnings = FALSE)
  dir.create(paste0("/home/lv71468/mfeigl/FSO_mHM/Results/flux_states/", run_name))
  basin_folders <- list.files(paste0(run_folder, "/FSO_mHM_major_basins/output"))
  for(basin in basin_folders){
    file.copy(paste0(run_folder, "/FSO_mHM_major_basins/output/", basin, "/mHM_Fluxes_States.nc"),
              paste0("/home/lv71468/mfeigl/FSO_mHM/Results/flux_states/", run_name, "/", basin, "_fulx_states.nc"),
              overwrite = TRUE)
  }  
  # save restarts for L1_soilMoistFC
  dir.create("/home/lv71468/mfeigl/FSO_mHM/Results/restarts", showWarnings = FALSE)
  dir.create(paste0("/home/lv71468/mfeigl/FSO_mHM/Results/restarts/", run_name))
  basin_folders <- list.files(paste0(run_folder, "/FSO_mHM_major_basins/output"))
  for(basin in basin_folders){
    file.copy(paste0(run_folder, "/FSO_mHM_major_basins/output/", basin, "/mHM_restart_001.nc"),
              paste0("/home/lv71468/mfeigl/FSO_mHM/Results/restarts/", run_name, "/", basin, "_restart.nc"),
              overwrite = TRUE)
  } 
}



