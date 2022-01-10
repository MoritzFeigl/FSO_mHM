# mHM standard runs
# FSO-mHM project
# Moritz Feigl, Apr 2021

# define run name to create environment automatically
cat("Start mhm standard parameter run", run_to_validate_name, "\n")
run_name <- "mhm_standard_parameter"
gcc <- ""
# 0. Setup -------------------------------------------------------------------------------
# Fodler structure
dir.create("/home/lv71468/mfeigl/FSO_mHM/Runs", showWarnings = FALSE)
run_folder <- paste0("/home/lv71468/mfeigl/FSO_mHM/Runs/", run_name)
result_folder <- paste0("/home/lv71468/mfeigl/FSO_mHM/Results/", run_name)
dir.create(run_folder, showWarnings = FALSE)
dir.create(result_folder, showWarnings = FALSE)
setwd(run_folder)

# mhm setup
# create run dir
if(!dir.exists(paste0(run_folder, "/FSO_mHM_major_basins"))){
  dir.create(paste0(run_folder, "/FSO_mHM_major_basins"))
  file.copy("/home/lv71468/mfeigl/FSO_mHM/mHM_setup/FSO_mHM_major_basins/config",
            paste0(run_folder, "/FSO_mHM_major_basins"), recursive=TRUE)
  file.copy("/home/lv71468/mfeigl/FSO_mHM/mHM_setup/FSO_mHM_major_basins/forcing",
            paste0(run_folder, "/FSO_mHM_major_basins"), recursive=TRUE)
  file.copy("/home/lv71468/mfeigl/FSO_mHM/mHM_setup/FSO_mHM_major_basins/output",
            paste0(run_folder, "/FSO_mHM_major_basins"), recursive=TRUE)
  file.copy("/home/lv71468/mfeigl/FSO_mHM/mHM_setup/FSO_mHM_major_basins/static",
            paste0(run_folder, "/FSO_mHM_major_basins"), recursive=TRUE)
  
  # adapt mhm paths
  config_paths <- list.files("FSO_mHM_major_basins/config", full.names = TRUE)
  for(cpath in config_paths){
    mhm_nml <- readLines(paste0(cpath, "/mhm.nml"))
    mhm_nml <- gsub("/home/lv71468/mfeigl/FSO_mHM/mHM_setup", run_folder, mhm_nml)
    writeLines(mhm_nml, paste0(cpath, "/mhm.nml"))
  }
  
  
}

if(!dir.exists(paste0(run_folder, "/mhm"))){
  file.copy("/home/lv71468/mfeigl/FSO_mHM/mHM_setup/mhm", 
            run_folder, recursive=TRUE, overwrite = TRUE)
  unlink(paste0(run_folder, "/mhm/build"), recursive = TRUE) 
  # Create mhm run scripts in project folder
  scripts <- list.files("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_scripts",
                        paste0("run_basin_[0-9]", gcc, ".sh"), full.names = TRUE)
  dir.create("mhm_run_scripts", showWarnings = FALSE)
  file.copy(scripts, "mhm_run_scripts")
  file.rename(list.files("mhm_run_scripts", full.names = TRUE),
              gsub("_gcc", "", list.files("mhm_run_scripts", full.names = TRUE)))
  all_scripts <- list.files("mhm_run_scripts",
                            "run_basin_[0-9].sh", full.names = TRUE)
  for(script in all_scripts){
    run_script <- suppressWarnings(readLines(script))
    run_script <- gsub("PROJECTPATH=/home/lv71468/mfeigl/FSO_mHM/mHM_setup",
                       paste0("PROJECTPATH=", run_folder),
                       run_script)
    writeLines(run_script, script)
  }
}

# copy and adapt first_run_preparation script
if(!dir.exists(paste0(run_folder, "/05_first_run_preparation"))){
  first_run_prep_script_name <- paste0("05_first_run_preparation", gcc, ".sh")
  file.copy(paste0("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/",
                   first_run_prep_script_name), run_folder, overwrite = TRUE)
  file.rename(list.files(pattern=first_run_prep_script_name), "05_first_run_preparation.sh")
  run_preparation <- readLines(paste0(run_folder, "/05_first_run_preparation.sh"))
  run_preparation <- gsub("PROJECTPATH=/home/lv71468/mfeigl/FSO_mHM/mHM_setup",
                          paste0("PROJECTPATH=", run_folder), run_preparation)
  writeLines(run_preparation, paste0(run_folder, "/05_first_run_preparation.sh"))
}
# Get last R codes
code_path <- "/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/"
file.copy(paste0(code_path, "03_vae_generators.R"),
          "03_vae_generators.R", overwrite = TRUE)
file.copy(paste0(code_path, "04_FSO_mhm_utils.R"),
          "04_FSO_mhm_utils.R", overwrite = TRUE)


# utils
source("04_FSO_mhm_utils.R")

# Change tfs and parameters to standard --------------------------------------------------
path <- paste0(run_folder, "/FSO_mHM_major_basins/config")
basin_folders <- list.files(path, pattern = "sub_")

# Adapt mpr.nml
mpr_nml <- readLines("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_master_namelists/mpr_fso.nml")
to_file_start <- grep("to_file(", mpr_nml, fixed = TRUE)
mpr_nml[to_file_start] # KSat: 28, 29, FieldCap: 37, 38
mpr_nml[to_file_start + 4] <- "                     .true., .true., .false., .false., .false., .false.,"
mpr_nml[to_file_start + 5] <- "                     .false., .false., .false., .true., .true., .false.,"
mpr_nml[to_file_start + 7] <- "                     .false., .false., .false., .false., .true., .true.,"
mpr_nml[to_file_start + 8] <- "                     .true., .false., .false., .false., , .true., .false."
mpr_nml[to_file_start + 9] <- "                     .false., .false., , , .true., .false., .false., , , , , ,"
mpr_nml[to_file_start + 11] <- "                     , , , .true., .true., , .false., , , , .false., .false.,"

# variables
variables <- c("bd", "sand", "clay", "slope", "aspect", "ThetaS", "KSat",
               "vGenu_n", "dem")
alltill_variables <- c("bd", "sand", "clay", "ThetaS", "KSat", "vGenu_n", "aspect", "dem", "slope")
till_variables <- c("bd_till", "sand_till", "clay_till", "ThetaS_till", "KSat_till",
                    "vGenu_n_till", "slope_till", "aspect_till", "dem_till")
notill_variables <- c("bd_notill", "sand_notill", "clay_notill", "ThetaS_notill",
                      "KSat_notill", "vGenu_n_notill", "slope_notill", "aspect_notill", "dem_notill")

# Change KSat
KSat_line <- grep('KSat_dummy', mpr_nml)
Ksat_fun <- "PTF_Ks_curveSlope * exp((PTF_Ks_constant + PTF_Ks_sand * sand - PTF_Ks_clay * clay) * log(Ks_c_base))"
mpr_nml[KSat_line] <- paste0("                           '", Ksat_fun, "',")
fun_vars <- unlist(strsplit(gsub(" ", "", Ksat_fun), "(?=[+-/*)()])", perl = TRUE))
horizon_vas <- c("aspect_horizon", "slope_horizon", "dem_horizon")
used_var <- unique(fun_vars[fun_vars %in% c(variables, horizon_vas)])
data_array_line <- grep("from_data_arrays(1:3,26)", mpr_nml, fixed = TRUE)
mpr_nml[data_array_line] <- paste0('    from_data_arrays(1:', length(used_var), ',26',
                                   ') = ',
                                   paste0("'", paste0(used_var, collapse = "', '"), "'"),
                                   collapse = "")
# change FieldCap tfs
FieldCap_till <- "ThetaS_till * exp(FieldCap_c1 * (FieldCap_c2 + log10(KSat_till)) * log(vGenu_n_till))"
FieldCap_notill <- "ThetaS_notill * exp(FieldCap_c1 * (FieldCap_c2 + log10(KSat_notill)) * log(vGenu_n_notill))"
FieldCap_till_line <- grep('FieldCap_till_dummy', mpr_nml)
FieldCap_notill_line <- grep('FieldCap_notill_dummy', mpr_nml)
mpr_nml[FieldCap_till_line] <- paste0("                           '", FieldCap_till, "',")
mpr_nml[FieldCap_notill_line] <- paste0("                           '", FieldCap_notill, "',")

# data array FieldCap_till
fun_vars <- unlist(strsplit(gsub(" ", "", FieldCap_till), "(?=[+-/*)()])", perl = TRUE))
used_var <- unique(fun_vars[fun_vars %in% till_variables])
data_array_line <- grep("from_data_arrays(1:3,37)", mpr_nml, fixed = TRUE)
mpr_nml[data_array_line] <- paste0('    from_data_arrays(1:', length(used_var), ',37',
                                   ') = ',
                                   paste0("'", paste0(used_var, collapse = "', '"), "'"),
                                   collapse = "")
# change FieldCap tfs
# data array FieldCap_notill
fun_vars <- unlist(strsplit(gsub(" ", "", FieldCap_notill), "(?=[+-/*)()])", perl = TRUE))
used_var <- unique(fun_vars[fun_vars %in% notill_variables])
data_array_line <- grep("from_data_arrays(1:3,38)", mpr_nml, fixed = TRUE)
mpr_nml[data_array_line] <- paste0('    from_data_arrays(1:', length(used_var), ',38',
                                   ') = ',
                                   paste0("'", paste0(used_var, collapse = "', '"), "'"),
                                   collapse = "")
# get standard numeric parameters
parameters <- suppressWarnings(
  read_nml("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_master_namelists/mhm_parameter.nml"))

# write new mpr.nml and mhm_parameter.nml for all basins
for(basin in basin_folders){
  basin_nml <- gsub("<BASIN>", basin, mpr_nml)
  writeLines(basin_nml, con = paste0(path, "/", basin, "/mpr.nml"))
  write_mhm_nml(parameters, paste0(path, "/", basin, "/mhm_parameter.nml"))
}


# change mhm_output.nml of all basins
path_conf <- paste0(run_folder, "/FSO_mHM_major_basins/config")
basin_conf_folders <- list.files(path_conf, pattern = "sub_", full.names = TRUE)
for(basin in basin_conf_folders){
  outputs <- readLines(paste0(basin, "/mhm_outputs.nml"))
  output_option <- grep("timeStep_model_outputs = ", outputs, fixed = TRUE)
  outputs[output_option] <- "timeStep_model_outputs = -1"
  writeLines(outputs, paste0(basin, "/mhm_outputs.nml"))
}


# Training time and parameters -----------------------------------------------------------
# Change dates to training dates
path_conf <- paste0(run_folder, "/FSO_mHM_major_basins/config")
config_paths <- list.files(path_conf, pattern = "sub_", full.names = TRUE)
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

system(paste0("chmod 775 ", run_folder, "/05_first_run_preparation.sh"))
system(paste0("bash ", run_folder, "/05_first_run_preparation.sh"))
for(basin_nr in 1:7){
  system(paste0("chmod 775 ", run_folder, "/mhm_run_scripts/run_basin_", basin_nr, ".sh"))
}
cat("running mhm\n")
system(paste0("bash ", run_folder, "/mhm_run_scripts/run_basin_1.sh &> ", run_folder, "/mhm_run_scripts/run_basin_1.log & ",
              "bash ", run_folder, "/mhm_run_scripts/run_basin_2.sh &> ", run_folder, "/mhm_run_scripts/run_basin_2.log & ",
              "bash ", run_folder, "/mhm_run_scripts/run_basin_3.sh &> ", run_folder, "/mhm_run_scripts/run_basin_3.log & ",
              "bash ", run_folder, "/mhm_run_scripts/run_basin_4.sh &> ", run_folder, "/mhm_run_scripts/run_basin_4.log & ",
              "bash ", run_folder, "/mhm_run_scripts/run_basin_5.sh &> ", run_folder, "/mhm_run_scripts/run_basin_5.log & ",
              "bash ", run_folder, "/mhm_run_scripts/run_basin_6.sh &> ", run_folder, "/mhm_run_scripts/run_basin_6.log & ",
              "process_id=$! bash ", run_folder, "/mhm_run_scripts/run_basin_7.sh &> ", run_folder,
              "/mhm_run_scripts/run_basin_7.log & wait $process_id"
))
cat("Finished running mhm!\n")
# read outputs and get KGE/NSE
mhm_results <- data.frame("Basin" = basin_folders, "KGE" = NA, "NSE" = NA, "lNSE" = NA)
for(basin in basin_folders){
  output <- readLines(paste0(path, "/", basin, "/output.txt"))
  lines <- grep("KGE of daily discharge", output)
  KGE <- as.numeric(substring(output[lines], 39))
  NSE <- as.numeric(substring(output[lines + 1], 39))
  lNSE <- as.numeric(substring(output[lines + 2], 39))
  try(mhm_results[mhm_results$Basin == basin, c("KGE", "NSE", "lNSE")] <- c(KGE, NSE, lNSE))
}

KGE <- mean(mhm_results$KGE)
NSE <- mean(mhm_results$NSE)
lNSE <- mean(mhm_results$lNSE)
multi_obj_nse_loss <- sign(mhm_results$NSE) * 0.5*abs(mhm_results$NSE)^6 +
  sign(mhm_results$lNSE) * 0.5*abs(mhm_results$lNSE)^6
multi_obj_nse_loss <- ifelse(multi_obj_nse_loss < 0,
                             -1*(multi_obj_nse_loss*-1)^(1/6),
                             multi_obj_nse_loss^(1/6))
wmulti_obj_loss <- weighted.mean(x = multi_obj_nse_loss, w = 1 - multi_obj_nse_loss)
domain_multi_nse_loss <- matrix(multi_obj_nse_loss, ncol = nrow(mhm_results))
colnames(domain_multi_nse_loss) <- sapply(mhm_results$Basin, function(x) gsub("sub", "basin", x))

# get numeric parameters
parameter_df <- do.call(rbind, parameters$mhm_parameters)
ind_of_num_paras <- which(parameter_df[, 4] == 1)
num_para_bounds <- parameter_df[ind_of_num_paras, c(1, 2)]
standard_parameters <- parameter_df[ind_of_num_paras, 3]

# save training results
nums_matrix <-  matrix(standard_parameters, ncol = length(standard_parameters))
colnames(nums_matrix) <- names(standard_parameters)
trainin_results_overview <- data.frame(KSat = Ksat_fun,
                                       FieldCap = FieldCap_till,
                                       wmulti_loss = wmulti_obj_loss,
                                       NSE = NSE,
                                       lNSE = lNSE,
                                       KGE = KGE,
                                       domain_multi_nse_loss,
                                       nums_matrix,
                                       stringsAsFactors = FALSE, check.names = FALSE)
write.csv(trainin_results_overview,
          paste0(result_folder, "/training_results_overview.csv"), row.names = FALSE, 
)
write.csv(mhm_results,
          paste0(result_folder, "/training_results_metrics.csv"), row.names = FALSE, 
)

# reset output to monthly variables
for(basin in config_paths){
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


# Validation time ------------------------------------------------------------------------
# Change dates in nml
for(cpath in config_paths){
  mhm_nml <- readLines(paste0(cpath, "/mhm.nml"))
  start_date <- as.POSIXct("1965-01-01") # period start 
  # end is always 31.12.1999
  mhm_nml[grep("eval_per(1)%yend = ", mhm_nml, fixed = TRUE)] <- 
    "    eval_per(1)%yend = 1999"
  # start depends on catchment data
  mhm_nml[grep("eval_per(1)%ystart = ", mhm_nml, fixed = TRUE)] <- 
    paste0("    eval_per(1)%ystart = ", format(start_date, "%Y"))
  mhm_nml[grep("eval_per(1)%dstart = ", mhm_nml, fixed = TRUE)] <- 
    paste0("    eval_per(1)%dstart = ", as.integer(format(start_date, "%d"))) 
  mhm_nml[grep("eval_per(1)%mstart = ", mhm_nml, fixed = TRUE)] <- 
    paste0("    eval_per(1)%mstart = ", as.integer(format(start_date, "%m")))
  writeLines(mhm_nml, paste0(cpath, "/mhm.nml"))
}

cat("running mhm\n")
system(paste0("bash ", run_folder, "/mhm_run_scripts/run_basin_1.sh &> ", run_folder, "/mhm_run_scripts/run_basin_1.log & ",
              "bash ", run_folder, "/mhm_run_scripts/run_basin_2.sh &> ", run_folder, "/mhm_run_scripts/run_basin_2.log & ",
              "bash ", run_folder, "/mhm_run_scripts/run_basin_3.sh &> ", run_folder, "/mhm_run_scripts/run_basin_3.log & ",
              "bash ", run_folder, "/mhm_run_scripts/run_basin_4.sh &> ", run_folder, "/mhm_run_scripts/run_basin_4.log & ",
              "bash ", run_folder, "/mhm_run_scripts/run_basin_5.sh &> ", run_folder, "/mhm_run_scripts/run_basin_5.log & ",
              "bash ", run_folder, "/mhm_run_scripts/run_basin_6.sh &> ", run_folder, "/mhm_run_scripts/run_basin_6.log & ",
              "process_id=$! bash ", run_folder, "/mhm_run_scripts/run_basin_7.sh &> ", run_folder,
              "/mhm_run_scripts/run_basin_7.log & wait $process_id"
))
cat("Finished running mhm!\n")


# get validation results -----------------------------------------------------------------
# read outputs and get KGE/NSE
mhm_results <- data.frame("Basin" = basin_folders, "KGE" = NA, "NSE" = NA, "lNSE" = NA)
for(basin in basin_folders){
  output <- readLines(paste0(path, "/", basin, "/output.txt"))
  lines <- grep("KGE of daily discharge", output)
  KGE <- as.numeric(substring(output[lines], 39))
  NSE <- as.numeric(substring(output[lines + 1], 39))
  lNSE <- as.numeric(substring(output[lines + 2], 39))
  try(mhm_results[mhm_results$Basin == basin, c("KGE", "NSE", "lNSE")] <- c(KGE, NSE, lNSE))
}


KGE <- mean(mhm_results$KGE)
NSE <- mean(mhm_results$NSE)
lNSE <- mean(mhm_results$lNSE)
multi_obj_nse_loss <- sign(mhm_results$NSE) * 0.5*abs(mhm_results$NSE)^6 +
  sign(mhm_results$lNSE) * 0.5*abs(mhm_results$lNSE)^6
multi_obj_nse_loss <- ifelse(multi_obj_nse_loss < 0,
                             -1*(multi_obj_nse_loss*-1)^(1/6),
                             multi_obj_nse_loss^(1/6))
wmulti_obj_loss <- weighted.mean(x = multi_obj_nse_loss, w = 1 - multi_obj_nse_loss)
domain_multi_nse_loss <- matrix(multi_obj_nse_loss, ncol = nrow(mhm_results))
colnames(domain_multi_nse_loss) <- sapply(mhm_results$Basin, function(x) gsub("sub", "basin", x))
# get numeric parameters
parameter_df <- do.call(rbind, parameters$mhm_parameters)
ind_of_num_paras <- which(parameter_df[, 4] == 1)
num_para_bounds <- parameter_df[ind_of_num_paras, c(1, 2)]
standard_parameters <- parameter_df[ind_of_num_paras, 3]

# save validation results
nums_matrix <-  matrix(standard_parameters, ncol = length(standard_parameters))
colnames(nums_matrix) <- names(standard_parameters)
validation_results <- data.frame(KSat = Ksat_fun,
                                 FieldCap = FieldCap_till,
                                 wmulti_loss = wmulti_obj_loss,
                                 NSE = NSE,
                                 lNSE = lNSE,
                                 KGE = KGE,
                                 domain_multi_nse_loss,
                                 nums_matrix,
                                 stringsAsFactors = FALSE, check.names = FALSE)
write.csv(validation_results,
          paste0(result_folder, "/validation_resultions_overview.csv"), row.names = FALSE, 
)
write.csv(mhm_results,
          paste0(result_folder, "/validation_resultions_metrics.csv"), row.names = FALSE, 
)

