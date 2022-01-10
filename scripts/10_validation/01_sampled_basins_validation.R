# mHM FSO optimization
# FSO-mHM project
# Moritz Feigl, Nov 2020

# define run name to create environment automatically
run_name <- "small_basins_validation"
gcc <- ""
runs <- c(paste0("FSO_SCE_NSE_run_", 1:4), "FSO_DDS_NSE_run_5")
#runs <- c("mhm_standard_parameter")

for(run_to_validate_name in runs){
  cat("Start sampled validation of run", run_to_validate_name, "\n")
  # 0. Setup -------------------------------------------------------------------------------
  # Fodler structure
  dir.create("/home/lv71468/mfeigl/FSO_mHM/Runs", showWarnings = FALSE)
  run_folder <- paste0("/home/lv71468/mfeigl/FSO_mHM/Runs/", run_name)
  result_folder <- paste0("/home/lv71468/mfeigl/FSO_mHM/Results/", run_name)
  run_to_validate <- paste0("/home/lv71468/mfeigl/FSO_mHM/Results/", run_to_validate_name)
  dir.create(run_folder, showWarnings = FALSE)
  dir.create(result_folder, showWarnings = FALSE)
  setwd(run_folder)
  
  # mhm setup
  if(!dir.exists(paste0(run_folder, "/mhm"))){
    file.copy("/home/lv71468/mfeigl/FSO_mHM/mHM_setup/mhm", 
              run_folder, recursive=TRUE, overwrite = TRUE)
    unlink(paste0(run_folder, "/mhm/build"), recursive = TRUE) 
  }
  dir.create(paste0(run_folder, "/Results"), showWarnings = FALSE)
  
  # copy and adapt run_mhm_preparation script
  run_preparation_script_name <- paste0("01_compile_for_small_basins", gcc, ".sh")
  file.copy(paste0("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/10_validation/",
                   run_preparation_script_name), 
            paste0(run_folder, "/01_compile_for_small_basins.sh"), overwrite = TRUE)
  # Get last R codes
  code_path <- "/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/"
  file.copy(paste0(code_path, "04_FSO_mhm_utils.R"),
            "04_FSO_mhm_utils.R", overwrite = TRUE)
  # basin infos 
  basin_info <- read.table("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/10_validation/LUT_german_basins_with_properties_and_water_balance.txt", 
                           sep = ";", header = TRUE)
  basins <- basin_info$Stat_ID
  basins[grep(579085, basins)] <- "0579085"
  
  # sampled validation basins
  sampled_basins <- read.csv("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/10_validation/sampled_validation_basins.csv")
  path <- "/home/lv71468/mfeigl/FSO_mHM/mHM_setup/FSO_mHM_all_basins/config"
  all_basin_folders <- list.files(path, pattern = "sub_")
  all_basin_folders <- all_basin_folders[- grep("backup", all_basin_folders)]
  basin_folders <- all_basin_folders[apply(sampled_basins, 1, function(x) grep(x, all_basin_folders))]
  basins <- basins[basins %in% sampled_basins$sampled_validation_basins]
  
  # Create mhm run scripts in project folder
  dir.create("mhm_run_scripts", showWarnings = FALSE)
  run_script <- readLines("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/10_validation/01_run_script_template.sh")
  count <- 0
  for(basin in basins){
    count <- count + 1
    basin_script <- gsub("I=6335304", paste0("I=", basin), run_script, fixed = TRUE)
    writeLines(basin_script, paste0("mhm_run_scripts/sampled_run_basin_", count, ".sh"))
  }
  
  # adapt mhm paths
  config_paths <- list.files("/home/lv71468/mfeigl/FSO_mHM/mHM_setup/FSO_mHM_all_basins/config", 
                             pattern = "sub_",
                             full.names = TRUE)
  zink <- read.table("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/10_validation/lut_gauges_de.mpeg", 
                     header = TRUE)
  
  for(cpath in config_paths){
    # change path
    try({
      mhm_nml <- readLines(paste0(cpath, "/mhm.nml"))
      mhm_nml <- gsub("/work/ottor/FSO_mHM_major_basins", 
                      "/home/lv71468/mfeigl/FSO_mHM/mHM_setup/FSO_mHM_all_basins", mhm_nml,
                      fixed = TRUE)
      cp_id <- tail(strsplit(cpath, "_")[[1]], 1)
      zink_start_yr <- zink[as.character(zink$Stat_ID) == cp_id, "from_year"] + 5
      zink_end_yr <- zink[as.character(zink$Stat_ID) == cp_id, "to_year"]
      # change dates
      cat("\n", cp_id, ":")
      if(length(zink_start_yr) == 1){
        start_date <- as.POSIXct(paste0(zink_start_yr, "-01-01"))
        cat(as.character(start_date))
      } else {
        start_date <- as.POSIXct("1965-01-01")
      }
      
      end_year <- ifelse(length(zink_end_yr) == 1, 
                         zink_end_yr,
                         "1999")
      # end is always Dec. 31
      mhm_nml[grep("eval_per(1)%yend = ", mhm_nml, fixed = TRUE)] <- paste0("    eval_per(1)%yend = ", end_year)
      # start depends on catchment data
      mhm_nml[grep("eval_per(1)%ystart = ", mhm_nml, fixed = TRUE)] <- 
        paste0("    eval_per(1)%ystart = ", format(start_date, "%Y"))
      mhm_nml[grep("eval_per(1)%dstart = ", mhm_nml, fixed = TRUE)] <- 
        paste0("    eval_per(1)%dstart = ", as.integer(format(start_date, "%d"))) 
      mhm_nml[grep("eval_per(1)%mstart = ", mhm_nml, fixed = TRUE)] <- 
        paste0("    eval_per(1)%mstart = ", as.integer(format(start_date, "%m")))
      writeLines(mhm_nml, paste0(cpath, "/mhm.nml"))
    })
  }
  
  # utils
  source("04_FSO_mhm_utils.R")
  library(magrittr)
  
  # prepare model parameters
  parameters <- suppressWarnings(
    read_nml("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_master_namelists/mhm_parameter.nml"))
  parameter_df <- do.call(rbind, parameters$mhm_parameters)
  ind_of_num_paras <- which(parameter_df[, 4] == 1)
  num_para_bounds <- parameter_df[ind_of_num_paras, c(1, 2)]
  standard_parameters <- parameter_df[ind_of_num_paras, 3]
  
  # get current best TFs and parameters
  training_results <- read.csv(paste0(run_to_validate, "/training_results_overview.csv"))
  start_col_num_paras <- which(colnames(training_results) == "canopyInterceptionFactor")
  numeric_parameters <- training_results[1, start_col_num_paras:ncol(training_results)]
  
  if(run_to_validate_name != "mhm_standard_parameter"){
    
    
    # numerics of TFs
    KSat_tf <- training_results$KSat
    FieldCap_tf <- training_results$FieldCap
    
    # run mhm from tf for sampled validation basins basins
    variables <- c("bd", "sand", "clay", "slope", "aspect", "ThetaS", "KSat",
                   "vGenu_n", "dem")
    alltill_variables <- c("bd", "sand", "clay", "ThetaS", "KSat", "vGenu_n", "aspect", "dem", "slope")
    till_variables <- c("bd_till", "sand_till", "clay_till", "ThetaS_till", "KSat_till",
                        "vGenu_n_till", "slope_till", "aspect_till", "dem_till")
    notill_variables <- c("bd_notill", "sand_notill", "clay_notill", "ThetaS_notill",
                          "KSat_notill", "vGenu_n_notill", "slope_notill", "aspect_notill", "dem_notill")
    scaling_bounds <- list("slope" = c(0, 90),
                           "aspect" = c(0, 360),
                           "bd" = c(0, 2.3),
                           "sand" = c(0, 100),
                           "clay" = c(0, 100),
                           "dem" = c(0, 4000),
                           "KSat" = c(9, 244),
                           "vGenu_n" = c(1,2),
                           "ThetaS" = c(0.24, 0.51))
    
    # Add scaling to given transfer functions
    # Ksat
    KSat_from <- as.numeric(
      read.csv("/home/lv71468/mfeigl/FSO_mHM/FSO_Data/VAE_data/KSat/distribution_cutoff_points.csv"))
    KSat_to <- c(1.1, 1000.0)
    a <- (KSat_from[2] - KSat_from[1]) %>% format(nsmall = 1)
    b <- (KSat_to[2] - KSat_to[1]) %>% format(nsmall = 1)
    KSat_from[1] <- KSat_from[1] %>% format(nsmall = 1)
    KSat_to[1] <- KSat_to[1] %>% format(nsmall = 1)
    Ksat_fun <- paste0(KSat_to[1], "+((", KSat_tf, "-", KSat_from[1], ")", "*", b, ")/", a)
    Ksat_fun <- gsub("--", "+", Ksat_fun, fixed = TRUE)
    # FieldCap
    FieldCap_from <- as.numeric(
      read.csv("/home/lv71468/mfeigl/FSO_mHM/FSO_Data/VAE_data/FieldCap/distribution_cutoff_points.csv"))
    FieldCap_to <- c(0.01, 0.55)
    a <- (FieldCap_from[2] - FieldCap_from[1]) %>% format(nsmall = 1)
    b <- (FieldCap_to[2] - FieldCap_to[1]) %>% format(nsmall = 1)
    FieldCap_from[1] <- FieldCap_from[1] %>% format(nsmall = 1)
    FieldCap_to[1] <- FieldCap_to[1] %>% format(nsmall = 1)
    FieldCap_fun <- paste0(FieldCap_to[1], "+((", FieldCap_tf, "-", FieldCap_from[1], ")", "*", b, ")/", a)
    FieldCap_fun <- gsub("--", "+", FieldCap_fun, fixed = TRUE)
    # prepare tfs with scaling bounds
    prepared_KSat_tf <- prepare_tf_for_mhm(Ksat_fun, scaling_bounds = scaling_bounds)
    prepared_FieldCap_tf <- prepare_tf_for_mhm(FieldCap_fun, scaling_bounds = scaling_bounds)
    # Ksat
    Ksat_fun <- prepared_KSat_tf$scaled_function
    KSat_nums <- prepared_KSat_tf$numerics
    # FieldCap
    FieldCap_fun <- prepared_FieldCap_tf$scaled_function
    FieldCap_nums <- prepared_FieldCap_tf$numerics
    # make numeric vectors with parameter names
    numeric_paras <- c(KSat_nums, FieldCap_nums) %>% unique()
    numeric_paras_names <- paste0("FSO", seq_along(numeric_paras))
    numeric_paras <- numeric_paras[order(nchar(numeric_paras), decreasing = TRUE)]
    # fill in numerics in functions
    for(k in seq_along(numeric_paras)){
      Ksat_fun <- gsub(numeric_paras[k], numeric_paras_names[k], Ksat_fun)
      FieldCap_fun <- gsub(numeric_paras[k], numeric_paras_names[k], FieldCap_fun)
    }
    
    # change KSat variables to _horizon if necessary
    for(var in c("aspect", "slope", "dem")){
      Ksat_fun <- gsub(var, paste0(var, "_horizon"), Ksat_fun)
    }
    
    # change all variables into _till and _notill
    FieldCap_till <- FieldCap_fun
    for(var in alltill_variables){
      FieldCap_till <- gsub(var, paste0(var, "_till"), FieldCap_till)
    }
    FieldCap_notill <- FieldCap_fun
    for(var2 in alltill_variables){
      FieldCap_notill <- gsub(var2, paste0(var2, "_notill"), FieldCap_notill)
    }
    
    # 1. Change transfer functions in nml
    mpr_nml <- readLines("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_master_namelists/mpr_fso.nml")
    # change KSat tf
    KSat_line <- grep('KSat_dummy', mpr_nml)
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
    # data array FieldCap_notill
    fun_vars <- unlist(strsplit(gsub(" ", "", FieldCap_notill), "(?=[+-/*)()])", perl = TRUE))
    used_var <- unique(fun_vars[fun_vars %in% notill_variables])
    data_array_line <- grep("from_data_arrays(1:3,38)", mpr_nml, fixed = TRUE)
    mpr_nml[data_array_line] <- paste0('    from_data_arrays(1:', length(used_var), ',38',
                                       ') = ',
                                       paste0("'", paste0(used_var, collapse = "', '"), "'"),
                                       collapse = "")
    
    # get parameter names, parameter_values
    parameter_name_line <- grep('parameter_names', mpr_nml)
    parameter_value_line <- grep('parameter_values', mpr_nml)
    parameters_end <- grep('&data_arrays', mpr_nml, fixed = TRUE) - 3
    parameter_names_initial <- mpr_nml[parameter_name_line:(parameter_value_line-1)]
    parameter_values_initial <- mpr_nml[parameter_value_line:parameters_end]
    # add new parameter names and values
    nlines <- length(parameter_names_initial)
    parameter_names_initial[nlines] <- paste0(parameter_names_initial[nlines], ",")
    parameter_names <- c(
      parameter_names_initial,
      paste0("                            ",
             paste0(paste0("'", numeric_paras_names, "'"), collapse = ", ")))
    nlines <- length(parameter_values_initial)
    parameter_values_initial[nlines] <- paste0(parameter_values_initial[nlines], ",")
    parameter_values <- c(
      parameter_values_initial,
      paste0("                            ", paste0(numeric_paras, collapse = ", ")))
    old_nr_of_numeric_paras <- mpr_nml[parameter_name_line] %>%
      strsplit(split = "(", fixed = TRUE) %>%
      unlist() %>%
      substr(., 3, 4) %>% paste0(., collapse = "") %>%
      as.numeric()
    total_nr_of_numeric_paras <- old_nr_of_numeric_paras + length(numeric_paras)
    parameter_names[1] <- paste0("    parameter_names(1:",total_nr_of_numeric_paras, ") = ",
                                 strsplit(parameter_names[1], " = ", fixed = TRUE)[[1]][2])
    parameter_values[1] <- paste0("    parameter_values(1:",total_nr_of_numeric_paras, ") = ",
                                  strsplit(parameter_values[1], " = ", fixed = TRUE)[[1]][2])
    # add new parmeter names and values
    new_mpr_nml <- c(mpr_nml[1:(parameter_name_line-1)],
                     parameter_names,
                     parameter_values,
                     mpr_nml[(parameters_end+1):length(mpr_nml)])
    # Change numeric parameters
    parameters <- suppressWarnings(
      read_nml("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_master_namelists/mhm_parameter.nml"))
    for(i_numerics in seq_along(numeric_parameters)){
      i_para <- ind_of_num_paras[i_numerics]
      parameters$mhm_parameters[[i_para]][3] <- numeric_parameters[i_numerics]
    }
    
    # write new mpr.nml and mhm_parameter.nml for all basins
    for(basin in all_basin_folders){
      basin_nml <- gsub("<BASIN>", basin, new_mpr_nml)
      writeLines(basin_nml, con = paste0(path, "/", basin, "/mpr.nml"))
      write_mhm_nml(parameters, paste0(path, "/", basin, "/mhm_parameter.nml"))
    }
  } else {
    # Standard parameter run
    # Adapt mpr.nml
    mpr_nml <- readLines("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_master_namelists/mpr_fso.nml")

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
    KSat_tf <- Ksat_fun
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
    FieldCap_tf <- FieldCap_till
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
    for(basin in all_basin_folders){
      basin_nml <- gsub("<BASIN>", basin, mpr_nml)
      writeLines(basin_nml, con = paste0(path, "/", basin, "/mpr.nml"))
      write_mhm_nml(parameters, paste0(path, "/", basin, "/mhm_parameter.nml"))
    }
    
    
  }
  # 2. Run python script and compile mhm
  system("chmod 775 01_compile_for_small_basins.sh")
  system("bash 01_compile_for_small_basins.sh")
  # copy mhm to all basin folders
  for(cp in config_paths){
    file.copy(paste0(run_folder, "/mhm/build/mhm"), 
              paste0(cp, "/mhm"), overwrite = TRUE)
    file.copy(paste0(run_folder, "/mhm/build/mhm-0.5.9"), 
              paste0(cp, "/mhm-0.5.9"), overwrite = TRUE)
  }
  for(basin_nr in 1:length(basin_folders)){
    system(paste0("chmod 775 ", run_folder, "/mhm_run_scripts/sampled_run_basin_", basin_nr, ".sh"))
  }
  
  
  cat("running mhm\n")
  steps <- 5
  basin_sets <- floor(length(basins) / steps)
  rest_basins <- length(basins) %% steps
  set_start <- 1
  # Run mhm in set steps
  for(set in 1:basin_sets){
    sets <- set_start:(set_start + (steps-1))
    bash_cmd <- NULL
    for(set_id in sets[-steps]){
      bash_cmd <- c(bash_cmd, 
                    paste0("bash ", run_folder, "/mhm_run_scripts/sampled_run_basin_",set_id, ".sh &> ",
                           run_folder, "/mhm_run_scripts/sampled_run_basin_", set_id, ".log & ")
      )
    }
    bash_cmd <- c(bash_cmd, 
                  paste0("process_id=$! bash ", run_folder, "/mhm_run_scripts/sampled_run_basin_", 
                         sets[steps], ".sh &> ", run_folder, 
                         "/mhm_run_scripts/sampled_run_basin_", sets[steps], ".log & wait $process_id"))
    cat("Started computing basins:", paste0(sets, collpase = ", "), "\n")
    system(paste0(bash_cmd, collapse = ""))
    set_start <- set_start + steps
  }
  
  cat("Finished running mhm!\n")
  
  
  # read outputs and get KGE/NSE
  basin_results <- data.frame("Basin" = basin_folders, "KGE" = NA, "NSE" = NA, "lNSE" = NA)
  for(basin in basin_folders){
    try({
      output <- readLines(paste0(path, "/", basin, "/output.txt"))
      lines <- grep("KGE of daily discharge", output)
      KGE <- as.numeric(substring(output[lines], 39))
      NSE <- as.numeric(substring(output[lines + 1], 39))
      lNSE <- as.numeric(substring(output[lines + 2], 39))
      basin_results[basin_results$Basin == basin, c("KGE", "NSE", "lNSE")] <- c(KGE, NSE, lNSE)
    })
  }
  
  KGE <- mean(basin_results$KGE)
  NSE <- mean(basin_results$NSE)
  lNSE <- mean(basin_results$lNSE)
  multi_obj_nse_loss <- sign(basin_results$NSE) * 0.5*abs(basin_results$NSE)^6 +
    sign(basin_results$lNSE) * 0.5*abs(basin_results$lNSE)^6
  multi_obj_nse_loss <- ifelse(multi_obj_nse_loss < 0,
                               -1*(multi_obj_nse_loss*-1)^(1/6),
                               multi_obj_nse_loss^(1/6))
  median_loss <- median(multi_obj_nse_loss)
  domain_multi_nse_loss <- matrix(multi_obj_nse_loss, ncol = nrow(basin_results))
  
  sampled_results <- data.frame(KSat = KSat_tf,
                                FieldCap = FieldCap_tf,
                                median_loss = median_loss,
                                NSE = NSE,
                                lNSE = lNSE,
                                KGE = KGE,
                                domain_multi_nse_loss,
                                numeric_parameters,
                                stringsAsFactors = FALSE, check.names = FALSE)
  write.csv(sampled_results,
            paste0(run_to_validate, "/sampled_validation_basins_overview.csv"),
            row.names = FALSE)
  write.csv(basin_results,
            paste0(run_to_validate, "/sampled_validation_basins_metrics.csv"),
            row.names = FALSE)
  
}








