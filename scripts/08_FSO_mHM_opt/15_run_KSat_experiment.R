# mHM FSO optimization
# FSO-mHM project
# Moritz Feigl, Nov 2020

# define run name to create environment automatically
run_name <- "KSat_experiment"
numIter <- 3000
iterations_per_run <- NA
gcc <- ""
dir.create("/home/lv71468/mfeigl/FSO_mHM/Runs", showWarnings = FALSE)
run_folder <- paste0("/home/lv71468/mfeigl/FSO_mHM/Runs/", run_name)
result_folder <- paste0("/home/lv71468/mfeigl/FSO_mHM/Results/", run_name)
dir.create(run_folder, showWarnings = FALSE)
setwd(run_folder)
.libPaths("/home/lv71468/mfeigl/R/x86_64-pc-linux-gnu-library/4.0")
library(magrittr)
library(ncdf4)
library(hydroGOF)
if(!dir.exists(paste0(run_folder, "/FSO_mHM_major_basins"))){
  # create run dir
  dir.create(paste0(run_folder, "/FSO_mHM_major_basins"))
  file.copy("/home/lv71468/mfeigl/FSO_mHM/mHM_setup/FSO_mHM_major_basins/config",
            paste0(run_folder, "/FSO_mHM_major_basins"), recursive=TRUE)
  file.copy("/home/lv71468/mfeigl/FSO_mHM/mHM_setup/FSO_mHM_major_basins/forcing",
            paste0(run_folder, "/FSO_mHM_major_basins"), recursive=TRUE)
  file.copy("/home/lv71468/mfeigl/FSO_mHM/mHM_setup/FSO_mHM_major_basins/output",
            paste0(run_folder, "/FSO_mHM_major_basins"), recursive=TRUE)
  file.copy("/home/lv71468/mfeigl/FSO_mHM/mHM_setup/FSO_mHM_major_basins/static",
            paste0(run_folder, "/FSO_mHM_major_basins"), recursive=TRUE)
  # mhm
  file.copy("/home/lv71468/mfeigl/FSO_mHM/mHM_setup/mhm",
            run_folder, recursive=TRUE, overwrite = TRUE)
  
  # remove build folder
  try(unlink(paste0(run_folder, "/mhm/build"), recursive = TRUE))
  dir.create(paste0(run_folder, "/Results"))
  
  # adapt mhm paths
  config_paths <- list.files("FSO_mHM_major_basins/config", full.names = TRUE)
  for(cpath in config_paths){
    mhm_nml <- readLines(paste0(cpath, "/mhm.nml"))
    mhm_nml <- gsub("/home/lv71468/mfeigl/FSO_mHM/mHM_setup", run_folder, mhm_nml)
    writeLines(mhm_nml, paste0(cpath, "/mhm.nml"))
  }
  
  # copy and adapt run_mhm_preparation script
  run_preparation_script_name <- paste0("05_run_mhm_preparation", gcc, ".sh")
  file.copy(paste0("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/",
                   run_preparation_script_name), run_folder, overwrite = TRUE)
  file.rename(list.files(pattern=run_preparation_script_name), "05_run_mhm_preparation.sh")
  
  run_preparation <- readLines(paste0(run_folder, "/05_run_mhm_preparation.sh"))
  run_preparation <- gsub("PROJECTPATH=/home/lv71468/mfeigl/FSO_mHM/mHM_setup",
                          paste0("PROJECTPATH=", run_folder), run_preparation)
  writeLines(run_preparation, paste0(run_folder, "/05_run_mhm_preparation.sh"))
  
  # copy and adapt first_run_preparation script
  first_run_prep_script_name <- paste0("05_first_run_preparation", gcc, ".sh")
  file.copy(paste0("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/",
                   first_run_prep_script_name), run_folder, overwrite = TRUE)
  file.rename(list.files(pattern=first_run_prep_script_name), "05_first_run_preparation.sh")
  run_preparation <- readLines(paste0(run_folder, "/05_first_run_preparation.sh"))
  run_preparation <- gsub("PROJECTPATH=/home/lv71468/mfeigl/FSO_mHM/mHM_setup",
                          paste0("PROJECTPATH=", run_folder), run_preparation)
  writeLines(run_preparation, paste0(run_folder, "/05_first_run_preparation.sh"))
  
  # get sbatch script
  file.copy("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/15_run_KSat_experiment.sh",
            "15_run_KSat_experiment.sh", overwrite = TRUE)
  start_fso_mhm <- readLines(paste0(run_folder, "/15_run_KSat_experiment.sh"))
  if(gcc == "_gcc"){
    start_fso_mhm <- start_fso_mhm[-grep("#SBATCH -L intel@vsc", start_fso_mhm, fixed = TRUE)]
  }
  start_fso_mhm <- gsub("cd /home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt",
                        paste0("cd ", run_folder), start_fso_mhm)
  writeLines(start_fso_mhm, paste0(run_folder, "/15_run_KSat_experiment.sh"))
  file.copy("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/15_run_KSat_experiment.R",
            "15_run_KSat_experiment.R", overwrite = TRUE)
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
  
  # Results folder
  dir.create(result_folder, showWarnings = FALSE)
  
  # Get last R codes
  code_path <- "/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/"
  file.copy(paste0(code_path, "03_vae_generators.R"),
            "03_vae_generators.R", overwrite = TRUE)
  file.copy(paste0(code_path, "04_FSO_mhm_utils.R"),
            "04_FSO_mhm_utils.R", overwrite = TRUE)
  
  
  dir.create("standard_parameter_discharge", showWarnings = FALSE)
  # adapt mhm paths
  stand_output_paths <- list.files("/home/lv71468/mfeigl/FSO_mHM/Runs/mhm_standard_parameter/FSO_mHM_major_basins/output", 
                                   full.names = TRUE)
  
  #install.packages("ncdf4", 
  #                 configure.args = "--with-nc-config=/opt/sw/x86_64/glibc-2.17/ivybridge-ep/netcdf_C/4.4.1.1/intel/18/intel-mpi/2018/hdf5/1.8.12-MPI/pnetcdf/1.10.0/bin/nc-config --with-netcdf-include=/opt/sw/x86_64/glibc-2.17/ivybridge-ep/netcdf_C/4.4.1.1/intel/18/intel-mpi/2018/hdf5/1.8.12-MPI/pnetcdf/1.10.0/include/ --with-netcdf-lib=/opt/sw/x86_64/glibc-2.17/ivybridge-ep/netcdf_C/4.4.1.1/intel/18/intel-mpi/2018/hdf5/1.8.12-MPI/pnetcdf/1.10.0/lib/")
  
  standard_q <- NULL
  standard_ksat <- list()
  for(sopath in stand_output_paths){
    basin_q <- read.table(paste0(sopath, "/daily_discharge.out"), header = TRUE)
    basin_q <- basin_q[, -c(1, 5)]
    colnames(basin_q)[4] <- colnames(basin_q)[4] %>% 
      gsub("Qsim_", "", ., fixed = TRUE) %>% 
      as.integer %>% 
      paste0("sub_", .)
    if(is.null(standard_q)) {
      standard_q <- basin_q
    } else {
      standard_q <- cbind(standard_q, basin_q[, 4])
      names(standard_q)[ncol(standard_q)] <- colnames(basin_q)[4]
    }
    
    # ksat from netcdf
    nc_data <- nc_open(paste0(sopath, "/mHM_parameters.nc"))
    ksat_nc <- ncvar_get(nc_data, "KSat_notill")
    standard_ksat[[colnames(basin_q)[4]]] <- ksat_nc
  }
  saveRDS(standard_q, "standard_discharge.rds")
  saveRDS(standard_ksat, "standard_ksat.rds")
}


# utils
source("04_FSO_mhm_utils.R")
# Load VAE generators
gen_try <- try(source("03_vae_generators.R"), silent = TRUE)
if(class(gen_try) == "try-error") source("03_vae_generators.R")

# prepare model parameters
parameters <- suppressWarnings(
  read_nml("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_master_namelists/mhm_parameter.nml"))
parameter_df <- do.call(rbind, parameters$mhm_parameters)
ind_of_num_paras <- which(parameter_df[, 4] == 1)
num_para_bounds <- parameter_df[ind_of_num_paras, c(1, 2)]
standard_parameters <- parameter_df[ind_of_num_paras, 3]

# Objective function
objective_function <- function(point){
  latent_dim <- 6
  # Make a counter for the optimization
  counter <- try(get("counter", envir = .GlobalEnv), silent = TRUE)
  if(class(counter) == "try-error") counter <- 0
  counter <- counter + 1
  assign("counter", counter, envir = .GlobalEnv)
  # numIter
  numIter <- try(get("numIter", envir = .GlobalEnv), silent = TRUE)
  if(class(numIter) == "try-error") numIter <- -999
  
  cat("\n****** Point nr.", paste0(counter, "/", numIter), "******\n")
  # calculate the transfer function of the point
  point_tf <- character(1)
  names(point_tf) <- "KSat"
  point_tf[1] <- tf_generator(point = point[1:6],
                              latent_dim = latent_dim,
                              generator = ksat_generator,
                              dictionary = dictionary_KSat,
                              variables = var_KSat)
  cat("KSat = ", point_tf[1], "\n")
  
  # sample tan/tanh
  point_tf <- sapply(point_tf, tan_switch)
  
  # check for only numeric functions
  numeric_test_tf <- gsub("(", "", point_tf, fixed = TRUE)
  numeric_test_tf <- gsub(")", "", numeric_test_tf, fixed = TRUE)
  numeric_test <- suppressWarnings(sapply(numeric_test_tf, as.numeric))
  if(sum(is.na(numeric_test)) != 1) point_tf <- NA
  if("" %in% point_tf | sum(is.na(point_tf)) > 0) {
    cat("No valid function found!\nResulting overall loss: NA\n")
    result_tracker <- try(get("result_tracker", envir = .GlobalEnv), silent = TRUE)
    if(class(result_tracker) != "try-error"){
      cat("\nThe best functions are:\n")
      cat("mean NSE:", result_tracker[nrow(result_tracker), "best_NSE"], "\n")
      cat("mean lNSE:", result_tracker[nrow(result_tracker), "best_lNSE"], "\n")
      cat("mean KGE:", result_tracker[nrow(result_tracker), "best_KGE"], "\n")
      cat("mean spaef:", result_tracker[nrow(result_tracker), "best_spaef"], "\n")
      cat("weighted mean multi-obj loss:", result_tracker[nrow(result_tracker), "best_wmulti_loss"], "\n")
      cat("KSat = ", result_tracker[nrow(result_tracker), "best_KSat"], "\n")
    }
    
    assign("counter", counter, envir = .GlobalEnv)
    return(NA)
  }
  
  # Run mHM
  run_folder <- try(get("run_folder", envir = .GlobalEnv), silent = TRUE)
  mhm_results <- run_mhm_from_ksat(KSat_tf = point_tf[1],
                                   run_folder = run_folder)
  if(sum(is.na(mhm_results)) == 0){
    # get results
    out_paths <- list.files("FSO_mHM_major_basins/output", full.names = TRUE)
    sim_q <- NULL
    sim_ksat <- list()
    for(sopath in out_paths){
      basin_q <- read.table(paste0(sopath, "/daily_discharge.out"), header = TRUE)
      basin_q <- basin_q[, -c(1, 5)]
      colnames(basin_q)[4] <- colnames(basin_q)[4] %>% 
        gsub("Qsim_", "", ., fixed = TRUE) %>% 
        as.integer %>% 
        paste0("sub_", .)
      if(is.null(sim_q)) {
        sim_q <- basin_q
      } else {
        sim_q <- cbind(sim_q, basin_q[, 4])
        names(sim_q)[ncol(sim_q)] <- colnames(basin_q)[4]
      }
      
      # ksat from netcdf
      try({
        nc_data <- nc_open(paste0(sopath, "/mHM_parameters.nc"))
        ksat_nc <- ncvar_get(nc_data, "KSat_notill")
        sim_ksat[[colnames(basin_q)[4]]] <- ksat_nc
      })
    }
    standard_q <- get("standard_q")
    standard_ksat <- get("standard_ksat")
    basin_names <- names(sim_q)[-c(1:3)]
    mhm_results <- data.frame("Basin" = basin_names, "KGE" = NA, "NSE" = NA, "lNSE" = NA, "SPAEF" = NA)
    for(basin in basin_names){
      # KGE, NSE, log NSE from simulated and synthetic discharge series
      KGE <- hydroGOF::KGE(sim_q[, basin], standard_q[, basin])
      NSE <- hydroGOF::NSE(sim_q[, basin], standard_q[, basin])
      lNSE <- hydroGOF::NSE(sim_q[, basin], standard_q[, basin], FUN = "log")
      spaef_layers <- numeric(3)
      for(layer in 1:3){
        # get ksat values of specific soil layer as vector
        sim_lay <- as.numeric(sim_ksat[[basin]][, , layer])
        stand_lay <- as.numeric(standard_ksat[[basin]][, , layer])
        # remove values that are NA in the simulated output
        stand_lay <- stand_lay[!is.na(sim_lay)]
        sim_lay <- sim_lay[!is.na(sim_lay)]
        # compute SPAEFF
        spaef_layers[layer] <- SPAEF(stand_lay, sim_lay)
      }
      spaef <- mean(spaef_layers)
      mhm_results[mhm_results$Basin == basin, c("KGE", "NSE", "lNSE", "SPAEF")] <- 
        c(KGE, NSE, lNSE, spaef)
    }
    
    KGE <- mean(mhm_results$KGE)
    NSE <- mean(mhm_results$NSE)
    lNSE <- mean(mhm_results$lNSE)
    spaef <- mean(mhm_results$SPAEF)
    functions_splitted <- lapply(point_tf, .function_splitter)
    size_loss <- length(unlist(functions_splitted)) * 0.001/2
    
    multi_obj_nse_loss <- sign(mhm_results$NSE) * 1/3*abs(mhm_results$NSE)^6 +
      sign(mhm_results$lNSE) * 1/3*abs(mhm_results$lNSE)^6 +
      sign(mhm_results$SPAEF) * 1/3*abs(mhm_results$SPAEF)^6
    multi_obj_nse_loss <- ifelse(multi_obj_nse_loss < 0,
                                 -1*(multi_obj_nse_loss*-1)^(1/6),
                                 multi_obj_nse_loss^(1/6))
    wmulti_obj_loss <- weighted.mean(x = multi_obj_nse_loss, w = 1 - multi_obj_nse_loss)
    loss <- wmulti_obj_loss - size_loss
  } else {
    KGE <- NSE <- lNSE <- spaef <- wmulti_obj_loss <- loss <- size_loss <- NA
  }
  
  
  inter_loss <- ifelse(is.na(loss), -1e12, loss)
  domain_multi_nse_loss <- matrix(wmulti_obj_loss, ncol = nrow(mhm_results))
  colnames(domain_multi_nse_loss) <- sapply(mhm_results$Basin, function(x) gsub("sub", "basin", x))
  current_point <- data.frame(matrix(point, ncol = length(point)))
  colnames(current_point) <- paste0("FSO", 1:6)
  
  # get result_tracker, or initialize
  result_tracker <- try(get("result_tracker", envir = .GlobalEnv), silent = TRUE)
  
  if(class(result_tracker) == "try-error"){
    result_tracker <- NULL
    old_min_loss <- -1e18
  } else {
    old_min_loss <- result_tracker$best_loss[nrow(result_tracker)]
    old_min_loss <- ifelse(is.na(old_min_loss) | length(old_min_loss) == 0, -1e18, old_min_loss)
  }
  
  if(inter_loss > old_min_loss){
    # full loss is < old_min_loss -> new best functions
    assign("result_tracker", rbind(result_tracker,
                                   data.frame(
                                     timestamp = as.character(Sys.time()),
                                     best_KSat = point_tf[1],
                                     best_FieldCap = "standard",
                                     best_loss = round(loss, 3),
                                     best_wmulti_loss = round(wmulti_obj_loss, 3),
                                     best_size_loss = round(size_loss, 3),
                                     best_NSE = round(NSE, 3),
                                     best_KGE = round(KGE, 3),
                                     best_lNSE = round(lNSE, 3),
                                     best_spaef = round(spaef, 3),
                                     KSat = point_tf[1],
                                     FieldCap = "standard",
                                     loss = loss,
                                     wmulti_loss = round(wmulti_obj_loss, 3),
                                     NSE = NSE,
                                     lNSE = lNSE,
                                     size_loss = size_loss,
                                     KGE = KGE,
                                     spaef = spaef,
                                     domain_multi_nse_loss,
                                     n_iteration_used = 1,
                                     n_iterations_since_improvement = 0,
                                     current_point,
                                     stringsAsFactors = FALSE, check.names = FALSE)),
           envir = .GlobalEnv)
  } else {
    # if full loss is > old_min_loss -> old best stays
    last_line <- result_tracker[nrow(result_tracker), ]
    assign("result_tracker",
           rbind(result_tracker,
                 data.frame(
                   timestamp = as.character(Sys.time()),
                   best_KSat = last_line[,"best_KSat"],
                   best_FieldCap = "standard",
                   best_loss = last_line[, "best_loss"],
                   best_wmulti_loss = last_line[, "best_wmulti_loss"],
                   best_size_loss =  last_line[, "best_size_loss"],
                   best_NSE =  last_line[, "best_NSE"],
                   best_KGE =  last_line[, "best_KGE"],
                   best_lNSE = last_line[, "best_lNSE"],
                   best_spaef = last_line[, "best_spaef"],
                   KSat = point_tf[1],
                   FieldCap = "standard",
                   loss = loss,
                   wmulti_loss = round(wmulti_obj_loss, 3),
                   NSE = NSE,
                   lNSE = lNSE,
                   size_loss = size_loss,
                   KGE = KGE,
                   spaef = spaef,
                   domain_multi_nse_loss,
                   n_iteration_used = 1,
                   n_iterations_since_improvement =  last_line[, "n_iterations_since_improvement"] + 1,
                   current_point,
                   stringsAsFactors = FALSE, check.names = FALSE)),
           envir = .GlobalEnv)
    
  }
  cat("\nSCE optimization results:\n")
  cat("mean NSE:", NSE, "\n")
  cat("mean lNSE:", lNSE, "\n")
  cat("mean KGE:", KGE, "\n")
  cat("mean spaef:", spaef, "\n")
  cat("weighted mean multi-obj loss:", wmulti_obj_loss, "\n")
  cat("overall loss:", loss, "\n")
  
  result_tracker <- try(get("result_tracker", envir = .GlobalEnv))
  if(class(result_tracker) != "try-error"){
    cat("\nThe best functions are:\n")
    cat("mean NSE:", result_tracker[nrow(result_tracker), "best_NSE"], "\n")
    cat("mean lNSE:", result_tracker[nrow(result_tracker), "best_lNSE"], "\n")
    cat("mean KGE:", result_tracker[nrow(result_tracker), "best_KGE"], "\n")
    cat("mean spaef:", result_tracker[nrow(result_tracker), "best_spaef"], "\n")
    cat("weighted mean multi-obj loss:", result_tracker[nrow(result_tracker), "best_wmulti_loss"], "\n")
    cat("overall loss:", result_tracker[nrow(result_tracker), "best_loss"], "\n")
    cat("KSat = ", result_tracker[nrow(result_tracker), "best_KSat"], "\n")
    
    write.csv(result_tracker,
              paste0(result_folder, "/result_tracker.csv"), row.names = FALSE,
    )
  }
  return(loss)
}


# SCE optimization -----------------------------------------------------------------------

# Function space dds
xBounds <- data.frame(lower = rep(-5, 6),
                      upper = rep(5, 6))
# start point with rnorm for VAE and standardvalue for paras
start_point <- rnorm(6)

# get previously computed optimizer states if available
if(!is.null(result_folder) & file.exists(paste0(result_folder, "/optimization_state.rds"))){
  opt_state <- readRDS(paste0(result_folder, "/optimization_state.rds"))
  counter <- readRDS(paste0(result_folder, "/counter.rds"))
  result_tracker <- read.csv(paste0(result_folder, "/result_tracker.csv"), check.names = FALSE)
  names(result_tracker)[29:34] <- paste0("FSO", 1:6)
} else {
  counter <- 0
  suppressWarnings(rm(result_tracker))
}
# first run preparation
system("chmod 775 05_first_run_preparation.sh")
if(counter == 0) system("bash 05_first_run_preparation.sh")

# load synthetic KSat and discharge
standard_q <- readRDS("standard_discharge.rds")
standard_ksat <- readRDS("standard_ksat.rds")



sceDefaults <- function(){
  list(ncomplex = 2,           ## number of complexes
       cce.iter = 5,          ## number of iteration in inner loop (CCE algorithm)
       fnscale = -1,            ## function scaling factor (set to -1 for maximisation)
       elitism = 1,            ## controls amount of weighting in sampling towards the better parameter sets
       initsample = "latin",   ## sampling scheme for initial values -- "latin" or "random"
       reltol = 1e-5,          ## convergence threshold: relative improvement factor required in an SCE iteration
       tolsteps = 500,           ## number of iterations within reltol to confirm convergence
       maxit = 3000,          ## maximum number of iterations
       maxeval = numIter,          ## maximum number of function evaluations
       maxtime = Inf,          ## maximum duration of optimization in seconds
       returnpop = FALSE,      ## whether to return populations from all iterations
       trace = 3,              ## level of user feedback
       FSO_dims = 6,
       REPORT = 1) ## number of iterations between reports when trace >= 1
}
if(numIter - counter > 0){
  
  obj <- SCEoptim(objective_function, start_point, lower = xBounds$lower, upper = xBounds$upper,
                  control = list(trace = 3,
                                 ncomplex = 2,
                                 fnscale = -1,
                                 maxeval = numIter,
                                 cce.iter = 5),
                  state_folder = result_folder)
  saveRDS(counter, paste0(result_folder, "/counter.rds"))
  
}
if(counter < numIter) system("sbatch 15_run_KSat_experiment.sh")



