# mHM FSO optimization
# FSO-mHM project
# Moritz Feigl, Nov 2020

# define run name to create environment automatically
run_name <- "FSO_dds_run_1"
numIter <- 4
iterations_per_run <- 2
gcc <- ""
dir.create("/home/lv71468/mfeigl/FSO_mHM/Runs", showWarnings = FALSE)
run_folder <- paste0("/home/lv71468/mfeigl/FSO_mHM/Runs/", run_name)
result_folder <- paste0("/home/lv71468/mfeigl/FSO_mHM/Results/", run_name)
dir.create(run_folder, showWarnings = FALSE)
setwd(run_folder)
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
  if(gcc == "_gcc"){
    file.copy("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/CMakeLists_gcc.txt", 
              paste0(run_folder, "/mhm/CMakeLists.txt"), overwrite = TRUE)
  }
  # remove build folder
  unlink(paste0(run_folder, "/mhm/build"), recursive = TRUE) 
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
  file.copy("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/06_start_fso_mhm.sh",
            "06_start_fso_mhm.sh", overwrite = TRUE)
  start_fso_mhm <- readLines(paste0(run_folder, "/06_start_fso_mhm.sh"))
  if(gcc == "_gcc"){
    start_fso_mhm <- start_fso_mhm[-grep("#SBATCH -L intel@vsc", start_fso_mhm, fixed = TRUE)]
  }
  start_fso_mhm <- gsub("#SBATCH -J FSO_mhm", 
                        paste0("#SBATCH -J ", run_name), start_fso_mhm)
  start_fso_mhm <- gsub("cd /home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt", 
                        paste0("cd ", run_folder), start_fso_mhm)
  writeLines(start_fso_mhm, paste0(run_folder, "/06_start_fso_mhm.sh"))
  
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
  point_tf <- character(2)
  names(point_tf) <- c("KSat", "FieldCap")
  point_tf[1] <- tf_generator(point = point[1:6],
                              latent_dim = latent_dim,
                              generator = ksat_generator,
                              dictionary = dictionary_KSat,
                              variables = var_KSat)
  cat("KSat = ", point_tf[1], "\n")
  point_tf[2] <- tf_generator(point = point[7:12],
                              latent_dim = latent_dim,
                              generator = fieldcap_generator,
                              dictionary = dictionary_FieldCap,
                              variables = var_FieldCap)
  cat("FieldCap = ", point_tf[2], "\n")
  
  # sample tan/tanh
  point_tf <- sapply(point_tf, tan_switch)
  
  # check for only numeric functions
  numeric_test_tf <- gsub("(", "", point_tf, fixed = TRUE)
  numeric_test_tf <- gsub(")", "", numeric_test_tf, fixed = TRUE)
  numeric_test <- suppressWarnings(sapply(numeric_test_tf, as.numeric))
  if(sum(is.na(numeric_test)) != 2) point_tf <- c(NA, NA)
  if("" %in% point_tf | sum(is.na(point_tf)) > 0) {
    cat("No valid function found!\nResulting overall loss: NA\n")
    result_tracker <- try(get("result_tracker", envir = .GlobalEnv), silent = TRUE)
    if(class(result_tracker) != "try-error"){
      cat("\nThe best functions are:\n")
      cat("mean KGE:", result_tracker[nrow(result_tracker), "best_KGE"], "\n")
      cat("weighted KGE:", result_tracker[nrow(result_tracker), "best_wKGE"], "\n")
      cat("KSat = ", result_tracker[nrow(result_tracker), "best_KSat"], "\n")
      cat("FieldCap = ", result_tracker[nrow(result_tracker), "best_FieldCap"], "\n")
    }
    return(NA)
  }
  
  # Run mHM
  run_folder <- try(get("run_folder", envir = .GlobalEnv), silent = TRUE)
  mhm_results <- run_mhm_from_tf(KSat_tf = point_tf[1], FieldCap_tf = point_tf[2],
                                 numeric_parameters = point[13:length(point)],
                                 run_folder = run_folder)
  KGE <- mean(mhm_results$KGE)
  NSE <- mean(mhm_results$NSE)
  wKGE <- weighted.mean(x = mhm_results$KGE, w = 1 - mhm_results$KGE)
  wNSE <- weighted.mean(x = mhm_results$NSE, w = 1 - mhm_results$NSE)
  functions_splitted <- lapply(point_tf, .function_splitter)
  size_loss <- length(unlist(functions_splitted)) * 0.001/2
  loss <- wNSE - size_loss
  inter_loss <- ifelse(is.na(loss), -999, loss)
  domain_KGEs <- matrix(mhm_results$KGE, ncol = nrow(mhm_results))
  colnames(domain_KGEs) <- sapply(mhm_results$Basin, function(x) gsub("sub", "basin", x))
  current_point <- data.frame(matrix(point, ncol = length(point)))
  colnames(current_point) <- c(paste0("FSO", 1:12), row.names(num_para_bounds))
  
  # get result_tracker, or initialize
  result_tracker <- try(get("result_tracker", envir = .GlobalEnv), silent = TRUE)
  
  if(class(result_tracker) == "try-error"){
    result_tracker <- NULL
    old_min_KGE <- -9999
  } else {
    old_min_KGE <- result_tracker$best_loss[nrow(result_tracker)]
    if(!("FSO1" %in% colnames(result_tracker))) {
      current_point_filler <- matrix(NA, ncol = length(point))
      colnames(current_point_filler) <- colnames(current_point)
      result_tracker <- cbind(result_tracker, current_point_filler, 
                              stringsAsFactors = FALSE)
    }
  }
  
  if(inter_loss > old_min_KGE){
    # full loss is < old_min_KGE -> new best functions
    assign("result_tracker", rbind(result_tracker,
                                   data.frame(
                                     timestamp = as.character(Sys.time()),
                                     best_KSat = point_tf[1],
                                     best_FieldCap = point_tf[2],
                                     best_loss = round(loss, 3),
                                     best_KGE = round(KGE, 3),
                                     best_wKGE = round(wKGE, 3),
                                     best_size_loss = round(size_loss, 3),
                                     KSat = point_tf[1],
                                     FieldCap = point_tf[2],
                                     loss = loss,
                                     KGE = KGE,
                                     wKGE = wKGE,
                                     size_loss = size_loss,
                                     NSE = NSE,
                                     domain_KGEs,
                                     n_iteration_used = 1,
                                     n_iterations_since_improvement = 0,
                                     current_point,
                                     stringsAsFactors = FALSE, check.names = FALSE)),
           envir = .GlobalEnv)
  } else {
    # if full loss is > old_min_KGE -> old best stays
    last_line <- result_tracker[nrow(result_tracker), ]
    assign("result_tracker",
           rbind(result_tracker,
                 data.frame(
                   timestamp = as.character(Sys.time()),
                   best_KSat = last_line[,"best_KSat"],
                   best_FieldCap = last_line[, "best_FieldCap"],
                   best_loss = last_line[, "best_loss"],
                   best_KGE = last_line[, "best_KGE"],
                   best_wKGE = last_line[, "best_wKGE"],
                   best_size_loss = last_line[, "best_size_loss"],
                   KSat = point_tf[1],
                   FieldCap = point_tf[2],
                   loss = loss,
                   KGE =  KGE,
                   wKGE = wKGE,
                   size_loss = size_loss,
                   NSE =  NSE,
                   domain_KGEs,
                   n_iteration_used = 1,
                   n_iterations_since_improvement =  last_line[, "n_iterations_since_improvement"] + 1,
                   current_point,
                   stringsAsFactors = FALSE, check.names = FALSE)),
           envir = .GlobalEnv)
    
  }
  cat("\nGA optimization results:\n")
  cat("mean KGE:", KGE, "\n")
  cat("weighted KGE:", wKGE, "\n")
  cat("overall loss:", loss, "\n")
  
  result_tracker <- try(get("result_tracker", envir = .GlobalEnv))
  if(class(result_tracker) != "try-error"){
    cat("\nThe best functions are:\n")
    cat("mean KGE:", result_tracker[nrow(result_tracker), "best_KGE"], "\n")
    cat("weighted KGE:", result_tracker[nrow(result_tracker), "best_wKGE"], "\n")
    cat("KSat = ", result_tracker[nrow(result_tracker), "best_KSat"], "\n")
    cat("FieldCap = ", result_tracker[nrow(result_tracker), "best_FieldCap"], "\n")
    
    write.csv(result_tracker,
              paste0(result_folder, "/result_tracker.csv"), row.names = FALSE, 
    )
  }
  # write optimization position
  write.csv(data.frame("Iteration" = counter,
                       "KGE" = KGE,
                       "wKGE" = wKGE,
                       "loss" = loss,
                       "current_point" = matrix(point, ncol = length(point))),
            paste0(result_folder, "/current_runoptimization_results.csv"), 
            row.names = FALSE)
  return(loss = loss)
}


# GA optimization -----------------------------------------------------------------------

# Function space dds
lower_bounds <- c(rep(-5, 12), num_para_bounds[, 1])
upper_bounds <- c(rep(5, 12), num_para_bounds[, 2])
# start point with rnorm for VAE and standardvalue for paras
start_point <- c(rnorm(12), standard_parameters)

# get previously computed optimizer states if available
if(!is.null(result_folder) & file.exists(paste0(result_folder, "/result_tracker.csv"))){
  result_tracker <- read.csv(paste0(result_folder, "/result_tracker.csv"), check.names = FALSE)
  names(result_tracker)[24:94] <- c(paste0("FSO", 1:12), row.names(num_para_bounds))
  training_results <- result_tracker[which.max(result_tracker$wKGE), ]
  start_point <- training_results[1, 24:94]
  counter <- nrow(training_results)
} else {
  counter <- 0
  suppressWarnings(rm(result_tracker))
}
# first run preparation
system("chmod 775 05_first_run_preparation.sh")
if(counter == 0) system("bash 05_first_run_preparation.sh")
# Run optimization
if(numIter - counter < iterations_per_run) iterations_per_run <- numIter - counter

iterations <- iterations_per_run
pop_size <- 20
n_iter <- iterations/pop_size
GA_results <- GA::ga(type = "real-valued", fitness = objective_function,
                     lower = lower_bounds,
                     upper = upper_bounds,
                     popSize = pop_size,
                     maxiter = n_iter,
                     run = iterations,
                     suggestions = start_point,
                     keepBest = TRUE
)

if(counter < numIter) system("sbatch 06_start_fso_mhm.sh")




















