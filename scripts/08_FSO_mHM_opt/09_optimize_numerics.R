# mHM FSO optimization
# FSO-mHM project
# Moritz Feigl, Nov 2020

# define run name to create environment automatically
run_name <- "FSO_best_num_improvement"
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
            paste0(run_folder, "/FSO_mHM_major_basins"), recursive=TRUE, overwrite = TRUE)
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

# get current best TFs and parameters
runs <- list.files("/home/lv71468/mfeigl/FSO_mHM/Results", full.names = TRUE)
results <- NULL
for(run in runs){
  try({run_results <- read.csv(paste0(run, "/result_tracker.csv"))
  results <- rbind(results,
                   run_results)
  }, silent = TRUE)
}
training_results <- results[which.max(results$wmulti_loss), ]
start_col_num_paras <- which(colnames(training_results) == "canopyInterceptionFactor")
opt_numeric_parameters <- training_results[1, start_col_num_paras:ncol(training_results)]

# define function to optimize numerical values
get_numerics <- function(string){
  numerics <- strsplit(string, c("[/\\^\\(\\)\\*\\+\\-]")) %>%
    unlist() %>%
    sapply(as.numeric)
  numerics[!is.na(numerics)]
}
# numerics of TFs
KSat_tf <- training_results$best_KSat
FieldCap_tf <- training_results$best_FieldCap
nums_KSat <- get_numerics(KSat_tf)
nums_FieldCap <- get_numerics(FieldCap_tf)
numerics <- c(nums_KSat, nums_FieldCap)
numeric_function <- c(rep("KSat", length(nums_KSat)),
                      rep("FieldCap", length(nums_FieldCap)))
nr_nums <- length(numerics)



numerics_objective_function <- function(optimization_numerics){
  # Make a counter for the optimization
  counter <- try(get("counter", envir = .GlobalEnv), silent = TRUE)
  if(class(counter) == "try-error") counter <- 0
  counter <- counter + 1
  assign("counter", counter, envir = .GlobalEnv)
  # numIter
  cat("\n****** Point nr.", counter, "******\n")
  
  # put in new numerics
  KSat_tf_tmp <- KSat_tf
  FieldCap_tf_tmp <- FieldCap_tf
  for(i in seq_along(optimization_numerics[1:nr_nums])){
    if(numeric_function[i] == "KSat"){
      KSat_tf_tmp <- sub(format(numerics[i], nsmall = 2), 
                         round(optimization_numerics[i], 3), 
                         KSat_tf_tmp, fixed = TRUE)
    } else {
      FieldCap_tf_tmp <- sub(format(numerics[i], nsmall = 2), 
                             round(optimization_numerics[i], 3), 
                             FieldCap_tf_tmp, fixed = TRUE)
    }
  }
  
  # run mhm basins
  mhm_results <- run_mhm_from_tf(KSat_tf = KSat_tf_tmp, FieldCap_tf = FieldCap_tf_tmp,
                                 numeric_parameters = optimization_numerics[-c(1:nr_nums)],
                                 run_folder = run_folder)
  
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
  
  
  # get numerics_tracker
  result_tracker <- try(get("numerics_tracker", envir = .GlobalEnv), silent = TRUE)
  
  # add new results to numerics_tracker
  nums_matrix <-  matrix(optimization_numerics, ncol = length(optimization_numerics))
  colnames(nums_matrix) <- c(paste0("KSat", 1:length(nums_KSat)), paste0("FieldCap", 1:length(nums_FieldCap)),
                             names(opt_numeric_parameters))
  
  assign("numerics_tracker", rbind(numerics_tracker,
                                   data.frame(KSat = KSat_tf_tmp,
                                              FieldCap = FieldCap_tf_tmp,
                                              wmulti_loss = wmulti_obj_loss,
                                              NSE = NSE,
                                              lNSE = lNSE,
                                              KGE = KGE,
                                              domain_multi_nse_loss,
                                              nums_matrix,
                                              stringsAsFactors = FALSE)),
         envir = .GlobalEnv)
  cat("\nNumerics optimization results:\n")
  cat("mean NSE:", NSE, "\n")
  cat("mean lNSE:", lNSE, "\n")
  cat("mean KGE:", KGE, "\n")
  cat("weighted NSE/KGE loss:", wmulti_obj_loss, "\n")
  
  write.csv(numerics_tracker,
            paste0(result_folder, "/numerics_tracker.csv"), row.names = FALSE, 
  )
  # return numeric optimization loss
  loss <- ifelse(is.na(wmulti_obj_loss), NA, wmulti_obj_loss)
  return(loss)
}

# Optimization start point
init_nums <- as.numeric(c(numerics, opt_numeric_parameters))
nums_matrix <-  matrix(init_nums, ncol = length(init_nums))
colnames(nums_matrix) <- c(paste0("KSat", 1:length(nums_KSat)), paste0("FieldCap", 1:length(nums_FieldCap)),
                           names(opt_numeric_parameters))
assign("numerics_tracker", data.frame(KSat = training_results$best_KSat,
                                      FieldCap = training_results$best_FieldCap,
                                      wmulti_loss = training_results$best_wmulti_loss,
                                      NSE = training_results$best_NSE,
                                      lNSE =  training_results$best_lNSE,
                                      KGE =  training_results$best_KGE,
                                      training_results[, grep("basin_", names(training_results), fixed = TRUE)],
                                      nums_matrix,
                                      stringsAsFactors = FALSE), envir = .GlobalEnv)

# allow +- 5% of the overall range for optimization of mhm numerics
numerics_bound_range <- (num_para_bounds[, 2] - num_para_bounds[, 1]) * 0.05
lower_bounds_mhm_numerics <- opt_numeric_parameters -  numerics_bound_range
upper_bounds_mhm_numerics <- opt_numeric_parameters +  numerics_bound_range
too_low <- lower_bounds_mhm_numerics < num_para_bounds[, 1]
lower_bounds_mhm_numerics[too_low] <- num_para_bounds[too_low, 1]
too_high <- upper_bounds_mhm_numerics > num_para_bounds[, 2]
upper_bounds_mhm_numerics[too_high] <- num_para_bounds[too_high, 2]

lower_bounds <- c(numerics-0.5, as.numeric(lower_bounds_mhm_numerics))
upper_bounds <- c(numerics+0.5, as.numeric(upper_bounds_mhm_numerics))

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


# initial compiling and function checking
#system("chmod 775 05_first_run_preparation.sh")
system("bash 05_first_run_preparation.sh")

iterations <- 100
pop_size <- 10
n_iter <- iterations/pop_size
GA_results <- GA::ga(type = "real-valued", fitness = numerics_objective_function,
                     lower = lower_bounds,
                     upper = upper_bounds,
                     popSize = pop_size,
                     maxiter = n_iter,
                     run = iterations,
                     keepBest = TRUE
)





# Run best model again in validation time ------------------------------------------------
results <- read.csv("/home/lv71468/mfeigl/FSO_mHM/Results/FSO_best_num_improvement/numerics_tracker.csv")
training_results <- results[which.max(results$wmulti_loss), ]
start_col_num_paras <- which(colnames(training_results) == "canopyInterceptionFactor")
opt_numeric_parameters <- training_results[1, start_col_num_paras:ncol(training_results)]
# numerics of TFs
KSat_tf <- training_results$KSat
FieldCap_tf <- training_results$FieldCap


# Change dates in nml
path <- paste0(run_folder, "/FSO_mHM_major_basins/config")
config_paths <- list.files(path, pattern = "sub_", full.names = TRUE)
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




mhm_results <- run_mhm_from_tf(KSat_tf = KSat_tf, FieldCap_tf = FieldCap_tf,
                               numeric_parameters = opt_numeric_parameters,
                               run_folder = run_folder)

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


# add new results to numerics_tracker
nums_matrix <-  matrix(optimization_numerics, ncol = length(optimization_numerics))
colnames(nums_matrix) <- c(paste0("KSat", 1:length(nums_KSat)), paste0("FieldCap", 1:length(nums_FieldCap)),
                           names(opt_numeric_parameters))


data.frame(KSat = KSat_tf_tmp,
           FieldCap = FieldCap_tf_tmp,
           wmulti_loss = wmulti_obj_loss,
           NSE = NSE,
           lNSE = lNSE,
           KGE = KGE,
           domain_multi_nse_loss,
           nums_matrix,
           stringsAsFactors = FALSE)


system(paste0("bash ", run_folder, "/mhm_run_scripts/run_basin_3.sh &> ", run_folder, "/mhm_run_scripts/run_basin_3.log"))# & ",
system(paste0("bash ", run_folder, "/mhm_run_scripts/run_basin_7.sh &> ", run_folder, "/mhm_run_scripts/run_basin_7.log"))# & ",

# to run: 1, 3, 5, 7
system(paste0("bash ", run_folder, "/mhm_run_scripts/run_basin_1.sh &> ", run_folder, "/mhm_run_scripts/run_basin_1.log & ",
              "bash ", run_folder, "/mhm_run_scripts/run_basin_2.sh &> ", run_folder, "/mhm_run_scripts/run_basin_2.log & ",
              "bash ", run_folder, "/mhm_run_scripts/run_basin_4.sh &> ", run_folder, "/mhm_run_scripts/run_basin_4.log & ",
              "bash ", run_folder, "/mhm_run_scripts/run_basin_5.sh &> ", run_folder, "/mhm_run_scripts/run_basin_5.log & ",
              "process_id=$! bash ", run_folder, "/mhm_run_scripts/run_basin_6.sh &> ", run_folder,
              "/mhm_run_scripts/run_basin_6.log & wait $process_id"
))

# change dates back to training time periods






