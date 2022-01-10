# mHM FSO optimization
# FSO-mHM project
# Moritz Feigl, Nov 2020

# define run name to create environment automatically
run_name <- "FSO_SCE_NSE_run_1"
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

# get current best TFs and parameters
results <- read.csv(paste0(result_folder, "/result_tracker.csv"))
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
# Loss functions
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
  # compute loss
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
  if(length(nums_FieldCap) == 0){
    colnames(nums_matrix) <- c(paste0("KSat", 1:length(nums_KSat)), 
                               names(opt_numeric_parameters))
  }
  if(length(nums_KSat) == 0){
    colnames(nums_matrix) <- c(paste0("FieldCap", 1:length(nums_FieldCap)), 
                               names(opt_numeric_parameters))
  }
  if(length(nums_FieldCap) != 0 & length(nums_KSat) != 0){
    colnames(nums_matrix) <- c(paste0("KSat", 1:length(nums_KSat)), paste0("FieldCap", 1:length(nums_FieldCap)),
                               names(opt_numeric_parameters))
  }
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
if(length(nums_KSat) + length(nums_FieldCap) != 0){
  if(length(nums_FieldCap) == 0){
    colnames(nums_matrix) <- c(paste0("KSat", 1:length(nums_KSat)), 
                               names(opt_numeric_parameters))
  }
  if(length(nums_KSat) == 0){
    colnames(nums_matrix) <- c(paste0("FieldCap", 1:length(nums_FieldCap)), 
                               names(opt_numeric_parameters))
  }
  if(length(nums_FieldCap) != 0 & length(nums_KSat) != 0){
    colnames(nums_matrix) <- c(paste0("KSat", 1:length(nums_KSat)), paste0("FieldCap", 1:length(nums_FieldCap)),
                               names(opt_numeric_parameters))
  }
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
  
  if(!("numerics_tracker.csv" %in% list.files(result_folder))){
    # initial compiling and function checking
    system("chmod 775 05_first_run_preparation.sh")
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
  }
}

# Run best model again in training and validation time -----------------------------------
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

# Training time
# run model an dcompute losses
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

# save training results
nums_matrix <-  matrix(opt_numeric_parameters, ncol = length(opt_numeric_parameters))
colnames(nums_matrix) <- names(opt_numeric_parameters)
trainin_results_overview <- data.frame(KSat = KSat_tf,
                                       FieldCap = FieldCap_tf,
                                       wmulti_loss = wmulti_obj_loss,
                                       NSE = NSE,
                                       lNSE = lNSE,
                                       KGE = KGE,
                                       domain_multi_nse_loss,
                                       opt_numeric_parameters,
                                       stringsAsFactors = FALSE, check.names = FALSE)
write.csv(trainin_results_overview,
          paste0(result_folder, "/training_results_overview.csv"), row.names = FALSE, 
)
write.csv(mhm_results,
          paste0(result_folder, "/training_results_metrics.csv"), row.names = FALSE, 
)

# Validation time
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
# run model and compute losses
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
# save validation results
nums_matrix <-  matrix(opt_numeric_parameters, ncol = length(opt_numeric_parameters))
colnames(nums_matrix) <- names(opt_numeric_parameters)
validation_results <- data.frame(KSat = KSat_tf,
                                 FieldCap = FieldCap_tf,
                                 wmulti_loss = wmulti_obj_loss,
                                 NSE = NSE,
                                 lNSE = lNSE,
                                 KGE = KGE,
                                 domain_multi_nse_loss,
                                 opt_numeric_parameters,
                                 stringsAsFactors = FALSE, check.names = FALSE)
write.csv(validation_results,
          paste0(result_folder, "/validation_resultions_overview.csv"), row.names = FALSE, 
)
write.csv(mhm_results,
          paste0(result_folder, "/validation_resultions_metrics.csv"), row.names = FALSE, 
)





