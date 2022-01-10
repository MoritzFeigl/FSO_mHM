# Estimate TF distributions for FSO
# FSO-mHM project
# Moritz Feigl, Sep 2020

.libPaths("/home/lv71468/mfeigl/R/x86_64-pc-linux-gnu-library/4.0")
library(FSO)
setwd("/home/lv71468/mfeigl/FSO_mHM/FSO_Data/prepared_TFs/prepared_functions")

# Start script
variable_df_sampled <- feather::read_feather("/home/lv71468/mfeigl/FSO_mHM/mHM_setup/FSO_mHM_major_basins/FSO_sp_dataframe/spatial_predictors_sampled.feather")
variable_df_sampled <- variable_df_sampled[, -which(names(variable_df_sampled) == "basin")]
# Use only till variables for distribution estimation
variable_df_sampled$KSat_notill <- NULL
variable_df_sampled$vGenu_n_notill <- NULL
variable_df_sampled$ThetaS_notill <- NULL
names(variable_df_sampled)[7:9] <- c("KSat", "vGenu_n", "ThetaS")
# Compute function distributions
functions_for_dist <- "dummy_file"
# scaling_bounds
scaling_bounds <- list("slope" = c(0, 90),
                       "aspect" = c(0, 360),
                       "bd" = c(0, 2.3), #max(variable_df$bd)), # is actually 2.2769
                       "sand" = c(0, 100),
                       "clay" = c(0, 100),
                       "dem" = c(0, 4000),
                       "KSat" = c(9, 244),
                       "vGenu_n" = c(1,2),
                       "ThetaS" = c(0.24, 0.51)
)
# calc dist for all given files
functions_para <- fst::read_fst(functions_for_dist)
dist_name <- gsub("simplified", "distribution", functions_for_dist)
distribution_sampler(functions = functions_para,
                     variable_df = variable_df_sampled,
                     scaling_bounds = scaling_bounds,
                     file_name = dist_name,
                     no_cores = 32)


