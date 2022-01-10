# Create spatial predictors data frame from netcdf
# FSO-mHM project
# Moritz Feigl, Aug 2020

if("RNetCDF" %in% rownames(installed.packages()) == FALSE) {
  install.packages("RNetCDF")
}
if("feather" %in% rownames(installed.packages()) == FALSE) {
  install.packages("feather")
}
library(RNetCDF)
setwd("/home/cfgrammar/Dropbox/Diss/FSO_mHM/mHM_setup/FSO_mHM_major_basins/static")

if(!dir.exists("../FSO_sp_dataframe")) dir.create("../FSO_sp_dataframe")
basins <- list.files()
for (basin in basins){
  # load mpr input netcdf, aggregate layers and combine in data.frame
  cat("******", basin, "******\n")
  variables <- c("mpr/bd", "mpr/sand", "mpr/clay", "mpr/slope", "mpr/aspect", "routing/dem")
  for(variable in variables){
    cat(variable, "\n")
    nc <- open.nc(paste0(basin, "/", variable, ".nc"))
    var <- var.get.nc(nc, unlist(strsplit(variable, "/"))[2])
    if(length(dim(var)) == 3) var <- apply(var, c(2, 3), function(x) mean(x, na.rm = TRUE))
    var_vector <- as.data.frame(as.vector(var))
    names(var_vector) <- variable
    if(!exists("spatial_predictors")) {
      spatial_predictors <- var_vector
    } else {
      spatial_predictors <- cbind(spatial_predictors, var_vector)
    }
    close.nc(nc)
  }
  
  # mHM specific variables
  variables <- c("KSat_till", "KSat_notill", "vGenu_n_till", "vGenu_n_notill", "ThetaS_till", "ThetaS_notill")
  nc <- open.nc(paste0("../FSO_spatial_predictors/", basin, "_mHM_parameters.nc"))
  for(variable in variables){
    cat(variable, "\n")
    var <- var.get.nc(nc, variable)
    if(length(dim(var)) > 2) var2 <- apply(var, c(2, 3), function(x) mean(x, na.rm = TRUE))
    var_vector <- as.data.frame(as.vector(var))
    names(var_vector) <- variable
    spatial_predictors <- cbind(spatial_predictors, var_vector)
  }
  if(length(nas) != 0){ # added after it failed by simon schitz
    nas <- which(is.na(spatial_predictors), arr.ind = TRUE)
  }
  spatial_predictors <- spatial_predictors[-nas[, 1], ]
  feather::write_feather(spatial_predictors,
                         paste0("../FSO_sp_dataframe/",
                                basin, "_spatial_predictors.feather"))
  rm(spatial_predictors)
}

all_sp_files <- list.files("../FSO_sp_dataframe/", full.names = TRUE)
cat("Aggregating spatial predictor data frames\n")
for(basin in basins){
  cat(basin, "\n")
  sp_curr <- feather::read_feather(all_sp_files[grep(basin, all_sp_files)])
  sp_curr$basin <- basin
  if(!exists("spatial_predictors")) {
    spatial_predictors <- sp_curr
  } else {
    spatial_predictors <- rbind(spatial_predictors, sp_curr)
  }
  rm(sp_curr)
}
names(spatial_predictors) <- c("bd", "sand", "clay", "slope", "aspect", "dem", 
                               "KSat_till", "KSat_notill", "vGenu_n_till", 
                               "vGenu_n_notill", "ThetaS_till", "ThetaS_notill", "basin")
feather::write_feather(spatial_predictors,"../FSO_sp_dataframe/spatial_predictors.feather")


# Create sampled sp data.frame
spatial_predictors <- feather::read_feather("../FSO_sp_dataframe/spatial_predictors.feather")
for(basin in basins){
  basin_data <- spatial_predictors[spatial_predictors$basin == basin, ]
  # sample always 0.1 % of the basin cells
  samples <- ceiling(nrow(basin_data)*0.01)
  if(exists("variable_df_sampled")){
    sp_df_sampled <- rbind(sp_df_sampled,
                                 basin_data[sample(nrow(basin_data), samples), ])
  } else {
    sp_df_sampled <- basin_data[sample(nrow(basin_data), samples), ]
  }
}

feather::write_feather(sp_df_sampled,
                       "../FSO_sp_dataframe/spatial_predictors_sampled.feather")

