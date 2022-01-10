setwd("C:/Users/morit/Dropbox/Diss/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_master_namelists")

file.copy("mpr_fso.nml", 
          "mpr_fso_backup.nml", overwrite = TRUE)

nml <- readLines("mpr_fso.nml")

head(nml, 40)

# add slope, aspect, dem notill and till variables

# data arrays
variables_to_add <- c(paste0(c("aspect", "dem", "slope"), "_till"),
                      paste0(c("aspect", "dem", "slope"), "_notill"))
new_linie <- paste0("                  ", 
                    paste0(variables_to_add, collapse = ", "), ", ")
data_array_line <- grep("'KSat', 'KSat_till_temp', ", nml)

nml <- c(nml[1:(data_array_line-1)],
         new_linie,
         nml[data_array_line:length(nml)])
data_array_number <- grep("name(1:111) = 'land_cover'", nml, fixed = TRUE)
gsub("name\(1:[0-9]\)", )

# 1:19 alt
# 20:25 neu
# 26:117 alt


names(21) = "dem_till"
from_data_arrays(1,21) = "dem"
target_dim_names(1:4,21) = "land_cover_period", "lat", "lon", "horizon_till"
upscale_ops(1:4,21) = "1.0", "1.0", "1.0", "1.0"
to_file(21) = .false.






