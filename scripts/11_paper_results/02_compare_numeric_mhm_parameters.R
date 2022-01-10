

library(ggplot2)
#devtools::install_github("hrbrmstr/ggalt")
library(ggplot2)
library(ggalt)
library(budyko)
setwd("C:/Users/morit/Dropbox/Diss/FSO_mHM/2020_fso_mhm/scripts/11_paper_results")
path <- "C:/Users/morit/Dropbox/Diss/FSO_mHM/Results/"
# color palettes
pal <- c(rev(as.character(yarrr::piratepal("basel", length.out = 7)))[1:6], "#FEC10BFF")
pal6 <- pal[1:6]
pal5 <- pal[1:5]
pal2 <- pal[c(2, 7)]
# 1. compile runs ------------------------------------------------------------------------
run_names <- c(paste0("FSO_SCE_NSE_run_", 1:4), "FSO_DDS_NSE_run_5")

# load TF estimates, performance results
training_overview <- NULL
validation_overview <- NULL
metrics <- NULL
for(run in run_names){
  # overviews basins
  train <- cbind("Run" = run, "Basins" = "Training", "time_period" = "Training",
                 read.csv(paste0(path, run, "/training_results_overview.csv")))
  training_overview <- rbind(training_overview, train)
}


paras <- training_overview[, c(1, 17:75)]
# get para bounds and standard parameters
source("C:/Users/morit/Dropbox/Diss/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/04_FSO_mhm_utils.R")
parameters <- suppressWarnings(
  read_nml("C:/Users/morit/Dropbox/Diss/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_master_namelists/mhm_parameter.nml"))
parameter_df <- do.call(rbind, parameters$mhm_parameters)
ind_of_num_paras <- which(parameter_df[, 4] == 1)
num_para_bounds <- parameter_df[ind_of_num_paras, c(1, 2)]
standard_parameters <- data.frame("Run" = "standard", 
                             matrix(parameter_df[ind_of_num_paras, 3], ncol = ncol(paras) - 1))
colnames(standard_parameters) <- colnames(paras)
paras <- rbind(paras, standard_parameters)

.range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

paras_scaled <- paras
for(par in 1:nrow(num_para_bounds)){
  paras_scaled[, par+1] <- .range01(paras_scaled[, par+1], 
                                    num_para_bounds[par, 1],
                                    num_para_bounds[par, 2])
                                    
}


paras_plot <- reshape2::melt(paras_scaled, id = "Run")
paras_plot$Run <- factor(paras_plot$Run,
                         levels = c("standard", paste0("FSO_SCE_NSE_run_", 1:4),"FSO_DDS_NSE_run_5"),
                         labels = c("Standard Parameters", paste0("Run ", 1:5)))
paras_plot$variable <- factor(paras_plot$variable,
                              levels = rev(names(paras_scaled)[-1]))

ggplot(paras_plot, aes(x = value, y = variable, col = Run, shape = Run)) +
  geom_point(size = 2) +
  xlim(0, 1) + 
  scale_shape_manual( values = c(3, 15, 15, 15, 15, 15)) +
  scale_color_manual(values = c("black", pal)) +
  labs(color = "", y = "", x = "scaled parameter values") + guides(shape = FALSE) +
  guides(color = guide_legend(override.aes = list(shape = c(3, 15, 15, 15, 15, 15)))) +
  theme_light()+
  ggsave("plots/08_numeric_parameter_distribution.png", width = 8, height = 7, units = "in")


# KSat specific parameters
ksat_paras <- c("infiltrationShapeFactor", "imperviousStorageCapacity",
                "interflowStorageCapacityFactor", "interflowRecession_slope", 
                "fastInterflowRecession_forest", "slowInterflowRecession_Ks", 
                "exponentSlowInterflow", "rechargeCoefficient", "rechargeFactor_karstic")
ggplot(paras_plot[paras_plot$variable %in% ksat_paras & 
                    paras_plot$Run %in% c("Standard Parameters", paste0("Run ", 2:4)), ], 
       aes(x = value, y = variable, col = Run, shape = Run)) +
  geom_point(size = 4) +
  xlim(0, 1) + 
  scale_shape_manual( values = c(3, 15, 15, 15)) +
  scale_color_manual(values = c("black", pal5[c(2, 3, 4)])) +
  labs(color = "", y = "", x = "scaled parameter values") + guides(shape = FALSE) +
  guides(color = guide_legend(override.aes = list(shape = c(3, 15, 15, 15)))) +
  theme_light()+
  theme(axis.text = element_text(size=14), 
        axis.title=element_text(size=14),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  ggsave("plots/08_numeric_parameter_distribution_KSat.png", width = 9, height = 4, units = "in")


# FieldCap specific parameters
FieldCap_paras <- as.character(unique(paras_plot$variable)[10:24])
ggplot(paras_plot[paras_plot$variable %in% FieldCap_paras & 
                    paras_plot$Run %in% c("Standard Parameters", paste0("Run ", 2:4)), ], 
       aes(x = value, y = variable, col = Run, shape = Run)) +
  geom_point(size = 4) +
  xlim(0, 1) + 
  scale_shape_manual( values = c(3, 15, 15, 15)) +
  scale_color_manual(values = c("black", pal5[c(2, 3, 4)])) +
  labs(color = "", y = "", x = "scaled parameter values") + guides(shape = FALSE) +
  guides(color = guide_legend(override.aes = list(shape = c(3, 15, 15, 15)))) +
  theme_light()+
  theme(axis.text = element_text(size=14), 
        axis.title=element_text(size=14),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  ggsave("plots/08_numeric_parameter_distribution_FieldCap.png", width = 9, height = 6, units = "in")









