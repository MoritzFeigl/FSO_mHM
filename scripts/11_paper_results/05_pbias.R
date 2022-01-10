# mHM FSO results
# FSO-mHM project
# Moritz Feigl, Apr 2021


library(ggplot2)
library(budyko)
if(Sys.info()[1] == "Linux"){
  setwd("/mnt/Data/Dropbox/Diss/FSO_mHM/2020_fso_mhm/scripts/11_paper_results")
  path <- "/mnt/Data/Dropbox/Diss/FSO_mHM/Results/"
} else {
  setwd("C:/Users/morit/Dropbox/Diss/FSO_mHM/2020_fso_mhm/scripts/11_paper_results")
  path <- "C:/Users/morit/Dropbox/Diss/FSO_mHM/Results/"
}

pal <- c(rev(as.character(yarrr::piratepal("basel", length.out = 7)))[1:6], "#FEC10BFF")
pal6 <- pal[1:6]
pal5 <- pal[1:5]
pal2 <- pal[c(2, 7)]

# 1. compute bias ------------------------------------------------------------------------
run_names <- c(paste0("FSO_SCE_NSE_run_", 1:4), "FSO_DDS_NSE_run_5")
mean_pbias <- data.frame(Run = paste0("Run ", 1:5), mean_pbias = NA, check.names = FALSE)
all_pbias <- NULL
# compute pbias from Q obs and Q pred for all runs
for(run in seq_along(run_names)){
  # basin wise metrics
  q_val_date <- read.csv(paste0(path, run_names[run], "/small_basins_predictions.csv"))
  q_val <- q_val_date[, -which(names(q_val_date) == "date")]
  
  pbiases <- NULL
  for(i in 1:ncol(q_val)){
    if(i %% 2 == 0) next
    pbiases <- rbind(pbiases, 
                     data.frame(Basin = as.character(as.integer(gsub("Qobs_", "", names(q_val)[i], fixed = TRUE))),
                                pbias = hydroGOF::pbias(q_val[, i+1], q_val[, i], na.rm = TRUE)))
  }
  write.csv(pbiases, paste0(path, run_names[run], "/run", run, "_pbias_validation_basins.csv"))
  all_pbias <- rbind(all_pbias,
                     cbind("Run" = paste0("Run ", run), pbiases))
  mean_pbias$mean_pbias[run] <- mean(pbiases$pbias)
}  
medians <- aggregate(pbias ~ Run, all_pbias, 
                     function(x) round(median(x), 3))
medians <- medians[order(as.integer(substring(medians$Run, 17))), ]
medians <- format(median(nsmall = 2))
all_pbias$Run <- factor(all_pbias$Run, 
                                levels = paste0("Run ", 1:5),
                                #"mhm_standard_parameter"),
                                #labels = c(paste0("Run ", 1:5), "Standard\nparameter"))
                                labels = paste0("Run ", 1:5, "\n", medians$pbias))

ggplot(all_pbias[all_pbias$pbias < 400, ], aes(Run, pbias, fill = Run)) +
  geom_boxplot() +
  theme_bw() + 
  scale_fill_manual(values = pal5, name = "", labels = paste0("Run ", 1:5)) +
  labs(x = "", y = "Bias (%)") +
  ggsave("plots/10_pbias.png", width = 5, height = 4, units = "in")





# 1. compute bias for sampled validation basins ------------------------------------------
sval <- read.csv("../10_validation/sampled_validation_basins.csv")

run_names <- c(paste0("FSO_SCE_NSE_run_", 1:4), "FSO_DDS_NSE_run_5")
mean_sampled_pbias <- data.frame(Run = paste0("Run ", 1:6), mean_pbias = NA, check.names = FALSE)
all_sampled_pbias <- NULL
# compute pbias from Q obs and Q pred for all runs
for(run in seq_along(run_names)){
  # basin wise metrics
  q_val_date <- read.csv(paste0(path, run_names[run], "/small_basins_predictions.csv"))
  q_val <- q_val_date[, -which(names(q_val_date) == "date")]
  
  ids <- NULL
  for(sval_basin in sval$sampled_validation_basins){
    ids <- c(ids, grep(sval_basin, names(q_val)))
  }
  pbiases <- NULL
  for(i in ids){
    if(i %% 2 == 0) next
    pbiases <- rbind(pbiases, 
                     data.frame(Basin = as.character(as.integer(gsub("Qobs_", "", names(q_val)[i], fixed = TRUE))),
                                pbias = hydroGOF::pbias(q_val[, i+1], q_val[, i], na.rm = TRUE)))
  }
  all_sampled_pbias <- rbind(all_sampled_pbias,
                     cbind("Run" = paste0("Run ", run), pbiases))
  mean_sampled_pbias$mean_pbias[run] <- mean(pbiases$pbias)
}  

# median pbias
medians <- aggregate(pbias ~ Run, all_sampled_pbias, 
                     function(x) round(median(x), 3))
medians <- medians[order(as.integer(substring(medians$Run, 17))), ]
all_sampled_pbias$Run <- factor(all_sampled_pbias$Run, 
                         levels = paste0("Run ", 1:5),
                         #"mhm_standard_parameter"),
                         #labels = c(paste0("Run ", 1:5), "Standard\nparameter"))
                         labels = paste0("Run ", 1:5, "\n", medians$pbias))

ggplot(all_sampled_pbias, aes(Run, pbias, fill = Run)) +
  geom_boxplot() +
  scale_fill_manual(values = pal6, name = "", labels = paste0("Run ", 1:5)) +
  labs(x = "", y = "Bias (%)") +
  theme_bw() + 
  ggsave("plots/11_sampled_pbias.png", width = 5, height = 4, units = "in")


  eins <- read.csv(paste0(path, run_names[1], "/small_basins_predictions.csv"))
  zwei <- read.csv(paste0(path, run_names[2], "/small_basins_predictions.csv"))
testthat::is_equivalent_to() eins == zwei
