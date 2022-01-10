# mHM FSO results
# FSO-mHM project
# Moritz Feigl, Mar 2021


library(ggplot2)
library(budyko)
if(Sys.info()[1] == "Linux"){
  setwd("/mnt/Data/Dropbox/Diss/FSO_mHM/2020_fso_mhm/scripts/11_paper_results")
  path <- "/mnt/Data/Dropbox/Diss/FSO_mHM/Results/"
} else {
  setwd("C:/Users/morit/Dropbox/Diss/FSO_mHM/2020_fso_mhm/scripts/11_paper_results")
  path <- "C:/Users/morit/Dropbox/Diss/FSO_mHM/Results/"
}

# color palettes
pal <- c(rev(as.character(yarrr::piratepal("basel", length.out = 7)))[1:6], "#FEC10BFF")
pal6 <- pal[1:6]
pal5 <- pal[1:5]
pal2 <- pal[c(2, 7)]


# 1. compile runs ------------------------------------------------------------------------
run_names <- c(paste0("FSO_SCE_NSE_run_", 1:4), "FSO_DDS_NSE_run_5", "mhm_standard_parameter")
sampled_basins <- read.csv("../10_validation/sampled_validation_basins.csv")
zink <- read.table("../10_validation/lut_gauges_de.mpeg", 
                   header = TRUE)

# load TF estimates, performance results
training_overview <- NULL
validation_overview <- NULL
metrics <- NULL
sval_metrics <- NULL
for(run in run_names){
  # basin wise metrics
  train_metrics <- cbind("Run" = run, "Basins" = "Training", "time_period" = "Training",
                         read.csv(paste0(path, run, "/training_results_metrics.csv")))
  train_val_metrics <- cbind("Run" = run, "Basins" = "Training", "time_period" = "Validation",
                             read.csv(paste0(path, run, "/validation_resultions_metrics.csv")))
  val_metrics <- cbind("Run" = run, "Basins" = "Validation", "time_period" = "Validation",
                       read.csv(paste0(path, run, "/all_validation_basins_results.csv")))
  if(names(val_metrics)[4] == "X") val_metrics <- val_metrics[, -4]
  sval_metric <- cbind("Run" = run, "Basins" = "sampled_validation", "time_period" = "Validation",
                       read.csv(paste0(path, run, "/sampled_validation_basins_metrics.csv")))
  metrics <- rbind(metrics, train_metrics, train_val_metrics, val_metrics)
  sval_metrics <- rbind(sval_metrics, sval_metric)
  # overviews basins
  train <- cbind("Run" = run, "Basins" = "Training", "time_period" = "Training",
                 read.csv(paste0(path, run, "/training_results_overview.csv")))
  val <- cbind("Run" = run, "Basins" = "Training", "time_period" = "Validation",
               read.csv(paste0(path, run, "/validation_resultions_overview.csv")))
  training_overview <- rbind(training_overview, train, val)
  sval <- cbind("Run" = run, "Basins" = "sampled_validation", "time_period" = "Validation",
                read.csv(paste0(path, run, "/sampled_validation_basins_overview.csv")))
  names(sval)[grep("X[0-9]", names(sval))] <- sampled_basins$sampled_validation_basins
  all_val <- sval
  # metrics for all validation basins
  KGE <- mean(val_metrics$KGE)
  NSE <- mean(val_metrics$NSE)
  lNSE <- mean(val_metrics$lNSE)
  multi_obj_nse_loss <- sign(val_metrics$NSE) * 0.5*abs(val_metrics$NSE)^6 +
    sign(val_metrics$lNSE) * 0.5*abs(val_metrics$lNSE)^6
  multi_obj_nse_loss <- ifelse(multi_obj_nse_loss < 0,
                               -1*(multi_obj_nse_loss*-1)^(1/6),
                               multi_obj_nse_loss^(1/6))
  median_loss <- median(multi_obj_nse_loss)
  domain_multi_nse_loss <- matrix(rep(NA, 20), ncol = 20)
  all_val[, 6:29] <-  data.frame(median_loss = median_loss,
                                 NSE = NSE,
                                 lNSE = lNSE,
                                 KGE = KGE,
                                 domain_multi_nse_loss)
  all_val$Basins <- "Validation"
  validation_overview <- rbind(validation_overview, sval, all_val)
}
# add more info to metrics
metrics$Basin <- sapply(metrics$Basin, function(x) gsub("sub_", "", x, fixed = TRUE))
sval_metrics$Basin <- sapply(sval_metrics$Basin, function(x) gsub("sub_", "", x, fixed = TRUE))

multi_obj_nse_loss <- function(NSE, lNSE){
  multi_obj_nse <- sign(NSE)*0.5*abs(NSE)^6+sign(lNSE)*0.5*abs(lNSE)^6
  loss <- ifelse(multi_obj_nse < 0,
                 -1*(multi_obj_nse*-1)^(1/6),
                 multi_obj_nse^(1/6))
  return(loss)
}
metrics$loss <- multi_obj_nse_loss(metrics$NSE, metrics$lNSE)
sval_metrics$loss <- multi_obj_nse_loss(sval_metrics$NSE, sval_metrics$lNSE)
# save as csv
write.csv(training_overview, "training_overview.csv", row.names = FALSE)
write.csv(validation_overview, "validation_overview.csv", row.names = FALSE)
write.csv(metrics, "metrics.csv", row.names = FALSE)
write.csv(sval_metrics, "sval_metrics", row.names = FALSE)

#remove all basins from zink that were not part of the zink et al paper
zink <- zink[zink$Stat_ID %in% metrics$Basin, ]
# 2. Compare Runs ------------------------------------------------------------------------
# melt data for plotting
plot_data <- reshape2::melt(metrics, id.vars = c("Run", "Basin", "Basins", "time_period"))
plot_data$variable <- as.character(plot_data$variable)
plot_data$criteria <- "KGE"
plot_data$criteria[grep("loss", plot_data$variable)] <- "loss"
plot_data$criteria[grep("lNSE", plot_data$variable)] <- "log NSE"
plot_data$criteria[grepl("NSE", plot_data$variable) & !grepl("lNSE", plot_data$variable)] <- "NSE"
plot_data$variable <- NULL
plot_data$criteria <- factor(plot_data$criteria, levels = c("loss", "NSE", "log NSE", "KGE"),
                             labels = c("Loss", "NSE", "log NSE", "KGE"))
plot_data <- plot_data[plot_data$Run != "mhm_standard_parameter",]
plot_data$Run <- factor(plot_data$Run, 
                        levels = c("FSO_SCE_NSE_run_1", 
                                   "FSO_SCE_NSE_run_2", 
                                   "FSO_SCE_NSE_run_3", 
                                   "FSO_SCE_NSE_run_4", 
                                   "FSO_DDS_NSE_run_5"),
                        #"mhm_standard_parameter"),
                        labels = paste0("Run ", 1:5))#, "mHM\ndefault"))

# Training Baisins Plot
summary(plot_data[plot_data$Basins == "Training" & 
                    plot_data$time_period == "Training" &
                    plot_data$criteria == "log NSE", 5] -
          plot_data[plot_data$Basins == "Training" & 
                      plot_data$time_period == "Validation"&
                      plot_data$criteria == "log NSE", 5])

train_performance <- aggregate(value ~ criteria + Run, 
                               plot_data[plot_data$Basins == "Training" & 
                                           plot_data$time_period == "Training"&
                                           plot_data$criteria == "KGE", ], median)
val_performance <- aggregate(value ~ criteria + Run, 
                             plot_data[plot_data$Basins == "Training" & 
                                         plot_data$time_period == "Validation", ], median)
perf_diffs <- cbind(train_performance, "validation" = val_performance[, 3], 
                    diff = val_performance$value - train_performance$value)
mean(perf_diffs[perf_diffs$criteria == "NSE", "diff"])
mean(perf_diffs[perf_diffs$criteria == "log NSE", "diff"])
mean(perf_diffs[perf_diffs$criteria == "KGE", "diff"])

aggregate(value ~ criteria + Run, plot_data[plot_data$Basins == "Training", ], median)
vars <- aggregate(value ~ criteria + Run, plot_data[plot_data$Basins == "Training" & 
                                                      plot_data$time_period == "Training" & 
                                                      plot_data$criteria == "NSE", ], var)
mean(vars[c(1, 5), 3])
mean(vars[2:4, 3])


ggplot(plot_data[plot_data$Basins == "Training", ], 
       aes(Run, value, fill = time_period)) + # , fill = Run
  geom_boxplot() +
  #scale_fill_manual(values = pal6, name = "") +
  #scale_color_manual(values = c("grey35", "black"), name = "Time period", labels = c("Calibration", "Validation")) +
  scale_fill_manual(values = c("grey45", "grey80"), name = "Time period", labels = c("Calibration", "Validation")) +
  theme_bw() + 
  #guides(fill = "none") +
  labs(x = "", y = "") + 
  theme(strip.text = element_text(size=18),
        axis.text = element_text(size=16), 
        axis.title=element_text(size=18),
        legend.text = element_text(size=16), 
        legend.title=element_text(size=18)) + 
  facet_wrap(~criteria)
ggsave("plots/01_training_basins_runs.png", width = 12, height = 5, units = "in", dpi = 600) 


# sampled validation Baisins Plot
sval_plot <- reshape2::melt(sval_metrics[, -c(2, 3)], id.vars = c("Run", "Basin"))
sval_plot$variable <- factor(sval_plot$variable, levels = c("loss", "NSE", "lNSE", "KGE"),
                             labels = c("Loss", "NSE", "log NSE", "KGE"))
sval_plot <- sval_plot[sval_plot$Run != "mhm_standard_parameter", ]
sval_plot$Run <- factor(sval_plot$Run, 
                        levels = c("FSO_SCE_NSE_run_1", 
                                   "FSO_SCE_NSE_run_2", 
                                   "FSO_SCE_NSE_run_3", 
                                   "FSO_SCE_NSE_run_4", 
                                   "FSO_DDS_NSE_run_5"),
                        #"mhm_standard_parameter"),
                        #labels = c(paste0("Run ", 1:5), "Standard\nparameter"))
                        labels = paste0("Run ", 1:5))

# compute bias for sampled validation basins 
sval <- read.csv("../10_validation/sampled_validation_basins.csv")

run_names <- c(paste0("FSO_SCE_NSE_run_", 1:4), "FSO_DDS_NSE_run_5")
all_sampled_pbias <- NULL
# compute pbias from Q obs and Q pred for sampled validation
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
                     data.frame("Run" = paste0("Run ", run),
                                "Basin" = as.character(as.integer(gsub("Qobs_", "", names(q_val)[i], fixed = TRUE))),
                                "variable" = "PBIAS",
                                "value" = hydroGOF::pbias(q_val[, i+1], q_val[, i], na.rm = TRUE)/100))
  }
  all_sampled_pbias <- rbind(all_sampled_pbias, pbiases)
  print(mean(pbiases$value))
}  
all_sampled_pbias$Run <- factor(all_sampled_pbias$Run, 
                                levels = paste0("Run ", 1:5))
aggregate(value ~ Run + variable, sval_plot_all, median)


# add PBIaS to sampled validation results and plot
sval_plot_all <- cbind(rbind(all_sampled_pbias, sval_plot), time_period = "validation")

sval_plot_all$variable <- factor(sval_plot_all$variable, 
                                 levels = c("PBIAS", "NSE", "log NSE", "KGE", "Loss"))
ggplot(sval_plot_all[sval_plot_all$variable != "Loss", ], 
       aes(Run, value)) + #, fill = Run, col = time_period)) +
  geom_boxplot() +
  scale_colour_manual(values = "black") +
  scale_fill_manual(values = pal5, name = "") +
  theme_bw() + 
  facet_wrap(~variable) +
  labs(x = "", y = "") + 
  theme(strip.text = element_text(size=18),
        axis.text = element_text(size=16), 
        axis.title=element_text(size=18),
        legend.text = element_text(size=16), 
        legend.title=element_text(size=18))
#legend.position = "none")
ggsave("plots/02_sampled_validation_basins_performance.png", width = 12, 
       height = 5, units = "in", dpi = 600) 


# Compute Budyko- Curve and plot
budyko_sval <- merge(sval_metrics, zink, by.x = "Basin", by.y = "Stat_ID")
budyko_sval$NSE <- round(budyko_sval$NSE, 2)
budyko_sval$'PET.P' <- budyko_sval$PET.mm_a.1./budyko_sval$Precip.mm_a.1.
budyko_sval$'AET.P' <- budyko_sval$AET.mm_a.1./budyko_sval$Precip.mm_a.1.
NSE_limit <- -0.5

# Budyko cureve of sampled basins
ogbudyko <-  budyko_sim()
blankBC +
  geom_line(data = ogbudyko)+
  coord_cartesian(xlim = c(0, 1.5)) + 
  geom_point(aes(size = catArea), fill = "grey65",
             data = budyko_sval[budyko_sval$Run == "FSO_SCE_NSE_run_2", ], 
             alpha = 0.8, pch = 21, colour = "black") + 
  scale_size(range = c(2, 10)) + 
  labs(size = "Basin area in km") + 
  theme_minimal() + 
  theme(axis.text = element_text(size=12), 
        axis.title=element_text(size=14),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14)) 
ggsave("plots/03_sampled_validation_basins_budyko.png", 
       width = 7, height = 4, units = "in") 



# Budyko cureve vs FSO Run 2
ogbudyko <-  budyko_sim()
blankBC +
  geom_line(data = ogbudyko)+
  coord_cartesian(xlim = c(0, 1.5)) + 
  geom_point(aes(fill = NSE, #cut(NSE, breaks= c(0.49, 0.7, 0.8, 1)), 
                 size = catArea), 
             data = budyko_sval[budyko_sval$Run == "FSO_SCE_NSE_run_2", ], 
             alpha = 0.8, pch = 21, colour = "black") + 
  scale_fill_gradientn(colours =c("red", "yellow", "green", "blue"), values = scales::rescale(c(0, 0.35, 0.65, 1)), 
                       limits = c(0, 1)) +  
  scale_size(range = c(2, 10)) + 
  labs(size = "Basin area in km?",
       fill = "Run 2 NSE") + 
  theme_minimal() + 
  theme(axis.text = element_text(size=12), 
        axis.title=element_text(size=14),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14)) 
ggsave("plots/03_sampled_validation_basins_budyko_run2performance.png", 
       width = 7, height = 4, units = "in") 



# Check for validation dependencie on catchment size
val_metrics_with_prop <- merge(metrics[metrics$Basins == "Validation" & 
                                         metrics$time_period == "Validation", c("Run", "Basin", "NSE")], 
                               zink, by.x = "Basin", by.y = "Stat_ID")

stand_par <- val_metrics_with_prop[val_metrics_with_prop$Run == "mhm_standard_parameter", ]
summary(lm(NSE ~ catArea, stand_par))
plot(stand_par$catArea, stand_par$NSE)
summary(lm(NSE ~ . , stand_par[, c(3, 4, 6:9, 33:37)]))

# 3. Compare with Zink -------------------------------------------------------------------
# prepare dfs for plotting
zink_plot_data <- data.frame("Run" = "Zink et al. (2017)", 
                             "Basin" = zink$Stat_ID, 
                             NSE = c(zink$NSEd_median, zink$NSEd_min, zink$NSEd_max, zink$NSEd_p05, zink$NSEd_p95))
plot_data2 <- reshape2::melt(metrics[metrics$Basins == "Validation" & 
                                       metrics$time_period == "Validation", c("Run", "Basin", "NSE")], 
                             id.vars = c("Run", "Basin"))[, -3]
plot_data2 <- plot_data2[plot_data2$Run != "mhm_standard_parameter", ]
names(plot_data2)[3] <- "NSE"

# median KGE, lNSE
val_data_all <- metrics[metrics$Basins == "Validation" & 
                          metrics$time_period == "Validation", ]
aggregate(KGE ~ Run, val_data_all, function(x) round(median(x), 3))
aggregate(lNSE ~ Run, val_data_all, function(x) round(median(x), 3))
aggregate(NSE ~ Run, val_data_all, function(x) round(median(x), 3))


# FSO runs vs Zink median
compare_plot2 <- rbind(plot_data2, zink_plot_data)
medians <- aggregate(NSE ~ Run, compare_plot2, 
                     function(x) round(median(x), 3))
medians$NSE <- format(medians$NSE, nsmall = 3)
medians <- medians[order(as.integer(substring(medians$Run, 17))), ]
compare_plot2_medians <- compare_plot2
compare_plot2$Run <- factor(compare_plot2$Run, 
                            levels = c("FSO_SCE_NSE_run_1", 
                                       "FSO_SCE_NSE_run_2", 
                                       "FSO_SCE_NSE_run_3", 
                                       "FSO_SCE_NSE_run_4", 
                                       "FSO_DDS_NSE_run_5",
                                       "Zink et al. (2017)"),
                            labels = c(paste0("Run ", 1:5), "Zink et al. (2017)"))
# test differences of run 2 vs zink
test_data <- compare_plot2[compare_plot2$Run %in% c("Run 2", "Zink et al. (2017)"), ]
kruskal.test(test_data$NSE, test_data$Run)
# plot
compare_plot2_medians$Run <- factor(compare_plot2_medians$Run, 
                                    levels = c("FSO_SCE_NSE_run_1", 
                                               "FSO_SCE_NSE_run_2", 
                                               "FSO_SCE_NSE_run_3", 
                                               "FSO_SCE_NSE_run_4", 
                                               "FSO_DDS_NSE_run_5",
                                               "Zink et al. (2017)"),
                                    #"mhm_standard_parameter"),
                                    #labels = c(paste0("Run ", 1:5), "Standard\nparameter"))
                                    labels = c(paste0("Run ", 1:5, "\n", medians$NSE[1:5]),
                                               paste0("Zink et al. (2017)\n", medians$NSE[6])))


ggplot(compare_plot2[compare_plot2$Run %in% c(paste0("Run ", 2), "Zink et al. (2017)"), ], 
       aes(Run, NSE)) + 
  geom_violin(fill = "grey80") + ylim(0, 1) +
  geom_boxplot(width=.2) +
  theme_bw() + xlab("") +
  scale_fill_manual(values = pal[c(2, 7)], 
                    name = "", 
                    labels = c(paste0("Run ", 2), "Zink et al. (2017)")) + 
  theme(legend.position = "none", text = element_text(size=20)) 
ggsave("plots/03_validation_basins_NSE.png", width = 7, height = 4, units = "in", dpi = 600) 

# scatterplot FSO vs. Zink
# create df with baisins (row) and results (columns)
plot_data2$Run <- factor(plot_data2$Run, 
                         levels = c("FSO_SCE_NSE_run_1", 
                                    "FSO_SCE_NSE_run_2", 
                                    "FSO_SCE_NSE_run_3", 
                                    "FSO_SCE_NSE_run_4", 
                                    "FSO_DDS_NSE_run_5"),
                         labels = paste0("Run ", 1:5))
for(run in 1:5){
  basin_df <- merge(plot_data2[plot_data2$Run == paste0("Run ", run), ], 
                    zink_plot_data, 
                    by = "Basin", suffixes = c(".fso", ".zink"))
  basin_df$diff <- basin_df$NSE.fso - basin_df$NSE.zink
  basin_df$status <- "FSO > Zink et al."
  basin_df$status[basin_df$diff < 0] <- "FSO < Zink et al."
  ggplot(basin_df, aes(x = NSE.fso, y = NSE.zink, col = factor(status))) + geom_point() +
    labs(col = "", x = paste0("NSE FSO Run ", run) , y = "median NSE Zink et al. 2017") +
    theme_bw()  + 
    scale_color_manual(values = c("firebrick2", "green3")) +
    xlim(0, 1) + ylim(0, 1) +
    annotate("text", x = 0.2, y = 0.1, 
             label = paste0("FSO > Zink: ", sum(basin_df$status == "FSO > Zink et al."), " \n",
                            "FSO < Zink: ", sum(basin_df$status == "FSO < Zink et al.")))+
    ggsave(paste0("plots/04_FSO_vs_zink_scatter_run_", run, ".png"), 
           width = 5, height = 4, units = "in") 
}

# NSE densities
ggplot(compare_plot2[compare_plot2$Run %in% c(paste0("Run ", 2), "Zink et al. (2017)"),], 
       aes(x = NSE, fill = Run)) + 
  geom_density(alpha = 0.7, lwd = 0.6) + 
  scale_fill_manual(values = pal2) +
  theme_minimal() + 
  theme(axis.text = element_text(size=14), 
        axis.title=element_text(size=16),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16)) + 
  labs(fill = "", x = "NSE", y = "f(x)") + 
  theme(legend.title = element_blank()) +
  xlim(NA, 1) 
ggsave("plots/05_NSE_densities.png", width = 8, height = 5, units = "in", dpi = 600) 

# cum distribution
ggplot(compare_plot2[compare_plot2$Run %in% c(paste0("Run ", 2), "Zink et al. (2017)"),], 
       aes(NSE, colour = Run)) + 
  stat_ecdf(lwd = 1.2) + 
  labs(x = "NSE", y = "F(x)",
       colour = "", linetype = "") +
  scale_colour_manual(values = pal2) +
  xlim(0, 1) + 
  theme_minimal() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=12),
        legend.text = element_text(size=12),
        legend.key.size = unit(1,"cm")) +
  ggsave("plots/06_NSE_cum_densities.png", width = 6, height = 4, units = "in") 

# simple run 2 cum dist
ggplot(compare_plot2[compare_plot2$Run =="Run 2", ], 
       aes(NSE, colour = Run)) + 
  stat_ecdf(lwd = 1.2) + 
  labs(x = "NSE", y = "F(x)",
       colour = "", linetype = "") +
  scale_colour_manual(values = "forestgreen") +
  xlim(0, 1) + 
  theme_minimal() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14),
        legend.text = element_text(size=14),
        legend.key.size = unit(1,"cm")) +
  ggsave("plots/06_minimal_run2_cum_density.png", width = 6, height = 4, units = "in") 



# Budyko cureve of FSO Run 2 validation
budyko_val_metrics <- metrics[metrics$Basins == "Validation", ]
budyko_val_metrics$Basin <- as.integer(gsub("sub_", "", budyko_val_metrics$Basin))
budyko_val <- merge(budyko_val_metrics, zink, by.x = "Basin", by.y = "Stat_ID")
budyko_val$NSE <- round(budyko_val$NSE, 2)
budyko_val$'PET.P' <- budyko_val$PET.mm_a.1./budyko_val$Precip.mm_a.1.
budyko_val$'AET.P' <- budyko_val$AET.mm_a.1./budyko_val$Precip.mm_a.1.
NSE_limit <- 0.55
# linear model to test budyko vs NSE
lm_data <- budyko_val[budyko_val$Run == "FSO_SCE_NSE_run_2", ]
lm_data$Altitude[lm_data$Altitude == -999] <- NA
lm_data$Altitude[lm_data$Altitude == -99] <- NA
summary(lm(NSE ~ PET.P + catArea + Altitude, lm_data))
plot(lm_data$NSE, log(lm_data$catArea))
plot(lm_data$NSE, lm_data$Altitude)
cor(lm_data[, c("Altitude", "PET.P",  "catArea", "Easting_X", "Norting_Y")], 
    use = "pairwise.complete.obs")
car::vif(lm(NSE ~ PET.P + catArea + Altitude, lm_data))

ogbudyko <-  budyko_sim()
blankBC +
  geom_line(data = ogbudyko)+
  coord_cartesian(xlim = c(0, 1.5)) + 
  geom_point(aes(fill = NSE, 
                 size = catArea), 
             data = budyko_val[budyko_val$Run == "FSO_SCE_NSE_run_2", ], 
             alpha = 0.8, pch = 21, colour = "black") + 
  #scale_fill_gradient(low = "red", high = "green", limits = c(0.4, 0.9)) + 
  scale_fill_gradientn(colours =c("red", "yellow", "green", "blue"), 
                       values = scales::rescale(c(0.4, 0.55, 0.65, 0.8, 0.9)), 
                       limits = c(0.4, 0.9)) +  
  scale_size(range = c(2, 15)) + 
  labs(size = "Basin area in km?",
       fill = "Run 2 NSE",
       y = "Evaporative Index (AET/P)",
       x = "Aridity Index (PET/P)") + 
  theme_minimal() + 
  theme(axis.text = element_text(size=16), 
        axis.title=element_text(size=18),
        legend.text = element_text(size=16),
        legend.title = element_text(size=18)) 
ggsave("plots/09_validation_basins_budyko_run2performance.png", 
       width = 7, height = 5, units = "in", dpi = 600) 

# NSE vs long/lat
side_plot <- budyko_val[, c("NSEd_median", "NSE", "Easting_X", "Norting_Y")]
side_plot$diff <- side_plot$NSE - side_plot$NSEd_median
summary(lm(NSE ~ -1 + Easting_X + Norting_Y + Easting_X*Norting_Y, side_plot))

ggplot(side_plot, aes(x = Easting_X, y = diff)) + geom_point() + ylim(-1, 1)
ggplot(side_plot, aes(x = Norting_Y, y = diff)) + geom_point() + ylim(-1, 1)

# inside vs outside
dbf <- foreign::read.dbf("validation_basins_position_info.dbf", as.is = FALSE)
dbf <- dbf[, c("basin", "inBasin")]
names(dbf)[1] <- c("Basin")
in_out <- merge(compare_plot2, dbf, all.x = TRUE)
aggregate(NSE ~ Run + inBasin, in_out, median)
# test if inside/outside has the same NSE distribution
kruskal.test(in_out[in_out$Run == "Run 2", "NSE"], 
       in_out[in_out$Run == "Run 2", "inBasin"])
kruskal.test(in_out[in_out$Run == "Zink et al. (2017)", "NSE"], 
             in_out[in_out$Run == "Zink et al. (2017)", "inBasin"])

ks.test(in_out[in_out$Run == "Run 2" & in_out$inBasin == 0, "NSE"], "pnorm")
ks.test(in_out[in_out$Run == "Run 2" & in_out$inBasin == 1, "NSE"], "pnorm")


t.test(in_out[in_out$Run == "Run 2" & in_out$inBasin == 0, "NSE"], 
       in_out[in_out$Run == "Run 2" & in_out$inBasin == 1, "NSE"])


in_out <- in_out[in_out$Run %in% c("Run 2", "Zink et al. (2017)"), ]
count <- 0
for(basin in unique(in_out$Basin)){
  best <- in_out[in_out$Basin == basin, "Run"][which.max(in_out[in_out$Basin == basin, "NSE"])]
  if(best == "Run 2") count <- count +1
}

count/length(unique(in_out$Basin))
