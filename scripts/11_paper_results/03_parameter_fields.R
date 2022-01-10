## Code to calculate empirical probability of snowfall for every day of year based on 
# mean temperature and rainfall threshold

# load libraries
library(raster)
library(ggplot2)
library(rasterVis)
library(magrittr)
library(RColorBrewer)
library(viridis)  # better colors for everyone
library(ggthemes) # theme_map()

if(Sys.info()[1] == "Linux"){
  # path to Tiff
  tiffpath <- "/mnt/Data/Dropbox/Diss/FSO_mHM/Results/parameter_maps" #"D:/FSO/WRR_Data/parameter_maps/"
  # path to outputs
  outpath <- "/mnt/Data/Dropbox/Diss/FSO_mHM/2020_fso_mhm/scripts/11_paper_results/maps/"
} else {
  # path to Tiff
  tiffpath <- "C:/Users/morit/Dropbox/Diss/FSO_mHM/Results/parameter_maps" #"D:/FSO/WRR_Data/parameter_maps/"
  # path to outputs
  outpath <- "C:/Users/morit/Dropbox/Diss/FSO_mHM/2020_fso_mhm/scripts/11_paper_results/maps/"
  
}
# get file names
Tiff_files <- list.files(path = tiffpath, pattern = '\\.tiff$', full.names = TRUE)

# color palette
pal <- c(rev(as.character(yarrr::piratepal("basel", length.out = 7)))[1:6], "#FEC10BFF")
pal <- RColorBrewer::brewer.pal(n = 7, "Dark2")
pal6 <- pal[1:6]
pal5 <- pal[1:5]
pal2 <- pal[c(2, 7)]

# 1. assess variable cascades in each run ------------------------------------------------

# Get limits
ksat_vars <- c("KSat notill", "VarKSat horizontal relative", 
               "VarKSat vertical relative", "SoilMoistureSaturationDeficit", "L1 Alpha", "L1 Kperco")
fc_vars <- c("FieldCap notill", "L1 SatSoilMoisture", "L1 FieldCap")
ranges <- data.frame(var = c(ksat_vars, fc_vars), low = NA, high = NA)
tiff_vars <- Tiff_files[grep(paste0("run_", 1), Tiff_files, fixed = TRUE)] %>% 
  sapply(., function(x) tail(unlist(strsplit(x, "/")), 1)) %>% 
  gsub(".tiff", "", ., fixed = TRUE) %>% 
  sapply(., function(x) paste(unlist(strsplit(x, "_"))[-c(1:5)], collapse = " ")) 
names(tiff_vars) <- NULL

for(run in c(2:4, 6)){
  if(run == 6) {
    run_files <- Tiff_files[grep("mhm_standard_parameter", Tiff_files, fixed = TRUE)]
  } else {
    run_files <- Tiff_files[grep(paste0("run_", run), Tiff_files, fixed = TRUE)]
  }
  for(i in seq_along(run_files)){
    map_values <- values(raster(run_files[i]))
    ranges$low[ranges$var == tiff_vars[i]] <- min(
      ranges$low[ranges$var == tiff_vars[i]],
      min(map_values, na.rm = TRUE), na.rm = TRUE)
    ranges$high[ranges$var == tiff_vars[i]] <- max(
      ranges$high[ranges$var == tiff_vars[i]],
      max(map_values, na.rm = TRUE), na.rm = TRUE)
  }
}

ranges$low <- floor(ranges$low)
ranges$high <- ceiling(ranges$high)
#ranges$high[ranges$var == "KSat notill"] <- 250
#ranges$high[ranges$var == "VarKSat horizontal relative"] <- 5
#ranges$high[ranges$var == "VarKSat vertical relative"] <- 5

# 
# for(run in 1:5){
#   run_files <- Tiff_files[grep(paste0("run_", run), Tiff_files, fixed = TRUE)]
#   # load tiffs (and aggregate if necessary)
#   tiff_list <- list()
#   for(i in seq_along(run_files)){
#     map <- raster(run_files[i])
#     # remove outliers from ksat 
#     if(grepl("KSat_notill", run_files[i], fixed = TRUE)) map[map>500] <- NA
#     if(grepl("VarKSat_vertical", run_files[i], fixed = TRUE)) {
#       map[map>quantile(values(map), 0.9, na.rm = TRUE)] <- NA
#     }
#     if(grepl("VarKSat_horizontal", run_files[i], fixed = TRUE)) {
#       map[map>quantile(values(map), 0.9, na.rm = TRUE)] <- NA
#     }
#     tiff_list[i] <- map
#   }
#   # variable names
#   for(i in 1:length(tiff_list)){
#     tiff_list[[i]]@data@names <- tiff_vars[i]
#   }
# 
#   # Plot Ksat cascade
#   mfgs <- list(c(1, 2), c(2, 1), c(2, 2), c(2, 3), c(3, 1), c(3, 3))
#   png(file = paste0(outpath, "KSat_cascade_run_", run, ".png"),
#       units = "cm", res = 700, pointsize = 20, width=48, height=48)
#   par(mar = c(2, 2, 3, 5))
#   par(mfrow = c(3, 3))
#   plot.new( )
#   pal <- RColorBrewer::brewer.pal('BrBG', n=11)
# 
#   
#   
#   for(k in seq_along(ksat_vars)){
#     # color limits and breaks
#     range_sub <- ranges[ranges$var == ksat_vars[k], c("low", "high")]
#     if(range_sub$high > 10){
#     breaks <- as.integer(seq(range_sub$low, range_sub$high, length.out = 11))
#     } else {
#       breaks <- round(seq(range_sub$low, range_sub$high, length.out = 11), 1)
#     }
#     label_set <- breaks[c(1, 3, 5, 7, 9, 11)]
#     if(ksat_vars[k] == "L1 Alpha"){
#       breaks <- c(round(seq(range_sub$low, range_sub$high/2, length.out = 9), 1),
#                   round(seq(range_sub$high/2, range_sub$high, length.out = 3), 1)[-1])
#       label_set <- breaks[c(1, 3, 5, 7, 9, 10, 11)]
#     }
#     arg <- list(at=label_set, labels=label_set)
#     # plot position
#     par(mfg = mfgs[[k]])
#     # plot
#     map_id <- which(tiff_vars == ksat_vars[k])
#     plot_map <- flip(tiff_list[[map_id]], direction = "y")
#     plot_map[plot_map > range_sub$high] <- NA
#     plot(plot_map, main = names(tiff_list[[map_id]]@data),
#          legend.args=list(text = "", side=4, font=2, line=2.5, 
#                           cex=2), 
#          breaks = breaks, 
#          col = pal,
#          axis.args=arg)
#     
#   }
#   dev.off()
# 
#   
#   # Plot FieldCap cascade
#   mfgs <- list(c(1, 2), c(2, 1), c(2, 2))
#   png(file = paste0(outpath, "FieldCap_cascade_run_", run, ".png"),
#       units = "cm", res = 700, pointsize = 15, width=30, height=34)
#   par(mar = c(2, 3, 3, 5))
#   par(mfrow = c(2, 2))
#   for(k in seq_along(fc_vars)){
#     # color limits and breaks
#     range_sub <- ranges[ranges$var == fc_vars[k], c("low", "high")]
#     if(range_sub$high > 10){
#       breaks <- as.integer(seq(range_sub$low, range_sub$high, length.out = 11))
#     } else {
#       breaks <- round(seq(range_sub$low, range_sub$high, length.out = 11), 1)
#     }    
#     label_set <- breaks[c(1, 3, 5, 7, 9, 11)]
#     arg <- list(at=label_set, labels=label_set)
#     
#     par(mfg = mfgs[[k]])
#     map_id <- which(tiff_vars == fc_vars[k])
#     plot(flip(tiff_list[[map_id]], direction = "y"), main = names(tiff_list[[map_id]]@data),
#          legend.args=list(text = "", side=4, font=2, line=2.5, cex=2),
#          breaks = breaks, 
#          col = pal,
#          axis.args=arg)
#     
#   }
#   dev.off()
# }
# 

# 2. Compute parameter distributions -----------------------------------------------------
f1 <- Tiff_files[grep("FieldCap_notill", Tiff_files, fixed = TRUE)]
f2 <- Tiff_files[grep("KSat_notill", Tiff_files, fixed = TRUE)]
f3 <- Tiff_files[grep("L1_FieldCap", Tiff_files, fixed = TRUE)]
f4 <- Tiff_files[grep("L1_Kperco", Tiff_files, fixed = TRUE)]
f5 <- Tiff_files[grep("L1_Alpha", Tiff_files, fixed = TRUE)]
# sort correctly and select only Run 2-4 and standard parameters
f1 <- f1[order(sapply(f1, function(x) tail(unlist(strsplit(x, "_")), 3)[1]))]#[c(2:4, 6)]
f2 <- f2[order(sapply(f2, function(x) tail(unlist(strsplit(x, "_")), 3)[1]))]#[c(2:4, 6)]
f3 <- f3[order(sapply(f3, function(x) tail(unlist(strsplit(x, "_")), 3)[1]))]#[c(2:4, 6)]
f4 <- f4[order(sapply(f4, function(x) tail(unlist(strsplit(x, "_")), 3)[1]))]#[c(2:4, 6)]
f5 <- f5[order(sapply(f5, function(x) tail(unlist(strsplit(x, "_")), 3)[1]))]#[c(2:4, 6)]

fieldcap <- flip(stack(f1), direction = "y")
ksat <- flip(stack(f2), direction = "y")
L1_FieldCap_4km <- flip(stack(f3), direction = "y")
L1_Kperco_4km <- flip(stack(f4), direction = "y")
L1_Alpha_4km <- flip(stack(f5), direction = "y")
# remove outliers from KSat of run5
#aggregate from 100x100 resolution to 4000x4000 (factor = 40)
fieldcap_4km <- fieldcap#aggregate(fieldcap, fact = 40)
ksat_4km <- ksat#aggregate(ksat, fact = 40)


run_values <- function(variable){
  map_values <- NULL
  for(run in 1:6){
    map <- subset(get(paste0(variable, "_4km")), run)
    data <- values(map)
    data <- data[!is.na(data)]
    if(run != 6){
      map_values <- rbind(map_values,
                          data.frame("Run" = paste0("Run ", run), 
                                     values = data))
    } else {
      map_values <- rbind(map_values,
                          data.frame("Run" = "default mHM", 
                                     values = data))
    }
  }
  return(map_values)
}
ksat_values <- run_values("ksat")
aggregate(values ~Run, ksat_values, summary)
fc_values <- run_values("fieldcap")
aggregate(values ~Run, fc_values, summary)



density_plotter <- function(variable, var_text, ylim, xlim = c(0, NA)){
  map_values <- NULL
  for(run in c(2:4, 6)){
    map <- subset(get(paste0(variable, "_4km")), run)
    data <- values(map)
    data <- data[!is.na(data)]
    if(run != 6){
      map_values <- rbind(map_values,
                          data.frame("Run" = paste0("Run ", run), 
                                     values = data))
    } else {
      map_values <- rbind(map_values,
                          data.frame("Run" = "default mHM", 
                                     values = data))
    }
  }
  if(Sys.info()[1] == "Linux"){
    path <- "/mnt/Data/Dropbox/Diss/FSO_mHM/2020_fso_mhm/scripts/11_paper_results/plots/"
  } else {
    path <- "C:/Users/morit/Dropbox/Diss/FSO_mHM/2020_fso_mhm/scripts/11_paper_results/plots/"
  }
  map_values$Run <- factor(map_values$Run, levels = c("Run 2", "Run 3", "Run 4", "default mHM"))
  # NSE densities
  p <- ggplot(map_values, 
              aes(x = values, fill = Run)) + 
    geom_density(alpha = 0.7) + 
    scale_fill_manual(values = RColorBrewer::brewer.pal(n = 4, "Dark2"))+#pal[c(2:4, 7)]) +
    theme_minimal() + 
    theme(axis.text = element_text(size=18), 
          axis.title=element_text(size=20),
          legend.text = element_text(size=18),
          legend.title = element_text(size=20)) +
    labs(fill = "", x = var_text, y = "f(x)") + 
    theme(legend.title = element_blank()) +
    coord_cartesian(ylim=c(0, ylim)) +
    xlim(xlim) 
  ggsave(paste0(path, "07_", variable, "_notill_densities.png"), 
         width = 10, height = 5, units = "in", dpi = 600) 
  
  # grey single densities
  for(run in unique(map_values$Run)){
  p <- ggplot(map_values[map_values$Run == run, ], 
              aes(x = values)) + 
    geom_density(alpha = 0.7, fill = "grey") + 
    #scale_fill_manual(values = RColorBrewer::brewer.pal(n = 4, "Dark2"))+#pal[c(2:4, 7)]) +
    theme_minimal() + 
    theme(axis.text = element_text(size=18), 
          axis.title=element_text(size=20),
          legend.text = element_text(size=18),
          legend.title = element_text(size=20)) +
    labs(fill = "", x = var_text, y = "f(x)") + 
    theme(legend.title = element_blank()) +
    coord_cartesian(ylim=c(0, ylim)) +
    xlim(xlim) 
  ggsave(paste0(path, "07_", variable, run, "_notill_densities.png"), 
         width = 6, height = 5, units = "in", dpi = 600) 
  
  }
  
  
}

density_plotter("ksat", "Saturated hydraulic conductivity (cm/day)", 0.075, xlim = c(0, 250))
density_plotter("fieldcap", "Field capacity (-)", 200, xlim = c(0, 0.3))
density_plotter("L1_FieldCap", "L1 Field capacity (mm)", 0.1)
density_plotter("L1_Kperco", "L1 Kperco ()", 0.5)
density_plotter("L1_Alpha", "L1 Alpha ()", 30, xlim = c(0, 1))

# 3. Plot Parameter maps and distribution of runs ----------------------------------------

var_maps_plotter <- function(var, limits){
  #var_map <- get(paste0(var, "_4km"))
  var_map <- subset(get(paste0(var, "_4km")), c(2:4, 6))
  
  if(class(var_map) != "RasterBrick") var_map <- brick(var_map)
  var_map@data@names <- c(paste0("Run ", 2:4), "Standard parameter")
  myPal <- RColorBrewer::brewer.pal('BrBG', n=9)
  myTheme <- rasterTheme(region = myPal)
  
  cat("limits:", limits)
  # plot var for all runs
  mains <- paste0(var, " run 2 to 4 and standard parameter run")
  flab <- paste0(outpath, "01_parameters_single_runs_", var, ".png")
  png(file = flab, units = "cm", res = 450, pointsize = 12,
      width=20.5, height=12.7)
  print(levelplot(var_map,
                  scales = list(x = list(cex =.7), y = list(cex =.7), xlab=list(cex=.7)),
                  ylab=list(cex=.25), main=list(label=mains, cex = 1.1), sub = list(cex = 0.57),
                  colorkey = list(labels=list(cex =.6)),
                  par.settings = myTheme, at = seq(0, limits, length.out=100)))
  dev.off()
  
  # plot var distribution
  var_map_sub <- stack(subset(var_map, 1, drop = FALSE), 
                       subset(var_map, 2, drop = FALSE), 
                       subset(var_map, 3, drop = FALSE))
  var_map_agg <- stack(c("mean" = mean(var_map_sub), 
                         "max" = max(var_map_sub), 
                         "min" = min(var_map_sub), 
                         "range" = max(var_map_sub) - min(var_map_sub)))
  mains <- paste0(var, " run 2 to 4")
  flab <- paste0(outpath, "02_parameters_var_dist_", var, ".png")
  png(file = flab, units = "cm", res = 450, pointsize = 12,
      width=20.5, height=12.7)
  print(levelplot(var_map_agg,
                  scales = list(x = list(cex =.7), y = list(cex =.7), xlab=list(cex=.7)),
                  ylab=list(cex=.25), main=list(label=mains, cex = 1.1), sub = list(cex = 0.57),
                  colorkey = list(labels=list(cex =.6)),
                  par.settings = myTheme, at = seq(0, limits, length.out=100)))
  dev.off()
  
  # # plot var distribution run 2-4
  # best_runs <- stack(subset(var_map, 2, drop = FALSE), 
  #                    subset(var_map, 3, drop = FALSE), 
  #                    subset(var_map, 4, drop = FALSE))
  # mains <- paste0(var, " run 2 to 4")
  # best_map_agg <- stack(c("mean" = mean(best_runs),  
  #                         "max" = max(best_runs), 
  #                         "min" = min(best_runs),
  #                         "range" = max(best_runs) - min(best_runs)))
  # 
  # flab <- paste(outpath, "03_parameters_var_dist_run_1_to_4_", var, ".png", sep = "")
  # png(file = flab, units = "cm", res = 450, pointsize = 12,
  #     width=20.5, height=12.7)
  # print(levelplot(best_map_agg,
  #                 scales = list(x = list(cex =.7), y = list(cex =.7), xlab=list(cex=.7)),
  #                 ylab=list(cex=.25), main=list(label=mains, cex = 1.1), sub = list(cex = 0.57),
  #                 colorkey = list(labels=list(cex =.6)),
  #                 par.settings = myTheme, at = seq(0, limits, length.out=100)))
  # dev.off()
}

var_maps_plotter("ksat", ranges[ranges$var == "KSat notill", "high"])
var_maps_plotter("fieldcap", 0.3)#ranges[ranges$var == "FieldCap notill", "high"])
var_maps_plotter("L1_FieldCap", ranges[ranges$var == "L1 FieldCap", "high"])
var_maps_plotter("L1_Kperco", ranges[ranges$var == "L1 Kperco", "high"])
var_maps_plotter("L1_Alpha", 1)#ranges[ranges$var == "L1 Alpha", "high"])



tiff_plotter <- function(var, limits){
  #var_map <- get(paste0(var, "_4km"))
  var_map <- subset(get(paste0(var, "_4km")), c(2:4, 6))
  
  if(class(var_map) != "RasterBrick") var_map <- brick(var_map)
  var_map@data@names <- c(paste0("Run ", 2:4), "Standard parameter")
  for(plot_id in c(1:4)){
    test <- subset(var_map, plot_id)
    test_spdf <- as(test, "SpatialPixelsDataFrame")
    test_df <- as.data.frame(test_spdf)
    colnames(test_df) <- c("value", "x", "y")
    tiff(paste0("C:/Users/morit/Dropbox/Diss/FSO_mHM/Manuscript/Figures/figure basis",
                "/para_map_", var, "_", plot_id, ".tiff"), compression = "lzw",
         width = 10, height = 5, units = "in", res = 200)
    p1 <- ggplot() +  
      geom_tile(data=test_df, aes(x=x, y=y, fill=value), alpha=0.8) + 
      scale_fill_fermenter(palette = "Spectral") +
      #scale_fill_viridis() +
      coord_equal() +
      theme_map() +
      theme(legend.position="none") +
      theme(legend.key.width=unit(2, "cm"))
    dev.off()
    ggsave(paste0("C:/Users/morit/Dropbox/Diss/FSO_mHM/Manuscript/Figures/figure basis",
                  "/para_map_", var, "_", plot_id, ".svg"), p1, compression = "lzw")
  }
}
tiff_plotter("ksat", ranges[ranges$var == "KSat notill", "high"])
tiff_plotter("fieldcap", 0.5)#ranges[ranges$var == "FieldCap notill", "high"])
tiff_plotter("L1_FieldCap", ranges[ranges$var == "L1 FieldCap", "high"])
tiff_plotter("L1_Kperco", ranges[ranges$var == "L1 Kperco", "high"])
tiff_plotter("L1_Alpha", 1)#ranges[ranges$var == "L1 Alpha", "high"])




# 4. Plot runoff components --------------------------------------------------------------
run_names <- c(paste0("FSO_SCE_NSE_run_", 1:4), "FSO_DDS_NSE_run_5", "mhm_standard_parameter")
tiffpath2 <- "C:/Users/morit/Dropbox/Diss/FSO_mHM/Results/flux_states/" 

# Plot runoff components
runoff_vars <- c("QB", "QD", "QIf", "QIs")
for(run in 1:6){
  run_files <- paste0(tiffpath2, run_names[run], "_", runoff_vars, ".tiff")
  # load tiffs (and aggregate if necessary)
  tiff_list <- list()
  for(i in seq_along(run_files)){
    tiff_list[i] <- raster(run_files[i])
  }
  # variable names
  tiff_vars <- sapply(run_files, function(x) tail(unlist(strsplit(x, "/")), 1)) %>% 
    gsub(".tiff", "", ., fixed = TRUE) %>% 
    sapply(., function(x) paste(unlist(strsplit(x, "_"))[-c(1:5)], collapse = " ")) 
  names(tiff_vars) <- NULL
  for(i in 1:length(tiff_list)){
    tiff_list[[i]]@data@names <- tiff_vars[i]
  }
  # Plot runoff components
  var_map <- stack(c("QB" = flip(tiff_list[[1]], direction = "y"), 
                     "QD" = flip(tiff_list[[2]], direction = "y"), 
                     "QIf" = flip(tiff_list[[3]], direction = "y"), 
                     "QIs" = flip(tiff_list[[4]], direction = "y")))
  
  myPal <- RColorBrewer::brewer.pal('BrBG', n=9)
  myTheme <- rasterTheme(region = myPal)
  mains <- paste0("Runoff components of run ", run)
  flab <- paste(outpath, "04_runoff_components_run", run, ".png", sep = "")
  if(run == 6) {
    mains <- "Runoff components of run with standard parameter"
    flab <- paste(outpath, "04_runoff_components_standard_parameter.png", sep = "")
  }
  png(file = flab, units = "cm", res = 450, pointsize = 12,
      width=20.5, height=12.7)
  print(levelplot(var_map,
                  scales = list(x = list(rot = 90, cex =.7), y = list(cex =.7), xlab=list(cex=.7)),
                  ylab=list(cex=.25), 
                  main=list(label=mains, cex = 0.7), 
                  sub = list(cex = 0.57),
                  colorkey = list(labels=list(cex =.6)),
                  par.settings = myTheme, 
                  at = seq(0, 1, length.out=100),
                  labels = TRUE,
                  margin = FALSE))
  dev.off()
  
}

# 5. SWE fractions -----------------------------------------------------------------------
run_names <- c(paste0("FSO_SCE_NSE_run_", 1:4), "FSO_DDS_NSE_run_5", "mhm_standard_parameter")
tiffpath2 <- "C:/Users/morit/Dropbox/Diss/FSO_mHM/Results/flux_states/" 
swe_vars <- c("SWC_fraction_01", "SWC_fraction_02", "SWC_fraction_03")

# Plot runoff components
for(run in 1:6){
  run_files <- paste0(tiffpath2, run_names[run], "_", swe_vars, ".tiff")
  # load tiffs (and aggregate if necessary)
  tiff_list <- list()
  for(i in seq_along(run_files)){
    tiff_list[i] <- raster(run_files[i])
  }
  # variable names
  tiff_vars <- sapply(run_files, function(x) tail(unlist(strsplit(x, "/")), 1)) %>% 
    gsub(".tiff", "", ., fixed = TRUE) %>% 
    sapply(., function(x) paste(unlist(strsplit(x, "_"))[-c(1:5)], collapse = " ")) 
  names(tiff_vars) <- NULL
  for(i in 1:length(tiff_list)){
    tiff_list[[i]]@data@names <- tiff_vars[i]
  }
  # Plot runoff components
  var_map <- stack(c("Depth 0 - 0.05 m" = flip(tiff_list[[1]], direction = "y"), 
                     "Depth 0.05 - 0.25 m" = flip(tiff_list[[2]], direction = "y"), 
                     "Depth 0.05 - 0.25 m" = flip(tiff_list[[3]], direction = "y")))
  myPal <- rev(RColorBrewer::brewer.pal('BrBG', n=9))
  myTheme <- rasterTheme(region = myPal)
  mains <- paste0("Fraction of days were SWC < L1_soilMoistFC of run ", run)
  flab <- paste(outpath, "05_SWC_fractions_run", run, ".png", sep = "")
  if(run == 6) {
    mains <- "Fraction of days were SWC < L1_soilMoistFC of run with standard parameter"
    flab <- paste(outpath, "05_SWC_fractions_standard_parameter.png", sep = "")
  }
  png(file = flab, units = "cm", res = 450, pointsize = 12,
      width=20.5, height=12.7)
  print(levelplot(var_map,
                  scales = list(x = list(rot = 90, cex =.7), y = list(cex =.7), xlab=list(cex=.7)),
                  ylab=list(cex=.25), 
                  main=list(label=mains, cex = 0.7), 
                  sub = list(cex = 0.57),
                  colorkey = list(labels=list(cex =.6)),
                  par.settings = myTheme, 
                  at = seq(0, 1, length.out=100),
                  labels = TRUE,
                  margin = FALSE,
                  names.attr=c("Depth 0 - 0.05 m", "Depth 0.05 - 0.25 m", "Depth 0.25 - 2 m")))
  dev.off()
  
}
