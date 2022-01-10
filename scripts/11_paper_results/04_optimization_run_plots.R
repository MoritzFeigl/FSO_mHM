# Compare low iteration vs high iteration runs
library(ggplot2)
if(Sys.info()[1] == "Linux"){
  setwd("/mnt/Data/Dropbox/Diss/FSO_mHM/2020_fso_mhm/scripts/11_paper_results")
  path <- "/mnt/Data/Dropbox/Diss/FSO_mHM/Results/"
} else {
  setwd("C:/Users/morit/Dropbox/Diss/FSO_mHM/2020_fso_mhm/scripts/11_paper_results")
  path <- "C:/Users/morit/Dropbox/Diss/FSO_mHM/Results/"
}

# colors
pal <- c(rev(as.character(yarrr::piratepal("basel", length.out = 7)))[1:6], "#FEC10BFF")
pal <- RColorBrewer::brewer.pal(n = 7, "Dark2")
pal6 <- pal[1:6]
pal5 <- pal[1:5]
pal2 <- pal[c(2, 7)]


# 1. compile runs ------------------------------------------------------------------------
run_names <- c(paste0("FSO_SCE_NSE_run_", 1:4), "FSO_DDS_NSE_run_5")


# load TF estimates, performance results
compare <- NULL
results <- NULL
for(run in run_names){
  # basin wise metrics
  result_tracker <- cbind("Run" = run, "Basins" = "Training", "time_period" = "Training",
                          read.csv(paste0(path, run, "/result_tracker.csv")))
  cat(run, ":", nrow(result_tracker), "\n")
  result_tracker$iteration <- as.integer(rownames(result_tracker))
  # Compare Iteration 2000 with best Iteration
  dim(result_tracker)
  old_best <- which(round(result_tracker$loss, 3) == result_tracker[2000, "best_loss"])[1]
  compare_df <- rbind(cbind("iter" = old_best, result_tracker[old_best, ]),
                      cbind("iter" = which.max(result_tracker$loss), 
                            result_tracker[which.max(result_tracker$loss), ]))
  
  get_diffs <- function(col){
    if(class(col[1]) == "character"){
      return(col[1] == col[2])
    } else {
      return( round(((col[2] - col[1])/col[1]) * 100, 1))
    }
  }
  compare_df <- rbind(compare_df, 
                      sapply(compare_df, get_diffs))
  compare_df$iter[3] <- compare_df$iter[2] -compare_df$iter[1]
  compare_df$Run[3] <- paste0(run, "_diff")
  compare <- rbind(compare, compare_df)
  results <- rbind(results, result_tracker)
  cat(run, "\nIterations since improvement:", tail(result_tracker[, 29], 1), "\n\n")
}
write.csv(compare, "comparing_run2000_with_best_runs.csv", row.names = FALSE)
results$Run <- factor(results$Run,
                      levels = run_names,
                      labels = paste0("Run ", 1:5))
results <- results[results$n_iterations_since_improvement <= 1000 | results$iteration < 3000, ] 


# read numerics optimization results and create iterations vs loss data frame
results_full <- NULL
for(run_id in 1:5){
  run <- paste0("Run ", run_id)
  run_path <- run_names[run_id]
  numerics_tracker <- cbind("Run" = run, "Basins" = "Training", "time_period" = "Training",
                            read.csv(paste0(path, run_path, "/numerics_tracker.csv")))
  numerics_tracker <- numerics_tracker[, c("Run", "wmulti_loss")]
  names(numerics_tracker)[2] <- "best_loss"
  tmp <- results[results$Run == run, c("Run", "best_loss", "iteration")]
  max_it <- max(tmp$iteration)
  numerics_tracker <- cbind(numerics_tracker, 
                            "iteration" = (max_it+1):(max_it+nrow(numerics_tracker)))
  if(numerics_tracker$best_loss[1] < tail(tmp$best_loss, 1)){
    numerics_tracker$best_loss[1] <- tail(tmp$best_loss, 1)
  }
  for(i in 2:nrow(numerics_tracker)){
    numerics_tracker$best_loss[i] <- max(numerics_tracker$best_loss[i], numerics_tracker$best_loss[i-1])
  }
  
  results_full <- rbind(results_full,
                        rbind(cbind(tmp, "opt"="initial"), 
                              cbind(numerics_tracker, "opt"="numeric"))
  )
}

ggplot(results_full, aes(iteration, best_loss, col = Run)) + 
  geom_line(size = 1.2) + 
  ylim(0.25, 0.8) +
  labs(y = "Loss", x = "Iterations", col = "") + 
  theme_bw() + 
scale_color_manual(values = pal6, name = "") +
  theme(axis.text = element_text(size=20), 
        axis.title=element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))
ggsave("plots/runs_loss_developement.png", height = 7, width = 9, units = "in")


final_result <- aggregate(best_loss ~ Run + opt, results_full[results_full$opt == "numeric",], max)
final_result <- merge(final_result,
                                aggregate(iteration ~ Run, 
                                          results_full[results_full$opt == "initial",], max))
#final_result <- merge(final_result, results_full, all.x=TRUE, all.y=FALSE)
#final_result <- aggregate(iteration ~ Run + best_loss + opt, final_result, max)[, c(1, 2, 4, 3)]
ggplot(results_full[results_full$opt == "initial",], aes(iteration, best_loss, col = Run)) + 
  geom_line(size = 1.2, linetype="solid") +
  geom_point(data = final_result, aes(iteration, best_loss, col = Run, shape = "opt"), size = 2) +
  scale_shape_manual("", values = 16, labels = "Optimized \ncoefficients") +
  ylim(0.25, 0.8) +
  labs(y = "Loss", x = "Iterations", col = "") + 
  theme_bw() + 
  scale_color_manual(values = pal6, name = "") +
  theme(axis.text = element_text(size=20), 
        axis.title=element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20)) + 
  guides(color = guide_legend(override.aes = list(size = 1, shape = NA)),
         shape = guide_legend(override.aes = list(size = 3)))
ggsave("plots/runs_loss_developement_v2.png", height = 7, width = 9, units = "in")


best_runs <- compare[c(2, 5, 8, 11, 14), ]
# PCA of FSO dimensions
library(factoextra)
results_pos <- results[,  paste0("FSO", 1:12)]
fso_pos <- best_runs[, paste0("FSO", 1:12)]
rownames(fso_pos) <- paste0("Run ", 1:5)
rownames(fso_pos) <- apply(cbind(best_runs$best_KSat, best_runs$best_FieldCap), 1, 
                           function(x) paste0("KSat=", x[1], "\n", "FC=", x[2]))
pos <- rbind(fso_pos, results_pos)
res.pca <- prcomp(pos)
fviz_eig(res.pca)
fviz_pca_ind(res.pca,
             select.ind = list(name = rownames(fso_pos)), #paste0("Run ", 1:5)),
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             legend.title = "Quality of\nrepresentation"
) 
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)

