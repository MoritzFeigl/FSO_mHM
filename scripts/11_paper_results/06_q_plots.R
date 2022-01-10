setwd("C:/Users/morit/Dropbox/Diss/FSO_mHM/Results/mhm_standard_parameter")
library(ggplot2)
q <- read.csv("small_basins_predictions.csv")
q2 <- read.csv("../FSO_SCE_NSE_run_2/small_basins_predictions.csv")

head(q)
id <- "9316479" # 9316479
q_sub <- cbind(date = q$date, q[, grep(id, names(q))],
               q2[, grep(id, names(q2))])
names(q_sub)[2:5] <- c(paste0(paste0(c("Qobs_", "Qsim_"), id), "_default"),
                       paste0(paste0(c("Qobs_", "Qsim_"), id),  "_run2"))
q_sub <- q_sub[, -4]
q_plot <- reshape2::melt(q_sub, id.vars = "date")
q_plot <- q_plot[!is.na(q_plot$value), ]
q_plot$date <- as.POSIXct(q_plot$date)
start_date <- q_plot$date[nrow(q_plot) - 365*10]
q_plotq_plot <- q_plot[q_plot$date > start_date, ]
ggplot(q_plot, aes(date, value, col = variable, linetype = variable)) +
  geom_line() + theme_bw()

hydroGOF::N
ggplot(q_plot, aes(date, value, col = variable)) +
  geom_line() + theme_bw()



