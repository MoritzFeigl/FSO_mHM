# Train FSO VAE
# FSO-mHM project
# Moritz Feigl, Sep 2020


library(tensorflow)
if(Sys.info()["sysname"] == "Linux"){
  Sys.setenv(RETICULATE_PYTHON="/home/cfgrammar/.virtualenvs/r-tensorflow/bin/python")
} else {
  Sys.setenv(RETICULATE_PYTHON="C:/Users/morit/anaconda3/envs/tf-gpu/python.exe")
}
library(keras)
library(tfdatasets)
#library(tfautograph)
library(tfprobability)
library(glue)
library(magrittr)
library(purrr, warn.conflicts = FALSE)

tf$test$is_gpu_available()
## KSat or FieldCap
#for(var in c("KSat", "FieldCap")){
var <- "KSat"

if(Sys.info()[["nodename"]] == "cfgrammar"){
  setwd(paste0("/mnt/Data/Dropbox/Diss/FSO_mHM/FSO_Data/VAE_data/", var, "/batches"))
} else {
  if(Sys.info()["sysname"] == "Linux"){
    setwd(paste0("/home/lv71468/mfeigl/FSO_mHM/FSO_Data/VAE_data/", var, "/batches"))
  } else {
    setwd(paste0("C:/Users/morit/Dropbox/Diss/FSO_mHM/FSO_Data/VAE_data/", var, "/batches"))
  }  
}
#strategy <- tf$distribute$MirroredStrategy()

# Hyperparameters-------------------------------------------------------------------------
batch_size <- 250L #500L
# Encoder
embedding_length <- 10L
latent_dim <- 6L
lstm_units <- 512L
inter_dim <- 100L
KL_weight <- batch_size
dist_loss_weight <- 10
# Decoder
TCN_kernel_sizes <- c(4L, 4L, 4L, 4L, 4L)
TCN_dilation_rates <- c(1L, 2L, 4L, 8L, 8L)
TCN_units <- 100L

# Distribution encoder/decoder
TCN_kernel_sizes2 <- c(2L, 2L, 2L)
TCN_dilation_rates2 <- c(1L, 2L, 4L)
TCN_units2 <- 30L

# Load Datasets and VAE model ------------------------------------------------------------
# source("../../../../2020_fso_mhm/scripts/07_FSO_VAE/02_load_data.R")
# source("../../../../2020_fso_mhm/scripts/07_FSO_VAE/03_FSO_VAE.R")
source("/mnt/Data/Dropbox/Diss/FSO_mHM/2020_fso_mhm/scripts/07_FSO_VAE/02_load_data.R")
source("/mnt/Data/Dropbox/Diss/FSO_mHM/2020_fso_mhm/scripts/07_FSO_VAE/03_FSO_VAE.R")

# define checkpoints
if (!file.exists("../checkpoints")) dir.create("../checkpoints")
checkpoints <- callback_model_checkpoint(
  paste0("../checkpoints/10WeightedDist_weights_FINAL.{epoch:02d}-{val_loss:.2f}.hdf5"),
  save_weights_only = TRUE)

# get steps per epoch
example_file <- list.files("training", full.names = TRUE)[1]
length_check <- read.csv(example_file)
steps_per_epoch <- floor(nrow(length_check) * length(list.files("training"))/batch_size)
validation_steps  <- floor(nrow(length_check) * length(list.files("validation"))/batch_size)
rm(length_check)

# check model prediction function
# source("../../../../2020_fso_mhm/scripts/07_FSO_VAE/04_check_model_prediction.R")
source("/mnt/Data/Dropbox/Diss/FSO_mHM/2020_fso_mhm/scripts/07_FSO_VAE/04_check_model_prediction.R")
if (!file.exists("../training_performance")) dir.create("../training_performance")

epochs <- 12
tensorflow::tf$random$set_seed(678)
history <- vae_model %>%
  fit(train_dataset,
      epochs = epochs,
      validation_data = val_dataset, 
      steps_per_epoch = steps_per_epoch,
      validation_steps = validation_steps,
      callbacks = checkpoints
  )

training_results <- do.call(cbind, history$metrics) %>% round(digits = 3)
write.csv(training_results, 
          paste0("../training_history_", format(Sys.time(), "%Y-%m-%d"), ".csv"), 
          row.names = FALSE)

check_model(n = 30, 
            dict = dictionaries[[var]], 
            epoch = epochs, 
            distribution_scale_parameters = "../distribution_scale_parameters.csv", 
            save_in = "../training_performance")
#}




# Load model weights (optional)
list.files("../checkpoints", pattern = "weights")
vae_model <- vae_model %>% 
  load_model_weights_hdf5("../checkpoints/10WeightedDist_weights_2_.11-13070.37.hdf5")



