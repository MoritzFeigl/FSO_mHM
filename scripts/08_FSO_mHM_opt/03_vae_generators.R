# VAE generators for FSO
# FSO-mHM project
# Moritz Feigl, Sep 2020

# Libraries
.libPaths("/home/lv71468/mfeigl/R/x86_64-pc-linux-gnu-library/4.0")
library(tensorflow)
library(keras)
#library(tfdatasets)
library(tfprobability)
library(magrittr)
#library(purrr, warn.conflicts = FALSE)

# Hyperparameters
batch_size <- 250L #500L
embedding_length <- 10L
latent_dim <- 6L
lstm_units <- 512L
inter_dim <- 100L
KL_weight <- batch_size
dist_loss_weight <- 10
TCN_kernel_sizes <- c(4L, 4L, 4L, 4L, 4L)
TCN_dilation_rates <- c(1L, 2L, 4L, 8L, 8L)
TCN_units <- 100L
TCN_kernel_sizes2 <- c(2L, 2L, 2L)
TCN_dilation_rates2 <- c(1L, 2L, 4L)
TCN_units2 <- 30L

# dictionaries
dict_from_list <- function(...){
  dictionary <- c(...)
  return(dictionary)
}
f <- c('exp', 'log10', 'log', 'sin', 'cos', 'tan', 'abs', 'acos',
       'asin', 'atan', 'cosh', 'sinh', 'sqrt')
var_KSat <- c("bd", "sand", "clay", "slope", "aspect", "dem")
var_FieldCap <- c("bd", "sand", "clay", "slope", "aspect", "dem", "ThetaS", "KSat", "vGenu_n")
operators <- c('+', '-', '*', '/', "^", "(", ")")
numbers <- seq(0.05, 3, 0.05) %>% format(nsmall=2)
dictionary_KSat <- dict_from_list(var_KSat, f, operators, numbers)
dictionary_FieldCap <- dict_from_list(var_FieldCap, f, operators, numbers)
dictionaries <- list(KSat = dictionary_KSat,
                     FieldCap = dictionary_FieldCap)

# VAE parts
vae_loss <- function(y_true, y_pred){
  k_sum(k_sparse_categorical_crossentropy(y_true, y_pred))
}
vae_metric <- function(y_true, y_pred){
  k_mean(k_sum(k_sparse_categorical_crossentropy(y_true, y_pred), axis = 2))
}
dist_loss <- function(y_true, y_pred){
  k_sum(abs(y_true - y_pred)) * dist_loss_weight
}
# custom highway layer
highway_class <- R6::R6Class(
  "KerasLayer",
  inherit = KerasLayer,
  
  public = list(
    size = NULL,
    activation = NULL,
    W_T = NULL,
    b_T = NULL,
    W = NULL,
    b = NULL,
    H = NULL,
    Tr = NULL,
    initialize = function(size, activation) {
      self$size <- size
      self$activation <- activation
    },
    
    build = function(input_shape) {
      self$W_T <- self$add_weight(
        name = 'W_T',
        shape = list(self$size, self$size),
        initializer = tf$keras$initializers$TruncatedNormal(stddev = 0.1),
        trainable = TRUE
      )
      self$b_T <- self$add_weight(
        name = 'b_T',
        shape = list(self$size),
        initializer = tf$keras$initializers$Constant(0.1),
        trainable = TRUE
      )
      self$W <- self$add_weight(
        name = 'W',
        shape = list(self$size, self$size),
        initializer = tf$keras$initializers$TruncatedNormal(stddev = 0.1),
        trainable = TRUE
      )
      self$b <- self$add_weight(
        name = 'b',
        shape = list(self$size),
        initializer = tf$keras$initializers$Constant(0.1),
        trainable = TRUE
      )
      
    },
    
    call = function(x, mask = NULL) {
      Tr <- tf$keras$activations$sigmoid(k_dot(x, self$W_T) + self$b_T)
      H <- self$activation(tf$keras$backend$dot(x, self$W) + self$b)
      tf$add(tf$multiply(H, Tr), tf$multiply(x, 1-Tr))
    }
  )
)
# Create layer wrapper function
layer_highway <- function(object, size, activation, name = NULL, trainable = TRUE) {
  create_layer(highway_class, object,
               list(
                 size = as.integer(size),
                 activation = activation,
                 name = name,
                 trainable = trainable
               ))
}
# TCN layer
TCN_block <- function(n_outputs, kernel_size,
                      dilation_rate, strides = 1L, dropout = 0.2) {
  keras_model_custom(function(self) {
    self$conv1 <- tf$keras$layers$Conv1D(filters = n_outputs, kernel_size = kernel_size,
                                         dilation_rate = dilation_rate, strides = strides,
                                         padding = "causal", activation = "relu")
    self$bn1 <- tf$keras$layers$LayerNormalization()
    self$drop1 <- layer_dropout(rate = dropout, noise_shape = c(NULL, 1, n_outputs))
    self$conv2 <- tf$keras$layers$Conv1D(filters = n_outputs, kernel_size = kernel_size,
                                         dilation_rate = dilation_rate, strides = strides,
                                         padding = "causal")
    self$bn2 <- tf$keras$layers$LayerNormalization()
    self$drop2 <- layer_dropout(rate = dropout, noise_shape = c(NULL, 1, n_outputs))
    self$downsampler <- tf$keras$layers$Dense(units = as.integer(n_outputs))
    function(inputs, mask = NULL, training = TRUE) {
      x <- inputs %>%
        self$conv1() %>%
        self$bn1(training = training) %>%
        self$drop1(training = training) %>%
        self$conv2() %>%
        self$bn2(training = training) %>%
        self$drop2(training = training)
      if(dim(inputs)[3] != n_outputs){
        inputs <- inputs %>% self$downsampler()
      }
      
      return(tf$nn$relu(x + inputs))
    }
  })
}
index_reconstructor <- function(pred_matrix){
  index <- integer(nrow(pred_matrix))
  for(i in 1:nrow(pred_matrix)) index[i] <- which.max(pred_matrix[i, ])
  return(index)
}
# Structure parameters for TCNs
TCN_filters <- rep(TCN_units, length(TCN_kernel_sizes))
TCN_filters2 <- rep(TCN_units2, length(TCN_kernel_sizes2))

create_generator_vae <- function(dictionary, vae_folder, checkpoint){
  dictionary_length <- as.integer(length(dictionary)) + 1L
  padding_length <- read.csv(paste0(vae_folder, "/padding_length.csv"))[, 3]
  input_shape <- as.integer(padding_length) 
  # function token inputs
  x1 <- layer_input(shape = input_shape, name = "tf_input")
  embedding <- tf$keras$layers$Embedding(input_dim = as.integer(dictionary_length),
                                         output_dim = as.integer(embedding_length),
                                         input_length = as.integer(padding_length),
                                         trainable = TRUE, name = "embedding")
  # function encoder
  bilstm <- tf$keras$layers$Bidirectional(
    layer = layer_lstm(units = lstm_units,
                       return_sequences = TRUE,
                       return_state = FALSE, 
                       recurrent_activation = "sigmoid"), merge_mode = "sum")
  highway1 <- tf$keras$layers$Dense(units = inter_dim, activation = "selu")
  highway2 <- layer_highway(size = inter_dim, activation = tf$nn$selu)
  highway3 <- layer_highway(size = inter_dim, activation = tf$nn$selu)
  highway4 <- layer_highway(size = inter_dim, activation = tf$nn$selu)
  x1_processed <- x1 %>%
    layer_flatten(input_shape = input_shape) %>%
    embedding %>% 
    bilstm %>%
    layer_flatten() %>%
    highway1 %>%
    highway2 %>%
    highway3 %>%
    highway4
  # function distribution input
  x2 <- layer_input(shape = 9L, name = "dist_input")
  layer_norm <- tf$keras$layers$LayerNormalization()
  dist_encoder1 <- tf$keras$layers$Dense(units = 9 * 5,
                                         activation = "selu")
  dist_reshape1 <- tf$keras$layers$Reshape(target_shape = c(9L, 5L))
  # distribution encoder
  lstm2 <- tf$keras$layers$Bidirectional(
    layer = layer_lstm(units = 30,
                       return_sequences = TRUE,
                       return_state = FALSE, 
                       recurrent_activation = "sigmoid"), merge_mode = "sum")
  dense_inter <- tf$keras$layers$Dense(units = 100, activation = "selu")
  x2_processed <- x2 %>%
    dist_encoder1 %>%
    dist_reshape1 %>%
    lstm2 %>%
    layer_flatten() %>%
    dense_inter
  # Latent Space
  dist_par <- tf$keras$layers$Dense(units = params_size_multivariate_normal_tri_l(latent_dim))
  nv_dist <- layer_multivariate_normal_tri_l(event_size = latent_dim)
  kl_loss <- layer_kl_divergence_add_loss(
    distribution = tfd_independent(
      tfd_normal(loc = rep(0, latent_dim), scale = 1),
      reinterpreted_batch_ndims =1
    ),
    weight = KL_weight)
  z <- layer_concatenate(list(x1_processed, x2_processed)) %>%
    dist_par %>%
    nv_dist %>%
    kl_loss
  # function string decoder
  decoder1 <- layer_dense(units = padding_length * latent_dim * 2,
                          activation = "selu")
  reshape_z <- layer_reshape(target_shape = c(padding_length, 2* latent_dim))
  TCN_input <- z %>%
    decoder1 %>%
    reshape_z
  TCN_output <- TCN_input
  for(i in 1:length(TCN_kernel_sizes)){
    assign(paste0("tf_block", i),
           TCN_block(n_outputs = TCN_filters[i],
                     kernel_size = TCN_kernel_sizes[i], dilation_rate = TCN_dilation_rates[i]))
    tf_block_tmp <- get(paste0("tf_block", i))
    TCN_output <- TCN_output %>% tf_block_tmp
  }
  fc <- layer_dense(units = dictionary_length, activation = "softmax", name = "tf_output")
  y1 <- TCN_output %>% fc
  # function distribution decoder
  dist_decoder1 <- layer_dense(units = 9*2*latent_dim,
                               activation = "selu")
  dist_reshape_z <- layer_reshape(target_shape = c(9, 2*latent_dim))
  dist_TCN_input <- z %>%
    dist_decoder1 %>%
    dist_reshape_z
  dist_TCN_output <- dist_TCN_input
  for(i in 1:length(TCN_kernel_sizes2)){
    assign(paste0("dist_block", i),
           TCN_block(n_outputs = TCN_filters2[i],
                     kernel_size = TCN_kernel_sizes2[i],
                     dilation_rate = TCN_dilation_rates2[i]))
    dist_block_tmp <- get(paste0("dist_block", i))
    dist_TCN_output <- dist_TCN_output %>% dist_block_tmp
  }
  fc1 <- layer_dense(units = 30, activation = "relu")
  fc2 <- layer_dense(units = 9L, name = "dist_output")
  y2 <- dist_TCN_output %>%
    layer_flatten() %>%
    fc1 %>%
    fc2
  # keras functional model uniting them both
  vae_model <- keras_model(c(x1, x2), c(y1, y2))
  vae_model <- vae_model %>% 
    load_model_weights_hdf5(paste0(vae_folder, "/checkpoints/", checkpoint))
  encoder_model <- keras_model(c(x1, x2), z)
  decoder_input <- layer_input(shape = latent_dim)
  generator <- decoder_input %>%
    decoder1() %>%
    reshape_z() %>%
    tf_block1() %>%
    tf_block2() %>%
    tf_block3() %>%
    tf_block4() %>%
    tf_block5() %>%
    fc()
  decoder_model <- keras_model(decoder_input, generator)
  return(decoder_model)
}

# Load generators ------------------------------------------------------------------------
ksat_generator <- create_generator_vae(dictionary = dictionaries[["KSat"]],
                                       vae_folder = "/home/lv71468/mfeigl/FSO_mHM/FSO_Data/VAE_data/KSat/",
                                       checkpoint = "10WeightedDist_weights_FINAL.01-13060.11.hdf5")
fieldcap_generator <- create_generator_vae(dictionary = dictionaries[["FieldCap"]],
                                       vae_folder = "/home/lv71468/mfeigl/FSO_mHM/FSO_Data/VAE_data/FieldCap/",
                                       checkpoint = "10WeightedDist_weights_FINAL.12-13322.24.hdf5")

# Generator functions --------------------------------------------------------------------
index_reconstructor <- function(pred_matrix){
  index <- integer(nrow(pred_matrix))
  for(i in 1:nrow(pred_matrix)) index[i] <- which.max(pred_matrix[i, ])
  return(index)
}
index_sampler <- function(pred_matrix){
  index <- integer(nrow(pred_matrix))
  for(i in 1:nrow(pred_matrix)) index[i] <- sample(ncol(pred_matrix), 1, prob = pred_matrix[i, ])
  return(index)
}

# helper function that makes the function evalutaion
tf_evaluation <- function(predicted_tf, variables){
  predicted_tf_num <- predicted_tf
  for(i in variables){
    predicted_tf_num <- gsub(i, "1.0", predicted_tf_num)
  }
  tf_eval <- try({
    eval(parse(text = paste('f_test <- function() {' ,  predicted_tf_num , '}', sep='')))
    f_test()
  }, silent = TRUE)
  return(tf_eval)
}
# main function for index prediction
tf_prediction <- function(index_pred, dictionary, variables){
  index_prediction <- index_reconstructor(index_pred)
  NULL_dictionary <- c("", dictionary)
  predicted_tf <- NULL_dictionary[index_prediction] %>% paste0(collapse = "")
  tf_eval <- tf_evaluation(predicted_tf, variables = variables)
  fail_count <- 0
  # random sample until a valid function is generated
  while(class(tf_eval) == "try-error" & fail_count < 2000){
      fail_count <- fail_count + 1
      index_prediction <- index_sampler(index_pred)
      predicted_tf <- NULL_dictionary[index_prediction] %>% paste0(collapse = "")
      tf_eval <- tf_evaluation(predicted_tf, variables = variables)
      if(is.null(tf_eval)) {
        tf_eval <- "try-error"
        class(tf_eval) <- "try-error"
      }
  }
  
  tf_eval <- tf_evaluation(predicted_tf, variables = variables)
  if(class(tf_eval) == "try-error" | is.null(tf_eval)) return(NA)
  return(predicted_tf)
}
tf_generator <- function(point, latent_dim, generator,
                         dictionary, variables){
  point <- matrix(as.numeric(point), ncol = latent_dim)
  index_prob_prediction <- predict(generator, point, batch_size = 1)
  point_tf <- suppressWarnings(apply(index_prob_prediction, 1, tf_prediction,
                                     dictionary = dictionary,
                                     variables = variables)
  )
  for(i in 1:5) point_tf <- gsub("--", "-", point_tf, fixed = TRUE)
  for(i in 1:5) point_tf <- gsub("++", "-", point_tf, fixed = TRUE)
  return(point_tf)
}




