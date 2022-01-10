# FSO VAE
# FSO-mHM project
# Moritz Feigl, Sep 2020

#cat("Available GPUs:", strategy$num_replicas_in_sync, "\n")
# Loss functions
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

# Perception field for choosing TCN layers
percept_field <- function(kernel_sizes, dilation_rates){
  pf <- 1
  for(i in 1:length(kernel_sizes)){
    pf <- pf + (kernel_sizes[i]-1) * dilation_rates[i]
  }
  pf
}


index_reconstructor <- function(pred_matrix){
  index <- integer(nrow(pred_matrix))
  for(i in 1:nrow(pred_matrix)) index[i] <- which.max(pred_matrix[i, ])
  return(index)
}

# Structure parameters for TCNs
TCN_filters <- rep(TCN_units, length(TCN_kernel_sizes))
cat("TF decoder perception field:", percept_field(TCN_kernel_sizes, TCN_dilation_rates), "\n")
TCN_filters2 <- rep(TCN_units2, length(TCN_kernel_sizes2))
cat("Distribution decoder perception field:", percept_field(TCN_kernel_sizes2, TCN_dilation_rates2))

# Encoder ------------------------------------------------------------------------------
#with (strategy$scope(), {
input_shape <- as.integer(padding_length) 
# function token inputs
x1 <- layer_input(shape = input_shape, name = "tf_input")
embedding <- tf$keras$layers$Embedding(input_dim = as.integer(dictionary_length),
                                       output_dim = as.integer(embedding_length),
                                       input_length = as.integer(padding_length),
                                       trainable = TRUE, name = "embedding")
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
  dense_inter()

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
  kl_loss()

# Decoder --------------------------------------------------------------------------------
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
metric_tf <- tf$keras$metrics$SparseCategoricalCrossentropy
metric_dist <- tf$keras$metrics$MeanAbsoluteError
optimizer <- tf$optimizers$Adam()
vae_model %>% compile(optimizer = optimizer, loss = list("tf_output" = vae_loss,
                                                         "dist_output" = dist_loss),
                      loss_weights = list(1, 2),
                      metrics = list("tf_output" = vae_metric,
                                     "dist_output" = metric_dist()))
#})