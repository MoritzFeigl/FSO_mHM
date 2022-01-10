check_model <- function(n, dict = dictionary, epoch, save_in,
                        distribution_scale_parameters){
  # define encoder/decoder
  folder <- paste0(save_in, "/epoch", epoch, "-generated_results_", format(Sys.time(), "%Y-%m-%d"))
  if(!dir.exists(folder)) dir.create(folder)
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
  dist_input <- layer_input(shape = latent_dim)
  dist_generator <- dist_input %>%
    dist_decoder1 %>%
    dist_reshape_z %>%
    dist_block1 %>%
    dist_block2 %>%
    dist_block3 %>%
    layer_flatten() %>%
    fc1 %>%
    fc2
  dist_decoder_model <- keras_model(dist_input, dist_generator)
  
  # Sample, generate and compare functions
  example_file <- list.files("validation", full.names = TRUE)
  example_file <- sample(example_file, 1)
  test_data <- read.csv(example_file)
  sampled_functions <- as.matrix(test_data[sample(nrow(test_data), n), ])
  tfs_encoded <- encoder_model(
    list(sampled_functions[, -c(1:9)], sampled_functions[, 1:9]))
  tfs_decoded <- as.array(decoder_model(tfs_encoded))
  pred_ind <- vector(mode = "list", length = n)
  for(i in 1:n){
    pred_ind[[i]] <- index_reconstructor(tfs_decoded[i, , ])
  }
  pred_ind <- do.call(rbind, pred_ind)
  true_ind <- sampled_functions[,-c(1:9)] + 1
  NULL_dictionary <- c("", dictionary)
  true_fun <- character(n)
  pred_fun <- character(n)
  for(i in 1:n){
    pred_fun[i] <- NULL_dictionary[pred_ind[i, ]] %>% paste0(collapse = "")
    true_fun[i] <- NULL_dictionary[true_ind[i, ]] %>% paste0(collapse = "")
  }
  compare_df <- data.frame("true functions" = true_fun,
                           "predicted functions" = pred_fun,
                           check.names = FALSE)
  row.names(compare_df) <- NULL
  
  write.csv(compare_df, paste0(folder, "/epoch", epoch, "-generated_functions.csv"), 
            row.names = FALSE)
  
  # Plot embedding
  embedding_weights <- t(get_weights(embedding)[[1]])
  pca <- prcomp(embedding_weights)
  pca_coordinates <- pca$rotation
  pca_coordinates <- data.frame(pca_coordinates,
                                Name = c("NULL", dictionary))
  library(ggplot2)
  ggplot(pca_coordinates, aes(x = PC1, y = PC2, label=Name)) +
    geom_text(aes(label=c("NULL", dictionary))) +
    ggsave(paste0(folder, "/epoch", epoch, "-embedding_pca12.png"),
           width = 9, height = 9, units = "in")
  ggplot(pca_coordinates, aes(x = PC3, y = PC4, label=Name)) +
    geom_text(aes(label=c("NULL", dictionary))) +
    ggsave(paste0(folder, "/epoch", epoch, "-embedding_pca34.png"),
           width = 9, height = 9, units = "in")
  ggplot(pca_coordinates, aes(x = PC5, y = PC6, label=Name)) +
    geom_text(aes(label=c("NULL", dictionary))) +
    ggsave(paste0(folder, "/epoch", epoch, "-embedding_pca56.png"),
           width = 9, height = 9, units = "in")
  # Plot latent space
  dists <- t(sampled_functions[, 1:9])
  dist_pca <- prcomp(dists)
  dist_pca_coordinates <- dist_pca$rotation
  
  latent_representation <- t(as.array(tfs_encoded$sample()))
  tf_pca <- prcomp(latent_representation)
  tf_pca_coordinates <- tf_pca$rotation
  pca_coordinates <- data.frame(tf_pca_coordinates,
                                dist_pca_coordinates,
                                Name = 1:n)
  
  ggplot(pca_coordinates, aes(x = PC1, y = PC1.1)) +geom_point(shape = 1) +
    geom_text(aes(label = Name, hjust = 1.5))+
    geom_smooth(method='lm') + labs(y = "PC1 of parameter distribution",
                                    x = "PC1 of latent space") +
    ggsave(paste0(folder, "/epoch", epoch, "-FSvsDIST.png"),
           width = 9, height = 9, units = "in")
  
  # Walk through functionspace
  weights <- c(seq(0, 1, 0.1))
  tf_ind <- sample(n, 2)
  # get encoded start and end point
  encoded_start_end <- as.array(tfs_encoded$sample())[tf_ind, ]
  #generate functions in between
  new_tfs <- matrix(NA, ncol = ncol(encoded_start_end), nrow = length(weights))
  for(i in seq_along(weights)){
    new_tfs[i, ] <- encoded_start_end[1,] * (1-weights[i]) + encoded_start_end[2,] * weights[i]
  }
  # generate tfs and distributions for all
  tfs_in_fs <- as.array(decoder_model(new_tfs))
  dist_in_fs <- as.data.frame(as.array(dist_decoder_model(new_tfs)))
  
  tf_pred <- vector(mode = "list", length = 11)
  for(i in 1:11){
    tf_pred[[i]] <- index_reconstructor(tfs_in_fs[i, , ])
  }
  tf_pred <- do.call(rbind, tf_pred)
  tf_pred_strings <- character(11)
  for(i in 1:11){
    tf_pred_strings[i] <- NULL_dictionary[tf_pred[i, ]] %>% paste0(collapse = "")
  }
  
  dist_in_fs <- cbind(name = seq(0, 1, 0.1), dist_in_fs)
  colnames(dist_in_fs) <- c("name", 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  scale_para <- read.csv(distribution_scale_parameters)
  for(col in 2:10){
    dist_in_fs[, col] <- dist_in_fs[, col] * scale_para$sd[col-1] + scale_para$mean[col-1]
  }
  # Plot
  plot_data <- reshape2::melt(dist_in_fs, id.vars = "name")
  names <- c("start", "step 1", "step 2", "step 3", "step 4", "step 5",
             "step 6", "step 7", "step 8", "step 9", "end" )
  plot_data$name <- factor(plot_data$name, levels = seq(0, 1, 0.1), labels = names)
  plot_data$variable <- as.numeric(as.character(plot_data$variable))
  
  cc <- scales::seq_gradient_pal("green", "red", "Lab")(seq(0,1,0.1))
  
  ggplot(plot_data, aes(x = value, y = variable, col = name)) + geom_line() +
    scale_colour_manual(values = cc ) + xlab("parameter values") +
    ylim(0, 1) +
    ylab("cumulative probability") + theme_classic() +
    ggsave(paste0(folder, "/walk_through_FS.png"), width = 7, height = 7, units = "in")
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))
  
  write.table(data.frame(weight_f2 = weights, tf_pred_strings),
              paste0(folder, "/walk_through_function_space.txt"),
              row.names = FALSE)
  
  
  
}
