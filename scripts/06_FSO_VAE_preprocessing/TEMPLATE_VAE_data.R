# Prepare TFs for VAE template
# FSO-mHM project
# Moritz Feigl, Sep 2020

# setup
setwd("/home/lv71468/mfeigl/FSO_mHM/FSO_Data/prepared_TFs/prepared_functions")
.libPaths("/home/lv71468/mfeigl/R/x86_64-pc-linux-gnu-library/4.0")
library(tensorflow)
library(tokenizers)
library(magrittr)

# Dictionaries
dict_from_list <- function(...){
  dict <- c(...)
  dictionary <- as.list(1:length(dict))
  names(dictionary) <- dict
  return(dictionary)
}
f <- c('exp', 'log10', 'log', 'sin', 'cos', 'tan', 'abs', 'acos',
       'asin', 'atan', 'cosh', 'sinh', 'sqrt')
var_KSat <- c("bd", "sand", "clay", "slope", "aspect", "dem")
var_FieldCap <- c("bd", "sand", "clay", "slope", "aspect", "dem", "thetas", "ksat", "vgenu_n")
operators <- c('+', '-', '*', '/', "^", "(", ")")
numbers <- seq(0.05, 3, 0.05) %>% format(nsmall=2)
dictionary_KSat <- dict_from_list(var_KSat, f, operators, numbers)
dictionary_FieldCap <- dict_from_list(var_FieldCap, f, operators, numbers)
dictionaries <- list(KSat = dictionary_KSat,
                     FieldCap = dictionary_FieldCap)

# string preparation function
prepare_text_for_tokenize <- function(text, operators = "\\^\\+\\-\\/\\*\\(\\)"){
  # Function to prepare TFs for tokenizing
  # input: strings of transfer functions
  # output: strings with ";" seperator
  text <- stringr::str_replace_all(text, paste0('([', operators, '])'), " \\1 ")
  text <- tokenizers::tokenize_words(text, strip_punct = FALSE, lowercase = FALSE) %>%
    sapply(FUN = function(x) paste0(c("", paste0(x, collapse = ";"), ""), collapse = ";"))
  return(text)
}
# Parameter
var <- "Dummy"
file <- "dummy_file"
batch <- 0

# Preprocessing start
file_data <- fst::read_fst(file)

# use only TFs in coutoff bounds
cutoff <- read.csv(paste0("../../VAE_data/", var, "/distribution_cutoff_points.csv"))
tfs_in_bounds <- which(file_data$`10%` > cutoff$lower_bound & 
                         file_data$`90%`< cutoff$upper_bound)
file_data <- file_data[tfs_in_bounds, ]
sample_tfs <- file_data$`function`

# Preprocess strings
# keras tokenizer with pre-defined dictionary
tt <- tf$keras$preprocessing$text$Tokenizer(filters = "\\;", lower=FALSE)
tt$word_index <- dictionaries[[var]]
prep_funs <- parallel::mclapply(sample_tfs, prepare_text_for_tokenize, mc.cores = 32)
funs_tokenized <- keras::texts_to_sequences(tt, prep_funs)
max_length <- max(lengths(funs_tokenized))
padding_length <- read.csv(paste0("../../VAE_data/", var, "/padding_length.csv"))
padding_length <- padding_length$padding_length
cat("Files max seqence length is", max_length, "and will be padded to", padding_length, "\n")
# remove TFs with length > padding_length
too_long_functions <- which(lengths(funs_tokenized) > padding_length)
if(length(too_long_functions) > 0){
  funs_tokenized <- funs_tokenized[-too_long_functions]
  file_data <- file_data[-too_long_functions, ]
}
funs_padded <- keras::pad_sequences(funs_tokenized, maxlen = padding_length)

# Preprocess numerical data
scale_para <- read.csv(paste0("../../VAE_data/", var, "/distribution_scale_parameters.csv"))
for(col in 2:ncol(file_data)){
  file_data[, col] <- (file_data[, col] - scale_para$mean[col-1]) / scale_para$sd[col-1]
}

# combina and shuffle
n <- nrow(file_data)
shuffle_ids <- sample(n, n)
data <- cbind(file_data[shuffle_ids, ], funs_padded[shuffle_ids, ])

# save as csv
write.csv(data[, -1], paste0("../../VAE_data/", var, "/batches/vae_data_batch", batch, ".csv"), 
          row.names = FALSE)