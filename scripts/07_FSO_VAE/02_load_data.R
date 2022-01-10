# Load VAE datasets 
# FSO-mHM project
# Moritz Feigl, Sep 2020

# Dictionaries
# dict_from_list <- function(...){
#   dict <- c(...)
#   dictionary <- as.list(1:length(dict))
#   names(dictionary) <- dict
#   return(dictionary)
# }
dict_from_list <- function(...){
  dictionary <- c(...)
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
dictionary <- dictionaries[[var]]
dictionary_length <- as.integer(length(dictionary)) + 1L

# Load csv data as tfdataset
padding_length <- read.csv("../padding_length.csv")[, 3]


# define train (80%) and validation (20%) datasets
if(!dir.exists("validation")){
  files <- list.files()
  dir.create("validation")
  nr_val_files <- round(length(files) * 0.2)
  validation_files <- sample(files, nr_val_files)
  file.copy(validation_files, paste0("validation/", validation_files))
  file.remove(validation_files)
  dir.create("training")
  file.copy(list.files(pattern = ".csv"), paste0("training/", list.files(pattern = ".csv")))
  file.remove(list.files(pattern = ".csv"))
}

# Define spec
example_file <- list.files("training", full.names = TRUE)[1]
data_spec <- csv_record_spec(
  example_file = example_file, 
  types = paste0(c(rep("d", 9), rep("i", padding_length)), collapse = "")
)

# Training dataset
train_dataset <- read_files("training/*.csv", text_line_dataset, 
                      record_spec = data_spec,
                      parallel_files = 10, 
                      parallel_interleave = 16) %>% 
  dataset_shuffle(buffer_size = 100000) %>% 
  dataset_repeat(count = 50) %>% 
  dataset_map(function(record) {
    record <- unname(record)
    list(
      list("tf_input" = record[-c(1:9)], "dist_input" = record[1:9]),
      list("tf_output" = record[-c(1:9)], "dist_output" = record[1:9]))
  }) %>% 
  dataset_batch(batch_size = batch_size) %>% 
  dataset_prefetch(50000) 

# Validation dataset
val_dataset <- read_files("validation/*.csv", text_line_dataset, 
                            record_spec = data_spec,
                            parallel_files = 4, 
                            parallel_interleave = 16) %>% 
  dataset_shuffle(buffer_size = 10000) %>% 
  dataset_repeat(count = 50) %>% 
  dataset_map(function(record) {
    record <- unname(record)
    list(
      list("tf_input" = record[-c(1:9)], "dist_input" = record[1:9]),
      list("tf_output" = record[-c(1:9)], "dist_output" = record[1:9]))
  }) %>% 
  dataset_batch(batch_size = batch_size) %>% 
  dataset_prefetch(2) 


