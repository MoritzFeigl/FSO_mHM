# Check VAE distribution prediction
# FSO-mHM project
# Moritz Feigl, May 2021




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

check_model(n = 12, 
            dict = dictionaries[[var]], 
            epoch = epochs, 
            distribution_scale_parameters = "../distribution_scale_parameters.csv", 
            save_in = "../training_performance")