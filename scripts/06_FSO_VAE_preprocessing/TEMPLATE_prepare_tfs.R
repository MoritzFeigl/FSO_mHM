# Prepare transfer functions for FSO
# FSO-mHM project
# Moritz Feigl, Aug 2020


setwd("/home/lv71468/mfeigl/FSO_mHM/FSO_Data/prepared_TFs")
.libPaths("/home/lv71468/mfeigl/R/x86_64-pc-linux-gnu-library/4.0")
library(FSO)
# Sample variables and numerics
f <- c('exp', 'log10', 'log', 'sin', 'cos', 'tan', 'abs', 'acos',
       'asin', 'atan', 'cosh', 'sinh', 'sqrt', 'tanh')
var_KSat <- c("bd", "sand", "clay", "slope", "aspect", "dem")
var_FieldCap <- c(var_KSat, "ThetaS", "KSat", "vGenu_n")
numbers <- seq(0.1, 3, 0.01)
KSat_variable_list <- list("var" = var_KSat,
                           "numeric" = numbers,
                           "f" = f)
FieldCap_variable_list <- list("var" = var_FieldCap,
                               "numeric" = numbers,
                               "f" = f)
variable_input(functions = "dummy_file",
               variable_list = KSat_variable_list,
               single_var_to_remove = "numeric",
               necessary_var_to_be_included = "var",
               n_iter = 10, no_cores = 16, seed = 4323,
               file_name = "KSat_functions")

simplify_functions(files = "dummy_files",
                   no_cores = 32)
