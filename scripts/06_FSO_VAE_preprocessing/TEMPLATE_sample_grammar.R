# Sample Grammar for FSO TEMPLATE FOR VSC
# FSO-mHM project
# Moritz Feigl, Aug 2020


setwd("/home/lv71468/mfeigl/FSO_mHM/FSO_Data/Grammar_samples")
.libPaths("/home/lv71468/mfeigl/R/x86_64-pc-linux-gnu-library/4.0")
library(FSO)
# Define variables and operators of the CFG
operators <- c('+', '-', '*', '/')
# Define grammar
grammar <- create_grammar(eq = "<eq> <op> <eq>,
                          <eq> <op> numeric,
                          <eq> <op> var,
                          <eq> <op> (<eq>),
                          var,
                          f(var),
                          f(<eq>),
                          (<eq>)^(<pm>numeric),
                          numeric",
                          op = paste(operators, collapse = ","),
                          pm = "+, -")
# Sample grammar
funs <- grammar_sampler(n = 36000000, grammar = grammar, max_depth = 5,
                        no_cores = 32, seed = 333, save = TRUE,
                        file_name = "36mil_sampled_grammar")
