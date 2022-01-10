# FSO-mHM utility functions
# FSO-mHM project
# Moritz Feigl, Nov 2020

# function that adds correct scaling numerical values to tf
prepare_tf_for_mhm <- function(tf, scaling_bounds){
  for(i in 1:length(scaling_bounds)){
    scaling_bounds[[i]] <- scaling_bounds[[i]] %>% format(nsmall = 1) %>% gsub(" ", "", .)
  }

  # prepare
  tf <- gsub(" ", "", tf)
  tf <- gsub("^", "**", tf, fixed = TRUE)
  # split tf
  tf_splitted <- tf %>%
    strsplit(., c("[/^()*+-]")) %>%
    unlist()
  tf_splitted <- tf_splitted[tf_splitted != ""]
  var <- character()
  variables <- names(scaling_bounds)
  for(i in seq_along(tf_splitted)){
    if(tf_splitted[i] %in% variables)  var <- c(var, tf_splitted[i])
  }
  var <- unique(var)
  # change variables to scaled variables in tf
  for(variable in var){
    bounds <- scaling_bounds[[variable]]
    tf <- gsub(pattern = variable,
               replacement = paste0("((", variable, "-", bounds[1], ")/(",
                                    bounds[2], "-", bounds[1], "))"),
               tf)
  }
  # get numerical values
  tf_splitted <- tf %>%
    strsplit(., c("[/^()*+-]")) %>%
    unlist()
  tf_splitted <- tf_splitted[tf_splitted != ""]
  numerics <- character()
  suppressWarnings(
    for(i in seq_along(tf_splitted)){
      if(!is.na(as.numeric(tf_splitted[i]))) numerics <- c(numerics, tf_splitted[i])
    }
  )
  # make integers to float
  ints <- grep(".", numerics, invert = TRUE, fixed = TRUE)
  numerics[ints] <- paste0(numerics[ints], ".0")

  # put space into function
  tf <- tf %>%
    gsub("+", " + ", ., fixed = TRUE) %>%
    gsub("-", " - ", ., fixed = TRUE) %>%
    gsub("*", " * ", ., fixed = TRUE) %>%
    gsub("*  *", "**", ., fixed = TRUE) %>%
    gsub("/", " / ", ., fixed = TRUE)

  if(substr(tf, 1, 1) == " ") tf <- substring(tf, 2)
  return(list("scaled_function" = tf, "numerics" = numerics))
}

# Runs mHM with tfs given as strings
run_mhm_from_tf <- function(KSat_tf, FieldCap_tf, numeric_parameters, run_folder){
  path <- paste0(run_folder, "/FSO_mHM_major_basins/config")
  basin_folders <- list.files(path, pattern = "sub_")

  variables <- c("bd", "sand", "clay", "slope", "aspect", "ThetaS", "KSat",
                 "vGenu_n", "dem")
  alltill_variables <- c("bd", "sand", "clay", "ThetaS", "KSat", "vGenu_n", "aspect", "dem", "slope")
  till_variables <- c("bd_till", "sand_till", "clay_till", "ThetaS_till", "KSat_till",
                      "vGenu_n_till", "slope_till", "aspect_till", "dem_till")
  notill_variables <- c("bd_notill", "sand_notill", "clay_notill", "ThetaS_notill",
                        "KSat_notill", "vGenu_n_notill", "slope_notill", "aspect_notill", "dem_notill")
  scaling_bounds <- list("slope" = c(0, 90),
                         "aspect" = c(0, 360),
                         "bd" = c(0, 2.3),
                         "sand" = c(0, 100),
                         "clay" = c(0, 100),
                         "dem" = c(0, 4000),
                         "KSat" = c(9, 244),
                         "vGenu_n" = c(1,2),
                         "ThetaS" = c(0.24, 0.51))

  # Add scaling to given transfer functions
  # Ksat
  KSat_from <- as.numeric(
    read.csv("/home/lv71468/mfeigl/FSO_mHM/FSO_Data/VAE_data/KSat/distribution_cutoff_points.csv"))
  KSat_to <- c(1.1, 1000.0)
  a <- (KSat_from[2] - KSat_from[1]) %>% format(nsmall = 1)
  b <- (KSat_to[2] - KSat_to[1]) %>% format(nsmall = 1)
  KSat_from[1] <- KSat_from[1] %>% format(nsmall = 1)
  KSat_to[1] <- KSat_to[1] %>% format(nsmall = 1)
  Ksat_fun <- paste0(KSat_to[1], "+((", KSat_tf, "-", KSat_from[1], ")", "*", b, ")/", a)
  Ksat_fun <- gsub("--", "+", Ksat_fun, fixed = TRUE)
  # FieldCap
  FieldCap_from <- as.numeric(
    read.csv("/home/lv71468/mfeigl/FSO_mHM/FSO_Data/VAE_data/FieldCap/distribution_cutoff_points.csv"))
  FieldCap_to <- c(0.01, 0.55)
  a <- (FieldCap_from[2] - FieldCap_from[1]) %>% format(nsmall = 1)
  b <- (FieldCap_to[2] - FieldCap_to[1]) %>% format(nsmall = 1)
  FieldCap_from[1] <- FieldCap_from[1] %>% format(nsmall = 1)
  FieldCap_to[1] <- FieldCap_to[1] %>% format(nsmall = 1)
  FieldCap_fun <- paste0(FieldCap_to[1], "+((", FieldCap_tf, "-", FieldCap_from[1], ")", "*", b, ")/", a)
  FieldCap_fun <- gsub("--", "+", FieldCap_fun, fixed = TRUE)
  # prepare tfs with scaling bounds
  prepared_KSat_tf <- prepare_tf_for_mhm(Ksat_fun, scaling_bounds = scaling_bounds)
  prepared_FieldCap_tf <- prepare_tf_for_mhm(FieldCap_fun, scaling_bounds = scaling_bounds)
  # Ksat
  Ksat_fun <- prepared_KSat_tf$scaled_function
  KSat_nums <- prepared_KSat_tf$numerics
  # FieldCap
  FieldCap_fun <- prepared_FieldCap_tf$scaled_function
  FieldCap_nums <- prepared_FieldCap_tf$numerics
  # make numeric vectors with parameter names
  numeric_paras <- c(KSat_nums, FieldCap_nums) %>% unique()
  numeric_paras_names <- paste0("FSO", seq_along(numeric_paras))
  numeric_paras <- numeric_paras[order(nchar(numeric_paras), decreasing = TRUE)]
  # fill in numerics in functions
  for(k in seq_along(numeric_paras)){
    Ksat_fun <- gsub(numeric_paras[k], numeric_paras_names[k], Ksat_fun)
    FieldCap_fun <- gsub(numeric_paras[k], numeric_paras_names[k], FieldCap_fun)
  }

  # change KSat variables to _horizon if necessary
  for(var in c("aspect", "slope", "dem")){
    Ksat_fun <- gsub(var, paste0(var, "_horizon"), Ksat_fun)
  }

  # change all variables into _till and _notill
  FieldCap_till <- FieldCap_fun
  for(var in alltill_variables){
    FieldCap_till <- gsub(var, paste0(var, "_till"), FieldCap_till)
  }
  FieldCap_notill <- FieldCap_fun
  for(var2 in alltill_variables){
    FieldCap_notill <- gsub(var2, paste0(var2, "_notill"), FieldCap_notill)
  }

  # 1. Change transfer functions in nml
  mpr_nml <- readLines("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_master_namelists/mpr_fso.nml")
  # change KSat tf
  KSat_line <- grep('KSat_dummy', mpr_nml)
  mpr_nml[KSat_line] <- paste0("                           '", Ksat_fun, "',")
  fun_vars <- unlist(strsplit(gsub(" ", "", Ksat_fun), "(?=[+-/*)()])", perl = TRUE))
  horizon_vas <- c("aspect_horizon", "slope_horizon", "dem_horizon")
  used_var <- unique(fun_vars[fun_vars %in% c(variables, horizon_vas)])
  data_array_line <- grep("from_data_arrays(1:3,26)", mpr_nml, fixed = TRUE)
  mpr_nml[data_array_line] <- paste0('    from_data_arrays(1:', length(used_var), ',26',
                                     ') = ',
                                     paste0("'", paste0(used_var, collapse = "', '"), "'"),
                                     collapse = "")

  # change FieldCap tfs
  FieldCap_till_line <- grep('FieldCap_till_dummy', mpr_nml)
  FieldCap_notill_line <- grep('FieldCap_notill_dummy', mpr_nml)
  mpr_nml[FieldCap_till_line] <- paste0("                           '", FieldCap_till, "',")
  mpr_nml[FieldCap_notill_line] <- paste0("                           '", FieldCap_notill, "',")

  # data array FieldCap_till
  fun_vars <- unlist(strsplit(gsub(" ", "", FieldCap_till), "(?=[+-/*)()])", perl = TRUE))
  used_var <- unique(fun_vars[fun_vars %in% till_variables])
  data_array_line <- grep("from_data_arrays(1:3,37)", mpr_nml, fixed = TRUE)
  mpr_nml[data_array_line] <- paste0('    from_data_arrays(1:', length(used_var), ',37',
                                     ') = ',
                                     paste0("'", paste0(used_var, collapse = "', '"), "'"),
                                     collapse = "")
  # data array FieldCap_notill
  fun_vars <- unlist(strsplit(gsub(" ", "", FieldCap_notill), "(?=[+-/*)()])", perl = TRUE))
  used_var <- unique(fun_vars[fun_vars %in% notill_variables])
  data_array_line <- grep("from_data_arrays(1:3,38)", mpr_nml, fixed = TRUE)
  mpr_nml[data_array_line] <- paste0('    from_data_arrays(1:', length(used_var), ',38',
                                     ') = ',
                                     paste0("'", paste0(used_var, collapse = "', '"), "'"),
                                     collapse = "")

  # get parameter names, parameter_values
  parameter_name_line <- grep('parameter_names', mpr_nml)
  parameter_value_line <- grep('parameter_values', mpr_nml)
  parameters_end <- grep('&data_arrays', mpr_nml, fixed = TRUE) - 3
  parameter_names_initial <- mpr_nml[parameter_name_line:(parameter_value_line-1)]
  parameter_values_initial <- mpr_nml[parameter_value_line:parameters_end]
  # add new parameter names and values
  nlines <- length(parameter_names_initial)
  parameter_names_initial[nlines] <- paste0(parameter_names_initial[nlines], ",")
  parameter_names <- c(
    parameter_names_initial,
    paste0("                            ",
           paste0(paste0("'", numeric_paras_names, "'"), collapse = ", ")))
  nlines <- length(parameter_values_initial)
  parameter_values_initial[nlines] <- paste0(parameter_values_initial[nlines], ",")
  parameter_values <- c(
    parameter_values_initial,
    paste0("                            ", paste0(numeric_paras, collapse = ", ")))
  old_nr_of_numeric_paras <- mpr_nml[parameter_name_line] %>%
    strsplit(split = "(", fixed = TRUE) %>%
    unlist() %>%
    substr(., 3, 4) %>% paste0(., collapse = "") %>%
    as.numeric()
  total_nr_of_numeric_paras <- old_nr_of_numeric_paras + length(numeric_paras)
  parameter_names[1] <- paste0("    parameter_names(1:",total_nr_of_numeric_paras, ") = ",
                               strsplit(parameter_names[1], " = ", fixed = TRUE)[[1]][2])
  parameter_values[1] <- paste0("    parameter_values(1:",total_nr_of_numeric_paras, ") = ",
                                strsplit(parameter_values[1], " = ", fixed = TRUE)[[1]][2])

  # add new parmeter names and values
  new_mpr_nml <- c(mpr_nml[1:(parameter_name_line-1)],
                   parameter_names,
                   parameter_values,
                   mpr_nml[(parameters_end+1):length(mpr_nml)])
  # Change numeric parameters
  parameters <- suppressWarnings(
    read_nml("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_master_namelists/mhm_parameter.nml"))
  for(i_numerics in seq_along(numeric_parameters)){
    i_para <- ind_of_num_paras[i_numerics]
    parameters$mhm_parameters[[i_para]][3] <- numeric_parameters[i_numerics]
  }

  # write new mpr.nml and mhm_parameter.nml for all basins
  for(basin in basin_folders){
    basin_nml <- gsub("<BASIN>", basin, new_mpr_nml)
    writeLines(basin_nml, con = paste0(path, "/", basin, "/mpr.nml"))
    write_mhm_nml(parameters, paste0(path, "/", basin, "/mhm_parameter.nml"))
  }

  # 2. Run python script and compile mhm
  output_error <- try({
    system(paste0("chmod 775 ", run_folder, "/05_run_mhm_preparation.sh"))
    system(paste0("bash ", run_folder, "/05_run_mhm_preparation.sh"))
    for(basin_nr in 1:7){
      system(paste0("chmod 775 ", run_folder, "/mhm_run_scripts/run_basin_", basin_nr, ".sh"))
    }
    cat("running mhm\n")
    system(paste0("bash ", run_folder, "/mhm_run_scripts/run_basin_1.sh &> ", run_folder, "/mhm_run_scripts/run_basin_1.log & ",
                  "bash ", run_folder, "/mhm_run_scripts/run_basin_2.sh &> ", run_folder, "/mhm_run_scripts/run_basin_2.log & ",
                  "bash ", run_folder, "/mhm_run_scripts/run_basin_3.sh &> ", run_folder, "/mhm_run_scripts/run_basin_3.log & ",
                  "bash ", run_folder, "/mhm_run_scripts/run_basin_4.sh &> ", run_folder, "/mhm_run_scripts/run_basin_4.log & ",
                  "bash ", run_folder, "/mhm_run_scripts/run_basin_5.sh &> ", run_folder, "/mhm_run_scripts/run_basin_5.log & ",
                  "bash ", run_folder, "/mhm_run_scripts/run_basin_6.sh &> ", run_folder, "/mhm_run_scripts/run_basin_6.log & ",
                  "process_id=$! bash ", run_folder, "/mhm_run_scripts/run_basin_7.sh &> ", run_folder,
                  "/mhm_run_scripts/run_basin_7.log & wait $process_id"
    ))
    cat("Finished running mhm!\n")
    # read outputs and get KGE/NSE
    basin_results <- data.frame("Basin" = basin_folders, "KGE" = NA, "NSE" = NA, "lNSE" = NA)
    for(basin in basin_folders){
      output <- readLines(paste0(path, "/", basin, "/output.txt"))
      lines <- grep("KGE of daily discharge", output)
      KGE <- as.numeric(substring(output[lines], 39))
      NSE <- as.numeric(substring(output[lines + 1], 39))
      lNSE <- as.numeric(substring(output[lines + 2], 39))
      try(basin_results[basin_results$Basin == basin, c("KGE", "NSE", "lNSE")] <- c(KGE, NSE, lNSE))
    }
  })
  # If mHM  produces no output return NA
  if(class(output_error) == "try-error"){
    basin_results <- data.frame("Basin" = basin_folders, "KGE" = NA, "NSE" = NA)
  }

  return(basin_results)
}


.function_splitter <- function(point_tf){
  function_splitted <- unlist(strsplit(point_tf, c("[/^()*+-]")))
  function_splitted <- gsub(" ", "", function_splitted)
  function_splitted <- function_splitted[function_splitted != ""]
  return(function_splitted)
}

.dds <- function(xBounds.df, numIter, OBJFUN, search_dim, start_point = NULL,
                 state_folder = NULL, iterations_per_run = numIter){
  # INPUTS:
  # xBounds.df must be a dataframe with 1st column as minimum, 2nd column as maximum
  # numIter is an integer
  # OBJFUN is a function which returns a scalar value, for which we are trying to minimize.
  #
  # OUTPUTS:
  # outputs.df is a two entry list, containing x_best and y_best, as they evolve over numIter iterations.

  colnames(xBounds.df) <- c("min", "max")

  if(!is.null(state_folder) & file.exists(paste0(state_folder, "/optimization_state.rds"))){
    opt_state <- readRDS(paste0(state_folder, "/optimization_state.rds"))
    peturbIdx <- opt_state[["peturbIdx"]]
    x_best <- opt_state[["x_best"]]
    y_best <- opt_state[["y_best"]]
    start_index <- opt_state[["i"]]
    if(start_index == numIter) return(cat("Optimization already complete."))
  } else {
    # Format xBounds.df colnames
    if(is.null(start_point)){
      # Generate initial first guess that is not NA
      x_init <- apply(xBounds.df, 1, function(x) runif(1, x[1], x[2]))
    } else {
      cat("Using pre-defined start point for optimization\n")
      x_init <- start_point
    }
    # Evaluate first cost function
    x_evaluated <- OBJFUN(x_init)
    while(is.na(x_evaluated$loss)){
      for(row in 1:search_dim){
        x_init[row] <- rnorm(1,
                             mean = x_init[row],
                             sd = (xBounds.df[row, "max"] - xBounds.df[row, "min"])/100)
        # check for boundary violation
        if(x_init[row] > xBounds.df[row, "max"]) x_init[row] <- xBounds.df[row, "max"]
        if(x_init[row] < xBounds.df[row, "min"]) x_init[row] <- xBounds.df[row, "min"]
      }
      x_evaluated <- OBJFUN(x_init)
    }

    x_ini <- x_evaluated$`Current point in function space`
    x_best <- matrix(x_init, nrow = 1)

    if(!is.na( x_evaluated$loss)){
      y_init <- x_evaluated$loss
    } else {
      y_init <- -999
    }
    y_best <- y_init

    # Select which entry to peturb at each iteration
    peturbIdx <- .probPeturb(xBounds.df, numIter)
    start_index <- 2
    iterations_per_run <- iterations_per_run - 2
  }

  r = 0.2
  # Peturb each entry by N(0,1)*r(x_max - x_min) reflecting if @ boundaries
  sigma <- xBounds.df$max - xBounds.df$min
  for (i in start_index:(start_index + iterations_per_run)){
    if(i > numIter) break
    # Set up test x
    x_test <- x_best[i-1, ]
    # Get entries we will peturb
    idx <- peturbIdx[[i]]
    # Initialize vector of peturbations initially zeros with same length of x so we will add this vector to peturb x
    peturbVec <- rep(0, length(x_test))
    # Generate the required number of random normal variables
    N <- rnorm(length(x_test), mean=0, sd=1)
    # Set up vector of peturbations
    peturbVec[idx] <- r*N[idx]*sigma[idx]
    # Temporary resulting x value if we peturbed it
    testPeturb <- x_test + peturbVec
    # Find the values in testPeturb OBJFUN <- wrapper_ofthat have boundary violations.  Store the indices in boundaryViolationsIdx
    boundaryViolationIdx <- which(testPeturb < xBounds.df$min | testPeturb > xBounds.df$max)
    # Reset those violated indices to the opposite peturbation direction
    peturbVec[boundaryViolationIdx] <- (-1*r*N[boundaryViolationIdx]*sigma[boundaryViolationIdx])
    # Find values still at violations of min or max and set them to the minimum or maximum values
    testPeturb <- x_test + peturbVec
    minViolationIdx <- which(testPeturb < xBounds.df$min)
    maxViolationIdx <- which(testPeturb > xBounds.df$max)
    testPeturb[minViolationIdx] <- xBounds.df$min[minViolationIdx]
    testPeturb[maxViolationIdx] <- xBounds.df$max[maxViolationIdx]
    # Peturb the test vector
    x_test <- x_test + peturbVec
    # Evaluate objective function
    x_evaluated <- OBJFUN(x_test)
    x_test <- x_evaluated$`Current point in function space`
    y_test <- x_evaluated$loss
    if(!is.na(y_test)) {
      y_best[i] <- max(c(y_test, y_best[i-1]))
      bestIdx <- which.max(c(y_test, y_best[i-1]))
    } else {
      y_best[i] <- y_best[i-1]
      bestIdx <- 2
    }
    x_choices <- cbind(x_test, x_best[i-1, ])
    x_best <- rbind(x_best, x_choices[,bestIdx])
    # save optimizer state
    if(!is.null(state_folder)){
      # save
      opt_state <- list("peturbIdx" = peturbIdx, "x_best" = x_best,
                        "y_best" = y_best, "i" = i)
      saveRDS(opt_state, paste0(state_folder, "/optimization_state.rds"))
    }
  }
  output.list <- list(t(x_best), y_best)
  return(output.list)
}

.probPeturb <- function(x, numIter){
  # Input: x and the number of iterations numIter
  # Output: list with indices that will be perturbed with length = numIter
  xDims <- nrow(x)
  probabilityVector <- 1- (log(1:numIter)/log(numIter))
  peturbIdx <- probabilityVector %>%
    lapply(function(x) as.logical(rbinom(xDims, 1, x))) %>%
    unlist() %>%
    matrix(byrow = TRUE, ncol = xDims) %>%
    apply(MARGIN = 1, FUN = which)
  return(peturbIdx)
}

write_mhm_nml  <-	function(glm_nml,file){
  sink(file)
  cat("! Emacs: -*- mode: f90 -*-\n",
      "!global_parameters\n",
      "!PARAMETER   lower_bound  upper_bound  value  FLAG  SCALING\n",
      "!interception)\n")
  print(glm_nml)
  sink()
}

read_nml  <-	function(nml_file){
  # skip all commented lines, return all variables and associated values
  # requires NO return line variables (all variables must be completely defined on a single line)
  c <- file(nml_file,"r")
  fileLines <- readLines(c)
  close(c)
  lineStart	<-	substr(fileLines,1,1)
  # ignore comment lines or empty lines
  ignoreLn	<-	lineStart=='!' | fileLines==""
  lineStart	<-	lineStart[!ignoreLn]
  fileLines	<-	fileLines[!ignoreLn]
  # find all lines which start with "&" * requires FIRST char to be value

  lineIdx		<- seq(1,length(lineStart))
  blckOpen	<-	lineIdx[lineStart=="&"]
  blckClse	<-	lineIdx[lineStart=="/"]

  nml <- list()
  for (i in seq_len(length(blckOpen))){
    blckName   <-	substr(fileLines[blckOpen[i]],
                         2, nchar(fileLines[blckOpen[i]]))
    blckName   <- gsub("\\s", "", blckName)
    oldNms	   <-	names(nml)
    nml[[i]]   <-	list()
    names(nml) <-	c(oldNms,blckName)

    carryover <- ''

    for (j in (blckOpen[i]+1):(blckClse[i]-1)){

      textLine	<-	paste(carryover,
                        gsub("\t", "", gsub(" ", "", fileLines[j])), sep = '')

      if(substr(textLine, 1, 1) != '!'){
        # Add a check here, sometimes, if there is a hanging comma,
        #and only sometimes that means add next row
        if(substr(textLine, nchar(textLine), nchar(textLine)) == ',' &&
           j+1 <= length(fileLines) &&
           !any(grep("=", fileLines[j + 1])) &&
           !any(grep("/", fileLines[j + 1]))){

          carryover = textLine
          next
        }else{
          carryover = ''
        }
        # else, line is commented out
        lineVal	  <-	buildVal(textLine, lineNum = j, blckName)
        nml[[i]]	<-	c(nml[[i]], lineVal)
      }
    }
  }
  nml <- .nml(nml)
  return(nml)
}

print.nml <- function(x, ...){
  glm_nml <- x
  for (i in seq_len(length(names(glm_nml)))){ # these are the blocks
    blckNm  <-	names(glm_nml)[i]
    cat("&")
    cat(blckNm)
    cat('\n')
    blckList	<-	glm_nml[[i]]
    for (j in seq_len(length(names(blckList)))){
      cat('   ')
      cat(names(blckList)[j])
      cat(' = ')
      if (length(blckList[[j]])>1){
        if (is.logical(blckList[[j]])){
          charText <- to.glm_boolean(blckList[[j]])
        } else {
          charText <- c(blckList[[j]])
        }
        writer	<-	paste(charText,collapse=', ')
      } else if (is.character(blckList[[j]])) {
        charText <- strsplit(blckList[[j]],',')
        writer <- paste(c("'",paste(c(charText[[1]]),collapse="','"),"'"),collapse='')
      } else if (is.logical(blckList[[j]])){
        writer <- to.glm_boolean(blckList[[j]])
      } else {
        writer <- blckList[[j]]
      }
      cat(writer)
      cat('\n')
    }
    cat('/\n')
  }
}

buildVal	<-	function(textLine, lineNum, blckName){
  #-----function appends nml list with new values-----
  # remove all text after comment string
  textLine	<-	strsplit(textLine,'!')[[1]][1]

  if (!any(grep("=", textLine))){
    stop(c("no hanging lines allowed in .nml, used ",textLine,'.\nSee line number:',lineNum,' in "&',blckName,'" section.'))
  }
  params	<-	strsplit(textLine,"=") # break text at "="
  parNm	  <-	params[[1]][1]
  parVl	  <-	params[[1]][2]
  # figure out what parval is...if string, remove quotes and keep as string
  # ***for boolean text, use "indentical" so that 0!= FALSE
  # can be: string, number, comma-sep-numbers, or boolean

  # special case for date:
  if (is.na(parVl)){
    stop('Empty values after "', textLine, '" on line ', lineNum,
         '. \nPerhaps the values are on the next line?', call. = FALSE)
  }
  if (nchar(parVl>17) & substr(parVl,14,14)==':' & substr(parVl,17,17)==':'){
    parVl<-paste(c(substr(parVl,1,11),' ',substr(parVl,12,nchar(parVl))),collapse='')
  }
  if (any(grep("'",parVl))){

    parVl	<-	gsub("'","",parVl)
  }else if (any(grep("\"",parVl))){
    parVl  <-	gsub("\"","",parVl)
  }else if (isTRUE(grepl(".true.",parVl) || grepl(".false.",parVl))){
    logicals <- unlist(strsplit(parVl,","))
    parVl <- from.glm_boolean(logicals)
  }else if (any(grep(",",parVl))){	# comma-sep-nums
    parVl	<-	c(as.numeric(unlist(strsplit(parVl,","))))
  }else {	# test for number
    parVl	<-	as.numeric(parVl)
  }
  lineVal	<-	list(parVl)
  names(lineVal)	<-	parNm
  return(lineVal)
}

from.glm_boolean <- function(values){
  logicals <- sapply(values, FUN = function(x){
    if (!isTRUE(grepl(".true.", x) || grepl(".false.", x))){
      stop(x, ' is not a .true. or .false.; conversion to TRUE or FALSE failed.',
           call. = FALSE)
    }
    return(ifelse(isTRUE(grepl(".true.", x)), TRUE, FALSE))
  })
  return(as.logical(logicals))
}

to.glm_boolean <- function(values){
  val.logical <- values
  values[val.logical] <- '.true.'
  values[!val.logical] <- '.false.'
  return(values)
}

findBlck	<-	function(nml,argName){

  # test for argName being a string
  if (!is.character(argName)){stop(c("parameter name must be a string"))}
  fau <- " "
  fault.string <- rep(fau,1000) # names fault matrix, only returned when empty match
  blockNames	<-	names(nml)
  blckI	<-	c()
  for (i in seq_len(length(blockNames))){
    if (any(argName %in% names(nml[[i]]))){
      blckI	<- c(blckI,i)
    } else {
      one.i <- which(fault.string==fau)[1]
      fault.string[one.i:(one.i+length(names(nml[[i]]))-1)]=names(nml[[i]])
    }

  }
  fault.string <- fault.string[!fault.string==fau] # is empty if found
  # test to see if a block match was made
  if (is.null(blckI)){stop(c("parameter name ",argName," not found in nml. Possible names:",paste(fault.string,collapse=', ')))}
  return(blckI)
}

setnmlList <- function(glm_nml,arg_list){
  if (!is.list(arg_list)){stop("arg_list must be a list")}
  if (any(nchar(names(arg_list)) == 0) | length(names(arg_list)) == 0){
    stop('arg_list must be a named list')
  }
  arg_names  <-	names(arg_list)
  for (i in seq_len(length(arg_names))){
    glm_nml <- set_nml(glm_nml,arg_name=arg_names[i],arg_val=arg_list[[i]])
  }
  return(glm_nml)
}


is_nml_file <- function(nml_file){
  is_nml <- FALSE
  fl_ext <- tail(strsplit(nml_file, "\\.")[[1]],1)
  if (fl_ext == 'nml'){
    is_nml <- TRUE
  }
  return(is_nml)
}

what_ascii <- function(file){
  response <- capture.output(showNonASCIIfile(file))
  return(response)
}

ascii_only <- function(file){
  response <- what_ascii(file)
  if (length(response) > 0){
    return(FALSE)
  } else {
    return(TRUE)
  }
}


get_block <- function(glm_nml, arg_name, warn=TRUE){
  arg_split = strsplit(arg_name,'::')[[1]]
  if (length(arg_split) > 1){
    blck = arg_split[1]
    arg_name = get_arg_name(arg_name)
  } else{
    blck	<-	findBlck(glm_nml,arg_name)
  }
  if(length(blck) > 1){
    if(warn){
      warning(arg_name, " found in ", 
              paste(names(glm_nml[blck]), collapse=' & '), 
              ", returning the first. Try ",
              names(glm_nml[blck])[1],"::",arg_name, " for explicit match")
    }
    blck = blck[1]
  }
  return(blck)
}

get_arg_name <- function(arg_name){
  arg_split = strsplit(arg_name,'::')[[1]]

  if (length(arg_split) > 1){
    blck = arg_split[1]
    arg_name = arg_split[2]
  }
  return(arg_name)
}

.nml <- function(list_obj){
  nml <- list_obj
  class(nml) <- "nml"
  invisible(nml)
}

tan_switch <- function(tan_tf){
  tan_tf <- gsub("atan", "ata*n", tan_tf, fixed  = TRUE)
  nr_of_tans <- gregexpr("tan",  tan_tf) %>% unlist() %>% length()
  tan_switch <- rbinom(nr_of_tans, 1, 0.5)
  for(tan in tan_switch){
    if(tan == 1){
      tan_tf <- sub("tan", "ta*nh", tan_tf, fixed  = TRUE)
    }
  }
  tan_tf <- gsub("ta*n", "tan", tan_tf, fixed  = TRUE)
  return(tan_tf)
}


SCEoptim <- function(FUN, par,
                     lower = -Inf, upper = Inf,
                     control = list(), state_folder = NULL, ...) {
  
  # Initial setup (always) ---------------------------------------------------------------
  FUN <- match.fun(FUN)
  stopifnot(is.numeric(par))
  stopifnot(length(par) > 0)
  stopifnot(is.numeric(lower))
  stopifnot(is.numeric(upper))
  ## allow `lower` or `upper` to apply to all parameters
  if (length(lower) == 1)
    lower <- rep(lower, length = length(par))
  if (length(upper) == 1)
    upper <- rep(upper, length = length(par))
  stopifnot(length(lower) == length(par))
  stopifnot(length(upper) == length(par))
  ## determine number of variables to be optimized
  NDIM <- length(par)
  ## update default options with supplied options
  stopifnot(is.list(control))
  control <- modifyList(sceDefaults(), control)
  isValid <- names(control) %in% names(sceDefaults())
  if (any(!isValid))
    stop("unrecognised options: ",
         toString(names(control)[!isValid]))
  returnpop <- control$returnpop
  trace <- control$trace
  nCOMPLEXES <- control$ncomplex
  CCEITER <- control$cce.iter
  MAXIT <- control$maxit
  MAXEVAL <- control$maxeval
  fnscale <- control$fnscale
  FSO_dims <- control$FSO_dims
  ## recommended number of CCE steps in Duan et al 1994:
  if (is.na(CCEITER))
    CCEITER <- 2 * NDIM + 1
  if (is.finite(MAXEVAL)) {
    ## upper bound on number of iterations to reach MAXEVAL
    MAXIT <- min(MAXIT, ceiling(MAXEVAL / (nCOMPLEXES * CCEITER)))
  }
  ## define number of points in each complex
  nPOINTS_COMPLEX <- 50
  ## define number of points in each simplex
  nPOINTS_SIMPLEX <- 10# NDIM+1
  ## define total number of points
  nPOINTS <- nCOMPLEXES * nPOINTS_COMPLEX
  
  
  
  # function definitions
  costFunction <- function(FUN, par, fnscale) {
    ## check lower and upper bounds
    i <- which(par < lower)
    if (any(i)) {
      i <- i[1]
      return( 1e12 + (lower[i] - par[i]) * 1e6 )
    }
    i <- which(par > upper)
    if (any(i)) {
      i <- i[1]
      return( 1e12 + (par[i] - upper[i]) * 1e6   )
    }
    funevals <<- funevals + 1
    result <- FUN(par) * fnscale
    if (is.na(result))
      result <- 1e12
    result
  }
  simplexStep <- function(P, FAC) {
    ## Extrapolates by a factor FAC through the face of the simplex across from
    ## the highest (i.e. worst) point.
    worst <- nPOINTS_SIMPLEX
    centr <- apply(P[-worst,,drop=FALSE], 2, mean)
    newpar <- centr*(1-FAC) + P[worst,]*FAC
    newpar
  }
  
  
  if(!file.exists(paste0(state_folder, "/optimization_state.rds"))){
    ## initialize population matrix ------------------------------------------------------
    POPULATION <- matrix(as.numeric(NA), nrow = nPOINTS, ncol = NDIM)
    if (!is.null(names(par)))
      colnames(POPULATION) <- names(par)
    POP.FITNESS <- numeric(length = nPOINTS)
    POPULATION[1,] <- par
    ## generate initial parameter values by random uniform sampling
    finitelower <- ifelse(is.infinite(lower), -(abs(par)+2)*5, lower)
    finiteupper <- ifelse(is.infinite(upper), +(abs(par)+2)*5, upper)
    
    if (control$initsample == "latin") {
      for (i in 1:NDIM) {
        tmp <- seq(finitelower[i], finiteupper[i], length = nPOINTS-1)
        tmp <- jitter(tmp, factor = 2)
        tmp <- pmax(finitelower[i], pmin(finiteupper[i], tmp))
        POPULATION[-1,i] <- sample(tmp)
      }
    } else {
      for (i in 1:NDIM)
        POPULATION[-1,i] <- runif(nPOINTS-1, finitelower[i], finiteupper[i])
    }
    
    ## only store all iterations if requested -- could be big!
    if (!is.finite(MAXIT)) {
      MAXIT <- 10000
      warning("setting maximum iterations to 10000")
    }
    if (returnpop) {
      POP.ALL <- array(as.numeric(NA), dim = c(nPOINTS, NDIM, MAXIT))
      if (!is.null(names(par)))
        dimnames(POP.ALL)[[2]] <- names(par)
    } else {
      POP.ALL <- NA
    }
    POP.FIT.ALL <- matrix(as.numeric(NA), ncol = nPOINTS, nrow = MAXIT)
    BESTMEM.ALL <- matrix(as.numeric(NA), ncol = NDIM, nrow = MAXIT)
    if (!is.null(names(par)))
      colnames(BESTMEM.ALL) <- names(par)
    ## the output object
    obj <- list()
    class(obj) <- c("SCEoptim", class(obj))
    obj$call <- match.call()
    obj$control <- control
    ## initialize timer
    tic <- as.numeric(Sys.time())
    toc <- 0
    ## initialize counters
    funevals <- 0
    ## calculate cost for each point in initial population
    for (i in 1:nPOINTS){
      POP.FITNESS[i] <- costFunction(FUN, POPULATION[i,], fnscale)
      if(POP.FITNESS[i] == 1e12){
        POPULATION[i, 1:FSO_dims] <- runif(FSO_dims, lower, upper)
        POP.FITNESS[i] <- costFunction(FUN, POPULATION[i,], fnscale)
      }
      
     }
    ## sort the population in order of increasing function values
    idx <- order(POP.FITNESS)
    POP.FITNESS <- POP.FITNESS[idx]
    POPULATION <- POPULATION[idx,,drop=FALSE]
    ## store one previous iteration only
    POP.PREV <- POPULATION
    POP.FIT.PREV <- POP.FITNESS
    if (returnpop) {
      POP.ALL[,,1] <- POPULATION
    }
    POP.FIT.ALL[1,] <- POP.FITNESS
    BESTMEM.ALL[1,] <- POPULATION[1,]

    ## initiliaize loop variable
    i <- 0
    # save optimizer states
    optimizer_state <- list(
      "tic" = tic,
      "i" = i,
      "POPULATION" = POPULATION,
      "POP.FITNESS" = POP.FITNESS,
      "POP.PREV" = POP.PREV,
      "POP.FIT.PREV" = POP.FIT.PREV,
      "POP.ALL" = POP.ALL,
      "POP.FIT.ALL" = POP.FIT.ALL,
      "BESTMEM.ALL" = BESTMEM.ALL,
      "funevals" = funevals
    )
    saveRDS(optimizer_state, paste0(state_folder, "/optimization_state.rds"))
    saveRDS(obj, paste0(state_folder, "/last_sce_obj.rds"))
    message("Finished computing loss for initial Population")
    return()
  } else {
    # Load states ------------------------------------------------------------------------
    obj <- readRDS(paste0(state_folder, "/last_sce_obj.rds"))
    optimizer_state <- readRDS(paste0(state_folder, "/optimization_state.rds"))
    i <- optimizer_state[["i"]]
    tic <- optimizer_state[["tic"]]
    POPULATION <- optimizer_state[["POPULATION"]]
    POP.FITNESS <- optimizer_state[["POP.FITNESS"]]
    POP.PREV <-  optimizer_state[["POP.PREV"]]
    POP.FIT.PREV <- optimizer_state[["POP.FIT.PREV"]]
    POP.ALL <- optimizer_state[["POP.ALL"]]
    POP.FIT.ALL <- optimizer_state[["POP.FIT.ALL"]]
    BESTMEM.ALL <- optimizer_state[["BESTMEM.ALL"]]
    funevals <- optimizer_state[["funevals"]]
  }
  
  # set flag and message
  EXITFLAG <- NA
  EXITMSG <- NULL
  ## store best solution from last two iterations
  prevBestVals <- rep(Inf, control$tolsteps)
  prevBestVals[1] <- POP.FITNESS[1]
  
  # optimization loop
  while (i < MAXIT) {
    
    i <- i + 1
    
    ## The population matrix POPULATION will now be rearranged into complexes.
    
    ## For each complex ...
    for (j in 1:nCOMPLEXES) {
      
      ## construct j-th complex from POPULATION
      
      k1 <- 1:nPOINTS_COMPLEX
      k2 <- (k1-1) * nCOMPLEXES + j
      
      COMPLEX <- POP.PREV[k2,,drop=FALSE]
      COMPLEX_FITNESS <- POP.FIT.PREV[k2]
      
      ## Each complex evolves a number of steps according to the competitive
      ## complex evolution (CCE) algorithm as described in Duan et al. (1992).
      ## Therefore, a number of 'parents' are selected from each complex which
      ## form a simplex. The selection of the parents is done so that the better
      ## points in the complex have a higher probability to be selected as a
      ## parent. The paper of Duan et al. (1992) describes how a trapezoidal
      ## probability distribution can be used for this purpose.
      
      for (k in 1:CCEITER) {
        
        ## select simplex by sampling the complex
        
        ## sample points with "trapezoidal" i.e. linear probability
        weights <- rev(ppoints(nPOINTS_COMPLEX))
        ## 'elitism' parameter can give more weight to the better results:
        weights <- weights ^ control$elitism
        LOCATION <- sample(seq(1,nPOINTS_COMPLEX), size = nPOINTS_SIMPLEX,
                           prob = weights)
        
        LOCATION <- sort(LOCATION)
        
        ## construct the simplex
        SIMPLEX <- COMPLEX[LOCATION,,drop=FALSE]
        SIMPLEX_FITNESS <- COMPLEX_FITNESS[LOCATION]
        
        worst <- nPOINTS_SIMPLEX
        
        ## generate new point for simplex
        
        ## first extrapolate by a factor -1 through the face of the simplex
        ## across from the high point,i.e.,reflect the simplex from the high point
        parRef <- simplexStep(SIMPLEX, FAC = -1)
        fitRef <- costFunction(FUN, parRef, fnscale)
        
        ## check the result
        if (fitRef <= SIMPLEX_FITNESS[1]) {
          ## gives a result better than the best point,so try an additional
          ## extrapolation by a factor 2
          parRefEx <- simplexStep(SIMPLEX, FAC = -2)
          fitRefEx <- costFunction(FUN, parRefEx, fnscale)
          if (fitRefEx < fitRef) {
            SIMPLEX[worst,] <- parRefEx
            SIMPLEX_FITNESS[worst] <- fitRefEx
            ALGOSTEP <- 'reflection and expansion'
          } else {
            SIMPLEX[worst,] <- parRef
            SIMPLEX_FITNESS[worst] <- fitRef
            ALGOSTEP <- 'reflection'
          }
        } else if (fitRef >= SIMPLEX_FITNESS[worst-1]) {
          ## the reflected point is worse than the second-highest, so look
          ## for an intermediate lower point, i.e., do a one-dimensional
          ## contraction
          parCon <- simplexStep(SIMPLEX, FAC = -0.5)
          fitCon <- costFunction(FUN, parCon, fnscale)
          if (fitCon < SIMPLEX_FITNESS[worst]) {
            SIMPLEX[worst,] <- parCon
            SIMPLEX_FITNESS[worst] <- fitCon
            ALGOSTEP <- 'one dimensional contraction'
          } else {
            ## can't seem to get rid of that high point, so better contract
            ## around the lowest (best) point
            SIMPLEX <- (SIMPLEX + rep(SIMPLEX[1,], each=nPOINTS_SIMPLEX)) / 2
            for (k in 2:nrow(SIMPLEX))
              SIMPLEX_FITNESS[k] <- costFunction(FUN, SIMPLEX[k,], fnscale)
            ALGOSTEP <- 'multiple contraction'
          }
        } else {
          ## if better than second-highest point, use this point
          SIMPLEX[worst,] <- parRef
          SIMPLEX_FITNESS[worst] <- fitRef
          ALGOSTEP <- 'reflection'
        }

          message(ALGOSTEP)
        
        ## replace the simplex into the complex
        COMPLEX[LOCATION,] <- SIMPLEX
        COMPLEX_FITNESS[LOCATION] <- SIMPLEX_FITNESS
        
        ## sort the complex
        idx <- order(COMPLEX_FITNESS)
        COMPLEX_FITNESS <- COMPLEX_FITNESS[idx]
        COMPLEX <- COMPLEX[idx,,drop=FALSE]
      }
      
      ## replace the complex back into the population
      POPULATION[k2,] <- COMPLEX
      POP.FITNESS[k2] <- COMPLEX_FITNESS
    }
    
    
    ## At this point, the population was divided in several complexes, each of which
    ## underwent a number of iteration of the simplex (Metropolis) algorithm. Now,
    ## the points in the population are sorted, the termination criteria are checked
    ## and output is given on the screen if requested.
    
    ## sort the population
    idx <- order(POP.FITNESS)
    POP.FITNESS <- POP.FITNESS[idx]
    POPULATION <- POPULATION[idx,,drop=FALSE]
    if (returnpop) {
      POP.ALL[,,i] <- POPULATION
    }
    POP.FIT.ALL[i,] <- POP.FITNESS
    BESTMEM.ALL[i,] <- POPULATION[1,]
    
    curBest <- POP.FITNESS[1]
    
    ## end the optimization if one of the stopping criteria is met
    
    prevBestVals <- c(curBest, head(prevBestVals, -1))
    reltol <- control$reltol
    if (all(abs(diff(prevBestVals)) <= reltol * (abs(curBest)+reltol))) {
      EXITMSG <- 'Change in solution over [tolsteps] less than specified tolerance (reltol).'
      EXITFLAG <- 0
    }
    
    ## give user feedback on screen if requested
    if (trace >= 1) {
      if (i == 1) {
        message(' Nr Iter  Nr Fun Eval    Current best function    Current worst function')
      }
      if ((i %% control$REPORT == 0) || (!is.na(EXITFLAG)))
      {
        message(sprintf(' %5.0f     %5.0f             %12.6g              %12.6g',
                        i, funevals, min(POP.FITNESS), max(POP.FITNESS)))
        if (trace >= 2)
          message("parameters: ", toString(signif(POPULATION[1,], 3)))
      }
    }
    
    if (!is.na(EXITFLAG))
      break
    
    if ((i >= control$maxit) || (funevals >= control$maxeval)) {
      EXITMSG <- 'Maximum number of function evaluations or iterations reached.'
      EXITFLAG <- 1
      break
    }
    
    toc <- as.numeric(Sys.time()) - tic
    if (toc > control$maxtime) {
      EXITMSG <- 'Exceeded maximum time.'
      EXITFLAG <- 2
      break
    }
    
    ## go to next iteration
    POP.PREV <- POPULATION
    POP.FIT.PREV <- POP.FITNESS
    
    # save optimizer states
    optimizer_state <- list(
      "tic" = tic,
      "i" = i,
      "POPULATION" = POPULATION,
      "POP.FITNESS" = POP.FITNESS,
      "POP.PREV" = POP.PREV,
      "POP.FIT.PREV" = POP.FIT.PREV,
      "POP.ALL" = POP.ALL,
      "POP.FIT.ALL" = POP.FIT.ALL,
      "BESTMEM.ALL" = BESTMEM.ALL,
      "funevals" = funevals
    )
    ## return solution
    obj$par <- POPULATION[1,]
    obj$value <- POP.FITNESS[1]
    obj$convergence <- EXITFLAG
    obj$message <- EXITMSG
    
    ## store number of function evaluations
    obj$counts <- funevals
    ## store number of iterations
    obj$iterations <- i
    ## store the amount of time taken
    obj$time <- toc
    
    saveRDS(optimizer_state, paste0(state_folder, "/optimization_state.rds"))
    saveRDS(obj, paste0(state_folder, "/last_sce_obj.rds"))
    if(i < MAXIT) return()
  }
  
  
  if (trace >= 1)
    message(EXITMSG)
  
  ## return solution
  obj$par <- POPULATION[1,]
  obj$value <- POP.FITNESS[1]
  obj$convergence <- EXITFLAG
  obj$message <- EXITMSG
  
  ## store number of function evaluations
  obj$counts <- funevals
  ## store number of iterations
  obj$iterations <- i
  ## store the amount of time taken
  obj$time <- toc
  
  if (returnpop) {
    ## store information on the population at each iteration
    obj$POP.ALL <- POP.ALL[,,1:i]
    dimnames(obj$POP.ALL)[[3]] <- paste("iteration", 1:i)
  }
  obj$POP.FIT.ALL <- POP.FIT.ALL[1:i,]
  obj$BESTMEM.ALL <- BESTMEM.ALL[1:i,]
  saveRDS(obj, paste0(state_folder, "/last_sce_obj.rds"))
  
  return(obj)
}

# Runs mHM with tfs given as strings
run_mhm_from_ksat <- function(KSat_tf, run_folder){
  path <- paste0(run_folder, "/FSO_mHM_major_basins/config")
  basin_folders <- list.files(path, pattern = "sub_")
  
  variables <- c("bd", "sand", "clay", "slope", "aspect", "ThetaS", "KSat",
                 "vGenu_n", "dem")
  alltill_variables <- c("bd", "sand", "clay", "ThetaS", "KSat", "vGenu_n", "aspect", "dem", "slope")
  till_variables <- c("bd_till", "sand_till", "clay_till", "ThetaS_till", "KSat_till",
                      "vGenu_n_till", "slope_till", "aspect_till", "dem_till")
  notill_variables <- c("bd_notill", "sand_notill", "clay_notill", "ThetaS_notill",
                        "KSat_notill", "vGenu_n_notill", "slope_notill", "aspect_notill", "dem_notill")
  scaling_bounds <- list("slope" = c(0, 90),
                         "aspect" = c(0, 360),
                         "bd" = c(0, 2.3),
                         "sand" = c(0, 100),
                         "clay" = c(0, 100),
                         "dem" = c(0, 4000),
                         "KSat" = c(9, 244),
                         "vGenu_n" = c(1,2),
                         "ThetaS" = c(0.24, 0.51))
  
  # Add scaling to given transfer functions
  # Ksat
  KSat_from <- as.numeric(
    read.csv("/home/lv71468/mfeigl/FSO_mHM/FSO_Data/VAE_data/KSat/distribution_cutoff_points.csv"))
  KSat_to <- c(1.1, 1000.0)
  a <- (KSat_from[2] - KSat_from[1]) %>% format(nsmall = 1)
  b <- (KSat_to[2] - KSat_to[1]) %>% format(nsmall = 1)
  KSat_from[1] <- KSat_from[1] %>% format(nsmall = 1)
  KSat_to[1] <- KSat_to[1] %>% format(nsmall = 1)
  Ksat_fun <- paste0(KSat_to[1], "+((", KSat_tf, "-", KSat_from[1], ")", "*", b, ")/", a)
  Ksat_fun <- gsub("--", "+", Ksat_fun, fixed = TRUE)
  prepared_KSat_tf <- prepare_tf_for_mhm(Ksat_fun, scaling_bounds = scaling_bounds)
  # Ksat
  Ksat_fun <- prepared_KSat_tf$scaled_function
  KSat_nums <- prepared_KSat_tf$numerics
  # make numeric vectors with parameter names
  numeric_paras <- KSat_nums %>% unique()
  numeric_paras_names <- paste0("FSO", seq_along(numeric_paras))
  numeric_paras <- numeric_paras[order(nchar(numeric_paras), decreasing = TRUE)]
  # fill in numerics in functions
  for(k in seq_along(numeric_paras)){
    Ksat_fun <- gsub(numeric_paras[k], numeric_paras_names[k], Ksat_fun)
  }
  
  # change KSat variables to _horizon if necessary
  for(var in c("aspect", "slope", "dem")){
    Ksat_fun <- gsub(var, paste0(var, "_horizon"), Ksat_fun)
  }
  
  # 1. Change transfer functions in nml
  mpr_nml <- readLines("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_master_namelists/mpr_fso.nml")
  # change KSat tf
  KSat_line <- grep('KSat_dummy', mpr_nml)
  mpr_nml[KSat_line] <- paste0("                           '", Ksat_fun, "',")
  fun_vars <- unlist(strsplit(gsub(" ", "", Ksat_fun), "(?=[+-/*)()])", perl = TRUE))
  horizon_vas <- c("aspect_horizon", "slope_horizon", "dem_horizon")
  used_var <- unique(fun_vars[fun_vars %in% c(variables, horizon_vas)])
  data_array_line <- grep("from_data_arrays(1:3,26)", mpr_nml, fixed = TRUE)
  mpr_nml[data_array_line] <- paste0('    from_data_arrays(1:', length(used_var), ',26',
                                     ') = ',
                                     paste0("'", paste0(used_var, collapse = "', '"), "'"),
                                     collapse = "")
  
  # change FieldCap tfs
  FieldCap_till <- "ThetaS_till * exp(FieldCap_c1 * (FieldCap_c2 + log10(KSat_till)) * log(vGenu_n_till))"
  FieldCap_notill <- "ThetaS_notill * exp(FieldCap_c1 * (FieldCap_c2 + log10(KSat_notill)) * log(vGenu_n_notill))"
  FieldCap_till_line <- grep('FieldCap_till_dummy', mpr_nml)
  FieldCap_notill_line <- grep('FieldCap_notill_dummy', mpr_nml)
  mpr_nml[FieldCap_till_line] <- paste0("                           '", FieldCap_till, "',")
  mpr_nml[FieldCap_notill_line] <- paste0("                           '", FieldCap_notill, "',")
  
  # data array FieldCap_till
  fun_vars <- unlist(strsplit(gsub(" ", "", FieldCap_till), "(?=[+-/*)()])", perl = TRUE))
  used_var <- unique(fun_vars[fun_vars %in% till_variables])
  data_array_line <- grep("from_data_arrays(1:3,37)", mpr_nml, fixed = TRUE)
  mpr_nml[data_array_line] <- paste0('    from_data_arrays(1:', length(used_var), ',37',
                                     ') = ',
                                     paste0("'", paste0(used_var, collapse = "', '"), "'"),
                                     collapse = "")
  # change FieldCap tfs
  # data array FieldCap_notill
  fun_vars <- unlist(strsplit(gsub(" ", "", FieldCap_notill), "(?=[+-/*)()])", perl = TRUE))
  used_var <- unique(fun_vars[fun_vars %in% notill_variables])
  data_array_line <- grep("from_data_arrays(1:3,38)", mpr_nml, fixed = TRUE)
  mpr_nml[data_array_line] <- paste0('    from_data_arrays(1:', length(used_var), ',38',
                                     ') = ',
                                     paste0("'", paste0(used_var, collapse = "', '"), "'"),
                                     collapse = "")
  
  # get parameter names, parameter_values
  parameter_name_line <- grep('parameter_names', mpr_nml)
  parameter_value_line <- grep('parameter_values', mpr_nml)
  parameters_end <- grep('&data_arrays', mpr_nml, fixed = TRUE) - 3
  parameter_names_initial <- mpr_nml[parameter_name_line:(parameter_value_line-1)]
  parameter_values_initial <- mpr_nml[parameter_value_line:parameters_end]
  # add new parameter names and values
  nlines <- length(parameter_names_initial)
  parameter_names_initial[nlines] <- paste0(parameter_names_initial[nlines], ",")
  parameter_names <- c(
    parameter_names_initial,
    paste0("                            ",
           paste0(paste0("'", numeric_paras_names, "'"), collapse = ", ")))
  nlines <- length(parameter_values_initial)
  parameter_values_initial[nlines] <- paste0(parameter_values_initial[nlines], ",")
  parameter_values <- c(
    parameter_values_initial,
    paste0("                            ", paste0(numeric_paras, collapse = ", ")))
  old_nr_of_numeric_paras <- mpr_nml[parameter_name_line] %>%
    strsplit(split = "(", fixed = TRUE) %>%
    unlist() %>%
    substr(., 3, 4) %>% paste0(., collapse = "") %>%
    as.numeric()
  total_nr_of_numeric_paras <- old_nr_of_numeric_paras + length(numeric_paras)
  parameter_names[1] <- paste0("    parameter_names(1:",total_nr_of_numeric_paras, ") = ",
                               strsplit(parameter_names[1], " = ", fixed = TRUE)[[1]][2])
  parameter_values[1] <- paste0("    parameter_values(1:",total_nr_of_numeric_paras, ") = ",
                                strsplit(parameter_values[1], " = ", fixed = TRUE)[[1]][2])
  
  # add new parmeter names and values
  new_mpr_nml <- c(mpr_nml[1:(parameter_name_line-1)],
                   parameter_names,
                   parameter_values,
                   mpr_nml[(parameters_end+1):length(mpr_nml)])
  
  # save KSat netcf
  to_file_start <- grep("to_file(", mpr_nml, fixed = TRUE)
  mpr_nml[to_file_start] # KSat: 28, 29
  mpr_nml[to_file_start + 4] <- "                     .true., .true., .false., .false., .false., .false.,"
  
  # Change numeric parameters
  parameters <- suppressWarnings(
    read_nml("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_master_namelists/mhm_parameter.nml"))

  # write new mpr.nml and mhm_parameter.nml for all basins
  for(basin in basin_folders){
    basin_nml <- gsub("<BASIN>", basin, new_mpr_nml)
    writeLines(basin_nml, con = paste0(path, "/", basin, "/mpr.nml"))
    write_mhm_nml(parameters, paste0(path, "/", basin, "/mhm_parameter.nml"))
  }
  
  # 2. Run python script and compile mhm
  output_error <- try({
    system(paste0("chmod 775 ", run_folder, "/05_run_mhm_preparation.sh"))
    system(paste0("bash ", run_folder, "/05_run_mhm_preparation.sh"))
    for(basin_nr in 1:7){
      system(paste0("chmod 775 ", run_folder, "/mhm_run_scripts/run_basin_", basin_nr, ".sh"))
    }
    cat("running mhm\n")
    system(paste0("bash ", run_folder, "/mhm_run_scripts/run_basin_1.sh &> ", run_folder, "/mhm_run_scripts/run_basin_1.log & ",
                  "bash ", run_folder, "/mhm_run_scripts/run_basin_2.sh &> ", run_folder, "/mhm_run_scripts/run_basin_2.log & ",
                  "bash ", run_folder, "/mhm_run_scripts/run_basin_3.sh &> ", run_folder, "/mhm_run_scripts/run_basin_3.log & ",
                  "bash ", run_folder, "/mhm_run_scripts/run_basin_4.sh &> ", run_folder, "/mhm_run_scripts/run_basin_4.log & ",
                  "bash ", run_folder, "/mhm_run_scripts/run_basin_5.sh &> ", run_folder, "/mhm_run_scripts/run_basin_5.log & ",
                  "bash ", run_folder, "/mhm_run_scripts/run_basin_6.sh &> ", run_folder, "/mhm_run_scripts/run_basin_6.log & ",
                  "process_id=$! bash ", run_folder, "/mhm_run_scripts/run_basin_7.sh &> ", run_folder,
                  "/mhm_run_scripts/run_basin_7.log & wait $process_id"
    ))
    cat("Finished running mhm!\n")
    # read outputs and get KGE/NSE
    basin_results <- data.frame("Basin" = basin_folders, "KGE" = NA, "NSE" = NA, "lNSE" = NA)
    for(basin in basin_folders){
      output <- readLines(paste0(path, "/", basin, "/output.txt"))
      lines <- grep("KGE of daily discharge", output)
      KGE <- as.numeric(substring(output[lines], 39))
      NSE <- as.numeric(substring(output[lines + 1], 39))
      lNSE <- as.numeric(substring(output[lines + 2], 39))
      try(basin_results[basin_results$Basin == basin, c("KGE", "NSE", "lNSE")] <- c(KGE, NSE, lNSE))
    }
  })
  # If mHM  produces no output return NA
  if(class(output_error) == "try-error"){
    basin_results <- data.frame("Basin" = basin_folders, "KGE" = NA, "NSE" = NA)
  }
  
  return(basin_results)
}

SPAEF <- function(observations, predictions){
  spaef_try <- try({
    #observations <- observations[is.finite(observations)]
    #predictions <- predictions[is.finite(predictions)]
    alpha <- cor(predictions, observations)
    beta <- (sd(observations)/mean(observations)) / (sd(predictions)/mean(predictions))
    
    # scale for histogram distance
    observations <- scale(observations)
    predictions <- scale(predictions)
    
    range_true <- max(observations) - min(observations)
    breaks <- seq(min(observations), max(observations), range_true/100)
    c1 <- as.integer(table(cut(observations, breaks = breaks)))
    c2 <- as.integer(table(cut(predictions, breaks = breaks)))
    
    c_min <- numeric()
    for(i in seq_along(c1)) c_min <- c(c_min, min(c1[i], c2[i]))
    gamma <- sum(c_min)/sum(c1)
    spaef <- 1 - sqrt((alpha-1)^2 + (beta-1)^2 + (gamma-1)^2)
  }, silent = TRUE)
  if(class(spaef_try) == "try-error") spaef <- NA
  return(spaef)
}





