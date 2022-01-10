# FSO-mHM start optimization runs
# FSO-mHM project
# Moritz Feigl, Nov 2020

setwd("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt")
dir.create("/home/lv71468/mfeigl/FSO_mHM/Runs", showWarnings = FALSE)
# Create run scripts for FSO-mHM runs
runs <-6:10
numIter <- 3000
iterations_per_run <- 70# floor(3*24*2.8) # DDS
gcc <- "_gcc" # "" for intel
optimizer <- "" # "" for DDS


for(run in runs){
  run_name <- paste0("FSO_DDS_NSE_run_", run)
  if(optimizer == "_ga") run_name <- paste0("FSO_GA_NSE_run_", run)
  if(optimizer == "_sce") run_name <- paste0("FSO_SCE_NSE_run_", run)
  run_folder <- paste0("/home/lv71468/mfeigl/FSO_mHM/Runs/", run_name)
  # create FSO opt script
  dir.create(run_folder, showWarnings = FALSE)
  if(optimizer != ""){
    file.copy(
      paste0("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/02_mHM_FSO",
             optimizer, ".R"),
      paste0(run_folder, "/02_mHM_FSO.R"), overwrite = TRUE)

  } else {
    file.copy("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/02_mHM_FSO.R",
              run_folder, overwrite = TRUE)
  }
  script <- readLines(paste0(run_folder, "/02_mHM_FSO.R"))
  script <- gsub('run_name <- "dummy"',
                 paste0('run_name <- "', run_name, '"'),
                 script, fixed = TRUE)

  script <- sub('iterations_per_run <- 2',
                paste0("iterations_per_run <- ", iterations_per_run),
                script, fixed = TRUE)

  script <- sub('numIter <- 4',
                paste0("numIter <- ", numIter),
                script, fixed = TRUE)

  script <- gsub('gcc <- ""',
                 paste0('gcc <- "', gcc, '"'),
                 script, fixed = TRUE)

  writeLines(script, paste0(run_folder, "/02_mHM_FSO.R"))
  setwd(run_folder)
  # get sbatch script
  file.copy("/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/06_start_fso_mhm.sh",
            "06_start_fso_mhm.sh", overwrite = TRUE)
  start_fso_mhm <- readLines(paste0(run_folder, "/06_start_fso_mhm.sh"))
  if(gcc == "_gcc"){
    start_fso_mhm <- start_fso_mhm[-grep("#SBATCH -L intel@vsc", start_fso_mhm, fixed = TRUE)]
  }
  start_fso_mhm <- gsub("#SBATCH -J FSO_mhm",
                        paste0("#SBATCH -J ", run_name), start_fso_mhm)
  start_fso_mhm <- gsub("cd /home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt",
                        paste0("cd ", run_folder), start_fso_mhm)
  writeLines(start_fso_mhm, paste0(run_folder, "/06_start_fso_mhm.sh"))
  # Get last R codes
  code_path <- "/home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/"
  file.copy(paste0(code_path, "03_vae_generators.R"),
            "03_vae_generators.R", overwrite = TRUE)
  file.copy(paste0(code_path, "04_FSO_mhm_utils.R"),
            "04_FSO_mhm_utils.R", overwrite = TRUE)
  
  # start optimization
  system("chmod 775 06_start_fso_mhm.sh")
  system("sbatch 06_start_fso_mhm.sh")
}


for(run in runs){
  run_name <- paste0("FSO_SCE_NSE_run_", run)
  run_folder <- paste0("/home/lv71468/mfeigl/FSO_mHM/Runs/", run_name)
  setwd(run_folder)
  system("sbatch 06_start_fso_mhm.sh")
}
