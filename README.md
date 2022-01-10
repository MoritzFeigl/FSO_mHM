[![DOI](https://zenodo.org/badge/446423033.svg)](https://zenodo.org/badge/latestdoi/446423033)

# Automatic regionalization of model parameters for hydrological models

<p align="center">
  <img width="460" src="https://github.com/MoritzFeigl/FSO_mHM/blob/master/fso.png">
</p>

Accompanying code for the publication "Automatic regionalization of model parameters for hydrological models" by Moritz Feigl, Stephan Thober, Robert Schweppe, Mathew Herrnegger, Luis Samaniego and Karsten Schulz.

The code in this repository was used to produce all results and figures in our manuscript. The data can be aquired from http://www.ufz.de/index.php?en=41160.

### Content of the repository


### Scripts
- `scripts` All scripts used to generate the paper results. This includes code for the mHM setup, FSO VAE training and generator functions, optimization routines, additional helper functions, model optimization and results analysis.
- `tests` Test for the mHM data preparation.
- `utils` utils functions for the nHN data preparation.

### mHM environment
The mHM environment is based on the setups of German basins based on Zink et al., 2017.
We need convert this environment to one that runs with a [version of mHM](https://git.ufz.de/ottor/mhm/-/tree/22-couple-mhm-with-mpr-1-0) with MPR (version 0.6.5) implemented.
For that the script [prepare_env.py](scripts/02_prepare_env_mHM/prepare_env.py) is run. There are some settings required:
- mHM compilation:
  - a working compiler and the netcdf-fortran library need to be loaded
  - this can be saved to a script, which can later be `source`d again: `deps/mpr/moduleLoadScripts/eve.gfortran83MPI`
  - then do:
      ```
      mkdir build
      cd build
      cmake ..
      make -j 8
      ```
- after that edit the prepare_env.py script (Global variables at top):
  - PROCESSES sets the parts of the script to be executed (configuration files for mHM and MPR based on ROOT_REPO_FOLDER)
  - BASIN_SELECTION sets the basins to be processed
  - USE_SYMLINKS sets whether symlinks are used or copies
  - ROOT_TARGET_FOLDER, ROOT_SOURCE_FOLDER and ROOT_REPO_FOLDER set the folders to work on
  - EXECUTABLE_PATH sets the path to the executable relative to ROOT_REPO_FOLDER
  - there are many more specifications, that can be tweaked

  ## License

  [Apache License 2.0](https://github.com/MoritzFeigl/FSO_mHM/blob/master/LICENSE)
