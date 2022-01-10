# FSO coupling to mHM
Project by Moritz Feigl, Robert Schweppe and Stephan Thober on the application of FSO to mHM. 

## Overview
See the [project plan](docs/project_plan/FSO_mHM.pdf) for details.
[Initial results](docs/presentations/Feigl_etal_2020_EGU2020.pdf) were presented during the EGU 2020 meeting.

## Documentation
See the [protocols](docs/protocols) for a history of project meetings.

## Scripts
See the [scripts](scripts) folder for an overview of processing scripts.

### mHM environment
The mHM environment is based on the setups of German basins based on Zink et al., 2017.
We need convert this environment to one that runs with a [version of mHM](https://git.ufz.de/ottor/mhm/-/tree/22-couple-mhm-with-mpr-1-0) with MPR (version 0.6.5) implemented (commit d90998f6aff96a5a25a05a3d5f0ac12921700d35).
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
  
### MPR environment
The script [concat_static_files.py](scripts/03_prepare_env_MPR/concat_static_files.py) is run.
It combines all static files found in the static folders into one big dataset for each variable.
These can then be used for MPR parallelization.

#### MPR subdomains

The script [split_domain.py](utils/split_domain.py) can cut domains in subdomains. It uses the libraries xarray and pyresample. It can be called on the command line.
Its usage is as follows:
- Provide the filename of file to be split ("--filename") and optionally the variable name "--variable" if only 1 variable shall be processed.
- Provide the target directory, where the files are saved to "--target_dir".
- If the coordinate information and projection information are not provided in the file according to the CF conventions, you need to provide them seperately ("--grid_kwargs"). The string "--grid_kwargs area_id=EPSG:4326 projection=EPSG:4326 shape=1800,3600 area_extent=-180,-90,180,90" is used parsed as following (keyword=value pairs):
    - area_id: arbitrary name of projection
    - projection: EPSG code of projection
    - shape: shape of grid
    - area_extent: area covered (xmin,ymin,xmax,ymax)
  
  See the [documentation of pyresample](https://pyresample.readthedocs.io/en/latest/geometry_utils.html#areadefinition-creation) for details.
- Provide the names of the x and y dimensions if there are other dimensions or the data are not (y, x)-ordered with the args "--x" and "--y".
- Subdomains can be created in 3 ways:
  1. by mask
     - provide "--maskfile" pointing to the netcdf file that has a variable "ids" with unique identifiers for all subdomains. Its shape must match that of the data variable.
  2. by a target number of subdomains
     - provide "--n_subdomains" for the requested number of subdomains to be created and "--split_dims" for the dimensions to be split over (`'x'`, `'y'`, `'x,y'`). `'x,y'` splits the domain in blocks.
  3. by a given subdomain size (and optionally a lower left corner)
     - provide "--dxdy" for the subdomain sizes and optionally "--ll" for the lower left corner coordinates and "--split_dims" for the dimensions to be split over (`'x'`, `'y'`, `'x,y'`). `'x,y'` splits the domain in blocks. 
     _Note: You must always provide a pair of values, even if one is not used when splitting over only one dimension._
  - if "--trim_to_mask" is provided, each subdomain is trimmed to the bounding box of valid data, if there are none, no file is written.

See details with `python split_domain.py --help`.

The script [merge_domain.py](utils/merge_domain.py) can reaggregate all files split previously. It can be called on the command line.
Its usage is as follows:
- Provide the source directory ("source_dir"), where the files with the pattern ("--pattern") exist.
- Provide the filename of target file to be written ("--target_file").
- If you used a mask during subdomain creation, you need to indicate that by providing "--used_mask".
  