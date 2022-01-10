#!/bin/bash
cd /home/lv71468/mfeigl/FSO_mHM/2020_fso_mhm/scripts/08_FSO_mHM_opt/mhm_scripts
parallel -j0 bash :::: <(ls run_basin_{1..7}.sh)
