#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l mem=100gb
#PBS -l walltime=1:59:00
#PBS -M zzfanyi@gmail.com
#PBS -m ae

#################### BacArena #################### 
#@2022-06-15
module load gcc/6.2.0
module load parallel/20190522
module load openmpi/4.0.1-intel
module load R/3.6.3 # with packageVersion('BacArena') # 1.8.2 (installed from CRAN)
# module load glpk/4.65
module load cplex/12.9.0-academic

cd /srv/scratch/z5095298/sponge_modeling/BacArena/BacArena_STY_v20210520_3autotrophs_OTU2n6n8_400grids_mineral_sw/Katana/R_3.6.3__BacArena_1.8.2_CRAN/mineral_sw_210716_8taxa_BacArena_simplify_ExImf_false_fba/Rscript_BacArena__Best_models_STY_8taxa_20220615__mtf_LxLyLMA_difspeed2_ChemotxNOB_Newdiet4_noH1/
