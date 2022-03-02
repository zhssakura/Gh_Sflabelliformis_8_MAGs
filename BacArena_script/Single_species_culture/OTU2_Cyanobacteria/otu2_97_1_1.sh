#!/bin/bash

#PBS -l nodes=1:ppn=4
#PBS -l mem=200gb
#PBS -l walltime=11:59:00
#PBS -M zzfanyi@gmail.com
#PBS -m ae

#################### BacArena #################### 

module load intel/19.0.0.117
module load gcc/6.2.0
module load parallel/20190522
module load openmpi/4.0.1-intel
module load R/3.6.3 # with packageVersion('BacArena') # 1.8.2 (installed from CRAN)
module load glpk/4.65
# module load cplex/12.9.0
module load cplex/12.9.0-academic

cd /srv/scratch/z5095298/sponge_modeling/BacArena/BacArena_STY_v20210515_2autotrophs_OTU06,08_400grids_mineral_sw/Single_specie/Katana/OTU2_Inbatch_differ_nitrate_con_use_optparse_mtf/

Rscript Rscript_BacArena_OTU2_nitrate_400grids_200iter__simuDiet_optparse_mtf_ka_v210907.R --no3_con 3.00E-04 --iter 6 --arena_mn 20 --death_r 0 --inocc_no 0.05 --cl_no 4 
