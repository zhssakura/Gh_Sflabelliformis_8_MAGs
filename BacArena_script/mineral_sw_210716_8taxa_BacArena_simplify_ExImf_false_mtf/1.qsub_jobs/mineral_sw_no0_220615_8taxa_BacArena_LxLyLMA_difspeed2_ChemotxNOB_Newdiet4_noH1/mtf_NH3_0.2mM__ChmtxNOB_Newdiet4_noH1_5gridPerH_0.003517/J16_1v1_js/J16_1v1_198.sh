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

Rscript Rscript_BacArena__Best_models_STY_8taxa_20220615__mtf_rmRate_tp_ExInf_replenish_growlimit_parL_diffModel_addSubsF_auto2.R --auto_num 198 --otu_name STY8taxa --difspeed 0.003517 --keywd J16_1v1 --NH3_con 1 --HT_con 1 --photon_con 5000 --arena_mn 40 --death_r 5 --inocc_no 0.8 --cl_no 1 --iter 1 --tstep 1 --rm_rate 5 --setAllExInf_value TRUE  --NutConPercent 20 --infile_BacarenaRDS_diet /srv/scratch/z5095298/software/gapseq/Diets_Sponges_July_2019/Diet_Sponges_July_2019_dl20190807_Shan_v20210425/Diet_Sponges_July_2019_dl20190807_Shan_v20210425_brief_v2_IDcheck_mineral_sw_210514/Diet_for_BacArena/mineral_sw_210514_SOBs_BacArena_realistic_simplified_999/mineral_sw_220615b_8taxa_BacArena_realistic_simplified_999_HT_ph100_x1000_merge.RDS
cd /srv/scratch/z5095298/qsub_BA_auto/qsub_diffModels_addSubsF_auto/mineral_sw_no0_220615_8taxa_BacArena_LxLyLMA_difspeed2_ChemotxNOB_Newdiet4_noH1/mtf_NH3_0.2mM__ChmtxNOB_Newdiet4_noH1_5gridPerH_0.003517/J16_1v1_js
qsub /srv/scratch/z5095298/qsub_BA_auto/qsub_diffModels_addSubsF_auto/mineral_sw_no0_220615_8taxa_BacArena_LxLyLMA_difspeed2_ChemotxNOB_Newdiet4_noH1/mtf_NH3_0.2mM__ChmtxNOB_Newdiet4_noH1_5gridPerH_0.003517/J16_1v1_js//J16_1v1_199.sh
