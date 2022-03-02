#PBS -l nodes=1:ppn=1
#PBS -l mem=200gb
#PBS -l walltime=5:59:00
#PBS -M zzfanyi@gmail.com
#PBS -m ae

#################### BacArena #################### 

module load gcc/6.2.0
module load parallel/20190522
module load openmpi/4.0.1-intel
module load R/3.6.3 # with packageVersion('BacArena') # 1.8.2 (installed from CRAN)
module load glpk/4.65
module load cplex/12.9.0-academic

cd /srv/scratch/z5095298/sponge_modeling/BacArena/BacArena_STY_v20210515_2autotrophs_OTU06,08_400grids_mineral_sw/Single_specie/Katana/OTU6_Inbatch_differ_nitrite_con_use_optparse_mtf/
# For single species analysis: only NO3, no NH3
dir_diet=/srv/scratch/z5095298/software/gapseq/Diets_Sponges_July_2019/Diet_Sponges_July_2019_dl20190807_Shan_v20210425/Diet_Sponges_July_2019_dl20190807_Shan_v20210425_brief_v2_IDcheck_mineral_sw_210514/Diet_for_BacArena/mineral_sw_210908_autotroph_BacArena_realistic_simplified_999/
infile_diet=mineral_sw_210514_OTU08_3_BacArena_realistic_simplified_999_NH3only.RDS

dir_model=/srv/scratch/z5095298/sponge_modeling/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/OTU06_renamed/20210515_STY_Merged_OTU06_mineral_sw_3_NoH2S/
key_word=OTU06_20210515__TrainDiet_NO2_NoH2S__mtf_lmtGrow_999_photosymbiont_ExInfT_NO2_H2S

Rscript Rscript_BacArena_8taxa_sulfur_400grids_200iter__simuDiet_optparse_mtf_limitgrow_ExInfT_ka_OTU6.R --otu_name OTU06 --infile_my_model ${dir_model}STY_Merged_OTU06.RDS --no2_con 3.00E-04 --infile_BacarenaRDS_diet ${dir_diet}${infile_diet} --keywd ${key_word} --arena_mn 20 --death_r 0 --inocc_no 0.05 --cl_no 1 --iter 10
