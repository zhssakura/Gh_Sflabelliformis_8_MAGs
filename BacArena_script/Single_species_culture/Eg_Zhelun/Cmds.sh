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

cd # PATH of the Rscript_BacArena.R


# dir_model: full path of your model
dir_model=/srv/scratch/z5095298/sponge_modeling/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/OTU06_renamed/20210515_STY_Merged_OTU06_mineral_sw_3_NoH2S/

# infile_model: file name of model
infile_model = your_model_name.RDS

# dir_diet: full path of your medium/diet
dir_diet=/srv/scratch/z5095298/software/gapseq/Diets_Sponges_July_2019/Diet_Sponges_July_2019_dl20190807_Shan_v20210425/Diet_Sponges_July_2019_dl20190807_Shan_v20210425_brief_v2_IDcheck_mineral_sw_210514/Diet_for_BacArena/mineral_sw_210908_autotroph_BacArena_realistic_simplified_999/

# infile_diet: file name of medium
infile_diet = your_medium_name.RDS

# any keyword that can be used for description
key_word=OTU06_20210515__TrainDiet_NO2_NoH2S__mtf_lmtGrow_999_photosymbiont_ExInfT_NO2_H2S

# You can also edit parameters for --arena_mn,  --death_r, --inocc_no, --cl_no, --iter in the cmd line below:
# --arena_mn: An integer indicating the length of an arena, defalt: 20
# --death_r: A percentage of biomass reduce due to the nutrient limitation, defalt: 0
# --inocc_no: inocculum, defalt: 0.05
# --cl_no: Number of replicates, defalt: 4
# --iter: Number of iteration, defalt:200


# Cmd line:
Rscript Rscript_BacArena.R --otu_name OTU06 --infile_my_model ${dir_model}/${infile_model} --infile_BacarenaRDS_diet ${dir_diet}/${infile_diet} --keywd ${key_word} --arena_mn 20 --death_r 0 --inocc_no 0.05 --cl_no 1 --iter 10
