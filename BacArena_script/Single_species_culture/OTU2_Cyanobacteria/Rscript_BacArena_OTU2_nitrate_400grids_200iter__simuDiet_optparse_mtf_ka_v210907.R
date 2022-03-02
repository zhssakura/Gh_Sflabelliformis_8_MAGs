################################################################# katana ####################################################################################################
############ 
library(lattice)
library(ReacTran)
library(rootSolve)
library(deSolve)
library(shape)
library(sybil)
library(Matrix)
library(BacArena)
library(parallel)
library(cplexAPI)
SYBIL_SETTINGS("SOLVER", "cplexAPI")
######################## argument ########################

option_list = list(
  optparse::make_option(c("-n", "--no3_con"),   type="double",       default=0.5,       help="NH3 concentration, defalt: 0.5"),
  optparse::make_option(c("-a", "--arena_mn"),  type="double",       default=20,        help="integer indicating the length of an arena, defalt: 20"),
  optparse::make_option(c("-d", "--death_r"),   type="double",       default=0,         help="A percentage of biomass reduce due to the nutrient limitation, defalt: 0"),
  optparse::make_option(c("-i", "--inocc_no"),  type="double",       default=0.05,      help="inocculum, defalt: 0.05"),
  optparse::make_option(c("-c", "--cl_no"),     type="double",       default=4,         help="Number of replicates, defalt: 4"),
  optparse::make_option(c("-t", "--iter"),      type="double",       default=6,         help="Number of iteration, defalt:6"));

opt_parser = optparse::OptionParser(option_list=option_list);
opt = optparse::parse_args(opt_parser);

NO3_con =  opt$no3_con
grid_no =  opt$arena_mn
death_r =  opt$death_r
Inocc_no = opt$inocc_no
Cl_no =    opt$cl_no
iter_no =  opt$iter

#####################################################################################################################################################################
############ Katana ############ 
# diet <- readRDS('/srv/scratch/z5095298/software/gapseq/Diets_Sponges_July_2019/Diet_Sponges_July_2019_dl20190807_Shan_v20210425/Diet_Sponges_July_2019_dl20190807_Shan_v20210425_brief_v2_IDcheck_mineral_sw_210514/Diet_for_BacArena/mineral_sw_210514_OTU08_3_BacArena_realistic.RDS')
# diet <- readRDS('/srv/scratch/z5095298/software/gapseq/Diets_Sponges_July_2019/Diet_Sponges_July_2019_dl20190807_Shan_v20210425/Diet_Sponges_July_2019_dl20190807_Shan_v20210425_brief_v2_IDcheck_mineral_sw_210514/Diet_for_BacArena/mineral_sw_210514_OTU08_3_BacArena_realistic_full.RDS')
# diet <- readRDS('/srv/scratch/z5095298/software/gapseq/Diets_Sponges_July_2019/Diet_Sponges_July_2019_dl20190807_Shan_v20210425/Diet_Sponges_July_2019_dl20190807_Shan_v20210425_brief_v2_IDcheck_mineral_sw_210514/Diet_for_BacArena/mineral_sw_210514_OTU08_3_BacArena_realistic_simplified.RDS')
diet <- readRDS('/srv/scratch/z5095298/software/gapseq/Diets_Sponges_July_2019/Diet_Sponges_July_2019_dl20190807_Shan_v20210425/Diet_Sponges_July_2019_dl20190807_Shan_v20210425_brief_v2_IDcheck_mineral_sw_210514/Diet_for_BacArena/mineral_sw_210514_OTU08_3_BacArena_realistic_simplified_9999_NH3only.RDS')
#####################################################################################################################################################################
getwd <- getwd()
getwd

# create a new folder:
new_folder = paste('OTU2_NO3_',NO3_con,'_iter_',iter_no,'/',sep = '')
dir.create(paste(getwd,'/',new_folder,sep = ''))

dir_my_model = '/srv/scratch/z5095298/sponge_modeling/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/' # change this only.

# Synechococcus use nitrate -> NH3
STYtaxon_2='OTU02_renamed/20210505_diet_sw_no_oxygen_STY_Merged_OTU02-draft__12_NO3/STY_Merged_OTU02.RDS' ## 20210505_diet_sw_no_oxygen_STY_Merged_OTU02-draft__12_NO3 is used in single species test and 3-taxa co-culture.
## 20210509_diet_sw_no_oxygen_STY_Merged_OTU02-draft__12_NH3_n_NO3 is used in the 8-taxa co-culture!

model2 <- readRDS(paste(dir_my_model,STYtaxon_2,sep = "")) #d__Bacteria;p__Cyanobacteriota;c__Cyanobacteriia;o__Synechococcales_A;f__Cyanobiaceae;g__Synechococcus_C;s__GCA_001628295.1

replicates <- Cl_no
cores <- Cl_no
cl <- makeCluster(cores, type="PSOCK")
clusterExport(cl, c("diet","model2","getwd","new_folder","NO3_con","grid_no","death_r","Inocc_no","iter_no"))
clusterEvalQ(cl, sink(paste0(getwd,'/',new_folder, Sys.getpid(), ".txt")))

simlist <- parLapply(cl, 1:replicates, function(i){
  print("====================================================================")
  print("============================= OTU2 photosymbiont @2021-09-07 =======================================")
  print("====================================================================")
  bacterium2 <- BacArena::Bac(model2,lyse = F,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=FALSE,limit_growth = F)
  arena <- BacArena::Arena(n=grid_no,m=grid_no)
  inocc <- arena@n * arena@m * Inocc_no # 5% inocculum
  arena <- BacArena::addOrg(object = arena, specI = bacterium2, amount = inocc*1) # amount can be added by fraction.
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange, smax = diet$Input_mM, unit = "mM", add = T)
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[35], smax = 0, unit = "mM", add = F) #[1] EX_cpd00013_e0	N as ammonia
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[36], smax = 0, unit = "mM", add = F) #[1] EX_cpd00528_e0	N as Dinitrogen N2
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[37], smax = NO3_con, unit = "mM", add = F) #[1] EX_cpd00209_e0	N as nitrate
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[38], smax = 0, unit = "mM", add = F) #[1] EX_cpd00075_e0	N as nitrite
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[40], smax = 0, unit = "mM", add = F) #[1] EX_cpd00007_e0	Oxygen
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[42], smax = 1000, unit = "mM", add = F) #[1] EX_cpd11632_e0	Photon hn
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[44], smax = 1, unit = "mM", add = F)    #[1] EX_cpd00067_e0	Proton (H+)
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[46], smax = 0, unit = "mM", add = F) #[1] EX_cpd00074_e0	S as Element sulfur S
  # arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[47], smax = 0, unit = "mM", add = F) #[1] EX_cpd00048_e0	S as Sulfate
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[48], smax = 0, unit = "mM", add = F) #[1] EX_cpd00239_e0	S as Sulfide (H2S)
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[49], smax = 0, unit = "mM", add = F) #[1] EX_cpd00081_e0	S as Sulfite
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[50], smax = 0, unit = "mM", add = F) #[1] EX_cpd00210_e0	S as Taurine
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[51], smax = 0, unit = "mM", add = F) #[1] EX_cpd00268_e0	S as Thiosulfate
  arena@tstep <- 1
  simulation <- BacArena::simEnv(object = arena, time = iter_no, sec_obj = "mtf",continue = T) #minimize total flux
  # simulation <- BacArena::simEnv(object = arena, time = iter_no, continue = T) #FBA; iter 1
  
  # for (i in 2:(iter_no-1)) {
  #   simulation <- BacArena::simEnv(object = simulation, time = 1, continue = T) #FBA
  #   saveRDS(simulation, file = paste(getwd,'/',new_folder,'/BacArena_STY_v20210515_OTU2_mineral_sw_rl_',i,'.RDS',sep = ""))
  # }
  # simulation <- BacArena::simEnv(object = simulation, time = 1, continue = T) #FBA
})
stopCluster(cl)

saveRDS(simlist, file = paste(getwd,'/',new_folder,'/BacArena_STY_v20210515_OTU2_400grids_mineral_sw_rl_Nitrate.RDS',sep = ''))

# ############################ Download RDS file and run on MacOS ############################
getwd <- getwd()
getwd
simulation_loop <- readRDS(paste(getwd,'/',new_folder,'/BacArena_STY_v20210515_OTU2_400grids_mineral_sw_rl_Nitrate.RDS',sep = ''))

# # 00. plotGrowthCurve:
# # Plot growth curve for several simulation_loops.The function plotGrowthCurve takes a list of simulation_loops and plots the time course of species with standard deviation.
pdf(paste(getwd,'/',new_folder,'/001_plotGrowthCurve_biomass_STY_v20210515_OTU2_mineral_sw_rl.pdf',sep = ''),width = 5, height = 4)
plotGrowthCurve(simulation_loop,use_biomass = T)
dev.off()
use_biomass <- plotGrowthCurve(simulation_loop,use_biomass = T,ret_data = T)
write.csv(use_biomass, file = paste(getwd,'/',new_folder,'/001_plotGrowthCurve_biomass_STY_v20210515_OTU2_mineral_sw_rl.csv',sep = ''))
#
pdf(paste(getwd,'/',new_folder,'/002_plotGrowthCurve_number_STY_v20210515_OTU2_mineral_sw_rl.pdf',sep = ''),width = 5, height = 4)
plotGrowthCurve(simulation_loop,use_biomass = F)
dev.off()
Nouse_biomass <- plotGrowthCurve(simulation_loop,use_biomass = F,ret_data = T)
write.csv(Nouse_biomass, file = paste(getwd,'/',new_folder,'/002_plotGrowthCurve_number_STY_v20210515_OTU2_mineral_sw_rl.csv',sep = ''))
