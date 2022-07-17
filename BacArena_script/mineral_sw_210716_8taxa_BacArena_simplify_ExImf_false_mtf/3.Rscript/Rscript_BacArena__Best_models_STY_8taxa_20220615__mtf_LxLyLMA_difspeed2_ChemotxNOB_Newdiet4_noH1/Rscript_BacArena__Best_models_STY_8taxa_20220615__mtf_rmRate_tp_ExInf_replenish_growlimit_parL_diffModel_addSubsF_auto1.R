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
  optparse::make_option(c("-o", "--otu_name"),                 type="character",      default='8taxa',       help="OTU_name, defalt: 8taxa"),
  optparse::make_option(c("-D", "--infile_BacarenaRDS_diet"),  type="character",      default="/srv/scratch/z5095298/software/gapseq/Diets_Sponges_July_2019/Diet_Sponges_July_2019_dl20190807_Shan_v20210425/Diet_Sponges_July_2019_dl20190807_Shan_v20210425_brief_v2_IDcheck_mineral_sw_210514/Diet_for_BacArena/mineral_sw_210514_OTU08_3_BacArena_realistic_simplified_9999_NH3only.RDS",
                        help="Diet RDS file for bacarena with full path, defalt: /srv/scratch/z5095298/software/gapseq/Diets_Sponges_July_2019/Diet_Sponges_July_2019_dl20190807_Shan_v20210425/Diet_Sponges_July_2019_dl20190807_Shan_v20210425_brief_v2_IDcheck_mineral_sw_210514/Diet_for_BacArena/mineral_sw_210514_OTU08_3_BacArena_realistic_simplified_9999_NH3only.RDS"),
  optparse::make_option(c("-k", "--keywd"),                    type="character",      default="sw_v20210716",    help="(Nutirent/date) keyword used in the folder name of output files, defalt: sw_v20210716"),
  optparse::make_option(c("-r", "--rm_rate"),   type="double",       default=0,        help="removeM value, defalt: 0"),
  optparse::make_option(c("-p", "--tstep"),     type="double",       default=1,         help="Time step, defalt: 1h per iteration"),
  
  optparse::make_option(c("-a", "--arena_mn"),              type="double",       default=20,        help="integer indicating the length of an arena, defalt: 20"),
  optparse::make_option(c("-d", "--death_r"),               type="double",       default=0,         help="A percentage of biomass reduce due to the nutrient limitation, defalt: 0"),
  optparse::make_option(c("-i", "--inocc_no"),              type="double",       default=1,         help="inocculum, defalt: 1"),
  optparse::make_option(c("-c", "--cl_no"),                 type="double",       default=1,         help="Number of replicates, defalt: 1"),
  optparse::make_option(c("-s", "--setAllExInf_value"),     type="character",    default='FALSE',   help="setAllExInf, defalt: FALSE"),
  optparse::make_option(c("-g", "--NutConPercent"),         type="double",       default=100,       help="percentage of nutrients, defalt: 100 (unit: %)"),
  optparse::make_option(c("-n", "--auto_num"),              type="double",       default=1,         help="surfix in the simulation RDS file name."),
  optparse::make_option(c("-N", "--NH3_con"),               type="double",       default=0.5,       help="NH3 concentration, defalt: 0.5"),
  optparse::make_option(c("-H", "--HT_con"),                type="double",       default=0.5,       help="HT concentration, defalt: 0.5"),
  optparse::make_option(c("-S", "--difspeed"),              type="double",       default=1.7e-09,   help="A number indicating the diffusion speed (given by number of cells per iteration), defalt: 1.7e-09 cm2 h-1"),
  optparse::make_option(c("-P", "--photon_con"),            type="double",       default=1000,      help="A number indicating the photon concentration, defalt: 1000 mM"),
  optparse::make_option(c("-t", "--iter"),                  type="double",       default=5,         help="Number of iteration, defalt:5"));

opt_parser = optparse::OptionParser(option_list=option_list, add_help_option=FALSE);
opt = optparse::parse_args(opt_parser);

otu_name        = opt$otu_name
infile_diet     = opt$infile_BacarenaRDS_diet
keywd           = opt$keywd
grid_no   = opt$arena_mn
death_r   = opt$death_r
Inocc_no  = opt$inocc_no
Cl_no     = opt$cl_no
setAllExInf_value = opt$setAllExInf_value
NutConPercent   = opt$NutConPercent
iter_no   = opt$iter
rm_rate   = opt$rm_rate
time_step =  opt$tstep
auto_num  = opt$auto_num
NH3_con =  opt$NH3_con
HT_con =  opt$HT_con
difspeed = opt$difspeed
photon_con = opt$photon_con

#####################################################################################################################################################################

############ Katana ############ 
# diet <- readRDS('/srv/scratch/z5095298/software/gapseq/Diets_Sponges_July_2019/Diet_Sponges_July_2019_dl20190807_Shan_v20210425/Diet_Sponges_July_2019_dl20190807_Shan_v20210425_brief_v2_IDcheck_mineral_sw_210514/Diet_for_BacArena/mineral_sw_210514_OTU08_3_BacArena_realistic.RDS')
# diet <- readRDS('/srv/scratch/z5095298/software/gapseq/Diets_Sponges_July_2019/Diet_Sponges_July_2019_dl20190807_Shan_v20210425/Diet_Sponges_July_2019_dl20190807_Shan_v20210425_brief_v2_IDcheck_mineral_sw_210514/Diet_for_BacArena/mineral_sw_210514_OTU08_3_BacArena_realistic_full.RDS')
# diet <- readRDS('/srv/scratch/z5095298/software/gapseq/Diets_Sponges_July_2019/Diet_Sponges_July_2019_dl20190807_Shan_v20210425/Diet_Sponges_July_2019_dl20190807_Shan_v20210425_brief_v2_IDcheck_mineral_sw_210514/Diet_for_BacArena/mineral_sw_210514_OTU08_3_BacArena_realistic_simplified.RDS')
diet <- readRDS(infile_diet)
#####################################################################################################################################################################
getwd <- getwd()
getwd

# create a new folder:
new_folder = paste(otu_name,'_',keywd,'_death_',death_r,'_rmRate_',rm_rate,'_NutPer_',NutConPercent,'_Inoc_',Inocc_no,'_NH3',NH3_con,'_HT',HT_con, '_tstep_',time_step,'h_difspeed_',difspeed,'/',sep = '')
dir.create(paste(getwd,'/',new_folder,sep = ''))
############# User-define the following lines using your models on Katana ############# 
dir_gs_models='/srv/scratch/z5095298/sponge_modeling/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/'

# OTU1 v20210902:
STYtaxon_1='/OTU01_renamed/2021-07-01_ReactionPool_to_ObjectModel_use_bio1/07_adapt_addReact_R/04_Rscp/OTU01_MOM_20210825/20210831__OTU01_MOM_20210825__TrainDiet_O2_bothN_no0_sulfite_NoHT/STY_Merged_OTU01_MOM.RDS'
# OTU2 v20210902:
# STYtaxon_2='/OTU02_renamed/20210509_diet_sw_no_oxygen_STY_Merged_OTU02-draft__12_NH3_n_NO3/STY_Merged_OTU02.RDS'
# OTU2 v20210909 (Synechococcus use nitrate -> NH3. This was used in 3-taxa co-culture)
STYtaxon_2='OTU02_renamed/20210505_diet_sw_no_oxygen_STY_Merged_OTU02-draft__12_NO3/STY_Merged_OTU02.RDS'
# OTU3 v20210902:
STYtaxon_3='/OTU03_renamed/2021-07-01_ReactionPool_to_ObjectModel_use_bio1/07_adapt_addReact_R/04_Rscp/OTU03_MOM_20210802hypotau_sulfitetransport_nitritetransport_adh/20210802_diet_sw_STY_Merged_OTU03-MOM_20210802hypotau_sulfitetransport_nitritetransport_adh__12_hypotaurine_nitrite_NoO2_noVB_noNO3/STY_Merged_OTU03_MOM.RDS'
# OTU4 v20210902:
STYtaxon_4='/OTU04_renamed/2021-05-11_ReactionPool_to_ObjectModel/07_adapt_addReact_R/04_Rscp/OTU04_MOM_20210825/20210825__OTU04_MOM_20210825__TrainDiet_O2_bothN_no0_sulfite_NoHT/STY_Merged_OTU04_MOM.RDS'
# OTU5 v20210902:
STYtaxon_5='/OTU05_renamed/2021-07-01_ReactionPool_to_ObjectModel_use_bio1/07_adapt_addReact_R/04_Rscp/OTU05_MOM_20210813_sulftTrans/20210815_diet_sw_STY_Merged_OTU05-MOM_20210813_sulftTrans__12_SOB_NoHT_O2_bothN_Sulfite_no0/STY_Merged_OTU05_MOM.RDS'
# OTU6 v20210902:
STYtaxon_6='/OTU06_renamed/20210515_STY_Merged_OTU06_mineral_sw_3_NoH2S/STY_Merged_OTU06.RDS'
# OTU7 v20210902:
STYtaxon_7='/OTU07_renamed/2021-07-01_ReactionPool_to_ObjectModel_use_bio1/07_adapt_addReact_R/04_Rscp/OTU07_MOM_20210813/20210815_diet_sw_STY_Merged_OTU07-MOM_20210813__12_SOB_noHT_O2_bothN_sulfite_no0/STY_Merged_OTU07_MOM.RDS'
# OTU8 v20210902 (NO gas forming?):
# STYtaxon_8='/OTU08_renamed_arc/20210515_STY_Merged_OTU08_mineral_sw_3_unlimited_CO2/STY_Merged_OTU08.RDS'
# OTU8 v20210909 (AOA - NO forming deactivate!! This was used in 3-taxa co-culture):
STYtaxon_8='OTU08_renamed_arc/Step_protocals_using_R/07_gf_model/20210515_STY_Merged_OTU08_mineral_sw_3_unlimited_CO2/STY_Merged_OTU08_MgfM/STY_Merged_OTU08_MgfM.RDS'

model1 <- readRDS(paste(dir_gs_models,STYtaxon_1,sep = "")) #d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__UBA10353;f__LS-SOB;g__;s__
model2 <- readRDS(paste(dir_gs_models,STYtaxon_2,sep = "")) #d__Bacteria;p__Cyanobacteriota;c__Cyanobacteriia;o__Synechococcales_A;f__Cyanobiaceae;g__Synechococcus_C;s__GCA_001628295.1
model3 <- readRDS(paste(dir_gs_models,STYtaxon_3,sep = "")) #d__Bacteria;p__Myxococcota;c__UBA9160;o__UBA9160;f__UBA6930;g__;s__
model4 <- readRDS(paste(dir_gs_models,STYtaxon_4,sep = "")) #d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__UBA10353;f__LS-SOB;g__;s__
model5 <- readRDS(paste(dir_gs_models,STYtaxon_5,sep = "")) #d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__UBA6729;f__;g__;s__
model6 <- readRDS(paste(dir_gs_models,STYtaxon_6,sep = "")) #d__Bacteria;p__Nitrospirota;c__Nitrospiria;o__Nitrospirales;f__UBA8639;g__bin75;s__
model7 <- readRDS(paste(dir_gs_models,STYtaxon_7,sep = "")) #d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__;f__;g__;s__
model8 <- readRDS(paste(dir_gs_models,STYtaxon_8,sep = "")) #d__Archaea;p__Crenarchaeota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrosopumilaceae;g__Cenarchaeum;s__
############# ############# ############# ############# ############# ############# ############# ############# 


replicates <- Cl_no
cores <- Cl_no
cl <- makeCluster(cores, type="PSOCK")
clusterExport(cl, c("diet","model1","model2","model3","model4","model5","model6","model7","model8",
                    "getwd","new_folder","grid_no","death_r","Inocc_no","iter_no","replicates",
                    "cores","setAllExInf_value","rm_rate","time_step","NutConPercent","auto_num","NH3_con","HT_con","difspeed","photon_con"))
clusterEvalQ(cl, sink(paste0(getwd,'/',new_folder, Sys.time(), ".txt")))

simlist <- parLapply(cl, 1:replicates, function(i){
  sybil::SYBIL_SETTINGS("SOLVER", "cplexAPI") #https://github.com/euba/BacArena/issues/152
  
  print("====================================================================================================")
  print("============================= Best_models_STY_8taxa_20220615: auto step1 =======================================")
  print(paste("=============================", Sys.time(), '=======================================', sep = ' '))
  print(paste("==Synechococcus use nitrate -> NH3. This was used in 3-taxa co-culture ===================================", sep = ' '))
  print(paste("==AOA - NO forming deactivate!! This was used in 3-taxa co-culture =======================================", sep = ' '))
  print(paste("==OTU1 & 4: SOB and NOB =======================================", sep = ' '))
  print(paste("============================= auto_num=", auto_num, '=======================================', sep = ' '))
  print("====================================================================================================")
  # bacterium_SOB <- BacArena::Bac(my_model,lyse = F,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
  arena <- BacArena::Arena(n=grid_no,m=grid_no,Lx=0.012,Ly=0.012) # This grid length represents there are 100 cell distance between two microbial cells.
  ############# ############# User-define the following lines using your models on Katana ############# ############# 
  bacterium1        <- BacArena::Bac(model1,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
  bacterium2        <- BacArena::Bac(model2,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
  bacterium3        <- BacArena::Bac(model3,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
  bacterium4        <- BacArena::Bac(model4,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
  bacterium5        <- BacArena::Bac(model5,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
  bacterium6        <- BacArena::Bac(model6,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
  bacterium7        <- BacArena::Bac(model7,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
  bacterium8        <- BacArena::Bac(model8,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
  
  inocc <- arena@n * arena@m * Inocc_no # Inocc_no = 1
  # arena <- BacArena::addOrg(object = arena, specI = bacterium_SOB, amount = inocc*1) # amount can be added by fraction.
  arena <- BacArena::addOrg(object = arena, specI = bacterium1, amount = inocc*0.18220) # amount can be added by fraction.
  arena <- BacArena::addOrg(object = arena, specI = bacterium2, amount = inocc*0.02740)
  arena <- BacArena::addOrg(object = arena, specI = bacterium3, amount = inocc*0.11580)
  arena <- BacArena::addOrg(object = arena, specI = bacterium4, amount = inocc*0.20940)
  arena <- BacArena::addOrg(object = arena, specI = bacterium5, amount = inocc*0.05150)
  arena <- BacArena::addOrg(object = arena, specI = bacterium6, amount = inocc*0.03690)
  arena <- BacArena::addOrg(object = arena, specI = bacterium7, amount = inocc*0.09090)
  arena <- BacArena::addOrg(object = arena, specI = bacterium8, amount = inocc*0.28570)
  
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[1:20], difunc = "pde",
                             pde = "Diff2d", difspeed = difspeed, smax = diet$Input_mM[1:20]*NutConPercent/100, unit = "mM", add = F) #Replenish all nutrients by "replacing" (add=F meaning replacing while add=T meaning summing up) after part of population removed.
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[21:29], difunc = "pde",
                             pde = "Diff2d", difspeed = difspeed, smax = diet$Input_mM[21:29]*NutConPercent/100, unit = "mM", add = F) #Replenish all nutrients by "replacing" (add=F meaning replacing while add=T meaning summing up) after part of population removed.
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[22], difunc = "pde",
                             pde = "Diff2d", difspeed = difspeed, smax = NH3_con*NutConPercent/100, unit = "mM", add = F) #[1] replace EX_cpd00013_e0	N as ammonia
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[17], difunc = "pde",
                             pde = "Diff2d", difspeed = difspeed, smax = HT_con*NutConPercent/100, unit = "mM", add = F) #[1] replace EX_cpd00406_e0	S as HT
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[15], difunc = "pde",
                             pde = "Diff2d", difspeed = difspeed, smax = photon_con*NutConPercent/100, unit = "mM", add = F) #[1] replace EX_cpd11632_e0	photon (hn)
  # The following line should be added in the auto1 script only!:
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[30:60], difunc = "pde",# all the exchanged metabolites found in the 200th cycle @M27_7v1b
                             pde = "Diff2d", difspeed = difspeed, smax = diet$Input_mM[30:60]*NutConPercent/100, unit = "mM", add = F) #Replenish all nutrients by "replacing" (add=F meaning replacing while add=T meaning summing up) after part of population removed.

  arena@tstep <- time_step
  
  BacArena::chemotaxis(arena@specs[["STY_Merged_OTU06"]], arena, 1, chemo='EX_cpd00075_e0', arena@occupyM)#N as nitrite
  ############# ############# ############# ############# ############# ############# ############# ############# 
  
  simulation <- BacArena::simEnv(object = arena, time = 1, sec_obj = "mtf", continue = T) #pFBA; iter 1
  # simulation <- BacArena::simEnv(object = arena, time = 1, continue = T) #FBA; iter 1
  
})
stopCluster(cl)

saveRDS(simlist, file = paste(getwd,'/',new_folder,'/BacArena_',otu_name,'_400grids_mineral_sw_',auto_num,'.RDS',sep = ''))

# ############################ Download RDS file and run on MacOS ############################
getwd <- getwd()
getwd
simulation_loop <- readRDS(paste(getwd,'/',new_folder,'/BacArena_',otu_name,'_400grids_mineral_sw_',auto_num,'.RDS',sep = ''))

# 01. evalArena: 
# Check the spatial and temporal change of populations and/or concentrations.
# pdf(paste(getwd,'/',new_folder,'/011_evalArena_',auto_num,'.pdf',sep = ''),width = 10, height = 8)
# evalArena(simulation_loop[[1]], plot_items = "Population", phencol = F, retdata = T, time = (seq_along(simulation_loop[[1]]@simlist) - 1),show_legend = F)
# dev.off()

pdf(paste(getwd,'/',new_folder,'/012_evalArena_',auto_num,'.pdf',sep = ''),width = 10, height = 8)
evalArena(simulation_loop[[1]], plot_items = "Population", phencol = F, retdata = T, 
          time = (seq_along(simulation_loop[[1]]@simlist) - 1),show_legend = T,legend_pos = 'right')
dev.off()
# To get better plot, run 081: BacArena_STY_5taxa_400grids_spongeDiet_6iterations.R


######################################## exclude biomass ########################################
####time0
type_xy_md_time0 <- evalArena(simulation_loop[[1]], plot_items = "Population", phencol = F, retdata = T, 
                              time = (seq_along(simulation_loop[[1]]@simlist) - 1),show_legend = F)
write.table(type_xy_md_time0$Population$time0, file = paste(getwd,'/',new_folder,'081_ggplot_',auto_num,'.tsv',sep = ''))

dt_type_xy_md_time0 <- as.data.frame(read.table(file = paste(getwd,'/',new_folder,'081_ggplot_',auto_num,'.tsv',sep = '')))
# head(dt_type_xy_md_time0)
# str(dt_type_xy_md_time0)

# library('ggplot2')
# pdf(paste(getwd,'/',new_folder,'081_ggplot_',auto_num,'.pdf',sep = ''),width = 5, height = 4)
# ggplot(dt_type_xy_md_time0, aes(x=dt_type_xy_md_time0$x, y=dt_type_xy_md_time0$y)) + 
#   geom_point(color=dt_type_xy_md_time0$type,size=0.3)+
#   labs(title = "STY_8species_time0") + 
#   coord_fixed(ratio=1) + 
#   theme(plot.background = element_blank())+
#   theme(plot.background = element_blank(), 
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), 
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.x = element_blank(), 
#         axis.text.y = element_blank(),
#         axis.ticks = element_blank(),
#         legend.position=c(1,1), legend.justification=c(1,1))
# dev.off()
