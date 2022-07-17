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
time_step = opt$tstep
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
# dir.create(paste(getwd,'/',new_folder,sep = ''))

# auto_num = 2
auto_num_pre = as.character(as.numeric(auto_num) - 1)
simulation_loop <- readRDS(paste(getwd,'/',new_folder,'/BacArena_',otu_name,'_400grids_mineral_sw_',auto_num_pre,'.RDS',sep = ''))


replicates <- Cl_no
cores <- Cl_no
cl <- makeCluster(cores, type="PSOCK")
clusterExport(cl, c("simulation_loop","diet",
                    "getwd","new_folder","grid_no","death_r","Inocc_no","iter_no","replicates",
                    "cores","setAllExInf_value","rm_rate","time_step","NutConPercent","auto_num","NH3_con","HT_con","difspeed","photon_con"))
clusterEvalQ(cl, sink(paste0(getwd,'/',new_folder, Sys.time(), ".txt")))

simlist <- parLapply(cl, 1:replicates, function(i){
  sybil::SYBIL_SETTINGS("SOLVER", "cplexAPI")
  
  print("====================================================================================================")
  print("============================= Best_models_STY_8taxa_20220615: auto step2 =======================================")
  print(paste("=============================", Sys.time(), '=======================================', sep = ' '))
  print(paste("==Synechococcus use nitrate -> NH3. This was used in 3-taxa co-culture ===================================", sep = ' '))
  print(paste("==AOA - NO forming deactivate!! This was used in 3-taxa co-culture =======================================", sep = ' '))
  print(paste("==OTU1 & 4: SOB and NOB =======================================", sep = ' '))
  print(paste("============================= auto_num=", auto_num, '=======================================', sep = ' '))
  print("====================================================================================================")
  
  arena2 <- BacArena::getArena(simulation_loop[[1]], 1) # Add the arena in the simulist of 1st iteration.
  arena2@orgdat <- arena2@orgdat[-sample(nrow(arena2@orgdat), round(nrow(arena2@orgdat) * rm_rate/100)), ] #if rm_rate = 10, meaning to remove ~10% of individuals randomly
  
  ############# ############# User-define the following lines using your models on Katana ############# ############# 
  arena2 <- BacArena::addSubs(object = arena2, mediac = diet$Exchange[1:20], difunc = "pde",
                               pde = "Diff2d", difspeed = difspeed, smax = diet$Input_mM[1:20]*NutConPercent/100, unit = "mM", add = F) #Replenish all nutrients by "replacing" (add=F meaning replacing while add=T meaning summing up) after part of population removed.
  arena2 <- BacArena::addSubs(object = arena2, mediac = diet$Exchange[21:29], difunc = "pde",
                               pde = "Diff2d", difspeed = difspeed, smax = diet$Input_mM[21:29]*NutConPercent/100, unit = "mM", add = T) #Replenish all nutrients by "summing up" (add=F meaning replacing while add=T meaning summing up) after part of population removed.
  arena2 <- BacArena::addSubs(object = arena2, mediac = diet$Exchange[22], difunc = "pde",
                               pde = "Diff2d", difspeed = difspeed, smax = NH3_con*NutConPercent/100, unit = "mM", add = T) #[1] Sum up EX_cpd00013_e0	N as ammonia
  arena2 <- BacArena::addSubs(object = arena2, mediac = diet$Exchange[17], difunc = "pde",
                               pde = "Diff2d", difspeed = difspeed, smax = HT_con*NutConPercent/100, unit = "mM", add = F) #[1] replace EX_cpd00406_e0	S as HT
  arena2 <- BacArena::addSubs(object = arena2, mediac = diet$Exchange[15], difunc = "pde",
                               pde = "Diff2d", difspeed = difspeed, smax = photon_con*NutConPercent/100, unit = "mM", add = F) #[1] replace EX_cpd11632_e0	photon (hn)
  
  arena2@tstep <- time_step
  
  BacArena::chemotaxis(arena2@specs[["STY_Merged_OTU06"]], arena2, 1, chemo='EX_cpd00075_e0', arena2@occupyM)#N as nitrite
  ############# ############# ############# ############# ############# ############# ############# ############# 
  
  simulation <- BacArena::simEnv(object = arena2, time = 1, sec_obj = "mtf", continue = T) #pFBA for 1 iter
  # simulation <- BacArena::simEnv(object = arena2, time = 1, continue = T) #FBA for 1 iter
  
  
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
