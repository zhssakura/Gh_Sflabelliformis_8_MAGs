############ 
library(lattice)
library(ReacTran)
library(rootSolve)
library(deSolve)
library(shape)
library(sybil)
library(Matrix)
library(BacArena)

############ Read RDS files:############ 
# MacOS
# Assigning Directory as a Variable in R
getwd <- getwd()
setwd(getwd)

SCM1_model_Li_2018 = '/Input_files/SCM1_model_Li_2018.RDS'
model_SCM1_Li2018 <- readRDS(paste(getwd,SCM1_model_Li_2018,sep = ""))

# Assign name for output file:
outfile_043 = paste(paste(getwd,'/01_Get_ReactionPool/043_df_model_SCM1_Li2018_list_react.tsv',sep = ""))
outfile_045 = paste(paste(getwd,'/01_Get_ReactionPool/045_df_model_SCM1_Li2018_list_metabolites.tsv',sep = ""))
outfile_046 = paste(paste(getwd,'/01_Get_ReactionPool/046_model_SCM1_Li2018_Smat.csv',sep = ""))
############ ############ ############ 
# 043.What are those reactions in each model after remove/add reactions by SEED rxn IDs?
# Check number of reactions of ******model_SCM1_Li2018****: 
model_SCM1_Li2018@react_num
# [1] 765

model_SCM1_Li2018@react_attr
# This is in a strange format.

# Extract all reactions from the original model. (Here is SCM1_model_Li_2018.xml)
model_SCM1_Li2018_list <- list(model_SCM1_Li2018@react_id, model_SCM1_Li2018@react_name, model_SCM1_Li2018@react_rev, model_SCM1_Li2018@react_single, model_SCM1_Li2018@react_de,
                                           model_SCM1_Li2018@lowbnd, model_SCM1_Li2018@uppbnd, model_SCM1_Li2018@obj_coef, model_SCM1_Li2018@gprRules, model_SCM1_Li2018@gpr
                                           ,model_SCM1_Li2018@react_attr
                               )
df_model_SCM1_Li2018_list  <- as.data.frame(model_SCM1_Li2018_list, row.names = NULL)
colnames(df_model_SCM1_Li2018_list) <- c('react_id','react_name','react_rev','react_single','react_de','lowbnd','uppbnd','obj_coef','gprRules','gpr'
                                         ,'react_attr'
                                         )
                                                     # "rxn","name","ec","tc","qseqid","pident","evalue","bitscore","qcovs","stitle","sstart","send","pathway",
                                                     # "status","pathway.status","complex","exception","complex.status","seed","gs.origin")
write.table(df_model_SCM1_Li2018_list, file = outfile_043, quote = F)
# The output contains many KEGG reaction ids, which can be converted to SEED RXN ids by using seed_reactions.tsv.
### Important ###
# Strange lines were included in the ourput file, since this 'SCM1_model_Li_2018.RDS' was 
# converted from a .xml file.
# Need to manually correct those 'errors' to generate a new file '043_df_model_SCM1_Li2018_list_react_ECcorrect.csv'



############ ############ ############ 

# 045.What are the metabolites in this new model?

# number of reactions of ******model_SCM1_Li2018****: 
model_SCM1_Li2018@met_num
# [1] 825
model_SCM1_Li2018_met_name = model_SCM1_Li2018@met_name
model_SCM1_Li2018_met_id = model_SCM1_Li2018@met_id
# model_SCM1_Li2018: Generate a dataframe with metabolite info for each model
model_SCM1_Li2018_list_met <- list(model_SCM1_Li2018@met_id, model_SCM1_Li2018@met_name, model_SCM1_Li2018@met_comp, model_SCM1_Li2018@met_single, model_SCM1_Li2018@met_de)
df_model_SCM1_Li2018_list_met <- as.data.frame(model_SCM1_Li2018_list_met, row.names = NULL)
colnames(df_model_SCM1_Li2018_list_met) <- c('met_id','met_name','met_comp','met_single','met_de')
# R Sort a Data Frame using Order()
df_model_SCM1_Li2018_list_met_sort <- df_model_SCM1_Li2018_list_met[order(df_model_SCM1_Li2018_list_met$met_id),]
write.table(df_model_SCM1_Li2018_list_met_sort,file = outfile_045, quote = F)

############ ############ ############ 

# 046.What are the stiochiometric matrix (model_xxx@S) in this model?
nrow(model_SCM1_Li2018@S)
# [1] 825

ncol(model_SCM1_Li2018@S)
# [1] 765


model_SCM1_Li2018_Smat <- model_SCM1_Li2018@S
model_SCM1_Li2018_Smat_mx = as.matrix(model_SCM1_Li2018@S)
class(model_SCM1_Li2018_Smat_mx)
# [1] "matrix"

colnames(model_SCM1_Li2018_Smat_mx)
# NULL
model_SCM1_Li2018@react_num
# [1] 765
model_SCM1_Li2018@react_id
colnames(model_SCM1_Li2018_Smat_mx)<-cbind(model_SCM1_Li2018@react_id)

row.names(model_SCM1_Li2018_Smat_mx)
# NULL
model_SCM1_Li2018@met_num
# 825
model_SCM1_Li2018@met_id
row.names(model_SCM1_Li2018_Smat_mx)<-cbind(model_SCM1_Li2018@met_id)

write.csv(model_SCM1_Li2018_Smat_mx, file = outfile_046, quote = F)

# In this table, you can check the flux (values) of a compound (e.g. cpd16503[e0]) of this model based on the reaction ID (e.g. EX_cpd16503_e0).
# ############ ############ ############ 
# ####### Updates on 23/Apr/2020: example starts ####### 
# # 048. Johannes: this is a function which returns all active reactions which use or produce a certain compound. can start with this and extend it according to the needs.
# get_rxn <- function(mod, met){
#   sol <- sybil::optimizeProb(mod, retOptSol=F)
#   smatrix <- mod@S[,which(sol$fluxes!=0)]
#   
#   met.idx <- grep(met, mod@met_id)
#   
#   print(paste("Active reactions that uses the metabolite", met, mod@met_name[met.idx]))
#   
#   rxn.idx <- which(smatrix[met.idx,]!=0)
#   printReaction(mod, react=rxn.idx)
# }
# 
# aa <- get_rxn(model_SCM1_Li2018, 'cpd00001') # example: all reactions that uses 'cpd00001' will be gave. 
# aa <- get_rxn(model_SCM1_Li2018, 'cpd00002') # cpd00002: ATP
# aa <- get_rxn(model_SCM1_Li2018, 'cpd00004') # cpd00004: NADH
# 
# # [1] "Active reactions that uses the metabolite cpd00001 H2O" "Active reactions that uses the metabolite cpd00001 H2O"
# # abbreviation	equation
# # R00289	(1) cpd00062[c] + (1) cpd00089[c] --> (1) cpd00012[c] + (1) cpd00026[c]
# # R00374	(1) cpd00033[c] + (2) cpd11609[c] --> (1) cpd00011[c] + (1) cpd00150[c] + (2) cpd11610[c] + (1) cpd00067[c]
# # R01267	(1) cpd00132[c] <==> (1) cpd00001[c] + (1) cpd01651[c]
# # R01731	(1) cpd00032[c] + (1) cpd00616[c] <==> (1) cpd00041[c] + (1) cpd00219[c]
# # R05555	(1) cpd00350[c] + (1) cpd00113[c] <==> (1) cpd00012[c] + (1) cpd08211[c]
# # R07399	(1) cpd00001[c] + (1) cpd00020[c] + (1) cpd00022[c] <==> (1) cpd00010[c] + (1) cpd01700[c]
# # R00291	(1) cpd00026[c] <==> (1) cpd00043[c]
# ####### Updates on 23/Apr/2020: example ends ####### 
# # Check if there are active reaction uses the following nitrogen resources:
# nitrate <- get_rxn(model_SCM1_Li2018, 'cpd00209') #
# # None
# 
# 
# nitrogen <- get_rxn(model_SCM1_Li2018, 'cpd00528') #N2-e0
# #None
# 
# Nitrite <- get_rxn(model_SCM1_Li2018, 'cpd00075') #cpd00075: Nitrite
# # [1] "Active reactions that uses the metabolite cpd00075 Nitrite" "Active reactions that uses the metabolite cpd00075 Nitrite"
# # abbreviation	equation
# # R05046	(1) cpd00001[c] + (1) cpd03518[c] --> (1) cpd00047[c] + (1) cpd03519[c]
# # R05069	(1) cpd00498[c] <==> (1) cpd10162[c]
# # R08716	(1) cpd17535[c] --> (1) cpd08369[c]
# # R05812	(1) cpd00067[c] + (1) cpd00004[c] + (1) cpd08372[c] --> (1) cpd00003[c] + (1) cpd08373[c]
# # R05721	(1) cpd00007[c] + (1) cpd00099[c] <==> (1) cpd01054[c]
# 
# Nitrite_c0 <- get_rxn(model_SCM1_Li2018, 'cpd00075_e') #cpd00075_e: Nitrite-c0
# # [1] "Active reactions that uses the metabolite cpd00075_e Nitrite"
# # abbreviation	equation
# # R03632	(1) cpd00002[c] + (1) cpd11900[c] --> (1) cpd00008[c] + (1) cpd12167[c]
# # R03646	(1) cpd00002[c] + (1) cpd00051[c] + (1) cpd11907[c] --> (1) cpd00018[c] + (1) cpd00012[c] + (1) cpd12036[c]
# 
# 
# NH3 <- get_rxn(model_SCM1_Li2018, 'cpd00013') #NH3-c0
# # [1] "Active reactions that uses the metabolite cpd00013 NH3-c0" "Active reactions that uses the metabolite cpd00013 NH3-e0"
# # abbreviation	equation
# # rxn00048_c0	(1) cpd00067[c0] + (2) cpd02656[c0] --> (1) cpd00220[c0] + (1) cpd02882[c0]
# # rxn00086_c0	(2) cpd00042[c0] + (1) cpd00006[c0] <==> (1) cpd00067[c0] + (1) cpd00111[c0] + (1) cpd00005[c0]
# # ... ...
# 
# 
# 
# 
# ############ ############ ############ 
# # The following attribute shows the way reactions been added to the model. by Silvio on 23/Apr/2020 via Skype meeting.
# model_SCM1_Li2018@react_attr$gs.origin
# 
# ############ ############ ############ 
# # 050.What are the active reactions which use or produce a certain compound in ******model_SCM1_Li2018****: 
# # Updates @8-1-21: this will give a new column called 'sol$fluxes' with the same length to the 'reaction id', which shows fluxes of each reaction.
# mod = model_SCM1_Li2018
# dim(mod)
# # [1] 825 765
# 
# sol <- sybil::optimizeProb(mod, retOptSol=F) # find a list (of 7,sol), which include flux info of 765 reactions for this model.
# model_SCM1_Li2018_flux_list <- as.list(sol$fluxes)
# length(model_SCM1_Li2018_flux_list)
# # [1] 765
# 
# model_SCM1_Li2018_flux_list <- as.data.frame(sol$fluxes)
# dim(model_SCM1_Li2018_flux_list)
# # [1] 765   1
# 
# 
# model_SCM1_Li2018_Smat <- model_SCM1_Li2018@S
# 
# model_SCM1_Li2018_Smat <- read.csv(paste(getwd,'SCM1_model_Li_2018/046_model_SCM1_Li2018_Smat.csv',sep = ""),row.names = 'X')
# t_model_SCM1_Li2018_Smat <- as.data.frame(t(model_SCM1_Li2018_Smat))
# dim(t_model_SCM1_Li2018_Smat)
# # [1] 765 825
# 
# cbind_t_model_SCM1_Li2018_Smat<- as.data.frame(cbind(model_SCM1_Li2018_flux_list, t_model_SCM1_Li2018_Smat))
# dim(cbind_t_model_SCM1_Li2018_Smat)
# # [1] 765 825
# 
# write.csv(cbind_t_model_SCM1_Li2018_Smat, file = paste(getwd,'SCM1_model_Li_2018/050_model_SCM1_Li2018_Smat_fluxes.csv',sep = ""))

