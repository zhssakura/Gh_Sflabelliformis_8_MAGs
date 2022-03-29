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

# For the draft model:
STY_Merged_STY_Merged_OTU01_draft = '/Input_files/STY_Merged_OTU01-draft.RDS'
model_STY_Merged_OTU01_draft <- readRDS(paste(getwd,STY_Merged_STY_Merged_OTU01_draft,sep = ""))
# Assign name for output file:
outfile_043 = paste(paste(getwd,'/03_Get_ObjectModel/043_df_model_STY_Merged_OTU01_draft_list_react.tsv',sep = ""))
outfile_045 = paste(paste(getwd,'/03_Get_ObjectModel/045_df_model_STY_Merged_OTU01_draft_list_metabolites.tsv',sep = ""))
outfile_046 = paste(paste(getwd,'/03_Get_ObjectModel/046_model_STY_Merged_OTU01_draft_Smat.csv',sep = ""))

# For the final model(modified gap-filled models):
dir_nm = '/20210708_diet_sw_no_oxygen_STY_Merged_OTU01_MOM__12_sulfite/'
dir_nm = '/20210708_diet_sw_no_oxygen_STY_Merged_OTU01_MOM__12_taurine/'
model_STY_Merged_OTU01_draft <- readRDS(paste(getwd,'/Input_files/',dir_nm,'/STY_Merged_OTU01_MOM.RDS', sep = ""))
# create a new folder:
dir.create(paste(getwd,'/03_Get_ObjectModel/',dir_nm,sep = ''))
# Assign name for output file:
outfile_043 = paste(paste(getwd,'/03_Get_ObjectModel/',dir_nm,'/043_df_model_STY_Merged_OTU01_final_list_react.tsv',sep = ""))
outfile_045 = paste(paste(getwd,'/03_Get_ObjectModel/',dir_nm,'/045_df_model_STY_Merged_OTU01_final_list_metabolites.tsv',sep = ""))
outfile_046 = paste(paste(getwd,'/03_Get_ObjectModel/',dir_nm,'/046_model_STY_Merged_OTU01_final_Smat.csv',sep = ""))

############ ############ ############ 
# 043.What are those reactions in each model after remove/add reactions by SEED rxn IDs?
# Check number of reactions of ******model_STY_Merged_OTU01_draft****: 
model_STY_Merged_OTU01_draft@react_num
# [1] 765

# Extract all reactions from the original model. (Here is STY_Merged_STY_Merged_OTU01_draft.xml)
model_STY_Merged_OTU01_draft_list <- list(model_STY_Merged_OTU01_draft@react_id, model_STY_Merged_OTU01_draft@react_name, model_STY_Merged_OTU01_draft@react_rev, model_STY_Merged_OTU01_draft@react_single, model_STY_Merged_OTU01_draft@react_de,
                                           model_STY_Merged_OTU01_draft@lowbnd, model_STY_Merged_OTU01_draft@uppbnd, model_STY_Merged_OTU01_draft@obj_coef, model_STY_Merged_OTU01_draft@gprRules, model_STY_Merged_OTU01_draft@gpr,
                                           model_STY_Merged_OTU01_draft@react_attr)
df_model_STY_Merged_OTU01_draft_list  <- as.data.frame(model_STY_Merged_OTU01_draft_list, row.names = NULL)
colnames(df_model_STY_Merged_OTU01_draft_list) <- c('react_id','react_name','react_rev','react_single','react_de','lowbnd','uppbnd','obj_coef','gprRules','gpr'
                              ,'react_attr',"rxn","name","ec","tc","qseqid","pident","evalue","bitscore","qcovs","stitle","sstart"
                              ,"send","pathway","status","pathway.status","complex","exception","complex.status","gs.origin","annotation","MNX_ID","seedID","keggID","biggID","biocycID"
)

write.table(df_model_STY_Merged_OTU01_draft_list, file = outfile_043,quote = F,row.names = F,sep = '\t')

# The output contains many KEGG reaction ids, which can be converted to SEED RXN ids by using seed_reactions.tsv.



############ ############ ############ 

# 045.What are the metabolites in this new model?

# number of reactions of ******model_STY_Merged_OTU01_draft****: 
model_STY_Merged_OTU01_draft@met_num
# [1] 825
model_STY_Merged_OTU01_draft_met_name = model_STY_Merged_OTU01_draft@met_name
model_STY_Merged_OTU01_draft_met_id = model_STY_Merged_OTU01_draft@met_id
# model_STY_Merged_OTU01_draft: Generate a dataframe with metabolite info for each model
model_STY_Merged_OTU01_draft_list_met <- list(model_STY_Merged_OTU01_draft@met_id, model_STY_Merged_OTU01_draft@met_name, model_STY_Merged_OTU01_draft@met_comp, model_STY_Merged_OTU01_draft@met_single, model_STY_Merged_OTU01_draft@met_de)
df_model_STY_Merged_OTU01_draft_list_met <- as.data.frame(model_STY_Merged_OTU01_draft_list_met, row.names = NULL)
colnames(df_model_STY_Merged_OTU01_draft_list_met) <- c('met_id','met_name','met_comp','met_single','met_de')
# R Sort a Data Frame using Order()
df_model_STY_Merged_OTU01_draft_list_met_sort <- df_model_STY_Merged_OTU01_draft_list_met[order(df_model_STY_Merged_OTU01_draft_list_met$met_id),]
write.table(df_model_STY_Merged_OTU01_draft_list_met_sort,file = outfile_045,quote = F,sep = '\t',row.names = F)

############ ############ ############ 

# 046.What are the stiochiometric matrix (model_xxx@S) in this model?
nrow(model_STY_Merged_OTU01_draft@S)
# [1] 825

ncol(model_STY_Merged_OTU01_draft@S)
# [1] 765


model_STY_Merged_OTU01_draft_Smat <- model_STY_Merged_OTU01_draft@S
model_STY_Merged_OTU01_draft_Smat_mx = as.matrix(model_STY_Merged_OTU01_draft@S)
class(model_STY_Merged_OTU01_draft_Smat_mx)
# [1] "matrix"

colnames(model_STY_Merged_OTU01_draft_Smat_mx)
# NULL
model_STY_Merged_OTU01_draft@react_num
# [1] 765
model_STY_Merged_OTU01_draft@react_id
colnames(model_STY_Merged_OTU01_draft_Smat_mx)<-cbind(model_STY_Merged_OTU01_draft@react_id)

row.names(model_STY_Merged_OTU01_draft_Smat_mx)
# NULL
model_STY_Merged_OTU01_draft@met_num
# 825
model_STY_Merged_OTU01_draft@met_id
row.names(model_STY_Merged_OTU01_draft_Smat_mx)<-cbind(model_STY_Merged_OTU01_draft@met_id)

write.csv(model_STY_Merged_OTU01_draft_Smat_mx, file = outfile_046,quote = F,sep = ',')
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
# aa <- get_rxn(model_STY_Merged_OTU01_draft, 'cpd00001') # example: all reactions that uses 'cpd00001' will be gave. 
# aa <- get_rxn(model_STY_Merged_OTU01_draft, 'cpd00002') # cpd00002: ATP
# aa <- get_rxn(model_STY_Merged_OTU01_draft, 'cpd00004') # cpd00004: NADH
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
# nitrate <- get_rxn(model_STY_Merged_OTU01_draft, 'cpd00209') #
# # None
# 
# 
# nitrogen <- get_rxn(model_STY_Merged_OTU01_draft, 'cpd00528') #N2-e0
# #None
# 
# Nitrite <- get_rxn(model_STY_Merged_OTU01_draft, 'cpd00075') #cpd00075: Nitrite
# # [1] "Active reactions that uses the metabolite cpd00075 Nitrite" "Active reactions that uses the metabolite cpd00075 Nitrite"
# # abbreviation	equation
# # R05046	(1) cpd00001[c] + (1) cpd03518[c] --> (1) cpd00047[c] + (1) cpd03519[c]
# # R05069	(1) cpd00498[c] <==> (1) cpd10162[c]
# # R08716	(1) cpd17535[c] --> (1) cpd08369[c]
# # R05812	(1) cpd00067[c] + (1) cpd00004[c] + (1) cpd08372[c] --> (1) cpd00003[c] + (1) cpd08373[c]
# # R05721	(1) cpd00007[c] + (1) cpd00099[c] <==> (1) cpd01054[c]
# 
# Nitrite_c0 <- get_rxn(model_STY_Merged_OTU01_draft, 'cpd00075_e') #cpd00075_e: Nitrite-c0
# # [1] "Active reactions that uses the metabolite cpd00075_e Nitrite"
# # abbreviation	equation
# # R03632	(1) cpd00002[c] + (1) cpd11900[c] --> (1) cpd00008[c] + (1) cpd12167[c]
# # R03646	(1) cpd00002[c] + (1) cpd00051[c] + (1) cpd11907[c] --> (1) cpd00018[c] + (1) cpd00012[c] + (1) cpd12036[c]
# 
# 
# NH3 <- get_rxn(model_STY_Merged_OTU01_draft, 'cpd00013') #NH3-c0
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
# model_STY_Merged_OTU01_draft@react_attr$gs.origin
# 
# ############ ############ ############ 
# # 050.What are the active reactions which use or produce a certain compound in ******model_STY_Merged_OTU01_draft****: 
# # Updates @8-1-21: this will give a new column called 'sol$fluxes' with the same length to the 'reaction id', which shows fluxes of each reaction.
# mod = model_STY_Merged_OTU01_draft
# dim(mod)
# # [1] 825 765
# 
# sol <- sybil::optimizeProb(mod, retOptSol=F) # find a list (of 7,sol), which include flux info of 765 reactions for this model.
# model_STY_Merged_OTU01_draft_flux_list <- as.list(sol$fluxes)
# length(model_STY_Merged_OTU01_draft_flux_list)
# # [1] 765
# 
# model_STY_Merged_OTU01_draft_flux_list <- as.data.frame(sol$fluxes)
# dim(model_STY_Merged_OTU01_draft_flux_list)
# # [1] 765   1
# 
# 
# model_STY_Merged_OTU01_draft_Smat <- model_STY_Merged_OTU01_draft@S
# 
# model_STY_Merged_OTU01_draft_Smat <- read.csv(paste(getwd,'STY_Merged_STY_Merged_OTU01_draft/046_model_STY_Merged_OTU01_draft_Smat.csv',sep = ""),row.names = 'X')
# t_model_STY_Merged_OTU01_draft_Smat <- as.data.frame(t(model_STY_Merged_OTU01_draft_Smat))
# dim(t_model_STY_Merged_OTU01_draft_Smat)
# # [1] 765 825
# 
# cbind_t_model_STY_Merged_OTU01_draft_Smat<- as.data.frame(cbind(model_STY_Merged_OTU01_draft_flux_list, t_model_STY_Merged_OTU01_draft_Smat))
# dim(cbind_t_model_STY_Merged_OTU01_draft_Smat)
# # [1] 765 825
# 
# write.csv(cbind_t_model_STY_Merged_OTU01_draft_Smat, file = paste(getwd,'STY_Merged_STY_Merged_OTU01_draft/050_model_STY_Merged_OTU01_draft_Smat_fluxes.csv',sep = ""))

