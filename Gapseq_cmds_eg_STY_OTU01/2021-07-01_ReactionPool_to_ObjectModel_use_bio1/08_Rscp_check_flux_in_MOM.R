########## Compare get_ALL_rxn(mod, met) and function(mod, met) ########## 
library(stringr)
library(sybil)

getwd = getwd()
getwd

# Import the object model (OM): 
OTU08_draftMM = paste(getwd,'/07_adapt_addReact_R/04_Rscp/STY_Merged_OTU08-draft_MOM.RDS',sep = '')
model_OTU08_draftMM = readRDS(OTU08_draftMM)
model_OTU08_draftMM

met = 'cpd00001'
mod = model_OTU08_draftMM
dim(mod)
# [1] 1466 1429

infile_046 = paste(getwd,'/07_adapt_addReact_R/04_Rscp/046_model_OTU08_draft_Smat.csv',sep = "")
outfile_050 =paste(getwd,'/08_Check_flux_ModifiedModel/050_model_OTU08_draft_MOM_Smat_fluxes.csv',sep = "")
###################### 049. Shan modified: Get ALL reactions (both active and non-active reactions) which use or produce a certain compound in ******model2_con****: 
get_ALL_rxn <- function(mod, met)
{sol <- sybil::optimizeProb(mod, retOptSol=F) # find a list (of 7,sol), which include flux info of 2187 reactions for this model.
met.idx <- grep(met, mod@met_id)             # get the index of target metabolites from the original stoichiomatric matrix (1 and 1901).
for (i in as.list(met.idx)){
  print(i)
}

print(paste("All reactions (no matter avtive or not) that uses the metabolite", met, mod@met_name[met.idx])) #[1] "H2O-c0" "H2O-e0"

for (i in as.list(met.idx)){
  rxn.idx <- which(mod@S[i,]!=0)              # get the index of reaction(s) in which target metabolite is not =0 from the original stoichiomatric matrix 'mod@S'.
  printReaction(mod, react=rxn.idx)           # get the reactions in the model by providing index in the original toichiomatric matrix.
}
}



get_ALL_rxn(mod, met) # 736 reactions

# Following with a python script to detect reactions with flux not = 0.



############ 050. What are the active reactions which use or produce a certain compound in ******model_OTU08_draft****:
# this will give a new column called 'sol$fluxes' with the same length to the 'reaction id', which shows fluxes of each reaction.
?optimizeProb
# more functions please read the manual of R package sybil.
dim(mod)
# [1] 1466 1429

sol <- sybil::optimizeProb(mod, retOptSol=F) # find a list (of 7,sol), which include flux info of 765 reactions for this model.
model_OTU08_draft_flux_list <- as.list(sol$fluxes)
length(model_OTU08_draft_flux_list)
# [1] 1429

model_OTU08_draft_flux_list <- as.data.frame(sol$fluxes)
dim(model_OTU08_draft_flux_list)
# [1] 1429    1


model_OTU08_draft_Smat <- mod@S

model_OTU08_draft_Smat <- read.csv(file = infile_046,row.names = 'X')
t_model_OTU08_draft_Smat <- as.data.frame(t(model_OTU08_draft_Smat))
dim(t_model_OTU08_draft_Smat)
# [1] 1429 1466

cbind_t_model_OTU08_draft_Smat<- as.data.frame(cbind(model_OTU08_draft_flux_list, t_model_OTU08_draft_Smat))
dim(cbind_t_model_OTU08_draft_Smat)
# [1] 1429 1466

write.csv(cbind_t_model_OTU08_draft_Smat, file = outfile_050)
# IMPORTANT: check the flux value of ADD_BIOMASS, if it is negative, it is not right!


######################## Follow with python script of 'Get_active_rxn_by_met_from_final_model_argv.py'######################## 
# Get_active_rxn_by_met_from_final_model_argv.py

mycmd1 = 'python3 /Users/zzfanyi/PycharmProjects/EMP500/GapSeq/Get_active_rxn_by_met_from_final_model_argv.py'
mycmd2 = paste('-in ',getwd,'/08_Check_flux_ModifiedModel/',sep = '')
mycmd3 = '-file 050_model_OTU08_draft_MOM_Smat_fluxes.csv'
mycmd4 = '-met cpd00002'
system(paste(mycmd1,mycmd2,mycmd3,mycmd4,sep = ' '))
# Newly added reaction can not be exported with equations in output file.

