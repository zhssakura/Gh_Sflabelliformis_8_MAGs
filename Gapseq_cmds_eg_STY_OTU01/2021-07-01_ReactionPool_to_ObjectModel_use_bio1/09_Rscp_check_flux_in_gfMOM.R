# st6.4. Extract all reactions from the MOM_gf.

##################### st1. Get input model files and output docs ##################### 
library(stringr)
library(sybil)

getwd = getwd()
getwd

# Import the gap-filled modified object model (gfMOM): 
dir_MOM_gf = '/07_adapt_addReact_R/04_Rscp/OTU08_MOM_20210411_addReact/20210411_diet_sw_no_photo_STY_Merged_OTU08_arc__12_oxygen/'


OTU08_draft = paste(getwd, dir_MOM_gf,'STY_Merged_OTU08_MOM.RDS',sep = '')

# Assign output RDS file name:
outfile_add_043= paste(getwd, dir_MOM_gf, '/043_df_model_OTU08_draft_gfMOM_list_react.tsv',sep = '')
outfile_add_045= paste(getwd, dir_MOM_gf, '/045_df_model_OTU08_draft_gfMOM_list_metabolites.tsv',sep = '')
outfile_add_046= paste(getwd, dir_MOM_gf, '/046_model_OTU08_draft_gfMOM_Smat.csv',sep = '')
infile_046     = outfile_add_046
outfile_050    = paste(getwd, dir_MOM_gf, '/050_model_OTU08_draft_gfMOM_Smat_fluxes.csv',sep = "")

model_OTU08_draft = readRDS(OTU08_draft)
model_OTU08_draft
# model name:             STY_Merged_OTU08 
# number of compartments  3 
# c0 
# e0 
# p0 
# number of reactions:    1580 
# number of metabolites:  1532 
# number of unique genes: 1850 
# objective function:     +1 EX_cpd11416_c0 

########################## st4. Extract all reactions from the modified model (MM). #############################  
# 043.What are those reactions in the model?
model_OTU08_draft@react_num
# [1] 3118

# Extract all reactions from the addReact OM.
model_OTU08_draft_list <- list(model_OTU08_draft@react_id, model_OTU08_draft@react_name, model_OTU08_draft@react_rev, model_OTU08_draft@react_single, model_OTU08_draft@react_de,
                               model_OTU08_draft@lowbnd, model_OTU08_draft@uppbnd, model_OTU08_draft@obj_coef, model_OTU08_draft@gprRules, model_OTU08_draft@gpr,
                               model_OTU08_draft@react_attr)
df_model_OTU08_draft_list  <- as.data.frame(model_OTU08_draft_list, row.names = NULL)
colnames(df_model_OTU08_draft_list) <- c('react_id','react_name','react_rev','react_single','react_de','lowbnd','uppbnd','obj_coef','gprRules','gpr','react_attr')
# "rxn","name","ec","tc","qseqid","pident","evalue","bitscore","qcovs","stitle","sstart","send","pathway",
# "status","pathway.status","complex","exception","complex.status","seed","gs.origin")
write.table(df_model_OTU08_draft_list, file = outfile_add_043,quote = F,row.names = F,sep = '\t')
# The output contains many KEGG reaction ids, which can be converted to SEED RXN ids by using seed_reactions.tsv.


########################## st5. Extract all metabolites from the modified model (MM). #############################  
# 045.What are the metabolites in this new model?

# number of reactions of ******model_OTU08_draft****: 
model_OTU08_draft@met_num
# [1] 1481
model_OTU08_draft_met_name = model_OTU08_draft@met_name
model_OTU08_draft_met_id = model_OTU08_draft@met_id
# model_OTU08_draft: Generate a dataframe with metabolite info for each model
model_OTU08_draft_list_met <- list(model_OTU08_draft@met_id, model_OTU08_draft@met_name, model_OTU08_draft@met_comp, model_OTU08_draft@met_single, model_OTU08_draft@met_de)
df_model_OTU08_draft_list_met <- as.data.frame(model_OTU08_draft_list_met, row.names = NULL)
colnames(df_model_OTU08_draft_list_met) <- c('met_id','met_name','met_comp','met_single','met_de')
# R Sort a Data Frame using Order()
df_model_OTU08_draft_list_met_sort <- df_model_OTU08_draft_list_met[order(df_model_OTU08_draft_list_met$met_id),]
write.table(df_model_OTU08_draft_list_met_sort, file = outfile_add_045, quote = F,row.names = F,sep = '\t')


########################## st6. Extract Stoichiometric matrix from the modified model (MM). #############################  
# 046.What are the stiochiometric matrix (model_xxx@S) in this model?
nrow(model_OTU08_draft@S)
# [1] 1481

ncol(model_OTU08_draft@S)
# [1] 1447


model_OTU08_draft_Smat <- model_OTU08_draft@S
model_OTU08_draft_Smat_mx = as.matrix(model_OTU08_draft@S)
class(model_OTU08_draft_Smat_mx)
# [1] "matrix"

colnames(model_OTU08_draft_Smat_mx)
# NULL
model_OTU08_draft@react_num
# [1] 1447
model_OTU08_draft@react_id
colnames(model_OTU08_draft_Smat_mx)<-cbind(model_OTU08_draft@react_id)

row.names(model_OTU08_draft_Smat_mx)
# NULL
model_OTU08_draft@met_num
# 1481
model_OTU08_draft@met_id
row.names(model_OTU08_draft_Smat_mx)<-cbind(model_OTU08_draft@met_id)

write.csv(model_OTU08_draft_Smat_mx, file = outfile_add_046, quote = FALSE,row.names = T)
# In this table, you can check the flux (values) of a compound (e.g. cpd16503[e0]) of this model based on the reaction ID (e.g. EX_cpd16503_e0).


########################## st6. Check flux from the gfMOM. #############################  
############ 050. What are the active reactions which use or produce a certain compound in ******model_OTU08_draft****:
mod = model_OTU08_draft
dim(mod)
# [1] 2169 3118

# this will give a new column called 'sol$fluxes' with the same length to the 'reaction id', which shows fluxes of each reaction.
?optimizeProb
# more functions please read the manual of R package sybil.
dim(mod)
# [1] 2169 3118

sol <- sybil::optimizeProb(mod, retOptSol=F) # find a list (of 7,sol), which include flux info of 765 reactions for this model.
model_OTU08_draft_flux_list <- as.list(sol$fluxes)
length(model_OTU08_draft_flux_list)
# [1] 3118

model_OTU08_draft_flux_list <- as.data.frame(sol$fluxes)
dim(model_OTU08_draft_flux_list)
# [1] 3118    1


model_OTU08_draft_Smat <- mod@S

model_OTU08_draft_Smat <- read.csv(file = infile_046,row.names = 'X')
t_model_OTU08_draft_Smat <- as.data.frame(t(model_OTU08_draft_Smat))
dim(t_model_OTU08_draft_Smat)
# [1] 3118 2169

cbind_t_model_OTU08_draft_Smat<- as.data.frame(cbind(model_OTU08_draft_flux_list, t_model_OTU08_draft_Smat))
dim(cbind_t_model_OTU08_draft_Smat)
# [1] 3118 2170

write.csv(cbind_t_model_OTU08_draft_Smat, quote = FALSE,file = outfile_050)
# IMPORTANT: check the flux value of ADD_BIOMASS, if it is negative, it is not right!


######################## Follow with python script of 'Get_active_rxn_by_met_from_final_model_argv.py'######################## 
# Get_active_rxn_by_met_from_final_model_argv.py
#### IMPORTANT-2: this script used data in gapseq_v20200407, should be updated with the pipline version used for modeling.
#### update @2021-02-24: reaction direction in the output file is the same in file dat/seed_reactions_corrected.tsv
#### The true direction should combine with flux values (depending on either positive or negative)
##### The flux here does not equal to the fluxes in the gap-filling step!!!

# Create a list of compounds that you have intrests in their flux.
list_met = list('cpd00013', #NH3-c0	and NH3-e0
                'cpd00165', #Hydroxylamine
                # 'cpd00137', #citrate
                'cpd00011', #CO2
                'cpd00418', # NO
                'cpd00528', # N2
                'cpd00073', # Urea
                'cpd00075' #Nitrite
                # 'cpd00279', #Acetoacetyl-CoA
                # 'cpd03375', #3-Hydroxypropanoyl-CoA
                # 'cpd00002', #ATP-c0	1	NA	NA
                # 'cpd00022', #Acetyl-CoA
                # 'cpd00048', #Sulfate
                # 'cpd00193',#APS
                # 'cpd00081', #Sulfite
                # 'cpd00239', #Hydrogen_sulfide(H2S)
                # 'cpd00268' #Thiosulfate_(S2O3)2-
                )


mycmd1 = 'python3 /Users/zzfanyi/PycharmProjects/EMP500/GapSeq/Get_active_rxn_by_met_from_final_model_argv.py'
mycmd2 = paste('-in ',getwd, dir_MOM_gf,sep = '')
mycmd3 = '-file 050_model_OTU08_draft_gfMOM_Smat_fluxes.csv'

for (met in list_met) {
  print(met)
  mycmd4 = paste('-met ',met,sep = '')
  system(paste(mycmd1,mycmd2,mycmd3,mycmd4,sep = ' '))
  
}

# Newly added reaction can not be exported with equations in output file.

######################## Follow with python script of 'Assign_names_2_ID_of_active_rxn_in_Smax_arg.py'######################## 
# New on 2021-mar-30
# Assign_names_2_ID_of_active_rxn_in_Smax_arg.py
#### IMPORTANT-2: this script used data in gapseq_v20200407, should be updated with the pipline version used for modeling.
#### update @2021-02-24: reaction direction in the output file is the same in file dat/seed_reactions_corrected.tsv
#### The true direction should combine with flux values (depending on either positive or negative)
######################  run script in command line.
getwd
input_dir = paste(getwd,'/07_adapt_addReact_R/04_Rscp/OTU08_MOM_20210411_addReact/20210411_diet_sw_no_photo_STY_Merged_OTU08_arc__12_oxygen/',sep = '')
input_file = '050_model_OTU08_draft_gfMOM_Smat_fluxes.csv'
dict_customized_rxn_file = paste(getwd,'/02_Get_RP_equations/Reaction_id_in_ReactionPool_Addreact_20210411.txt',sep = '')
outfile = '051_model_OTU08_draft_gfMOM_Smat_fluxes.annotate.tsv'


mycmd1 = 'python3 /Users/zzfanyi/PycharmProjects/EMP500/GapSeq/Assign_names_2_ID_of_active_rxn_in_Smax_arg.py'

system(paste(mycmd1,'-in',input_dir,'-file',input_file,'-cus_rxn_file',dict_customized_rxn_file,'-outfile',outfile,sep = ' '))


