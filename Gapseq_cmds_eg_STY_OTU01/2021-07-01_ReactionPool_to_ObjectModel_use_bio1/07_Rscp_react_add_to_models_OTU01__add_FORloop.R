# st6.1. Add reactions to a gs model with FOR-loop
# st6.2. Use For-loop to add all RP to OM
# st6.3. Sets/changes the Objective Function to OM
# st6.4. Extract all reactions from the MOM.

##################### st1. Add reactions to a gs model with FOR loop ##################### 
library(stringr)
library(sybil)

getwd = getwd()
getwd

# Import the object model (OM): 
# addReact:
STY_Merged_OTU01_draft = paste(getwd,'/Input_files/STY_Merged_OTU01-draft.RDS',sep = '')
# adapt-add + addReact:
# STY_Merged_OTU01_draft = paste(getwd,'/Input_files/STY_Merged_OTU01-adapt.RDS',sep = '')

model_STY_Merged_OTU01_draft = readRDS(STY_Merged_OTU01_draft)
model_STY_Merged_OTU01_draft
# model name:             STY_Merged_OTU01 
# number of compartments  3 
# c0 
# e0 
# p0 
# number of reactions:    1394 
# number of metabolites:  1462 
# number of unique genes: 1841 
# objective function:     +1 bio1 

# Define the date you work:
# workdate = '20210327pm1'
# workdate = '20210327pm2'
# workdate = '20210327pm3'
# workdate = '20210330pm1'
# workdate = '20210410'
# workdate = '20210411'
# workdate = '20210411pm2'
# workdate = '20210707'
# workdate = '20210823'
# workdate = '20210823pm2'
workdate = '20210825'

# Import a data frame that include all reactions from reaction pool (RP):
React_id = read.delim(paste(getwd,'/02_Get_RP_equations/Reaction_id_in_ReactionPool_Addreact_',workdate,'.txt',sep = ''))

# create a new folder:
new_folder = paste('OTU01_MOM_',workdate,sep = '')
dir.create(paste(getwd,'/07_adapt_addReact_R/04_Rscp/',new_folder,sep = ''))

# Assign output RDS file name:
outputRDS      = paste(getwd, '/07_adapt_addReact_R/04_Rscp/',new_folder,'/STY_Merged_OTU01_MOM.RDS',sep = '')
outfile_rmbio1 = paste(getwd, '/07_adapt_addReact_R/04_Rscp/',new_folder,'/043_df_model_STY_Merged_OTU01_rmReact_list_react.tsv',sep = '')
outfile_add_043= paste(getwd, '/07_adapt_addReact_R/04_Rscp/',new_folder,'/043_df_model_STY_Merged_OTU01_rmReact_addReact_list_react.tsv',sep = '')
outfile_add_045= paste(getwd, '/07_adapt_addReact_R/04_Rscp/',new_folder,'/045_df_model_STY_Merged_OTU01_rmReact_addReact_list_metabolites.tsv',sep = '')
outfile_add_046= paste(getwd, '/07_adapt_addReact_R/04_Rscp/',new_folder,'/046_model_STY_Merged_OTU01_draft_Smat.csv',sep = '')
########### st1. Remove original biomass equation in OM! (Cause we will add new bio1 to the OM) ########### 
# model_STY_Merged_OTU01_draft <- rmReact(model_STY_Merged_OTU01_draft, react = 'bio1', rm_met = T)
# # Extract all reactions from the rmReact OM.
# model_STY_Merged_OTU01_draft_list <- list(model_STY_Merged_OTU01_draft@react_id, model_STY_Merged_OTU01_draft@react_name, model_STY_Merged_OTU01_draft@react_rev, model_STY_Merged_OTU01_draft@react_single, model_STY_Merged_OTU01_draft@react_de,
#                                model_STY_Merged_OTU01_draft@lowbnd, model_STY_Merged_OTU01_draft@uppbnd, model_STY_Merged_OTU01_draft@obj_coef, model_STY_Merged_OTU01_draft@gprRules, model_STY_Merged_OTU01_draft@gpr,
#                                model_STY_Merged_OTU01_draft@react_attr)
# df_model_STY_Merged_OTU01_draft_list  <- as.data.frame(model_STY_Merged_OTU01_draft_list, row.names = NULL)
# colnames(df_model_STY_Merged_OTU01_draft_list) <- c('react_id','react_name','react_rev','react_single','react_de','lowbnd','uppbnd','obj_coef','gprRules','gpr','react_attr')
# # "rxn","name","ec","tc","qseqid","pident","evalue","bitscore","qcovs","stitle","sstart","send","pathway",
# # "status","pathway.status","complex","exception","complex.status","seed","gs.origin")
# write.table(df_model_STY_Merged_OTU01_draft_list, file = outfile_rmbio1, quote = F,row.names = F,sep = '\t')
# # The output contains many KEGG reaction ids, which can be converted to SEED RXN ids by using seed_reactions.tsv.

##########################  Start with a test #############################  
########### add the first RP reaction to the OM from data frame ########### 

# Here is the 1st reaction from reaction pool (RP) that we want to add to the object model (OM):
React_id[1,]

# Check the last reaction before adding new reaction(s):
model_STY_Merged_OTU01_draft@react_id[model_STY_Merged_OTU01_draft@react_num]
# [1] "DM_cpd01042_c0"

as.character(React_id[1,]$Rid)
# [1] "AMO_R00148_Shan_STY_Merged_OTU01"

str_split(as.character(React_id[1,]$met_id),',')[[1]]
# [1] "cpd00109[c0]" "cpd00075[c0]" "cpd00001[c0]" "cpd00110[c0]" "cpd00209[c0]" "cpd00067[c0]"

as.integer(str_split(as.character(React_id[1,]$met_Scoef),',')[[1]])
# [1] -2 -1 -1  2  1  2

as.character(React_id[1,]$Rid)
# [1] "R247-RXN"

React_id[1,]$react_rev
# [1] FALSE
str_split(as.character(React_id[1,]$met_name),',')[[1]]
# [1] "Cytochrome c3+-c0" "Nitrite-c0"        "H2O-c0"            "Cytochrome c2+-c0" "Nitrate-c0"        "H+-c0"            

as.integer(str_split(as.character(React_id[1,]$met_comp),',')[[1]])
# [1] 1 1 1 1 1 1

as.integer(React_id[1,]$lowbnd)
# [1] 0
as.integer(React_id[1,]$uppbnd)
# [1] 1000


# The following commond should equal to:
# model_name <- addReact(
#   model_name,
#   id='Trans_NH4', 
#   met=c('cpd00013[c0]','cpd00067[e0]','cpd00013[e0]'),
#   metName = c('NH3','H+[e]','NH4'),
#   Scoef=c(1,1,-1),
#   metComp = c(1,1,1,1)
#   reactName = 'Trans_NH4',
#   reversible= FALSE,
#   lb = 0,
#   up = 1000
# )

model_STY_Merged_OTU01_draft_new <- addReact(
  model_STY_Merged_OTU01_draft,
  id         = as.character(React_id[1,]$Rid),
  met        = str_split(as.character(React_id[1,]$met_id),',')[[1]],
  Scoef      = as.integer(str_split(as.character(React_id[1,]$met_Scoef),',')[[1]]),
  reactName  = as.character(React_id[1,]$Rid),
  reversible = React_id[1,]$react_rev,
  metName    = str_split(as.character(React_id[1,]$met_name),',')[[1]],
  metComp    = as.integer(str_split(as.character(React_id[1,]$met_comp),',')[[1]]),
  lb         = as.integer(React_id[1,]$lowbnd),
  ub         = as.integer(React_id[1,]$uppbnd)
  
)

# check the latest reaction added:
model_STY_Merged_OTU01_draft_new@react_id[model_STY_Merged_OTU01_draft_new@react_num]
# [1] "AMO_R00148_Shan_STY_Merged_OTU01"

##########################  Test ends #############################  


########################## st2. Use For loop to add all RP to OM #############################  
model_STY_Merged_OTU01_draft

for (i in 1:nrow(React_id)) {
  model_STY_Merged_OTU01_draft <- addReact(
    model_STY_Merged_OTU01_draft,
    id         = as.character(React_id[i,]$Rid),
    met        = str_split(as.character(React_id[i,]$met_id),',')[[1]],
    Scoef      = as.integer(str_split(as.character(React_id[i,]$met_Scoef),',')[[1]]),
    reactName  = as.character(React_id[i,]$Rid),
    reversible = React_id[i,]$react_rev,
    metName    = str_split(as.character(React_id[i,]$met_name),',')[[1]],
    metComp    = as.integer(str_split(as.character(React_id[i,]$met_comp),',')[[1]]),
    lb         = as.integer(React_id[i,]$lowbnd),
    ub         = as.integer(React_id[i,]$uppbnd)
  )
  print(i)
 }

# check the latest reaction added:
model_STY_Merged_OTU01_draft@react_id[model_STY_Merged_OTU01_draft@react_num]
# [1] "rxn13892"

model_STY_Merged_OTU01_draft
# model name:             STY_Merged_OTU01 
# number of compartments  3 
# c0 
# e0 
# p0 
# number of reactions:    1537 
# number of metabolites:  1547 
# number of unique genes: 1841 
# objective function:     +1 bio1 

########################## st3. Sets/changes the Objective Function to OM #############################  
# # changeObjFunc: Sets/changes the Objective Function
# model_STY_Merged_OTU01_draft <- changeObjFunc(model_STY_Merged_OTU01_draft, react='ADD_BIOMASS', obj_coef = rep(1, length(react)))
# 
# model_STY_Merged_OTU01_draft
# # model name:             STY_Merged_OTU01 
# # number of compartments  3 
# # c0 
# # e0 
# # p0 
# # number of reactions:    1429 
# # number of metabolites:  1466 
# # number of unique genes: 1845 
# # objective function:     +1 ADD_BIOMASS 
# ############ ############ ############ 

########################## st4. Extract all reactions from the modified model (MM). #############################  
# 043.What are those reactions in the model?
model_STY_Merged_OTU01_draft@react_num
# [1] 1427

# Extract all reactions from the addReact OM.
model_STY_Merged_OTU01_draft_list <- list(model_STY_Merged_OTU01_draft@react_id, model_STY_Merged_OTU01_draft@react_name, model_STY_Merged_OTU01_draft@react_rev, model_STY_Merged_OTU01_draft@react_single, model_STY_Merged_OTU01_draft@react_de,
                               model_STY_Merged_OTU01_draft@lowbnd, model_STY_Merged_OTU01_draft@uppbnd, model_STY_Merged_OTU01_draft@obj_coef, model_STY_Merged_OTU01_draft@gprRules, model_STY_Merged_OTU01_draft@gpr,
                               model_STY_Merged_OTU01_draft@react_attr)
df_model_STY_Merged_OTU01_draft_list  <- as.data.frame(model_STY_Merged_OTU01_draft_list, row.names = NULL)
colnames(df_model_STY_Merged_OTU01_draft_list) <- c('react_id','react_name','react_rev','react_single','react_de','lowbnd','uppbnd','obj_coef','gprRules','gpr','react_attr')
# "rxn","name","ec","tc","qseqid","pident","evalue","bitscore","qcovs","stitle","sstart","send","pathway",
# "status","pathway.status","complex","exception","complex.status","seed","gs.origin")
write.table(df_model_STY_Merged_OTU01_draft_list, file = outfile_add_043,quote = F,row.names = F,sep = '\t')
# The output contains many KEGG reaction ids, which can be converted to SEED RXN ids by using seed_reactions.tsv.


########################## st5. Extract all metabolites from the modified model (MM). #############################  
# 045.What are the metabolites in this new model?

# number of reactions of ******model_STY_Merged_OTU01_draft****: 
model_STY_Merged_OTU01_draft@met_num
# [1] 1467
model_STY_Merged_OTU01_draft_met_name = model_STY_Merged_OTU01_draft@met_name
model_STY_Merged_OTU01_draft_met_id = model_STY_Merged_OTU01_draft@met_id
# model_STY_Merged_OTU01_draft: Generate a dataframe with metabolite info for each model
model_STY_Merged_OTU01_draft_list_met <- list(model_STY_Merged_OTU01_draft@met_id, model_STY_Merged_OTU01_draft@met_name, model_STY_Merged_OTU01_draft@met_comp, model_STY_Merged_OTU01_draft@met_single, model_STY_Merged_OTU01_draft@met_de)
df_model_STY_Merged_OTU01_draft_list_met <- as.data.frame(model_STY_Merged_OTU01_draft_list_met, row.names = NULL)
colnames(df_model_STY_Merged_OTU01_draft_list_met) <- c('met_id','met_name','met_comp','met_single','met_de')
# R Sort a Data Frame using Order()
df_model_STY_Merged_OTU01_draft_list_met_sort <- df_model_STY_Merged_OTU01_draft_list_met[order(df_model_STY_Merged_OTU01_draft_list_met$met_id),]
write.table(df_model_STY_Merged_OTU01_draft_list_met_sort, file = outfile_add_045, quote = F,row.names = F,sep = '\t')


########################## st6. Extract Stoichiometric matrix from the modified model (MM). #############################  
# 046.What are the stiochiometric matrix (model_xxx@S) in this model?
# nrow(model_STY_Merged_OTU01_draft@S)
# # [1] 1465
# 
# ncol(model_STY_Merged_OTU01_draft@S)
# # [1] 1420
# 
# 
# model_STY_Merged_OTU01_draft_Smat <- model_STY_Merged_OTU01_draft@S
# model_STY_Merged_OTU01_draft_Smat_mx = as.matrix(model_STY_Merged_OTU01_draft@S)
# class(model_STY_Merged_OTU01_draft_Smat_mx)
# # [1] "matrix"
# 
# colnames(model_STY_Merged_OTU01_draft_Smat_mx)
# # NULL
# model_STY_Merged_OTU01_draft@react_num
# # [1] 1427
# model_STY_Merged_OTU01_draft@react_id
# colnames(model_STY_Merged_OTU01_draft_Smat_mx)<-cbind(model_STY_Merged_OTU01_draft@react_id)
# 
# row.names(model_STY_Merged_OTU01_draft_Smat_mx)
# # NULL
# model_STY_Merged_OTU01_draft@met_num
# # 1467
# model_STY_Merged_OTU01_draft@met_id
# row.names(model_STY_Merged_OTU01_draft_Smat_mx)<-cbind(model_STY_Merged_OTU01_draft@met_id)
# 
# write.csv(model_STY_Merged_OTU01_draft_Smat_mx, file = outfile_add_046, quote = F,row.names = T)
# In this table, you can check the flux (values) of a compound (e.g. cpd16503[e0]) of this model based on the reaction ID (e.g. EX_cpd16503_e0).

#################### st7. Save the modified model (MOM) ##################
saveRDS(model_STY_Merged_OTU01_draft,file = outputRDS)
# Reactions were removed, then added to MOM.
# Go to ./07_adapt_addReact_R/04_Rscp to get MOM model.

