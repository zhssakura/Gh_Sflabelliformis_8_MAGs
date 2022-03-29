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
## addReact:
# NmSCM1_draft = paste(getwd,'/Input_files/N.maritimusSCM1-draft.RDS',sep = '')
## adapt-add + addReact:
# NmSCM1_draft = paste(getwd,'/Input_files/NmSCM1-adapt.RDS',sep = '')
## addReact (new draft created using customized biomass_archaea.json):
# NmSCM1_draft = paste(getwd,'/Input_files/draft_biomass_archaea_20210417/N.maritimusSCM1-draft.RDS',sep = '')
# NmSCM1_draft = paste(getwd,'/Input_files/draft_biomass_archaea_20210418/N.maritimusSCM1-draft.RDS',sep = '')
# NmSCM1_draft = paste(getwd,'/Input_files/draft_biomass_archaea_original/N.maritimusSCM1-draft.RDS',sep = '')
# NmSCM1_draft = paste(getwd,'/Input_files/draft_biomass_archaea_original/adapt_rm/N.maritimusSCM1-adapt.RDS',sep = '')
NmSCM1_draft = paste(getwd,'/Input_files/draft_biomass_archaea_original/N.maritimusSCM1-draft.RDS',sep = '')

model_NmSCM1_draft = readRDS(NmSCM1_draft)
model_NmSCM1_draft
# model name:             N.maritimusSCM1 
# number of compartments  3 
# c0 
# e0 
# p0 
# number of reactions:    1826 
# number of metabolites:  1835 
# number of unique genes: 630 
# objective function:     +1 bio1 


# Define the date you work:
# workdate = '20210327pm1'
# workdate = '20210327pm2'
# workdate = '20210327pm3'
# workdate = '20210330pm1'
# workdate = '20210410'
# workdate = '20210411'
# workdate = '20210411pm2'
# workdate = '20210418'
# workdate = '20210418pm2'
# workdate = '20210419pm1'
# workdate = '20210419pm2'
# workdate = '20210419pm2_2'
# workdate = '20210419pm1_2'
# workdate = '20210419pm1_3'
workdate = '20210419pm1_4'

# Import a data frame that include all reactions from reaction pool (RP):
React_id = read.delim(paste(getwd,'/02_Get_RP_equations/Reaction_id_in_ReactionPool_Addreact_',workdate,'.txt',sep = ''))

# create a new folder:
new_folder = paste('NmSCM1_MOM_',workdate,'_bioArc_original',sep = '')
dir.create(paste(getwd,'/07_adapt_addReact_R/04_Rscp/',new_folder,sep = ''))

# Assign output RDS file name:
outputRDS      = paste(getwd, '/07_adapt_addReact_R/04_Rscp/',new_folder,'/NmSCM1_MOM.RDS',sep = '')
outfile_rmbio1 = paste(getwd, '/07_adapt_addReact_R/04_Rscp/',new_folder,'/043_model_NmSCM1_MOM_rmReact_list_react.tsv',sep = '')
outfile_add_043= paste(getwd, '/07_adapt_addReact_R/04_Rscp/',new_folder,'/043_model_NmSCM1_MOM_list_react.tsv',sep = '')
outfile_add_045= paste(getwd, '/07_adapt_addReact_R/04_Rscp/',new_folder,'/045_model_NmSCM1_MOM_list_metabolites.tsv',sep = '')
outfile_add_046= paste(getwd, '/07_adapt_addReact_R/04_Rscp/',new_folder,'/046_model_NmSCM1_MOM_draft_Smat.csv',sep = '')
########### st1. Remove original biomass equation in OM! (Cause we will add new bio1 to the OM) ########### 
# model_NmSCM1_draft <- rmReact(model_NmSCM1_draft, react = 'bio1', rm_met = T)
# # Extract all reactions from the rmReact OM.
# model_NmSCM1_draft_list <- list(model_NmSCM1_draft@react_id, model_NmSCM1_draft@react_name, model_NmSCM1_draft@react_rev, model_NmSCM1_draft@react_single, model_NmSCM1_draft@react_de,
#                                model_NmSCM1_draft@lowbnd, model_NmSCM1_draft@uppbnd, model_NmSCM1_draft@obj_coef, model_NmSCM1_draft@gprRules, model_NmSCM1_draft@gpr,
#                                model_NmSCM1_draft@react_attr)
# df_model_NmSCM1_draft_list  <- as.data.frame(model_NmSCM1_draft_list, row.names = NULL)
# colnames(df_model_NmSCM1_draft_list) <- c('react_id','react_name','react_rev','react_single','react_de','lowbnd','uppbnd','obj_coef','gprRules','gpr','react_attr')
# # "rxn","name","ec","tc","qseqid","pident","evalue","bitscore","qcovs","stitle","sstart","send","pathway",
# # "status","pathway.status","complex","exception","complex.status","seed","gs.origin")
# write.table(df_model_NmSCM1_draft_list, file = outfile_rmbio1, quote = F,row.names = F,sep = '\t')
# # The output contains many KEGG reaction ids, which can be converted to SEED RXN ids by using seed_reactions.tsv.

##########################  Start with a test #############################  
########### add the first RP reaction to the OM from data frame ########### 

# Here is the 1st reaction from reaction pool (RP) that we want to add to the object model (OM):
React_id[1,]

# Check the last reaction before adding new reaction(s):
model_NmSCM1_draft@react_id[model_NmSCM1_draft@react_num]
# [1] "DM_cpd01042_c0"

as.character(React_id[1,]$Rid)
# [1] "R00148_Shan_OTU08"

as.character(React_id[1,]$react_name)
# [1] "ammonia,ubiquinol:oxygen oxidoreductase"

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

model_NmSCM1_draft_new <- addReact(
  model_NmSCM1_draft,
  id         = as.character(React_id[1,]$Rid),
  met        = str_split(as.character(React_id[1,]$met_id),',')[[1]],
  Scoef      = as.integer(str_split(as.character(React_id[1,]$met_Scoef),',')[[1]]),
  reactName  = as.character(React_id[1,]$react_name),
  reversible = React_id[1,]$react_rev,
  metName    = str_split(as.character(React_id[1,]$met_name),',')[[1]],
  metComp    = as.integer(str_split(as.character(React_id[1,]$met_comp),',')[[1]]),
  lb         = as.integer(React_id[1,]$lowbnd),
  ub         = as.integer(React_id[1,]$uppbnd)
  
)

# check the latest reaction added:
model_NmSCM1_draft_new@react_id[model_NmSCM1_draft_new@react_num]
# [1] "R00148_Shan_OTU08"
model_NmSCM1_draft_new@react_name[model_NmSCM1_draft_new@react_num]

##########################  Test ends #############################  


########################## st2. Use For loop to add all RP to OM #############################  
model_NmSCM1_draft

for (i in 1:nrow(React_id)) {
  model_NmSCM1_draft <- addReact(
    model_NmSCM1_draft,
    id         = as.character(React_id[i,]$Rid),
    met        = str_split(as.character(React_id[i,]$met_id),',')[[1]],
    Scoef      = as.integer(str_split(as.character(React_id[i,]$met_Scoef),',')[[1]]),
    reactName  = as.character(React_id[i,]$react_name),
    reversible = React_id[i,]$react_rev,
    metName    = str_split(as.character(React_id[i,]$met_name),',')[[1]],
    metComp    = as.integer(str_split(as.character(React_id[i,]$met_comp),',')[[1]]),
    lb         = as.integer(React_id[i,]$lowbnd),
    ub         = as.integer(React_id[i,]$uppbnd)
  )
  print(i)
 }

# check the latest reaction added:
model_NmSCM1_draft@react_id[model_NmSCM1_draft@react_num]
# [1] "R00238_PWY.5789_18/19_Shan_OTU08"

model_NmSCM1_draft
# model name:             N.maritimusSCM1 
# number of compartments  3 
# c0 
# e0 
# p0 
# number of reactions:    1850 
# number of metabolites:  1836 
# number of unique genes: 630 
# objective function:     +1 bio1 
########################## st3. Sets/changes the Objective Function to OM #############################  
# # changeObjFunc: Sets/changes the Objective Function
# model_NmSCM1_draft <- changeObjFunc(model_NmSCM1_draft, react='ADD_BIOMASS', obj_coef = rep(1, length(react)))
# 
# model_NmSCM1_draft
# # model name:             NmSCM1 
# # number of compartments  3 
# # c0 
# # e0 
# # p0 
# # number of reactions:    1429 
# # number of metabolites:  1466 
# # number of unique genes: 1845 
# # objective function:     +1 ADD_BIOMASS 
# ############ ############ ############ 

########################## st4. Extract all reactions from the modified object model (MOM). #############################  
# 043.What are those reactions in the MOM model?
model_NmSCM1_draft@react_num
# [1] 1850

# Extract all reactions from the addReact MOM.
model_NmSCM1_draft_list <- list(model_NmSCM1_draft@react_id, model_NmSCM1_draft@react_name, model_NmSCM1_draft@react_rev, model_NmSCM1_draft@react_single, model_NmSCM1_draft@react_de,
                               model_NmSCM1_draft@lowbnd, model_NmSCM1_draft@uppbnd, model_NmSCM1_draft@obj_coef, model_NmSCM1_draft@gprRules, model_NmSCM1_draft@gpr,
                               model_NmSCM1_draft@react_attr)
df_model_NmSCM1_draft_list  <- as.data.frame(model_NmSCM1_draft_list, row.names = NULL)

colnames(df_model_NmSCM1_draft_list) <- c('react_id','react_name','react_rev','react_single','react_de','lowbnd','uppbnd','obj_coef','gprRules','gpr'
                                          ,'react_attr',"seed","name","ec","tc","qseqid","pident","evalue","bitscore","qcovs","stitle","sstart"
                                          ,"send","pathway","blast.status","pathway.status","complex","exception","complex.status","gs.origin","annotation","MNX_ID","seedID","keggID","biggID","biocycID"
)

write.table(df_model_NmSCM1_draft_list, file = outfile_add_043,quote = F,row.names = F,sep = '\t')
# The output contains many KEGG reaction ids, which can be converted to SEED RXN ids by using seed_reactions.tsv.


########################## st5. Extract all metabolites from the modified model (MM). #############################  
# 045.What are the metabolites in this new model?

# number of reactions of ******model_NmSCM1_draft****: 
model_NmSCM1_draft@met_num
# [1] 1836
model_NmSCM1_draft_met_name = model_NmSCM1_draft@met_name
model_NmSCM1_draft_met_id = model_NmSCM1_draft@met_id
# model_NmSCM1_draft: Generate a dataframe with metabolite info for each model
model_NmSCM1_draft_list_met <- list(model_NmSCM1_draft@met_id, model_NmSCM1_draft@met_name, model_NmSCM1_draft@met_comp, model_NmSCM1_draft@met_single, model_NmSCM1_draft@met_de)
df_model_NmSCM1_draft_list_met <- as.data.frame(model_NmSCM1_draft_list_met, row.names = NULL)
colnames(df_model_NmSCM1_draft_list_met) <- c('met_id','met_name','met_comp','met_single','met_de')
# R Sort a Data Frame using Order()
df_model_NmSCM1_draft_list_met_sort <- df_model_NmSCM1_draft_list_met[order(df_model_NmSCM1_draft_list_met$met_id),]
write.table(df_model_NmSCM1_draft_list_met_sort, file = outfile_add_045, quote = F,row.names = F,sep = '\t')


########################## st6. Extract Stoichiometric matrix from the modified model (MM). #############################  
# 046.What are the stiochiometric matrix (model_xxx@S) in this model?
nrow(model_NmSCM1_draft@S)
# [1] 1836

ncol(model_NmSCM1_draft@S)
# [1] 1850


model_NmSCM1_draft_Smat <- model_NmSCM1_draft@S
model_NmSCM1_draft_Smat_mx = as.matrix(model_NmSCM1_draft@S)
class(model_NmSCM1_draft_Smat_mx)
# [1] "matrix"

colnames(model_NmSCM1_draft_Smat_mx)
# NULL
model_NmSCM1_draft@react_num
# [1] 1427
model_NmSCM1_draft@react_id
colnames(model_NmSCM1_draft_Smat_mx)<-cbind(model_NmSCM1_draft@react_id)

row.names(model_NmSCM1_draft_Smat_mx)
# NULL
model_NmSCM1_draft@met_num
# 1467
model_NmSCM1_draft@met_id
row.names(model_NmSCM1_draft_Smat_mx)<-cbind(model_NmSCM1_draft@met_id)

write.csv(model_NmSCM1_draft_Smat_mx, file = outfile_add_046, quote = F,row.names = T)
# In this table, you can check the flux (values) of a compound (e.g. cpd16503[e0]) of this model based on the reaction ID (e.g. EX_cpd16503_e0).

#################### st7. Save the modified model (MOM) ##################
saveRDS(model_NmSCM1_draft,file = outputRDS)
# Reactions were removed, then added to MOM.
# Go to ./07_adapt_addReact_R/04_Rscp to get MOM model.

