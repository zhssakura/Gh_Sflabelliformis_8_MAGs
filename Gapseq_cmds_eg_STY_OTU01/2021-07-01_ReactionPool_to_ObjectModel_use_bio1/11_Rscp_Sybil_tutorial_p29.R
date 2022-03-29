library(sybil)
data(Ec_core)
# 4.13 Interacting with the optimization process
# 
# Method optimizeProb provides a basic mechanism to run commands before and or after solving the optimization problem. In order to retrieve the reduced costs after the optimization, use argument poCmd.

opt <- optimizeProb(Ec_core, poCmd = list("getRedCosts"))
postProc(opt)


##################### st1. Get input model files and output docs ##################### 
library(stringr)
library(sybil)

getwd = getwd()
getwd

# Import the gap-filled modified object model (gfMOM): 
dir_MOM_gf = '/07_adapt_addReact_R/04_Rscp/OTU08_MOM_20210330pm1/20210330pm1_04.medium_Jo_no_photo_STY_Merged_OTU08_arc__04_Rscp_NOcitrate_sulfate/'


OTU08_draft = paste(getwd, dir_MOM_gf,'STY_Merged_OTU08_MOM.RDS',sep = '')

# # Assign output RDS file name:
# outfile_add_043= paste(getwd, dir_MOM_gf, '/043_df_model_OTU08_draft_gfMOM_list_react.tsv',sep = '')
# outfile_add_045= paste(getwd, dir_MOM_gf, '/045_df_model_OTU08_draft_gfMOM_list_metabolites.tsv',sep = '')
# outfile_add_046= paste(getwd, dir_MOM_gf, '/046_model_OTU08_draft_gfMOM_Smat.csv',sep = '')
# infile_046     = outfile_add_046
# outfile_050    = paste(getwd, dir_MOM_gf, '/050_model_OTU08_draft_gfMOM_Smat_fluxes.csv',sep = "")
outfile_052    = paste(getwd, dir_MOM_gf, '/052_model_OTU08_draft_gfMOM_DeadEndMet.csv',sep = '')
outfile_053    = paste(getwd, dir_MOM_gf, '/053_model_OTU08_draft_gfMOM_findExchReact.csv',sep = '')

model_OTU08_draft = readRDS(OTU08_draft)
model_OTU08_draft


library('dplyr')
opt <- optimizeProb(model_OTU08_draft, poCmd = list("getRedCosts")) #getRedCosts: get reduced costs of all variables after optimization.
result <- postProc(opt)
result@pa[[1]]


df <- data.frame('reducecost'=result@pa[[1]], 'react_id'=model_OTU08_draft@react_id, 'react_name'=model_OTU08_draft@react_name)
df1 <- dplyr::filter(df, df$reducecost!=0)

# change boundary of the old model.
new_model <- changeBounds(model_OTU08_draft,react = 'EX_cpd00007_e0', lb = -20, ub = 1000)# do this bf gapfilling
opt <- optimizeProb(new_model,algorithm='fba', retOptSol='F')
opt
# check mtf (minimizing total flux), fba
df <- data.frame('flux'=opt[["fluxes"]], 'react_id'=model_OTU08_draft@react_id, 'react_name'=model_OTU08_draft@react_name)

opt <- optimizeProb(new_model,algorithm='fba', retOptSol='F')


result_new <- postProc(opt)




#4.5 Minimize total ﬂux