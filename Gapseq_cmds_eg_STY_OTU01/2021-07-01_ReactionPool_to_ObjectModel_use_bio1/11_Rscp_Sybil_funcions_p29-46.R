# st6.4. Extract all reactions from the MOM_gf.

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
# model name:             STY_Merged_OTU08 
# number of compartments  3 
# c0 
# e0 
# p0 
# number of reactions:    1659 
# number of metabolites:  1606 
# number of unique genes: 1841 
# objective function:     +1 EX_cpd11416_c0 

class(model_OTU08_draft)
# [1] "modelorg"
# attr(,"package")
# [1] "sybil"
########################## st2. Functions of sybil from p29-46. #############################  
# changeRowsBnds:Change Row Bounds in the Optimization Problem
## S4 method for signature 'optObj_clpAPI' changeRowsBnds(lp, i, lb, ub)
# lp An object extending class optObj. 
# i A numeric vector containing the row indices of the constraints to change. 
# lb A numeric vector of the same length as i containing the lower bounds of the constraints to change. 
# ub A numeric vector of the same length as i containing the upper bounds of the constraints to change.

changeRowsBnds()



# Test, if a given algorithm can has a certain purpose.
# Returns TRUE if successful, otherwise FALSE.
checkAlgorithm('sbalg')



## show all current parameters
SYBIL_SETTINGS()

# checkDefaultMethod
checkDefaultMethod('glpkAPI',"simplex","lp")

# delProb-methods:Free Memory Associated to the Pointer to the Problem Object
delProb(lp = 'optObj_clpAPI')

# checkOptSol-methods:Summarized Information About an Object of Class Optsol
# The function checkOptSol is used by functions performing a linear optimization (e.g. optimizeProb). In that case, the argument onlywarn is set to TRUE. 
# If the optimization ends unsuccesfull, a warning will be produced.

Ec_f <- optimizeProb(model_OTU08_draft, retOptSol = TRUE) 
Ec_check <- checkOptSol(Ec_f)
# Return code:
# Code    #       meaning
# 0       1       solution process was successful
# 
# Solution status:
#   Code    #       meaning
# 5       1       solution is optimal




checkAlgorithm('cplexAPI')
# [1] FALSE
checkAlgorithm('glpkAPI')
# [1] FALSE
checkAlgorithm('cplAPI')
# [1] FALSE

# checkReactId:The function checkReactId evaluates a vector of reaction id’s if they are unique and appear in a given model.
checkReactId(model_OTU08_draft,'rxn00433_Shan_OTU08_20210330')
# #   position   reaction id
# [1]  1415       rxn00433_Shan_OTU08_20210330
# 
# number of reactions: 1


checkReactId(model_OTU08_draft,1415)
#   position   reaction id
# [1]  1415       rxn00433_Shan_OTU08_20210330
# 
# number of reactions: 1




# checkVersion-methods:checks Version of modelorg
# Returns TRUE if the version is correct. Otherwise returns a character stating the reason.
checkVersion(model_OTU08_draft)
# [1] TRUE


# deadEndMetabolites-methods: Identify Dead End Metabolites
df_dead_end_met<- deadEndMetabolites(model_OTU08_draft)
class(df_dead_end_met)

df = as.data.frame(cbind(df_dead_end_met$dem, df_dead_end_met$der))
colnames(df) <- c("dem", "der")

write.csv(df, file = outfile_052, row.names = F)
# Can do some work to assign names to the met and rxn.


# doubleFluxDel:Double Flux Deletion Experiment
res<- doubleFluxDel(model_OTU08_draft)
# solver:                                   glpkAPI
# method:                                   simplex
# algorithm:                                fba
# number of variables:                      1659
# number of constraints:                    1606
# number of problems to solve:              1659
# number of successful solution processes:  1659

res<- doubleFluxDel(model_OTU08_draft, 'rxn16236_c0','rxn15972_c0')
# solver:                                   glpkAPI
# method:                                   simplex
# algorithm:                                fba
# number of variables:                      1659
# number of constraints:                    1606
# return value of solver:                   solution process was successful
# solution status:                          solution is optimal
# value of objective function (fba):        -0.000000
# value of objective function (model):      -0.000000

class(res)
# [1] "optsol_fluxdel"
# attr(,"package")
# [1] "sybil"



# doubleGeneDel: Double Gene Deletion Experiment
# Predict the metabolic phenotype of of double-gene knock out mutants.
res2<- doubleGeneDel(model_OTU08_draft)
# solver:                                   glpkAPI
# method:                                   simplex
# algorithm:                                fba
# number of variables:                      1659
# number of constraints:                    1606
# number of problems to solve:              1840
# number of successful solution processes:  1840


# doubleReact:Identiﬁes Identical Reactions
# If no identical reactions were found, the return value is FALSE. Otherwise a list is returned, ordered by the number of metabolites used in each reaction. 
# Each element is a numerical vector containing the indices (column number fo the stoichiometrix matrix) of identical reactions.
doubleReact(model_OTU08_draft)
# [[12]]
# [1] 1458 1500

# Check the identicals:
checkReactId(model_OTU08_draft,1458)
checkReactId(model_OTU08_draft,1500)




# editEnvir:Environment Editor for Metabolic Networks
mod <-editEnvir(model_OTU08_draft)
# Error in edit.data.frame(exfr, ...) : X11 dataentry cannot be loaded


# findExchReact: Find Exchange Reactions
df<- findExchReact(model_OTU08_draft)
class(df)
# [1] "reactId_Exch"
# attr(,"package")
# [1] "sybil"

df[2]

df = as.data.frame(unlist(df))
# Error in as.data.frame.default(df) : 
  # cannot coerce class ‘structure("reactId_Exch", package = "sybil")’ to a data.frame

# write.table(,outfile_053,row.names = F)

# Followings checks the flux of exchanges:
ex <- findExchReact(model_OTU08_draft)
# run FBA 
opt <- optimizeProb(model_OTU08_draft)
# get flux distribution of exchange reactions 
flux1<- getFluxDist(opt, ex)
class(flux1)
# [1] "numeric"
df=as.data.frame(flux1)
write.csv(df,outfile_053,row.names = T)
# Can be later assigned with exchange reaction names


# fluxDistribution-class:Class "fluxDistribution"
fluxDistribution(model_OTU08_draft)
showClass("fluxDistribution")
# Class "fluxDistribution" [package "sybil"]
# 
# Slots:
#   
#   Name:         fluxes num_of_fluxes
# Class:        Matrix       integer
showClass("checksol")
# Class "checksol" [package "sybil"]
# 
# Slots:
#   
#   Name:     num_of_prob      exit_code       exit_num   exit_meaning    status_code     status_num status_meaning
# Class:        integer        integer        integer      character        integer        integer      character

