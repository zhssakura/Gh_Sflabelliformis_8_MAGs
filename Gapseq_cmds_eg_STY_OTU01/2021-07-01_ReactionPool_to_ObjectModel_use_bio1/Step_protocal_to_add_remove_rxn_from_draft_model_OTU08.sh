# St0. Assess reactions in OM, and select [Pathways] that you want to compare. 
# Run 'st1.Get_reaction_status_by_pathwayID_in_gs_report_Willis.py' to assess reaction status in Object model. 
# Input file: qsub_GapSeq_v20201228_OTU08_renamed_diet_Jo_arc.sh.o925287
# Output file: Cellular_pathway_input_Shan_20210128_MetaCyc.tsv
# Manually prepare 'pathways_EC.tsv', which includes your interested pathway IDs (metacyc).


# St1. Extract reactions in RP.
# Run 01_Rscript_compare_RDS_files_SCM1_model_Li_2018.R
# Reformat '043_df_model_SCM1_Li2018_list_react.tsv' to manually get 043_df_model_SCM1_Li2018_list_react_ECcorrect.csv.

# St2. Compare reactions of interested [Pathways] from OM and RP.
# Run 'st2.Get_rxn_EC_via_pathwayID.py' to get pathways_EC_out2.tsv, which includes reaction status from OM and reaction IDs from RP.


# St3. Get equation of reactions from RP that you want to add/remove from OM.
# Manually prepare file 'Reaction_id_in_ReactionPool.txt', which include [reactions] of your interested [pathways] from RP.
## The standard of choose reaction ids to 'Reaction_id_in_ReactionPool.txt':
# 1. The reaction belongs to or related to an interested pathway;
# 2. The reaction can be found in RP but not in OM.

# Run '02_Get_dict_from_Smax_Willis.py' to get 'DictFile_Rids_adapt_add_Li2018SCM1.tsv'  
# Run '06.1_Get_target_react_eqtions.py' to get 'Reaction_id_in_ReactionPool_Addreact.txt' (can rename the output file with date created).


# St4. Extract reactions in OM.
# Run 02_Rscript_compare_RDS_files_OTU08-draft.R

# St5. Remove (with care!!) reactions from OM.
# Mannually get 'Things_need_remove.tsv'
# Login on Katana.
# Follow the comonds in '03_Cmds_atapt_remove_bio1.sh'.
## (Be cautious: if the entire PATHWAY was removed from the model, that means all relavent reactions will be removed.
## Many reactions are commonly present in different pathways.)
# If Modified OM (MOM) was created, download it to ./06_adapt_rmReact_R folder.

# st6. Add reactions (by equations) to OM (using FOR-loop in R)
# Run '04_Rscp_react_add_to_models_OTU08-adapt__add_rxn_Li2018SCM1_FORloop.R' to get MOM.
# This will:
# st6.1. Add reactions to a gs model with FOR-loop
# st6.2. Use For-loop to add all RP to OM
# st6.3. Sets/changes the Objective Function to OM
# st6.4. Extract all reactions from the MOM.

# OR: run '03_Cmds_atapt_remove__adapt_add.sh'


# st7.Check flux in MOM
# Run '05_Rscp_check_flux_in_MOM.R'

# st8. Summarise the names of metacyc pathway in for your genome only based on blast result (gapseq -find).
# Run 'Way1_Get_the_MetaCyc_PWY.py'

# st9. Get ModelSEED reaction IDs by providing keys (KO id, EC number, Metacyc RXN id or gene name)
# @2021-03-20.
# This py script can be used to get Model SEED reaction ids by providing keys
# (e.g. gene name, metacyc id, kegg id and EC number)
# Following with 'st2.Find_if_SEED_rxnID_in_ObjectModel_arg.py'.
# /Users/zzfanyi/PycharmProjects/EMP500/GapSeq/Get_reaction_equation_by_keys/st0.Get_SEED_rxnID_by_keys_Willis_arg.py
