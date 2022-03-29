
#st1. get Model SEED reaction ids by providing keys (e.g. gene name, metacyc id, kegg id and EC number)
getwd = getwd()
gapseq_version      = 'v20210504'
dir_gapseq_files    = paste('/Users/zzfanyi/Bioinfo/software/gapseq/gapseq_',gapseq_version,'/dat/',sep = '')

dir_python_code     = '/Users/zzfanyi/PycharmProjects/EMP500/GapSeq/Get_reaction_equation_by_keys/'
python_code         = paste(dir_python_code, 'st0.Get_SEED_rxnID_by_keys_Willis_arg.py',sep = '')

dir_in              = paste(getwd,'/10_get_pwy_react/',sep = '')
infile_key_list     = paste(dir_in,'Key_to_get_MS_ID.txt',sep = '')

dir_out             = paste(getwd,'/10_get_pwy_react/Key_to_get_MS_ID_out/',sep = '')


commands = paste('python3 ',python_code,' -in ',infile_key_list,' -db ',dir_gapseq_files,' -out ',dir_out, sep = '')
system(commands)



# st2.1 Find if the MS ids in the last step can be found in Object Model (OM)_draft model! 
# First, manually make file 'op_merge_all_Shan.txt'.

python_code         = paste(dir_python_code, 'st2.Find_if_SEED_rxnID_in_ObjectModel_arg.py',sep = '')
Taxa_ID = 'STY_Merged_OTU07'

infile_key_list     = paste(dir_out,'op_merge_all_Shan.txt',sep = '')


dir_infile_ObjectModel = paste(getwd,'/03_Get_ObjectModel/',sep = '')
infile_ObjectModel = paste(dir_infile_ObjectModel,'043_df_model_',Taxa_ID,'_draft_list_react.tsv',sep = '')

outfile_key_list_in_Taxa_ID = paste(dir_out,'op_merge_all_Shan.in_',Taxa_ID,'_draft.txt',sep = '')

commands = paste('python3 ',python_code,' -in ',infile_key_list,' -taxa ',Taxa_ID,' -om ',infile_ObjectModel,' -outfile ',outfile_key_list_in_Taxa_ID, sep = '')
system(commands)
# st1. There are 2265 MS IDs found in STY_Merged_OTU07


# st2.2 Find if the MS ids in the last step can be found in Object Model (OM)_final model! 
# First, manually make file 'op_merge_all_Shan.txt'.

dir_final_model = '20210630_diet_sw_no_oxygen_STY_Merged_OTU07-draft__12_H2S'
# Then create a new directory: 
dir.create(paste(dir_out,'/',dir_final_model,sep = ''))

python_code         = paste(dir_python_code, 'st2.Find_if_SEED_rxnID_in_ObjectModel_arg.py',sep = '')
Taxa_ID = 'STY_Merged_OTU07'

infile_key_list     = paste(dir_out,'op_merge_all_Shan.txt',sep = '')


dir_infile_ObjectModel = paste(getwd,'/03_Get_ObjectModel/',dir_final_model,sep = '')
infile_ObjectModel     = paste(dir_infile_ObjectModel,'/043_df_model_',Taxa_ID,'_final_list_react.tsv',sep = '')

outfile_key_list_in_Taxa_ID = paste(dir_out,dir_final_model,'/op_merge_all_Shan.in_',Taxa_ID,'_final.txt',sep = '')

commands = paste('python3 ',python_code,' -in ',infile_key_list,' -taxa ',Taxa_ID,' -om ',infile_ObjectModel,' -outfile ',outfile_key_list_in_Taxa_ID, sep = '')
system(commands)
# st1. There are 2423 MS IDs found in STY_Merged_OTU07

################################################################################################################################################

# st3. Get reaction equations by MS ID:. 
python_code         = paste(dir_python_code, 'st3.Get_react_equations_by_SEED_rxnID_arg.py',sep = '')


dict_seed_reactions_corrected_tsv   = paste('/Users/zzfanyi/Bioinfo/software/gapseq/gapseq_',gapseq_version,'/dat/seed_reactions_corrected.tsv',sep = '')
infile_key_list                     = paste(dir_out,'op_merge_all_Shan.txt',sep = '')
outfile_key_list_equation           = paste(dir_out,'op_merge_all_Shan.equation.txt',sep = '')

commands = paste('python3 ',python_code,' -in ',infile_key_list,' -db ',dict_seed_reactions_corrected_tsv,' -outfile ',outfile_key_list_equation, sep = '')
system(commands)
# st1. There are 34822 MS IDs with associated equations found in gapseq v20210318



# st5. Get formatted reaction equations by MS ID:. 

python_code         = paste(dir_python_code, 'st5.Get_formatted_seed_equation_arg.py',sep = '')
dict_seed_reactions_corrected_tsv   = paste('/Users/zzfanyi/Bioinfo/software/gapseq/gapseq_',gapseq_version,'/Shan_output_files/seed_reactions_corrected.formatted_reat_add.tsv',sep = '')


infile_key_list                     = paste(dir_out,'op_merge_all_Shan.txt',sep = '')
outfile_key_list_equation_formatted = paste(dir_out,'op_merge_all_Shan.equation_formatted.txt',sep = '')

commands = paste('python3 ',python_code,' -in ',infile_key_list,' -db ',dict_seed_reactions_corrected_tsv,' -outfile ',outfile_key_list_equation_formatted, sep = '')
system(commands)

