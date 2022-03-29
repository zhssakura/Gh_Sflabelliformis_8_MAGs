# Get formatted reaction equations and fluxes  by providing MS ID (e.g.rxnxxxx_c0): 
#### IMPORTANT: this script can be applied to any input file with ModelSEED id, e.g. rxnxxxx_c0.
#### IMPORTANT: ModelSEED reaction id must be in the 1st column with header start with 'react_id'.

# python3 /Users/zzfanyi/PycharmProjects/EMP500/GapSeq/Get_formatted_equations_n_flux_of_MS_react_in_gfModel/Get_formatted_equations_n_flux_of_MS_react_in_gfModel.py -infile_MSid -outfile -infile_flux -gs -db
dir_out = getwd()

gapseq_version      = 'gapseq_v20210504'
python_code         = '/Users/zzfanyi/PycharmProjects/EMP500/GapSeq/Get_formatted_equations_n_flux_of_MS_react_in_gfModel/Get_formatted_equations_n_flux_of_MS_react_in_gfModel.py'
dict_seed_reactions_corrected_tsv   = paste('/Users/zzfanyi/Bioinfo/software/gapseq/',gapseq_version,'/Shan_output_files/seed_reactions_corrected.formatted_reat_add.tsv',sep = '')

#Change here:
infile_flux = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/OTU01_renamed/2021-07-01_ReactionPool_to_ObjectModel_use_bio1/07_adapt_addReact_R/04_Rscp/OTU01_MOM_20210707/20210708_diet_sw_no_oxygen_STY_Merged_OTU01_MOM__12_taurine/050_model_OTU01_gfMOM_Smat_fluxes.csv'

infile_key_list                     = paste(dir_out,'/20210708_diet_sw_no_oxygen_STY_Merged_OTU01_MOM__12_taurine_completely_filled_in_reactions.txt',sep = '')
outfile_key_list_equation_formatted = paste(dir_out,'/20210708_diet_sw_no_oxygen_STY_Merged_OTU01_MOM__12_taurine_completely_filled_in_reactions.equation_formatted_n_fluxes.txt',sep = '')

commands = paste('python3',python_code,'-infile_MSid',infile_key_list,'-infile_flux', infile_flux, '-db',dict_seed_reactions_corrected_tsv,'-outfile',outfile_key_list_equation_formatted,'-gs',gapseq_version, sep = ' ')
system(commands)

