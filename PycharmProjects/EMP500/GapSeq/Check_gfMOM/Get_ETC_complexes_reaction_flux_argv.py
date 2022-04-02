#### IMPORTANT-2: this script used data in gapseq_v20200407, should be updated with the pipline version used for modeling.
#### update @2021-02-24: reaction direction in the output file is the same in file dat/seed_reactions_corrected.tsv
#### The true direction should combine with flux values (depending on either positive or negative)
#### update @2021-07-01:
#### Give results of flux values of ETC chain complexes reactions
##### Input key file: /Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/EC_num_of_ETC_complexes.txt
###### could be updates if reactions were missing, or EC num is wrong.
######################  run script in command line.
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-E.g',
                    required=False,
                    help='python3 /Users/zzfanyi/PycharmProjects/my_pys/EMP500/GapSeq/Get_ETC_complexes_reaction_flux_argv.py -in_dict /Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/OTU07_renamed/Step_protocals_using_R/07_gf_model/20210630_diet_sw_no_oxygen_STY_Merged_OTU07-draft__12_thiosulfate/050_model_OTU07_gfMOM_Smat_fluxes.csv -in_key /Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/EC_num_of_ETC_complexes.txt -gs gapseq_v20210409 -out /Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/OTU07_renamed/Step_protocals_using_R/07_gf_model/20210630_diet_sw_no_oxygen_STY_Merged_OTU07-draft__12_thiosulfate/050_model_OTU07_gfMOM_Smat_fluxes_ETC_complexes.tsv')

parser.add_argument('-in_dict',
                    required=True,
                    help='/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/OTU07_renamed/Step_protocals_using_R/07_gf_model/20210630_diet_sw_no_oxygen_STY_Merged_OTU07-draft__12_thiosulfate/050_model_OTU07_gfMOM_Smat_fluxes.csv')

parser.add_argument('-in_key',
                    required=True,
                    default='/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/EC_num_of_ETC_complexes.txt',
                    help='Full path + file name of reactions added to MOM.')

parser.add_argument('-gs',
                    required=True,
                    default='gapseq_v20210409',
                    help='gapseq version')

parser.add_argument('-out',
                    required=True,
                    default='/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/OTU07_renamed/Step_protocals_using_R/07_gf_model/20210630_diet_sw_no_oxygen_STY_Merged_OTU07-draft__12_thiosulfate/050_model_OTU07_gfMOM_Smat_fluxes_ETC_complexes.tsv',
                    help='model_id' )

args = vars(parser.parse_args())
Input_dictfile  = args['in_dict']
Input_keyfile   = args['in_key']
gs_version      = args['gs']
Outfile         = args['out']


##################### Start of running run script within pycharm #####################
# Input_dictfile = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/OTU07_renamed/Step_protocals_using_R/07_gf_model/20210630_diet_sw_no_oxygen_STY_Merged_OTU07-draft__12_thiosulfate/050_model_OTU07_gfMOM_Smat_fluxes.csv'
# Input_keyfile = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/EC_num_of_ETC_complexes.txt'
# Outfile = '%s/%s' % ('/'.join(Input_dictfile.strip().split('/')[0:-1]), '050_model_OTU07_gfMOM_Smat_fluxes_ETC_complexes.tsv')
# # print(Outfile)
# gs_version = 'gapseq_v20210409'
##################### End of running script within pycharm #####################
# st1. create a dict:
dict_reactID_2_reactFlux = {}
for each in open(Input_dictfile):
    each_split = each.strip().split(',')
    # print(each_split)
    if not each.startswith(','):
        react_id = '_'.join(each_split[0].split('_')[0:-1])
        react_flux = each_split[1]
        # print(react_id,react_flux)
        dict_reactID_2_reactFlux[react_id] = react_flux


# st2. get output:
n = 0
output_handle = open(Outfile,'w')
for each in open(Input_keyfile):
    each_split = each.strip().split('\t')
    for_print = '%s\t%s\t%s\n' % ('react_id', 'react_flux', '\t'.join(each_split))
    if n ==0:
        output_handle.write(for_print)
        n+=1

    else:
        # print(each_split)
        if not each.startswith('Complex_id'):
            react_id = each_split[1]
            react_flux = 'NA'
            if react_id in dict_reactID_2_reactFlux:
                react_flux = dict_reactID_2_reactFlux[react_id]
                # print(react_id,react_flux)
            for_print = '%s\t%s\t%s\n' % (react_id,react_flux,'\t'.join(each_split))
            output_handle.write(for_print)
output_handle.close()
print('Done! Check result file: 050_model_OTU07_gfMOM_Smat_fluxes_ETC_complexes.tsv')