# 2021-apr-15 get a file of formatted MS reactions with assigned reaction and pathway status for that genome/draft model.

import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-E.g',
                    required=False,
                    help='python3 /Users/zzfanyi/PycharmProjects/EMP500/GapSeq/Check_gfMOM/Get_a_file_of_formatted_MS_reactions_with_assigned_reaction_and_pathway_status.py -dir_gfMOM /Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/OTU06_renamed/Step_protocals_using_R/07_gf_model/ -gs gapseq_v20210504 -dict_pwy_tbl /Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/OTU06_renamed/STY_Merged_OTU06-all-Pathways.tbl -dict_043 043_model_OTU06_gfMOM_list_react.tsv')

parser.add_argument('-dir_gfMOM',
                    required=True,
                    help='Input directory of file_052(Smax_fluxes.flux_lt0.annotated.tsv).')

parser.add_argument('-gs',
                    required=True,
                    default='gapseq_v20210504',
                    help='gapseq version')

parser.add_argument('-dict_pwy_tbl',
                    required=True,
                    help='Full path + file name of "xxx-all-Pathways.tbl".')

parser.add_argument('-dict_043',
                    required=True,
                    default='043_model_OTU06_gfMOM_list_react.tsv',
                    help='File name of file_043 (043_model_xxx_gfMOM_list_react.tsv).')


args = vars(parser.parse_args())
dir_gfMOM           = args['dir_gfMOM']
gapseq_version      = args['gs']
dict_Pathway_tbl    = args['dict_pwy_tbl']
dict_043            = args['dict_043']

######################## User define starts: ########################
# dir_gfMOM           = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210324_N.maritimusSCM1_STYOTU08_gapseq_v0324_MKSdb_Diet_sw_lALL_b_uncon/ncbi-genomes-2021-03-23/2021-03-25_ReactionPool_to_ObjectModel_use_bio1/07_adapt_addReact_R/04_Rscp/NmSCM1_MOM_20210411/20210412_diet_sw_no_photo_NmSCM1_MOM__12_oxygen/'
# dir_gfMOM           = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210324_N.maritimusSCM1_STYOTU08_gapseq_v0324_MKSdb_Diet_sw_lALL_b_uncon/ncbi-genomes-2021-03-23/2021-03-25_ReactionPool_to_ObjectModel_use_bio1/07_adapt_addReact_R/04_Rscp/NmSCM1_MOM_20210411/20210417_diet_sw_no_photo_NmSCM1_MOM__12_oxygen_new_biomassArc_-q/'
# dir_gfMOM           = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210324_N.maritimusSCM1_STYOTU08_gapseq_v0324_MKSdb_Diet_sw_lALL_b_uncon/ncbi-genomes-2021-03-23/2021-03-25_ReactionPool_to_ObjectModel_use_bio1/07_adapt_addReact_R/04_Rscp/NmSCM1_MOM_20210418_bioArc_original/20210418_diet_sw_no_photo_NmSCM1_MOM__12_oxygen_bioArc_original/'
# dir_gfMOM           = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210324_N.maritimusSCM1_STYOTU08_gapseq_v0324_MKSdb_Diet_sw_lALL_b_uncon/ncbi-genomes-2021-03-23/2021-03-25_ReactionPool_to_ObjectModel_use_bio1/07_adapt_addReact_R/04_Rscp/NmSCM1_MOM_20210419pm1_bioArc_original/20210419pm1_diet_sw_no_photo_NmSCM1_MOM__12_oxygen_bioArc_original/'
# dir_gfMOM           = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210324_N.maritimusSCM1_STYOTU08_gapseq_v0324_MKSdb_Diet_sw_lALL_b_uncon/ncbi-genomes-2021-03-23/2021-03-25_ReactionPool_to_ObjectModel_use_bio1/07_adapt_addReact_R/04_Rscp/NmSCM1_MOM_20210418pm2_bioArc_original/20210418pm2_diet_sw_no_photo_NmSCM1_MOM__12_oxygen_bioArc_original/'
# dir_gfMOM           = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210324_N.maritimusSCM1_STYOTU08_gapseq_v0324_MKSdb_Diet_sw_lALL_b_uncon/ncbi-genomes-2021-03-23/2021-03-25_ReactionPool_to_ObjectModel_use_bio1/07_adapt_addReact_R/04_Rscp/NmSCM1_MOM_20210419pm1_3_bioArc_original/20210419pm1_3_diet_sw_no_photo_NmSCM1_MOM__12_oxygen_bioArc_original/'
# dir_gfMOM           = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210324_N.maritimusSCM1_STYOTU08_gapseq_v0324_MKSdb_Diet_sw_lALL_b_uncon/ncbi-genomes-2021-03-23/2021-03-25_ReactionPool_to_ObjectModel_use_bio1/07_adapt_addReact_R/04_Rscp/NmSCM1_MOM_20210419pm1_3_bioArc_original/20210419pm1_3_diet_sw_no_photo_NmSCM1_MOM__12_oxygen_bioArc_original/'
# dir_gfMOM           = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/OTU08_renamed_arc/Step_protocals_using_R/07_gf_model/'
# dir_gfMOM           = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/OTU06_renamed/Step_protocals_using_R/07_gf_model/'
# dir_gfMOM           = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon_real0429/OTU06_renamed/Step_protocals_using_R/07_gf_model/'
# dir_gfMOM           = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/OTU06_renamed/Step_protocals_using_R/07_gf_model/'
#
# gapseq_version      = 'gapseq_v20210409'
# gapseq_version      = 'gapseq_v20210504'
#
# dict_Pathway_tbl    = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210324_N.maritimusSCM1_STYOTU08_gapseq_v0324_MKSdb_Diet_sw_lALL_b_uncon/ncbi-genomes-2021-03-23/N.maritimusSCM1-all-Pathways.tbl'
# dict_Pathway_tbl    = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/OTU08_renamed_arc/STY_Merged_OTU08-all-Pathways.tbl'
# dict_Pathway_tbl    = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/OTU06_renamed/STY_Merged_OTU06-all-Pathways.tbl'
# dict_Pathway_tbl    = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/OTU06_renamed/STY_Merged_OTU06-all-Pathways.tbl'
#
# dict_043            = dir_gfMOM + '043_model_NmSCM1_gfMOM_list_react.tsv'
# dict_043            = dir_gfMOM + '043_model_OTU08_gfMOM_list_react.tsv'
# dict_043            = '043_model_OTU06_gfMOM_list_react.tsv'
######################## User define ends: ########################
dict_043            = dir_gfMOM + dict_043
outfile_for_gfMOM   = dir_gfMOM + 'seed_reactions_corrected.formatted_reat_add.status.tsv'
dict_gs_seed_rxn    = '/Users/zzfanyi/Bioinfo/software/gapseq/'+gapseq_version+'/Shan_output_files/seed_reactions_corrected.formatted_reat_add.tsv'


# st1. get a dict of pwy id to pwy des:
dict_pwyID_2_des = {}
for each in open(dict_Pathway_tbl):
    each_split = each.strip().split('\t')
    if each.startswith('ID\t'):
        header_2 = each_split
        # print(header_2)
    else:
        pwy_id = each_split[0]
        pwy_des = '\t'.join(each_split)
        dict_pwyID_2_des[pwy_id] = pwy_des
# print(len(dict_pwyID_2_des))

# st2. get a dict of rxn id to rxn des:
dict_rxnID_2_des = {}
for each in open(dict_043):
    each_split = each.strip().split('\t')
    if each.startswith('react_id'):
        header_1 = '%s\t%s\t%s\t%s\t%s\t%s\t%s' % (each_split[0], each_split[9],each_split[11],each_split[12],each_split[13], '\t'.join(each_split[15:]), '\t'.join(header_2))
        # print(header_1)
    else:
        # print(each)
        rxn_id_print    = each_split[0]
        rxn_id          = '_'.join(each_split[0].split('_')[:-1])
        gpr             = each_split[9]
        metacyc_rxn     = each_split[11]
        metacyc_rxn_nm  = each_split[12]
        EC_num          = each_split[13]
        all_the_rest    = each_split[15:]
        pathway         = each_split[23]
        pwy_des         = 'NA'


        if pathway in dict_pwyID_2_des:
            pwy_des = dict_pwyID_2_des[pathway]

        description = '%s\t%s\t%s\t%s\t%s\t%s\t%s' % (
            rxn_id_print, gpr, metacyc_rxn, metacyc_rxn_nm, EC_num, '\t'.join(all_the_rest), pwy_des)
        dict_rxnID_2_des[rxn_id] = description
        # print(pwy_des)
        # print(description)

print(len(dict_rxnID_2_des))

# st3. create a dict of rxn if to rxn des in 'seed_reactions_corrected.formatted_reat_add.tsv'
output_handle = open(outfile_for_gfMOM,'w')
for each in open(dict_gs_seed_rxn):
    each_split = each.strip().split('\t')
    if each.startswith('Rid\t'):
        header_3 = '%s\t%s\n' % ('\t'.join(each_split),header_1)
        output_handle.write(header_3)
    else:
        rxn_id = each_split[0].split('_')[0]
        description = 'NA'
        if rxn_id in dict_rxnID_2_des:
            description = dict_rxnID_2_des[rxn_id]
            for_print = ('%s\t%s\n' % ('\t'.join(each_split), description))
            output_handle.write(for_print)
        # print(for_print)
        # print(description)
output_handle.close()

