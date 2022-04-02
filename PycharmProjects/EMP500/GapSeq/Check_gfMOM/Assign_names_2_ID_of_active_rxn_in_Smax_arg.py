#### IMPORTANT-2: this script used data in gapseq_v20200407, should be updated with the pipline version used for modeling.
#### update @2021-02-24: reaction direction in the output file is the same in file dat/seed_reactions_corrected.tsv
#### The true direction should combine with flux values (depending on either positive or negative)

# merge two Python dictionaries:
# https://stackoverflow.com/questions/38987/how-do-i-merge-two-dictionaries-in-a-single-expression-taking-union-of-dictiona
def merge_two_dicts(x, y):
    """Given two dictionaries, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z
# z = merge_two_dicts(x, y)
######################  run script in command line.
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-E.g',
                    required=False,
                    help='python3 /Users/zzfanyi/PycharmProjects/my_pys/EMP500/GapSeq/Get_active_rxn_by_met_from_final_model_argv(old).py -in /Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_Stylissa_20200407_MergedSuperGenomes_MetaC_kegg_seed_GapSeqv0405_Diet_TT_sw_lALL_b_archaea_uncon_phototest/OTU02_renamed/20220423_diet_Jo_photo_nitrate/ -file 050_model2_con_Smat_fluxes.R.csv -met cpd00001')

parser.add_argument('-in_dir',
                    required=True,
                    help='input bin folder')

parser.add_argument('-RP_rxn_file',
                    required=True,
                    default='/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210324_N.maritimusSCM1_STYOTU08_gapseq_v0324_MKSdb_Diet_sw_lALL_b_uncon/OTU08_renamed_arc/2021-03-26_ReactionPool_to_ObjectModel_use_bio1/02_Get_RP_equations/Reaction_id_in_ReactionPool_Addreact_20210330pm1.txt',
                    help='full path + file name of a customized MS reaction table (reaction pool) with MS id(col0) and MS name(col6).')

parser.add_argument('-gs',
                    required=True,
                    default='gapseq_v20210324',
                    help='gapseq version' )

parser.add_argument('-modNm',
                    required=True,
                    default='NmSCM1',
                    help='model_id' )

args = vars(parser.parse_args())

input_dir                   = args['in_dir']
dict_customized_rxn_file    = args['RP_rxn_file']
gapseq_version              = args['gs']
model_id                    = args['modNm']


# Step00. Define a function: find non-zero value list:
def rm_zero(num_list):

    num_list_no_zero = []
    for num in num_list:
        if num != '0':
            num_list_no_zero.append(num)

    return num_list_no_zero


##################### Start of running run script within pycharm #####################
# input_dir = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210409_STY8taxa_gapseq_v0409_MKSdb_Diet_sw_lALL_b_uncon/OTU08_renamed_arc/2021-04-10_ReactionPool_to_ObjectModel_use_bio1/07_adapt_addReact_R/04_Rscp/OTU08_MOM_20210411_addReact/20210411_diet_sw_no_photo_STY_Merged_OTU08_arc__12_oxygen/'
# dict_customized_rxn_file = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210409_STY8taxa_gapseq_v0409_MKSdb_Diet_sw_lALL_b_uncon/OTU08_renamed_arc/2021-04-10_ReactionPool_to_ObjectModel_use_bio1/02_Get_RP_equations/Reaction_id_in_ReactionPool_Addreact_20210411.txt'
# gapseq_version = 'gapseq_v20210324'
# model_id = 'NmSCM1'
##################### End of running script within pycharm #####################

# forward to working directory
os.chdir(input_dir)

# # # no need /
# print(input_dir)
if input_dir[-1] != "/":
    input_dir += "/"
# print(input_dir)


# Dict input files:
dict_seed_reactions_corrected_tsv = '/Users/zzfanyi/Bioinfo/software/gapseq/'+ gapseq_version +'/dat/seed_reactions_corrected.tsv'
dict_seed_met_edited        = '/Users/zzfanyi/Bioinfo/software/gapseq/' + gapseq_version + '/dat/seed_metabolites_edited.tsv'
dict_seed_met_edited_custom = '/Users/zzfanyi/Bioinfo/software/gapseq/' + gapseq_version + '/Shan_output_files/seed_metabolites_edited_Shan.tsv'
# Input and output files:
infile_Smat_flux            = input_dir + '050_model_' + model_id + '_gfMOM_Smat_fluxes.csv'
outfile_flux_lt0            = input_dir + '051_model_' + model_id + '_gfMOM_Smat_fluxes.flux_lt0.csv'
outfile_flux_lt0_annot      = input_dir + '052_model_' + model_id + '_gfMOM_Smat_fluxes.flux_lt0.annotate.tsv'


# step01. get a table with reactions' fluxes NOT = 0.
n = 0
output_handle = open(outfile_flux_lt0, 'w')
for each in open(infile_Smat_flux):
    if n == 0:
        each_split = each.strip().split(',')
        # print(each_split)
        n =+1
        header = each
        output_handle.write(header)
    else:
        rxn_ID = each.strip().split(',')[0]
        flux_index = float(each.strip().split(',')[1])
        if flux_index != 0:
            # print(each)
            # print(flux_index)
            output_handle.write(each)
output_handle.close()


# step2. make a dict of met id to met names:
dict_metID_2_nm = {}
for each in open(dict_seed_met_edited):
    each_split = each.strip().split('\t')
    if not each.startswith('id'):
        met_id = each_split[0]
        met_name = each_split[3]
        dict_metID_2_nm[met_id] = met_name

# step2.1 make a dict of (customized) met id to met names:
dict_metID_2_nm_custom = {}
for each in open(dict_seed_met_edited_custom):
    each_split = each.strip().split('\t')
    if not each.startswith('id'):
        met_id = each_split[0]
        met_name = each_split[3]
        dict_metID_2_nm_custom[met_id] = met_name

# Given two dictionaries, merge them into a new dict as a shallow copy.
# https://stackoverflow.com/questions/38987/how-do-i-merge-two-dictionaries-in-a-single-expression-taking-union-of-dictiona
dict_metID_2_nm = merge_two_dicts(dict_metID_2_nm, dict_metID_2_nm_custom)

# step3. make a dict of rxn id to rxn names:
dict_rxnID_2_nm = {}
for each in open(dict_seed_reactions_corrected_tsv):
    each_split = each.strip().split('\t')
    if not each.startswith('id'):
        rxn_id = each_split[0]
        rxn_name = each_split[2]
        dict_rxnID_2_nm[rxn_id] = rxn_name

# step4. make a dict of rxn id to rxn names from customized file:
dict_rxnID_2_nm_SZ = {}
for each in open(dict_customized_rxn_file):
    each_split = each.strip().split('\t')
    if not each.startswith('Rid'):
        rxn_id = each_split[0]
        rxn_name = each_split[6]
        dict_rxnID_2_nm_SZ[rxn_id] = rxn_name
# print(dict_rxnID_2_nm_SZ['R00132_PWY.5789_1_Shan_OTU08'])
# print(dict_rxnID_2_nm_SZ)

#### Merge reaction dictionaries
dict_rxn_merge_all = {**dict_rxnID_2_nm, **dict_rxnID_2_nm_SZ}


#step04. assign names to ids:
n=0
list_compartment = []
list_met_nm = []
output_handle = open(outfile_flux_lt0_annot, 'w')
for each in open(outfile_flux_lt0):
    each_split = each.strip().split(',')
    if n==0:
        # print(each_split)
        list_met = each_split[2:]
        for i in list_met:
            met_id = i.strip().split('[')[0]
            met_comp = i.strip().split('[')[1]
            # print(met_comp)
            list_compartment.append(met_comp)
            met_name = '%s%s%s' % (dict_metID_2_nm[met_id], '[', met_comp)
            list_met_nm.append(met_name)

        # print(list_met)
        # print(list_met_nm)
        # print(list_compartment)
        for_print = '%s\t%s\t%s\t%s\n' % ('', '', 'sol$fluxes', '\t'.join(list_met))
        output_handle.write(for_print)
        for_print = '%s\t%s\t%s\t%s\n' % ('', '', 'met_name', '\t'.join(list_met_nm))
        output_handle.write(for_print)
    else:
        rxn_name = 'NA'
        if each_split[0].startswith('rxn') and (each_split[0] not in dict_rxnID_2_nm_SZ):
            rxn_ID = each_split[0].split('_')[0]
            if rxn_ID in dict_rxnID_2_nm:
                rxn_name = dict_rxnID_2_nm[rxn_ID]
                for_print = '%s\t%s\t%s\n' % (each_split[0], rxn_name, '\t'.join(each_split[1:]))
                output_handle.write(for_print)
                # print(for_print)
            else:
                rxn_name = 'Not found in MS id dict.'
                for_print = '%s\t%s\t%s\n' % (each_split[0], rxn_name, '\t'.join(each_split[1:]))
                output_handle.write(for_print)
                # print(for_print)
        if each_split[0].startswith('EX'):
            rxn_ID = each_split[0].split('_')[1]
            if rxn_ID in dict_metID_2_nm:
                rxn_name = dict_metID_2_nm[rxn_ID]
                for_print = '%s\t%s%s\t%s\n' % (each_split[0], 'Exchange of ',rxn_name, '\t'.join(each_split[1:]))
                output_handle.write(for_print)
                # print(for_print)
            else:
                rxn_name = 'Not found in metabolite dict.'
                for_print = '%s\t%s\t%s\n' % (each_split[0], rxn_name, '\t'.join(each_split[1:]))
                output_handle.write(for_print)
                # print(for_print)

        if each_split[0] in dict_rxnID_2_nm_SZ:
            rxn_ID = each_split[0]
            rxn_name = dict_rxnID_2_nm_SZ[rxn_ID]
            for_print = '%s\t%s\t%s\n' % (each_split[0], rxn_name, '\t'.join(each_split[1:]))
            output_handle.write(for_print)
            # print(for_print)

    n += 1
output_handle.close()
print('Done!')