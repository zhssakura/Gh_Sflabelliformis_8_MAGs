#### IMPORTANT-2: this script used data in gapseq_v20200407, should be updated with the pipline version used for modeling.
#### update @2021-02-24: reaction direction in the output file is the same in file dat/seed_reactions_corrected.tsv
#### The true direction should combine with flux values (depending on either positive or negative)
#### update @2021-08-04: add new function to add reactants and products based on flux value and met_coef.
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

# parser.add_argument('-file',
#                     required=True,
#                     help='Smax flux csv file of a model')

parser.add_argument('-met',
                    required=True,
                    default='cpd00001',
                    help='target metabolite, default: fatsa ')

parser.add_argument('-dict',
                    required=True,
                    default='/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210324_N.maritimusSCM1_STYOTU08_gapseq_v0324_MKSdb_Diet_sw_lALL_b_uncon/ncbi-genomes-2021-03-23/2021-03-25_ReactionPool_to_ObjectModel_use_bio1/02_Get_RP_equations/Reaction_id_in_ReactionPool_Addreact_20210411.txt',
                    help='Full path + file name of reactions added to MOM.')

parser.add_argument('-gs',
                    required=True,
                    default='gapseq_v20210409',
                    help='gapseq version')

parser.add_argument('-modNm',
                    required=True,
                    default='NmSCM1',
                    help='model_id' )

args = vars(parser.parse_args())
input_dir       = args['in_dir']
# input_file      = args['file']
target_met      = args['met']
dict_cus_react  = args['dict']
gs_version      = args['gs']
model_id        = args['modNm']

#
#
#
#
# Step00. Define a function: find non-zero value list:
def rm_zero(num_list):

    num_list_no_zero = []
    for num in num_list:
        if num != '0':
            num_list_no_zero.append(num)

    return num_list_no_zero
# e.g.
# print(rm_zero(['1','0']))

##################### Start of running run script within pycharm #####################
# input_dir = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210324_N.maritimusSCM1_STYOTU08_gapseq_v0324_MKSdb_Diet_sw_lALL_b_uncon/ncbi-genomes-2021-03-23/2021-03-25_ReactionPool_to_ObjectModel_use_bio1/07_adapt_addReact_R/04_Rscp/NmSCM1_MOM_20210411/20210412_diet_sw_no_photo_NmSCM1_MOM__12_oxygen/'
# dict_cus_react = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210324_N.maritimusSCM1_STYOTU08_gapseq_v0324_MKSdb_Diet_sw_lALL_b_uncon/ncbi-genomes-2021-03-23/2021-03-25_ReactionPool_to_ObjectModel_use_bio1/02_Get_RP_equations/Reaction_id_in_ReactionPool_Addreact_20210411.txt'
# target_met = 'cpd00011'
# gs_version = 'gapseq_v20210409'
# model_id = 'NmSCM1'
##################### End of running script within pycharm #####################

# forward to working directory
os.chdir(input_dir)

# # # no need /
# print(input_dir)
if input_dir[-1] != "/":
    input_dir += "/"
# print(input_dir)

# Dict file:
dict_seed_reactions_corrected_tsv = '/Users/zzfanyi/Bioinfo/software/gapseq/' + gs_version + '/dat/seed_reactions_corrected.tsv'
# Formatted dict:
# dict_seed_reactions_corrected_tsv_format = '/Users/zzfanyi/Bioinfo/software/gapseq/' + gs_version + '/Shan_output_files/seed_reactions_corrected.formatted_reat_add.tsv'
dict_seed_reactions_corrected_tsv_format   = input_dir + 'seed_reactions_corrected.formatted_reat_add.status.tsv'

# Input/output files:
input_file_Smat_flux            = input_dir + '052_model_' + model_id + '_gfMOM_Smat_fluxes.flux_lt0.annotate.tsv'
output_file_target_met          = input_dir + '053_model_' + model_id + '_gfMOM_Smat_fluxes.flux_lt0.annotate.'+target_met+'.tsv'
output_file_target_rxn_equation = input_dir + '053_model_' + model_id + '_gfMOM_Smat_fluxes.flux_lt0.annotate.'+target_met+'.Equ.tsv'

# # step01. get a table with reactions' fluxes NOT = 0.

# step02. get reactions with target metabolites:
n = 1
list_cpd = []
index_list = []
index_list_cpdID = []
index_list_cpdName = []
output_handle = open(output_file_target_met,'w')
for each in open(input_file_Smat_flux):
    each_split = each.strip().split('\t')
    # For the 1st header row:
    if n == 1:
        index = 0
        for cpd in each_split:
            # print(cpd)
            if target_met in cpd:
                index_list.append(index)
                index_list_cpdID.append(cpd)
            index += 1
        for_print = ('\t\t%s\t%s\n' % (each_split[0], '\t'.join(index_list_cpdID)))
        output_handle.write(for_print)

    # then print out cpd names which associated with target metabolites:
    elif n == 2:
        # print(each_split)
        cpdName_list = []
        for indx in index_list:
            cpdName_list.append(each_split[indx])
        # print(index_list)
        # print(cpdName_list)
        for_print = ('\t\t%s\t%s\n' % (each_split[0], '\t'.join(cpdName_list)))
        # print(for_print)
        output_handle.write(for_print)


    elif n>=3:
    # else:
        rxn_id = each_split[0]
        rxn_name = each_split[1]
        fluxes = each_split[2]

        coefficient_list = []
        for indx in index_list:
            coefficient_list.append(each_split[indx+2]) # PS: here indx +2: from 3rd row, len(each_split) is differs compared to row 1 or 2.
        # use funtion to filter out redundant reactions:
        #print(coefficient_list)
        coefficient_list_no_zero = rm_zero(coefficient_list)
        #print(coefficient_list_no_zero)

        if coefficient_list_no_zero != []:

            # print('%s,%s,%s' % (rxn_id, fluxes, ','.join(coefficient_list)))
            for_print = ('%s\t%s\t%s\t%s\n' % (rxn_id, rxn_name, fluxes, '\t'.join(coefficient_list)))
            # print(for_print)
            output_handle.write(for_print)
    n += 1
    # print(len(each_split))

output_handle.close()


# step03new. create a seed reaction dictionary using Shan_output_files/seed_reactions_corrected.formatted_reat_add.tsv
dict_seed_rxn_id_2_description = {}
list_seed_rxn_id = []
for each in open(dict_seed_reactions_corrected_tsv_format):
    each_split = each.strip().split('\t')
    if each.startswith('Rid\t'):
        # header_1 = each_split[0:] # exclude the first element
        header_1 = each_split

    else:
        # print(each)
        rxn_id_dict = each_split[0].split('_Shan')[0]
        rxn_description_dict = each_split[0:] # exclude the first element
        # print(rxn_description_dict)
        # rxn_description_dict = dict_rxn_id_2_description['\t'.join(rxn_id_dict)]
        dict_seed_rxn_id_2_description[rxn_id_dict] = '\t'.join(rxn_description_dict)
        if rxn_id_dict not in list_seed_rxn_id:
            list_seed_rxn_id.append(rxn_id_dict)
print('%s %s %s %s' % ('There are', len(dict_seed_rxn_id_2_description), 'reactions in', gs_version))

# step03-2:create a seed reaction dictionary using costumized file (e.g.Reaction_id_in_ReactionPool_Addreact_20210411.txt)
dict_cus_rxn_id_2_description = {}
for each in open(dict_cus_react):
    # print(each)
    each_split = each.strip().split('\t')
    if not each.startswith('Rid\t'):
        rxn_id_dict = each_split[0]
        rxn_description_dict = each_split[0:] # exclude the first element
        dict_cus_rxn_id_2_description[rxn_id_dict] = '\t'.join(rxn_description_dict)
print('%s %s %s %s' % ('There are', len(dict_cus_rxn_id_2_description), 'reactions in', gs_version))

# step03-3:Merge dictionaries above:
# Two ** if two dicts in {}.
dict_merge_all = {**dict_cus_rxn_id_2_description, **dict_seed_rxn_id_2_description}
# print(len(dict_merge_all))
print('%s %s %s' % ('There are in total', len(dict_merge_all), 'reactions in merged dicts.'))

# step04. get reaction descriptions by providing rxn IDs:
num_rxn = 0
output_handle = open(output_file_target_rxn_equation,'w')
output_handle.write('%s\t%s\t%s\t%s\n' % ('target_rxn_IDs','target_rxn_name','flux','\t'.join(header_1)))
for each in open(output_file_target_met):
    each_split = each.strip().split('\t')
    if not each.startswith('\t\tsol$fluxes') and not each.startswith('\t\tmet_name') :
        target_rxn_id      = each_split[0]
        target_rxn_name    = each_split[1]
        target_rxn_id_flux = each_split[2]
        target_rxn_description = 'NA' # In case this variable is not found in dict_merge_all.
        num_rxn+=1

        if target_rxn_id.endswith('_c0') or target_rxn_id.endswith('_e0'):
            target_rxn_id_search = '_'.join(target_rxn_id.split('_')[0:-1])
            # print(target_rxn_id_search)
            if target_rxn_id_search in dict_merge_all:
                target_rxn_description = dict_merge_all[target_rxn_id_search]
                # print(target_rxn_id,target_rxn_description)
                for_print = ('%s\t%s\t%s\t%s\n' % (target_rxn_id,target_rxn_name,target_rxn_id_flux,target_rxn_description)) #target_rxn_description is not list.
                output_handle.write(for_print)
        else:
            target_rxn_id_search = target_rxn_id
            # print(target_rxn_id_search)
            target_rxn_description = dict_merge_all[target_rxn_id_search]
            for_print = ('%s\t%s\t%s\t%s\n' % (target_rxn_id,target_rxn_name,target_rxn_id_flux,target_rxn_description))
            output_handle.write(for_print)
output_handle.close()
print('Done! '+str(num_rxn)+' ACTIVE reactions that utilize '+target_met+' are found!')