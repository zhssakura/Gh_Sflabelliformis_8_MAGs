import os

gapseq_version = 'gapseq_v20210409'
gapseq_version = 'gapseq_v20210504'

dir_dict_seed_reactions_corrected_tsv   = '/Users/zzfanyi/Bioinfo/software/gapseq/' + gapseq_version + '/'
dir_shan_in_gs                          = dir_dict_seed_reactions_corrected_tsv + 'Shan_output_files/'
seed_reactions_corrected_tsv            = dir_dict_seed_reactions_corrected_tsv + 'dat/seed_reactions_corrected.tsv'
output_file                             = dir_shan_in_gs + 'seed_reactions_corrected.formatted_reat_add.tsv'

# Creat a folder if not present:
if not os.path.isdir(dir_shan_in_gs):
   os.makedirs(dir_shan_in_gs)

# forward to working directory
os.chdir(dir_shan_in_gs)

# # # no need /
# print(input_dir)
if dir_shan_in_gs[-1] != "/":
    dir_shan_in_gs += "/"
# print(input_dir)


output_file_handle = open(output_file, 'w')
output_file_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
                         ('Rid','met_id','met_name','met_Scoef','met_comp','react_id','react_name','react_rev','lowbnd','uppbnd'))
for each_line in open(seed_reactions_corrected_tsv):
    if not each_line.startswith('id'):
        each_line_split     = each_line.strip().split('\t')
        rxn_id              = '%s_%s' % (each_line_split[0], 'Shan')
        react_id            = each_line_split[0]
        react_name          = each_line_split[2]
        stoichiometry       = each_line_split[4]
        reversibility       = each_line_split[8]
        stoichiometry_split = stoichiometry.split(';')

        # get react_rev, lowbnd and uppbnd
        react_rev = ''
        lowbnd = ''
        uppbnd = ''
        if reversibility == '=':
            react_rev = 'TRUE'
            lowbnd = '-1000'
            uppbnd = '1000'
        if reversibility == '>':
            react_rev = 'FALSE'
            lowbnd = '0'
            uppbnd = '1000'
        if reversibility == '<':
            react_rev = 'TRUE'
            lowbnd = '-1000'
            uppbnd = '0'


        met_id_list = []
        met_name_list = []
        met_scoef_list = []
        met_comp_list = []
        for each_cpd in stoichiometry_split:
            each_cpd_split = each_cpd.split(':')

            # get met_id_list and met_name_list
            # met_comp = 1 --> cpd00106[c0]
            # met_comp = 2 --> cpd00106[e0]
            # met_comp = 3 --> cpd00106[p0]
            # refer to file 045_df_model_OTU06_draft_list_metabolites.tsv and DictFile_Rids_adapt_add_Li2018SCM1.tsv

            if each_cpd_split[2] == '0':
                met_id_list.append('%s[c0]' % each_cpd_split[1])
                met_name_list.append('%s[c0]' % each_cpd_split[4][1:-1])
            elif each_cpd_split[2] == '1':
                met_id_list.append('%s[e0]' % each_cpd_split[1])
                met_name_list.append('%s[e0]' % each_cpd_split[4][1:-1])
            elif each_cpd_split[2] == '2':
                met_id_list.append('%s[p0]' % each_cpd_split[1])
                met_name_list.append('%s[p0]' % each_cpd_split[4][1:-1])
            else:
                print(stoichiometry)

            # get met_scoef_list list
            met_scoef_list.append(each_cpd_split[0])

            # get met_comp_list list
            met_comp_list.append((int(each_cpd_split[2]) + 1))

        met_comp_list_str = [str(i) for i in met_comp_list]

        for_write = '%s\t"%s"\t"%s"\t"%s"\t"%s"\t%s\t%s\t%s\t%s\t%s\n' % (rxn_id,
                                                                        ','.join(met_id_list),
                                                                        ','.join(met_name_list),
                                                                        ','.join(met_scoef_list),
                                                                        ','.join(met_comp_list_str),
                                                                        react_id,
                                                                        react_name,
                                                                        react_rev,
                                                                        lowbnd,
                                                                        uppbnd)
        output_file_handle.write(for_write)

output_file_handle.close()