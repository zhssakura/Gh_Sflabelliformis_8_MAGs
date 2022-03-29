#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l mem=20gb
#PBS -l walltime=11:59:00
#PBS -M zzfanyi@gmail.com
#PBS -m ae

########################### Test: gapseq adapt add  #############################
# @2021-apr-09 suggestion

module load perl/5.28.0
module load git/2.22.0
module load bedtools/2.27.1
module load blast+/2.9.0
module load hmmer/3.2.1
module load glpk/4.65
module load barrnap/0.9
module load gcc/7.3.0
module load exonerate/2.2.0
module load parallel/20190522
module add R/3.6.1
cd /srv/scratch/z5095298/sponge_modeling/GapSeq_v20210409_STY8taxa_gapseq_v0409_MKSdb_Diet_sw_lALL_b_uncon/OTU08_renamed_arc


keyword=PWY-5789,AMMOXID-PWY

./gapseq.sh adapt add ${keyword} ./STY_Merged_OTU08-draft.RDS > ./STY_Merged_OTU08_draft__adapt_add_${keyword}.txt