#!/bin/bash

# variables:
# -> species = species
# -> genome = specific reference genome
# -> experiment = experiment name (name of expdata directory)
# -> sublibrary = sublibrary of the experiment being analysed (or comb)
# -> parse_analysis_subdir = all-well/DGE_{filtered,unfiltered}/
# -> do_integration = TRUE or FALSE

# set directories, activate conda env
# base_dir=/Volumes/TracerX/working/VHL_GERMLINE/tidda/ ; wd=$base_dir/vhl/ ; cd $wd
cd /camp/project/tracerX/working/VHL_GERMLINE/tidda/vhl/
. src/config.sh

# get parse_pipeline runs to analyse
Rscript src/get_runs.R

# run generate_qc_report
cat out/runs.tsv | sed 1d | grep FALSE | grep -P '230210|221202' |
grep -P 'comb|SHE5052A9_S101' | grep -P '\thg38\t' | grep DGE_unfiltered |
{
  while read path experiment genome sublibrary parse_analysis_subdir do_integration sample_subset mem out_dir ; do
    . src/submit_generate_qc_report.sh
  done
}

# pilot run
cat out/runs.tsv | sed 1d |
grep -P '\t221202_A01366_0326_AHHTTWDMXY\thg38\t' |
grep DGE_filtered | grep FALSE |
{
  while read path experiment genome sublibrary parse_analysis_subdir do_integration sample_subset mem out_dir dge_mtx_dir ; do
    . src/submit_generate_qc_report.sh
  done
}

# full final run
cat out/runs.tsv | sed 1d |
grep -P '230210_A01366_0351_AHNHCFDSX5,221202_A01366_0326_AHHTTWDMXY' |
grep DGE_filtered | grep FALSE |
{
  while read path experiment genome sublibrary parse_analysis_subdir do_integration sample_subset mem out_dir dge_mtx_dir ; do
    . src/submit_generate_qc_report.sh
  done
}

# Daqi's PDOs
cat out/runs.tsv | sed 1d |
grep "PDOs" |
{
  while read path experiment genome sublibrary parse_analysis_subdir do_integration sample_subset mem out_dir dge_mtx_dir ; do
    . src/submit_generate_qc_report.sh
  done
}

