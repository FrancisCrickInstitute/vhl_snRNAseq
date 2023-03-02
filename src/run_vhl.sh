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

# run analyse_snRNAseq
cat out/runs.tsv | sed 1d | grep FALSE | grep -P '230210|221202' |
grep -P 'comb|SHE5052A9_S101' | grep -P '\thg38\t' | grep DGE_filtered |
{
  while read path experiment genome sublibrary parse_analysis_subdir do_integration ; do
    . src/submit_generate_qc_report.sh
  done
}

# # run manually:
# cat out/runs.tsv | sed 1d |
# {
#   while read path experiment genome sublibrary parse_analysis_subdir do_integration ; do
#     Rscript src/run_analyse_snRNAseq.R   \
#       $experiment   $genome   $sublibrary  $parse_analysis_subdir $do_integration
#   done
# }
