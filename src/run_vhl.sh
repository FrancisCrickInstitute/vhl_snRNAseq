#!/bin/bash

# variables:
# -> species = species
# -> genome = specific reference genome
# -> exp = experiment name (name of expdata directory)
# -> sublibrary = sublibrary of the experiment being analysed (or comb)
# -> parse_analysis_subdir = all-well/DGE_{filtered,unfiltered}/
# -> do_integration = TRUE or FALSE

# set directories, load R
# base_dir=/Volumes/TracerX/working/VHL_GERMLINE/tidda/ ; wd=$base_dir/vhl/
cd /camp/project/tracerX/working/VHL_GERMLINE/tidda/vhl/
. src/config.sh

# get parse_pipeline runs to analyse
Rscript src/get_runs.R

# run
cat out/runs.tsv | sed 1d |
{
  while read path experiment genome sublibrary parse_analysis_subdir do_integration ; do
    . src/submit_generate_qc_report.sh
  done
}



# # run manually:
# cat out/runs.tsv | sed 1d |
# {
#   while read path experiment genome sublibrary parse_analysis_subdir ; do
#     Rscript src/run_generate_QC_report.R   $experiment   $genome   $sublibrary   $parse_analysis_subdir
#   done
# }
