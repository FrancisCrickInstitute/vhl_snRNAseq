#!/bin/bash

# set directories, activate R
base_dir=/camp/project/tracerX/working/VHL_GERMLINE/tidda/ # base_dir=/Volumes/TracerX/working/VHL_GERMLINE/tidda/
wkdir=$base_dir/vhl/
parse_pipeline_dir=$base_dir/parse_pipeline/
cd $wkdir
ml R/4.2.0-foss-2021b

# run
cat out/runs.tsv | sed 1d |
{
  while read path experiment genome sublibrary parse_analysis_subdir ; do
    . src/submit_generate_QC_report.sh
  done
}

# run manually:
cat out/runs.tsv | sed 1d |  head -3 | tail -1 |
{
  while read path experiment genome sublibrary parse_analysis_subdir ; do
    Rscript src/run_generate_QC_report.R   $experiment   $genome   $sublibrary   $parse_analysis_subdir
  done
}
