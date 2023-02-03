args <- commandArgs(trailingOnly=TRUE)
experiment <- args[1]
genome <- args[2]
sublibrary <- args[3]
parse_analysis_subdir <- args[4]

# testing: # experiment="221202_A01366_0326_AHHTTWDMXY";genome="premRNAhg38";sublibrary="SHE5052A9_S101";parse_analysis_subdir="all-well/DGE_filtered"

library(devtools) ; load_all()

generate_QC_report(
  experiment = experiment,
  parse_pipeline_dir = "/camp/project/tracerX/working/VHL_GERMLINE/tidda/parse_pipeline/",
  genome = genome,
  sublibrary = sublibrary,
  parse_analysis_subdir = parse_analysis_subdir,
  out_dir = NULL,
  sample_subset = NULL,
  do_timestamp = F)
