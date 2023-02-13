args <- commandArgs(trailingOnly=TRUE)
experiment <- args[1]
genome <- args[2]
sublibrary <- args[3]
parse_analysis_subdir <- args[4]
do_integration <- args[5]

# testing: pilot # experiment="221202_A01366_0326_AHHTTWDMXY";genome="hg38";sublibrary="SHE5052A9_S101";parse_analysis_subdir="all-well/DGE_filtered";do_integration=F
# testing: full  # experiment="230127_A01366_0343_AHGNCVDMXY";genome="hg38";sublibrary="SHE5052A11_S164";parse_analysis_subdir="all-well/DGE_filtered";do_integration=F

library(devtools) ; load_all()

generate_qc_report(
  experiment = experiment,
  parse_pipeline_dir = "/camp/project/tracerX/working/VHL_GERMLINE/tidda/parse_pipeline/",
  genome = genome,
  sublibrary = sublibrary,
  parse_analysis_subdir = parse_analysis_subdir,
  out_dir = NULL,
  sample_subset = NULL,
  do_timestamp = F,
  do_integration = do_integration)
