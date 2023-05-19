# dir
base_dir <- ifelse(Sys.info()["nodename"] == "Alexs-MacBook-Air-2.local",
                   "/Volumes/TracerX/",
                   "/camp/project/tracerX/")
setwd(paste0(base_dir,"working/VHL_GERMLINE/tidda/vhl/"))

# load package
library(devtools)
load_all()

# run
generate_qc_report(
  parse_dir = "/nemo/lab/turajlics/home/users/dengd/scRNA_seq/proj_PDOT_Parse/output/Parse_lib_1/",
  experiment = "",
  genome = "",
  sublibrary = "",
  max_percent_mito = 30,
  max_nCount_RNA = 30000,
  do_add_sample_metadata = F,
  do_add_summary_stats = F,
  out_dir = "out/PDOs/"
)

