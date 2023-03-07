args <- commandArgs(trailingOnly=TRUE)
experiment <- strsplit(args[1], ",")[[1]]
genome <- args[2]
sublibrary <- strsplit(args[3], ",")[[1]]
parse_analysis_subdir <- args[4]
do_integration <- args[5]
sample_subset <- args[6]

# set sample subset arg
if (sample_subset == "NA") {
  sample_subset <- NULL
} else {
  sample_subset <- strsplit(sample_subset, ",")[[1]]
}

base_dir <- ifelse(Sys.info()["nodename"] == "Alexs-MacBook-Air-2.local",
                   "/Volumes/TracerX/",
                   "/camp/project/tracerX/")

# testing: setwd # base_dir=ifelse(Sys.info()["nodename"]=="Alexs-MacBook-Air-2.local","/Volumes/TracerX/","/camp/project/tracerX/");setwd(paste0(base_dir,"working/VHL_GERMLINE/tidda/vhl/"));library(devtools);load_all();parse_dir=paste0(base_dir,"working/VHL_GERMLINE/tidda/parse_pipeline/");genome="hg38";sublibrary="comb";parse_analysis_subdir="all-well/DGE_unfiltered/";do_filtering=T;remove_doublets=T;min_nuclei_per_gene=NULL;min_nFeature_RNA=NULL;max_nCount_RNA=NULL;max_nFeature_RNA=NULL;max_percent_mito=NULL;n_dims=NULL;clustering_resolutions = seq(0.1, 0.8, by = 0.1);final_clustering_resolution=NULL;out_dir = NULL;sample_subset = NULL;do_timestamp = F;do_integration = F;integration_col="sample";
# testing: pilot # experiment="221202_A01366_0326_AHHTTWDMXY";genome="hg38";sublibrary="SHE5052A9_S101";parse_analysis_subdir="all-well/DGE_filtered";do_integration=F
# testing: full  # experiment="230127_A01366_0343_AHGNCVDMXY";genome="hg38";sublibrary="SHE5052A11_S164";parse_analysis_subdir="all-well/DGE_filtered";do_integration=F

# filters taken from Li et al., 2023 (PMID: 36563681)
# >300 total transcripts -> 100
# >200 & <10000 genes expressed
# >1000 UMIs
# <10% mt genes expressed

library(devtools) ; load_all() ; generate_qc_report(
  parse_dir = paste0(base_dir, "working/VHL_GERMLINE/tidda/parse_pipeline/"),
  experiment = experiment,
  genome = genome,
  sublibrary = sublibrary,
  parse_analysis_subdir = parse_analysis_subdir,
  remove_doublets = T,
  do_filtering = T,
  min_nuclei_per_gene = 5,
  min_nCount_RNA = 100,
  max_nCount_RNA = 10000,
  min_nFeature_RNA = 100,
  max_nFeature_RNA = 10000,
  max_percent_mito = 12,
  final_clustering_resolution = 0.3,
  do_integration = do_integration,
  do_timestamp = T,
  sample_subset = sample_subset)

