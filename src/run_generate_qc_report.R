args <- commandArgs(trailingOnly=TRUE)
names(args) <- c("data_dir",
                 "output_dir",
                 "sample_metadata_file")

# parse comma-delimited args
args[grepl(",", args)] <- strsplit(args[grepl(",", args)], ",")

# set NULL args
args <- as.list(args)
args[args=="NA"] <- list(NULL)

# testing and cacheing params
testing = F
do_timestamp = F
rerun = F

# # dir
# base_dir <- ifelse(Sys.info()["nodename"] == "Alexs-MacBook-Air-2.local",
#                    "/Volumes/TracerX/",
#                    "/camp/project/tracerX/")
# setwd(paste0(base_dir,"working/VHL_GERMLINE/tidda/vhl/"))

# testing: setwd, default params, load_all # base_dir=ifelse(Sys.info()["nodename"]=="Alexs-MacBook-Air-2.local","/Volumes/TracerX/","/camp/project/tracerX/");setwd(paste0(base_dir,"working/VHL_GERMLINE/tidda/vhl/"));library(devtools);load_all();parse_dir=paste0(base_dir,"working/VHL_GERMLINE/tidda/parse_pipeline/analysis/");genome="hg38";sublibrary="comb";parse_analysis_subdir="all-well/DGE_filtered/";do_filtering=T;remove_doublets=T;min_nuclei_per_gene=5;min_nFeature_RNA=NULL;min_nCount_RNA=NULL;max_nCount_RNA=NULL;max_nFeature_RNA=NULL;max_percent_mito=NULL;vars_to_regress=c("percent_mito","nCount_RNA");n_dims=NULL;clustering_resolutions = seq(0.1, 0.8, by = 0.1);final_clustering_resolution=0.3;out_dir = NULL;sample_subset=NULL;cell_subset=NULL;do_timestamp = F;do_integration = F;integration_col="sample";final_annotations = NULL;final_annotations_lvl="partition";testing=F;rerun=F
# testing: pilot # args$experiment="221202_A01366_0326_AHHTTWDMXY";sublibrary="SHE5052A9_S101";final_annotations_lvl="cluster";final_annotations=final_annotations_list[[args$experiment]][[final_annotations_lvl]]
# testing: 8 SLs + pilot # args$experiment=c("221202_A01366_0326_AHHTTWDMXY","230210_A01366_0351_AHNHCFDSX5");sublibrary=c("SHE5052A9_S101","comb");sample_subset=strsplit("N045_V008C,N045_V010,N045_V003,N045_V004,N045_N001,N059_V001,N059_M001,N059_V102A,N059_N001,N059_V003,N059_V103,N088_V006,N088_V004,N088_V008,N088_V106,N088_V108,N090_V116,N090_V124D,N090_V126,N090_V127,N090_V124A,N090_V128,N090_N002,K891_V014", ",")[[1]]
# testing: tumour_cells # args$experiment=c("221202_A01366_0326_AHHTTWDMXY","230210_A01366_0351_AHNHCFDSX5");sublibrary=c("SHE5052A9_S101","comb");sample_subset=strsplit("N045_V008C,N045_V010,N045_V003,N045_V004,N045_N001,N059_V001,N059_M001,N059_V102A,N059_N001,N059_V003,N059_V103,N088_V006,N088_V004,N088_V008,N088_V106,N088_V108,N090_V116,N090_V124D,N090_V126,N090_V127,N090_V124A,N090_V128,N090_N002,K891_V014", ",")[[1]];cell_subset=dplyr::filter(readr::read_tsv("out/230210_A01366_0351_AHNHCFDSX5_x_221202_A01366_0326_AHHTTWDMXY/hg38/comb_x_SHE5052A9_S101/all-well/DGE_filtered/unintegrated/cell_annotations.tsv"), cluster_annot=="tumour")$cell

# filters from Li et al., 2023 (PMID: 36563681) (snRNAseq paper)
# min_nCount_RNA = 300
# min_nFeature_RNA = 200
# max_nFeature_RNA = 10000
# max_percent_mito = 10
# max n UMIs = 10000
# min n UMIs = 1000
# no doublet filtering

# filters from Braun et al., 2021
# min n UMI = 1000
# min_nFeature_RNA = 500
# min_percent_mito = 15
# doublet filtering: scran::doubletCells

library(devtools)
load_all()

# run qc report
load_all() ; generate_qc_report(
  data_dir = args$data_dir,
  output_dir = args$output_dir,
  sample_metadata_file = args$sample_metadata_file,
  remove_doublets = T,
  do_filtering = T,
  do_cell_cycle_scoring = T,
  min_nuclei_per_gene = 5,
  min_nCount_RNA = 300,
  max_nCount_RNA = 10000,
  min_nFeature_RNA = 200,
  max_nFeature_RNA = 10000,
  max_percent_mito = 15)



