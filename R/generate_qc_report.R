#' analyse snRNAseq output from Parse Biosciences split-pipe
#'
#' @param data_dir Path to the directory containing the count matrix.
#' @param output_dir Optional. Default is `./output/`.
#' @param sample Optional. Name of the sample to be matched with the `sample` column in the `sample_metadata_file`.
#' @param sample_metadata_file Optional. Full path to sample metadata file. If none given and `do_add_sample_metadata` is `TRUE`, will take it from the experimental subdirectory of the `expdata/` directory in the Parse split-pipe output.
#' @param do_filtering If `TRUE`, apply quality control filters to genes and nuclei.
#' @param remove_doublets If `TRUE`, removes suspected doublet nuclei, as detected by the scDblFinder package.
#' @param do_cell_cycle_scoring If `TRUE`, scores cells by cell cycle stage.
#' @param min_nuclei_per_gene \emph{Gene filter}. The minimum number of nuclei in which a gene must be present. Genes that are present in very few nuclei are uninformative and unlikely to help in differentiating groups of nuclei. In general, most genes removed by this filtering will be those not detected in any nucleus.
#' @param min_nCount_RNA \emph{Nucleus filter}. The minimum number of transcripts detected per nucleus.
#' @param max_nCount_RNA \emph{Nucleus filter}. The maximum number of transcripts detected per nucleus. An unusually high number of RNA molecules suggests that the nucleus is a multiplet.
#' @param min_nFeature_RNA \emph{Nucleus filter}. The minimum number of genes detected per nucleus. An unusually low number of genes suggests that the nucleus has a damaged membrane and is therefore low quality.
#' @param max_nFeature_RNA \emph{Nucleus filter}. The maximum number of genes detected per nucleus. An unusually high number of genes suggests that the nucleus is a multiplet.
#' @param max_percent_mito \emph{Nucleus filter}. The maximum percentage of mitochondrial transcripts per nucleus. Overrepresentation of mitochondrial transcripts suggests cell death, loss of cytoplasmic RNA, or heightened apoptosis.
#' @param n_dims Optional. The dimensionality of the dataset to use for downstream analysis. This is the number of principal components believed to capture the majority of true biological signal in the dataset. This can be decided by consulting the elbow plot. If no value given, dimensionality is calculated using the `intrinsicDimensions` package.
#' @param vars_to_regress Vector of variables to regress out of the SCTransform residuals. This prevents these variables from contributing much to the PCA during dimensionality reduction, and thus confounding the analysis. Default is c("percent_mito", "nCount_RNA").
#' @param out_dir Optional. Output directory. If no value given, the output will be saved to `out/{experiment}/{genome}/{sublibrary}/{parse_analysis_subdir}/{integrated,unintegrated}`.
#' @param out_subdir Optional. Name of a subdirectory within the automatically generated output path to save output to.
#' @param sample_subset Optional. Vector of sample IDs to subset to.
#' @param cell_subset Optional. Vector of cell IDs to subset to.
#' @param do_timestamp If `TRUE`, will save the output to a time-stamped subdirectory (e.g. `20230217_105554/`).
#' @param do_integration If `TRUE`, will integrate the dataset and save the output to `integrated/` subdirectory. If `FALSE` (default), output will save to `unintegrated/` subdirectory.
#' @param integration_col If `do_integration` is `TRUE`, Seurat metadata column upon which to integrate the dataset.
#' @param final_annotations Named vector of celltype annotations for the clusters/partitions in the dataset.
#' @param final_annotations_lvl The level of grouping (cluster or partition) at which the `final_annotations` should be applied. Default is "partition".
#' @param testing If `TRUE`, a test run of the function is performed using a subset of a small dataset. The output is saved to `out/test/`.
#' @param rerun If `TRUE`, caches will be invalidated and all analyses will be run from scratch. If `FALSE` (default) and the pipeline has been run previously in the same directory, caches will be used.
#' @param do_add_sample_metadata If `TRUE`, will add sample metadata.
#' @param do_add_summary_stats If `TRUE`, will add summary statistics.
#' @return A Seurat object.
#'
#' @export
generate_qc_report <-
  function(data_dir,
           output_dir = "output/",
           sample = NULL,
           genome = "hg38",
           sample_metadata_file = NULL,
           do_filtering = T,
           remove_doublets = T,
           do_cell_cycle_scoring = T,
           min_nuclei_per_gene = 5,
           min_nCount_RNA = NULL,
           max_nCount_RNA = NULL,
           min_nFeature_RNA = NULL,
           max_nFeature_RNA = NULL,
           max_percent_mito = NULL,
           vars_to_regress = c("percent_mito", "nCount_RNA"),
           n_dims = NULL,
           out_dir = NULL,
           out_subdir = NULL,
           sample_subset = NULL,
           cell_subset = NULL,
           do_timestamp = F,
           do_integration = F,
           integration_col = "sample",
           final_annotations = NULL,
           final_annotations_lvl = "partition",
           testing = F,
           rerun = F,
           do_add_sample_metadata = T,
           do_add_summary_stats = T) {

  # for internal test runs
  # testing: setwd, default params, load_all # base_dir=ifelse(Sys.info()["nodename"]=="Alexs-MacBook-Air-2.local","/Volumes/TracerX/","/camp/project/tracerX/");setwd(paste0(base_dir,"working/VHL_GERMLINE/tidda/vhl/"));library(devtools);load_all();data_dir=paste0(base_dir,"working/VHL_GERMLINE/tidda/parse_pipeline/analysis/");experiment=NULL;genome="hg38";sublibrary="comb";parse_analysis_subdir="all-well/DGE_filtered/";do_filtering=T;remove_doublets=T;min_nuclei_per_gene=5;min_nFeature_RNA=NULL;min_nCount_RNA=NULL;max_nCount_RNA=NULL;max_nFeature_RNA=NULL;max_percent_mito=NULL;vars_to_regress=c("percent_mito","nCount_RNA");n_dims=NULL;out_dir = NULL;out_subdir=NULL;sample_subset=NULL;cell_subset=NULL;do_timestamp = F;do_integration = F;integration_col="sample";final_annotations = NULL;final_annotations_lvl="partition";testing=F;rerun=F;do_add_sample_metadata=T;do_add_summary_stats=T
  # testing: Li filters # remove_doublets=T;min_nuclei_per_gene = 5;min_nCount_RNA = 300;max_nCount_RNA = 10000;min_nFeature_RNA = 200;max_nFeature_RNA = 10000;max_percent_mito = 15
  # testing: pilot # experiment="221202_A01366_0326_AHHTTWDMXY";sublibrary="SHE5052A9_S101";final_annotations_lvl="cluster";final_annotations=final_annotations_list[[experiment]][[final_annotations_lvl]]
  # testing: 2 SLs # experiment="230127_A01366_0343_AHGNCVDMXY";sublibrary="comb"
  # testing: 8 SLs # experiment="230210_A01366_0351_AHNHCFDSX5";sublibrary="comb"
  # testing: 8 SLs + pilot # experiment=c("221202_A01366_0326_AHHTTWDMXY","230210_A01366_0351_AHNHCFDSX5");sublibrary=c("SHE5052A9_S101","comb");sample_subset=strsplit("N045_V008C,N045_V010,N045_V003,N045_V004,N045_N001,N059_V001,N059_M001,N059_V102A,N059_N001,N059_V003,N059_V103,N088_V006,N088_V004,N088_V008,N088_V106,N088_V108,N090_V116,N090_V124D,N090_V126,N090_V127,N090_V124A,N090_V128,N090_N002,K891_V014", ",")[[1]]
  # testing: args  # library(devtools);load_all(); args <- dget("out/230127_A01366_0343_AHGNCVDMXY/hg38/comb/all-well/DGE_filtered/args_for_generate_qc_report.R") ; list2env(args,globalenv())
  # testing: args  # args <- dget("out/221202_A01366_0326_AHHTTWDMXY/hg38/SHE5052A9_S101/all-well/DGE_filtered/unintegrated/args_for_generate_qc_report.R") ; args$data_dir <- "/Volumes/TracerX/working/VHL_GERMLINE/tidda/parse_pipeline/analysis/"
  # testing: PDOs  # dge_mtx_dir="/nemo/lab/turajlics/home/users/dengd/scRNA_seq/proj_PDOT_Parse/output/Parse_lib_1/all-well/DGE_unfiltered/";max_percent_mito=30;max_nCount_RNA=30000;do_add_sample_metadata=F;do_add_summary_stats=F;out_dir="out/PDOs/"
  # testing: PDOs v2 # data_dir = "/nemo/lab/turajlics/home/users/dengd/scRNA_seq/proj_PDOT_Parse/output/Parse_lib_1/";experiment = "";genome = "";sublibrary = "";max_percent_mito = 30;max_nCount_RNA = 30000;do_add_sample_metadata = F;do_add_summary_stats = T;out_dir = "out/PDOs/"

  # capture arguments
  args <- list(data_dir = data_dir,
               output_dir = output_dir,
               sample_metadata_file = sample_metadata_file,
               do_filtering = do_filtering,
               remove_doublets = remove_doublets,
               min_nuclei_per_gene = min_nuclei_per_gene,
               min_nCount_RNA = min_nCount_RNA,
               max_nCount_RNA = max_nCount_RNA,
               min_nFeature_RNA = min_nFeature_RNA,
               max_nFeature_RNA = max_nFeature_RNA,
               max_percent_mito = max_percent_mito,
               vars_to_regress = vars_to_regress,
               n_dims = n_dims,
               sample_subset = sample_subset,
               cell_subset = cell_subset,
               do_integration = do_integration,
               integration_col = integration_col,
               final_annotations = final_annotations,
               final_annotations_lvl = final_annotations_lvl,
               testing = testing,
               rerun = rerun,
               do_add_sample_metadata = do_add_sample_metadata,
               do_add_summary_stats = do_add_summary_stats)

  # create output directory
  dir.create(output_dir, showWarnings = F, recursive = T)

  # save args to output directory
  dput(args, paste0(output_dir, "/args_for_generate_qc_report.R"))

  # pass args to rmd
  rmarkdown::render(system.file("rmd", "generate_qc_report.rmd", package = "vhl"),
                    knit_root_dir = rprojroot::find_rstudio_root_file(),
                    output_dir = output_dir,
                    output_file = "qc_report",
                    params = list(args, out = out,
                                  groupings = groupings, statistics = statistics))

  }


# devtools::load_all();rmarkdown::render(system.file("rmd", "generate_qc_report.rmd", package = "vhl"),
#                   knit_root_dir = rprojroot::find_rstudio_root_file(),
#                   output_dir = "out/test/",
#                   output_file = "qc_report",
#                   params = list(testing = T))
