#' analyse snRNAseq output from Parse Biosciences split-pipe
#'
#' @param parse_dir Path to Parse Biosciences split-pipe directory.
#' @param experiment The name of the experiment run via split-pipe. If you want to combine multiple runs, pass a vector.
#' @param genome The reference genome build used by split-pipe for alignment. Default is "hg38".
#' @param sublibrary The sublibrary to analyse. Default is "comb". If you are combining multiple runs, pass a vector whose order matches that of the vector passed to `experiment`.
#' @param parse_analysis_subdir Parse Biosciences split-pipe DGE subdirectory. This directory must contain `DGE.mtx`, `all_genes.csv`, and `cell_metadata.csv` files. This will indicate which wells to include and whether to use the `DGE_filtered/` or `DGE_unfiltered/` matrix. Default is "all-well/DGE_filtered/".
#' @param do_filtering If `TRUE`, apply quality control filters to genes and nuclei.
#' @param remove_doublets If `TRUE`, removes suspected doublet nuclei, as detected by the scDblFinder package.
#' @param min_nuclei_per_gene \emph{Gene filter}. The minimum number of nuclei in which a gene must be present. Genes that are present in very few nuclei are uninformative and unlikely to help in differentiating groups of nuclei. In general, most genes removed by this filtering will be those not detected in any nucleus.
#' @param min_nCount_RNA \emph{Nucleus filter}. The minimum number of transcripts detected per nucleus.
#' @param max_nCount_RNA \emph{Nucleus filter}. The maximum number of transcripts detected per nucleus. An unusually high number of RNA molecules suggests that the nucleus is a multiplet.
#' @param min_nFeature_RNA \emph{Nucleus filter}. The minimum number of genes detected per nucleus. An unusually low number of genes suggests that the nucleus has a damaged membrane and is therefore low quality.
#' @param max_nFeature_RNA \emph{Nucleus filter}. The maximum number of genes detected per nucleus. An unusually high number of genes suggests that the nucleus is a multiplet.
#' @param max_percent_mito \emph{Nucleus filter}. The maximum percentage of mitochondrial transcripts per nucleus. Overrepresentation of mitochondrial transcripts suggests cell death, loss of cytoplasmic RNA, or heightened apoptosis.
#' @param n_dims Optional. The dimensionality of the dataset to use for downstream analysis. This is the number of principal components believed to capture the majority of true biological signal in the dataset. This can be decided by consulting the elbow plot. If no value given, dimensionality is calculated using the `intrinsicDimensions` package.
#' @param vars_to_regress Vector of variables to regress out of the SCTransform residuals. This prevents these variables from contributing much to the PCA during dimensionality reduction, and thus confounding the analysis. Default is c("percent_mito", "nCount_RNA").
#' @param cluster_resolutions Optional. A vector of clustering resolutions to test. A higher resolution will result in a larger number of communities. Default is 0.1-0.8.
#' @param final_clustering_resolution Optional. The chosen clustering resolution of the dataset to use for downstream analysis. This is the resolution at which clusters appear to capture true biological groupings of interest in the dataset. This can be decided by consulting the clustering tree. If no value given, downstream analysis will not proceed.
#' @param out_dir Optional. Output directory. If no value given, the output will be saved to `out/{experiment}/{genome}/{sublibrary}/{parse_analysis_subdir}/{integrated,unintegrated}`.
#' @param sample_subset Optional. Vector of sample IDs to subset to.
#' @param do_timestamp If `TRUE`, will save the output to a time-stamped subdirectory (e.g. `20230217_105554/`).
#' @param do_integration If `TRUE`, will integrate the dataset and save the output to `integrated/` subdirectory. If `FALSE` (default), output will save to `unintegrated/` subdirectory.
#' @param integration_col If `do_integration` is `TRUE`, Seurat metadata column upon which to integrate the dataset.
#' @param final_annotations Named vector of celltype annotations for the clusters/partitions in the dataset.
#' @param final_annotations_lvl The level of grouping (cluster or partition) at which the `final_annotations` should be applied. Default is "partition".
#' @param testing If `testing` is `TRUE`, a test run of the function is performed using a subset of a small dataset. The output is saved to `out/test/`.
#' @param rerun If `rerun` is `TRUE`, caches will be invalidated and all analyses will be run from scratch. If `FALSE` (default) and the pipeline has been run previously in the same directory, caches will be used.
#' @return A Seurat object.
#'
#' @export
generate_qc_report <-
  function(parse_dir = "/camp/project/tracerX/working/VHL_GERMLINE/tidda/parse_pipeline/",
           experiment,
           genome = "hg38",
           sublibrary = "comb",
           parse_analysis_subdir = "all-well/DGE_filtered/",
           do_filtering = T,
           remove_doublets = T,
           min_nuclei_per_gene = 5,
           min_nCount_RNA = NULL,
           max_nCount_RNA = NULL,
           min_nFeature_RNA = NULL,
           max_nFeature_RNA = NULL,
           max_percent_mito = NULL,
           vars_to_regress = c("percent_mito", "nCount_RNA"),
           n_dims = NULL,
           clustering_resolutions = seq(0.1, 0.8, by = 0.1),
           final_clustering_resolution = NULL,
           out_dir = NULL,
           sample_subset = NULL,
           do_timestamp = F,
           do_integration = F,
           integration_col = "sample",
           final_annotations = NULL,
           final_annotations_lvl = "partition",
           testing = F,
           rerun = F) {

  # for internal test runs
  # testing: setwd, default params, load_all # base_dir=ifelse(Sys.info()["nodename"]=="Alexs-MacBook-Air-2.local","/Volumes/TracerX/","/camp/project/tracerX/");setwd(paste0(base_dir,"working/VHL_GERMLINE/tidda/vhl/"));library(devtools);load_all();parse_dir=paste0(base_dir,"working/VHL_GERMLINE/tidda/parse_pipeline/");genome="hg38";sublibrary="comb";parse_analysis_subdir="all-well/DGE_filtered/";do_filtering=T;remove_doublets=T;min_nuclei_per_gene=5;min_nFeature_RNA=NULL;min_nCount_RNA=NULL;max_nCount_RNA=NULL;max_nFeature_RNA=NULL;max_percent_mito=NULL;vars_to_regress=c("percent_mito","nCount_RNA");n_dims=NULL;clustering_resolutions = seq(0.1, 0.8, by = 0.1);final_clustering_resolution=0.3;out_dir = NULL;sample_subset = NULL;do_timestamp = F;do_integration = F;integration_col="sample";final_annotations = NULL;final_annotations_lvl="partition";testing=F;rerun=F
  # testing: Li filters # remove_doublets=T;min_nuclei_per_gene = 5;min_nCount_RNA = 300;max_nCount_RNA = 10000;min_nFeature_RNA = 200;max_nFeature_RNA = 10000;max_percent_mito = 15
  # testing: pilot # experiment="221202_A01366_0326_AHHTTWDMXY";sublibrary="SHE5052A9_S101";final_annotations_lvl="cluster";final_annotations=final_annotations_list[[experiment]][[final_annotations_lvl]]
  # testing: 2 SLs # experiment="230127_A01366_0343_AHGNCVDMXY";sublibrary="comb"
  # testing: 8 SLs # experiment="230210_A01366_0351_AHNHCFDSX5";sublibrary="comb"
  # testing: 8 SLs + pilot # experiment=c("221202_A01366_0326_AHHTTWDMXY","230210_A01366_0351_AHNHCFDSX5");sublibrary=c("SHE5052A9_S101","comb");sample_subset=strsplit("N045_V008C,N045_V010,N045_V003,N045_V004,N045_N001,N059_V001,N059_M001,N059_V102A,N059_N001,N059_V003,N059_V103,N088_V006,N088_V004,N088_V008,N088_V106,N088_V108,N090_V116,N090_V124D,N090_V126,N090_V127,N090_V124A,N090_V128,N090_N002,K891_V014", ",")[[1]]
  # testing: args  # library(devtools);load_all(); args <- dget("out/230127_A01366_0343_AHGNCVDMXY/hg38/comb/all-well/DGE_filtered/args_for_generate_qc_report.R") ; list2env(args,globalenv())
  # testing: args  # args <- dget("out/221202_A01366_0326_AHHTTWDMXY/hg38/SHE5052A9_S101/all-well/DGE_filtered/unintegrated/args_for_generate_qc_report.R") ; args$parse_dir <- "/Volumes/TracerX/working/VHL_GERMLINE/tidda/parse_pipeline/"

  # capture arguments
  args <- list(parse_dir = parse_dir,
               experiment = experiment,
               genome = genome,
               sublibrary = sublibrary,
               parse_analysis_subdir = parse_analysis_subdir,
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
               clustering_resolutions = clustering_resolutions,
               final_clustering_resolution = final_clustering_resolution,
               out_dir = out_dir,
               sample_subset = sample_subset,
               do_timestamp = do_timestamp,
               do_integration = do_integration,
               integration_col = integration_col,
               final_annotations = final_annotations,
               final_annotations_lvl = final_annotations_lvl,
               testing = testing,
               rerun = rerun)

  # check arguments
  check_args(args)

  # define output directory
  out <- get_out(out_dir,
                 experiment,
                 genome,
                 sublibrary,
                 parse_analysis_subdir,
                 do_integration,
                 do_timestamp)
  dir.create(out$base, showWarnings = F, recursive = T)

  # set sample groupings to check at the clustering stage
  groupings <-
    c("sample", "percent_mito", "percent_ribo", "percent_globin",
      "date_prep", "nih_pid", "rin", "lesion_type", "tumour_size", "fuhrman_grade") %>%
    # if genome is human, do cell cycle scoring (doesn't work with other genomes)
    { if (grepl("hg38", genome)) c(., "Phase") else . } %>%
    # if checking for doublets, add to the groupings
    { if (remove_doublets) c(., "doublet") else . } %>%
    # if multiple experiments, add to the groupings
    { if (length(experiment) > 1) c(., "experiment") else . }

  # statistics of interest in the Parse analysis
  statistics <-
    c("median_tscp_per_cell",
      "median_genes_per_cell",
      "fraction_tscp_in_cells",
      "cell_tscp_cutoff",
      "valid_barcode_fraction"
    )

  # save args to output directory
  dput(args, paste0(out$base, "/args_for_generate_qc_report.R"))

  # pass args to rmd
  rmarkdown::render(system.file("rmd", "generate_qc_report.rmd", package = "vhl"),
                    knit_root_dir = rprojroot::find_rstudio_root_file(),
                    output_dir = out$base,
                    output_file = "qc_report",
                    params = list(args, out = out,
                                  groupings = groupings, statistics = statistics))

  }


# devtools::load_all();rmarkdown::render(system.file("rmd", "generate_qc_report.rmd", package = "vhl"),
#                   knit_root_dir = rprojroot::find_rstudio_root_file(),
#                   output_dir = "out/test/",
#                   output_file = "qc_report",
#                   params = list(testing = T))
