#' Add together two numbers
#'
#' @param experiment The name of the experiment directory in the Parse pipeline.
#' @param parse_pipeline_dir The full path to the directory where the Parse split-pipe pipeline was run. This should contain genomes/, expdata/, and analysis/ directories.
#' @param genome The name of the genome build used as a reference. This will also be the immediate subdirectory of the experiment's analysis.
#' @param parse_analysis_subdir The Parse output directory to be used.
#' @param out_dir The output directory. If none given, output will save to out/`experiment`/`genome`/.
#' @param sample_subset Sample IDs to subset to.
#' @return A Seurat object
analyse_parse <- function(
    experiment,
    parse_pipeline_dir = "/camp/project/tracerX/working/VHL_GERMLINE/tidda/parse_pipeline/",
    genome = "hg38",
    parse_analysis_subdir = "/all-well/DGE_filtered/",
    out_dir = NULL,
    sample_subset = NULL,
    min_genes = 3,
    min_cells = 100,
    do_timestamp = F
) {

  # capture function arguments (do not run when testing internally)
  args <- as.list(environment())

  # SETUP ----

  # for testing:
  if (Sys.info()["nodename"]=="Alexs-MacBook-Air-2.local") {
    setwd("/Volumes/TracerX/working/VHL_GERMLINE/tidda/vhl/")
  } else {
    setwd("/camp/project/tracerX/working/VHL_GERMLINE/tidda/vhl/")
  }
  library(devtools) ; load_all()
  experiment = "221202_A01366_0326_AHHTTWDMXY"
  genome = "hg38"
  parse_dir = "/camp/project/tracerX/working/VHL_GERMLINE/tidda/parse_pipeline/analysis/"
  parse_analysis_subdir = "/all-well/DGE_filtered/"
  out_dir = NULL
  sample_subset = c("N090_V1024A", "N090_V1027")
  min_genes = 3
  min_cells = 100
  do_timestamp = F

  # get output paths
  out <- define_out(experiment, genome, out_dir, do_timestamp)

  # save arguments
  # dput(args, file = out$args)

  # load parse output to seurat object
  seu <- load_parse_to_seurat(
    experiment, genome, parse_dir, parse_analysis_subdir,
    min_genes, min_cells, sample_subset
  )

  # 1) QUALITY CONTROL PER GENE AND CELL ----
  cat("1) Performaing cell/gene quality control...\n")

  # cell quality control
  # -> unusually high transcript/gene counts indicate multiplets
  # -> unusually low transcript/gene counts indicate barcoding of cells with
  #    damaged membranes
  # -> high % mito genes indicates loss of cytoplasmic RNA / increased apoptosis
  #    (for scRNA-seq, not snRNA-seq)

  cat("Saving pre-QC Seurat object to", out$seu_pre_qc, "\n")
  saveRDS(seu, file = out$seu_pre_qc)
  cat("Performing quality control...\n")

  # check proportion of relevant transcript types
  purrr::pwalk(transcript_types, function(...) {
    tt <- tibble::tibble(...)
    cat("Calculating %", tt$name, "genes...\n")
    cat(tt$message, "\n\n")
    # use `<<` for global assignment
    seu <<- seu %>%
      Seurat::PercentageFeatureSet(
        pattern = tt$pattern,
        col.name = paste0("percent.", tt$name)
      )
  })

  # plots
  cat("Visualising QC per cell and gene...\n")
  cat("Cell QC plots are saved to", out$base, "\n")
  pdf(paste0(out$base, "cell_qc_plots.pdf"), onefile = T)
  Seurat::VlnPlot(seu, c("percent.mito", "percent.ribo", "percent.globin"),
                  group.by = "sample",
                  pt.size = 0.1, ncol = 3)
  Seurat::FeatureScatter(seu, "nCount_RNA", "nFeature_RNA")
  Seurat::FeatureScatter(seu, "nCount_RNA", "percent.mito")
  Seurat::FeatureScatter(seu, "percent.globin", "percent.ribo")
  dev.off()

  # 2) NORMALISATION AND SCALING

  # 3) DIMENSIONALITY REDUCTION

  # 4) CLUSTERING

  # 5) CELL ANNOTATION

  # 6) DIFFERENTIAL GENE EXPRESSION

  # 7) ENRICHMENT ANALYSIS

  # 8) TRAJECTORY ANALYSIS

}
