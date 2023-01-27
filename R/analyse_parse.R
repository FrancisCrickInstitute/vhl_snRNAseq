#' Implement pipeline to QC and analyse Parse snRNAseq data
#'
#' @param experiment The name of the experiment directory in the Parse pipeline.
#' @param parse_pipeline_dir The full path to the directory where the Parse split-pipe pipeline was run. This should contain genomes/, expdata/, and analysis/ directories.
#' @param genome The name of the genome build used as a reference. This will also be the immediate subdirectory of the experiment's analysis. Default is human genome build 38 ("hg38").
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
    base_dir <- "/Volumes/TracerX/working/VHL_GERMLINE/tidda/"
  } else {
    base_dir <- "/camp/project/tracerX/working/VHL_GERMLINE/tidda/"
  }
  setwd(paste0(base_dir, "/vhl/"))
  library(devtools) ; load_all()
  experiment = "221202_A01366_0326_AHHTTWDMXY"
  genome = "hg38"
  parse_dir = paste0(base_dir, "/parse_pipeline/analysis/")
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
  #    damaged membranes (uninformative)
  # -> high % mito genes indicates loss of cytoplasmic RNA / increased apoptosis
  #    (for scRNA-seq, not snRNA-seq)
  # -> high % globin and low % ribo suggests erythrocytes

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

  # plot transcript type abundances
  cat("Visualising QC per cell and gene...\n")
  cat("Cell QC plots are saved to", out$base, "\n")
  pdf(paste0(out$base, "cell_qc_plots.pdf"), onefile = T)
  Seurat::VlnPlot(seu, c("percent.mito", "percent.ribo", "percent.globin"),
                  group.by = "sample", pt.size = 0.1, ncol = 3)
  Seurat::FeatureScatter(seu, "nCount_RNA", "nFeature_RNA", group.by = "sample")
  Seurat::FeatureScatter(seu, "nCount_RNA", "percent.mito", group.by = "sample")
  Seurat::FeatureScatter(seu, "percent.globin", "percent.ribo", group.by = "sample")
  dev.off()

  # most highly expressed genes
  pdf(paste0(out$base, "top_genes_plot.pdf"))
  boxplot_top_genes(seu, 20)
  dev.off()

  # mild cell filtering
  # -> < 8% mitochondrial RNA counts per cell
  # -> > 200 features per cell
  # -> < 5000 detected features per cell
  seu <- seu %>%
    subset(subset = nFeature_RNA > 200 &
                    nFeature_RNA < 5000 &
                    percent.mito < 8)

  # visualise new distribution
  pdf(paste0(out$base, "cell_post_qc_plot.pdf"))
  Seurat::VlnPlot(seu, c("nFeature_RNA", "percent.mito"),
                  group.by = "sample")
  dev.off()

  # 2) NORMALISATION AND SCALING
  cat("2) Performaing normalisation and scaling...\n")

  # log normalise
  seu <- seu %>%
    Seurat::NormalizeData(normalization.method = "LogNormalize",
                          scale.factor = 10000)

  # find variable features
  seu <- seu %>%
    Seurat::FindVariableFeatures(selection.method = "vst",
                                 nfeatures = 2000)

  # plot top 10 variable genes
  top_10 <- seu %>% Seurat::VariableFeatures() %>% head(10)
  pdf(paste0(out$base, "top_10_variable_genes_plot.pdf"))
  seu %>%
    Seurat::VariableFeaturePlot() %>%
    Seurat::LabelPoints(points = top_10, repel = TRUE)
  dev.off()

  # scale
  seu <- Seurat::ScaleData(seu, rownames(seu))

  # assign cell cycle scores
  s_genes <- Seurat::cc.genes$s.genes
  g2m_genes <- Seurat::cc.genes$g2m.genes
  seu <- Seurat::CellCycleScoring(
    seu, s.features = s_genes, g2m.features = g2m_genes,
    set.indent = T
  )

  # perform PCA and colour by phase
  # if there are large differences due to cell cycle phase,
  # the we might regress out variation due to cell cycle
  seu <- Seurat::RunPCA(seu)
  pdf(paste0(out$base, "pca_vs_cell_cycle_phase_plot.pdf"))
  Seurat::DimPlot(seu,
                  reduction = "pca",
                  group.by = "Phase")
  dev.off()
  
  # regressing out sources of unwanted variation with SCTransform
  # SCTransform is a better alternative to log transform normalisation.
  # It normalises data, performs variance stabilsation, and allows for 
  # additional covariates to be regressed out. 
  
  
  
  # 3) DIMENSIONALITY REDUCTION

  # 4) CLUSTERING

  # 5) CELL ANNOTATION
  # -> HIF metagene

  # 6) DIFFERENTIAL GENE EXPRESSION

  # 7) ENRICHMENT ANALYSIS

  # 8) TRAJECTORY ANALYSIS

}
