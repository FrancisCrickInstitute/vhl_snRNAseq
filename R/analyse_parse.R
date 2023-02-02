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
    parse_analysis_subdir = "/all-well/DGE_unfiltered/",
    out_dir = NULL,
    sample_subset = NULL,
    min_genes = 3,
    min_cells = 100,
    max_percent_mito = 5,
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
  parse_dir = paste0(base_dir, "/parse_pipeline/analysis/archive_230130/")
  parse_analysis_subdir = "/all-well/DGE_filtered/"
  out_dir = NULL
  sample_subset = c("N090_V1024A", "N090_V1027")
  min_genes_per_cell = 3
  min_cells_per_gene = 100
  max_percent_mito = 5
  do_timestamp = F

  # get output paths
  out <- define_out(experiment, genome, out_dir, do_timestamp)

  # save arguments
  # dput(args, file = out$args)

  # load parse output to seurat object
  seu <- load_parse_to_seurat(
    experiment, genome, parse_dir, parse_analysis_subdir,
    min_genes_per_cell, min_cells_per_gene, sample_subset
  )

  # 1) QUALITY CONTROL PER GENE AND CELL ----
  cat("1) Performaing cell/gene quality control...\n")

  # cell quality control
  # -> unusually high transcript/gene counts indicate multiplets
  # -> unusually low transcript/gene counts indicate barcoding of cells with
  #    damaged membranes (uninformative)
  # -> high % mito genes indicates death / loss of cytoplasmic RNA / increased 
  #    apoptosis (for scRNA-seq, not snRNA-seq)
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

  # mild cell filtering for doublets, dead cells, and ambient RNA
  # -> > 200 & < 5000 detected molecules per cell
  # -> < 10,000 genes per cell
  # -> < 8% mitochondrial RNA counts per cell
  seu <- seu %>%
    subset(subset = nFeature_RNA > 200 &
             nFeature_RNA < 5000 &
             gene_count < 10000 &
             percent.mito < max_percent_mito)
  # TODO: change percent.mito for snRNAseq? more stringent threshold?

  cat("Saving post-filtering Seurat object to", out$seu_post_filtering, "\n")
  saveRDS(seu, file = out$seu_post_filtering)
  
  # visualise new distribution
  pdf(paste0(out$base, "cell_post_qc_plot.pdf"))
  Seurat::VlnPlot(seu, c("nFeature_RNA", "percent.mito"),
                  group.by = "sample")
  dev.off()
  
  do_norm_and_scale <- FALSE
  if(do_norm_and_scale == TRUE) {
    
    # 2) NORMALISATION AND SCALING ----
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
    
    # perform PCA and colour by phase to check for cell cycle-based variation
    # if there are large differences due to cell cycle phase,
    # the we might regress out variation due to cell cycle
    seu <- Seurat::RunPCA(seu)
    pdf(paste0(out$base, "pca_vs_cell_cycle_phase_plot.pdf"))
    Seurat::DimPlot(seu,
                    reduction = "pca",
                    group.by = "Phase")
    dev.off()
    
  }
  
  # regressing out sources of unwanted variation with SCTransform
  # SCTransform is a better alternative to log transform normalisation.
  # It normalises data, performs variance stabilsation, and allows for 
  # additional covariates to be regressed out. 
  seu <- Seurat::SCTransform(seu, vars.to.regress = c("nCount_RNA", "percent.mito"))
  
  cat("Saving SCTransform-ed Seurat object to", out$seu_transformed, "\n")
  saveRDS(seu, file = out$seu_transformed)
  
  do_integration <- FALSE
  if(do_intergation == TRUE) {
    
    # TODO: integration????
    # Condition-specific clustering of the cells indicates that we need to 
    # integrate the cells across conditions to ensure that cells of the same 
    # cell type cluster together.
    
    # split cells into samples
    split_seu <- Seurat::SplitObject(seu, split.by = "sample")
    
    # iterate SCTransform across cycles
    split_seu <- split_seu %>%
      purrr::map(Seurat::SCTransform, vars.to.regress = c("nCount_RNA", "percent.mito"))
    
    cat("Saving SCTransform-ed split Seurat object to", out$seu_split_transformed, "\n")
    saveRDS(seu, file = out$seu_split_transformed)
    
  }
  
  # 3) DIMENSIONALITY REDUCTION ----
  
  # linear dimensionality reduction: perform PCA
  seu <- Seurat::RunPCA(seu, features = Seurat::VariableFeatures(seu))
  Seurat::VizDimLoadings(seu, dims = 1:2, reduction = "pca")
  
  # plot PCA
  Seurat::DimPlot(seu, reduction = "pca", group.by = "sample")
  Seurat::DimHeatmap(seu, dims = 1, cells = 500, balanced = TRUE)
  Seurat::DimHeatmap(seu, dims = 1:15, cells = 500, balanced = TRUE)
  
  # determine the dimensionality of the data
  Seurat::ElbowPlot(seu)
  n_dims <- 20
  
  # 4) CLUSTERING ----
  
  Seurat::DefaultAssay(seu) <- "SCT"
  seu <- Seurat::FindNeighbors(seu, dims = 1:n_dims)
  seu <- Seurat::FindClusters(seu, resolution = 0.5)
  head(Seurat::Idents(seu))
  
  # run non-linear dimensionality reduction: UMAP
  seu <- Seurat::RunUMAP(seu, dims = 1:n_dims)
  Seurat::DimPlot(seu, reduction = "umap")
  Seurat::DimPlot(seu, reduction = "umap", group.by = "sample")
  
  cat("Saving UMAP-ed Seurat object to", out$seu_umap, "\n")
  saveRDS(seu, file = out$seu_umap)

  # 5) FINDING MARKERS - DIFFERENTIAL GENE EXPRESSION ----

  # genes that define clusters
  seu_markers <- Seurat::FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  # top 10 genes per cluster
  top_markers_per_cluster <- seu_markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(
      rank = rank(-avg_log2FC)
    ) %>%
    dplyr::filter(rank <= 100) %>% 
    dplyr::arrange(rank, .by_group = T)
  
  # top gene per cluster
  seu %>% Seurat::FeaturePlot(
    features = dplyr::filter(top_markers_per_cluster, rank == 1)$gene
  )
  
  # top 10 genes per cluster
  seu %>%
    Seurat::DoHeatmap(
      features = dplyr::filter(top_markers_per_cluster, rank <= 10)$gene) + 
    Seurat::NoLegend()
  seu %>%
    dittoSeq::dittoDotPlot(
      vars = unique(dplyr::filter(top_markers_per_cluster, rank <= 3)$gene),
      group.by = "SCT_snn_res.0.5"
    )
  
  # 5) CELL ANNOTATION ----
  
  # celldex references from primary cell atlas
  ref <- celldex::HumanPrimaryCellAtlasData()
  seu_int_SingleR <- SingleR::SingleR(
    test = Seurat::GetAssayData(seu, slot = "data"),
    ref = ref,
    labels = ref$label.main)
  SingleR::plotScoreHeatmap(seu_int_SingleR)
  SingleR::plotDeltaDistribution(seu_int_SingleR)
  
  # add labels, cull labels with < 10 cells
  new_annots <- seu_int_SingleR$labels %>%
    table() %>%
    dplyr::as_tibble() %>%
    dplyr::rename(celltype = ".") %>%
    dplyr::mutate(SingleR_annot = replace(celltype, n < 10, "none")) %>%
    dplyr::right_join(dplyr::tibble(celltype = seu_int_SingleR$labels))
  seu$SingleR_annot <- new_annots$SingleR_annot
  dittoSeq::dittoDimPlot(seu, "SingleR_annot", size = 0.7)
  dittoSeq::dittoDimPlot(seu, "sample", size = 0.7)
  dittoSeq::dittoBarPlot(seu, var = "SingleR_annot", group.by = "sample")
  
  cat("Saving celldex-annotated Seurat object to", out$seu_celldex_annot, "\n")
  saveRDS(seu, file = out$seu_celldex_annot)
  
  # -> HIF metagene
  
  # -> literature markers
  
  # checking that all markers have been found
  unmatched_markers <- unlist(literature_markers)[
    !(unlist(literature_markers) %in% unique(rownames(seu@assays$RNA)))
  ]
  if(length(unmatched_markers) > 0) { 
    warning("Marker(s) ", paste(unmatched_markers, collapse = ", "), " not found!")
  }
  
  # generate module scores
  seu %>%
    Seurat::AddModuleScore(features = literature_markers)
  
  # 7) ENRICHMENT ANALYSIS ----

  # 8) TRAJECTORY ANALYSIS ----

}
