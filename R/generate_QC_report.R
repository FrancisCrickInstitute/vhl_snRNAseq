generate_QC_report <- function(
    experiment,
    parse_pipeline_dir = "/camp/project/tracerX/working/VHL_GERMLINE/tidda/parse_pipeline/",
    genome = "hg38",
    sublibrary = "comb",
    parse_analysis_subdir = "/all-well/DGE_unfiltered/",
    n_dims = 20,
    out_dir = NULL,
    sample_subset = NULL,
    do_timestamp = F,
    do_integration = F) {

  # testing: # library(devtools);load_all();experiment="221202_A01366_0326_AHHTTWDMXY";genome="hg38_mm10";sublibrary="SHE5052A9_S101";parse_analysis_subdir="all-well/DGE_filtered";parse_pipeline_dir = paste0(get_base_dir(), "/parse_pipeline/");n_dims = 20;out_dir = NULL;sample_subset = NULL;do_timestamp = F;do_integration = T

  # cell quality control
  # -> unusually high transcript/gene counts indicate multiplets
  # -> unusually low transcript/gene counts indicate barcoding of cells with
  #    damaged membranes (uninformative)
  # -> high % mito genes indicates death / loss of cytoplasmic RNA / increased
  #    apoptosis (for scRNA-seq, not snRNA-seq)
  # -> high % globin and low % ribo suggests erythrocytes

  # set directories
  data_dir <- paste(parse_pipeline_dir, "analysis", experiment, genome, sublibrary, sep = "/")
  dge_dir <- paste(data_dir, parse_analysis_subdir, sep = "/")
  out <- define_out(experiment, genome, sublibrary, parse_analysis_subdir,
                    out_dir, do_timestamp, do_integration)
  plots <- list()

  # plot summary stats
  cat("Plotting summary stats...\n")
  plots[["run_summary_stats"]] <- readr::read_csv(paste0(data_dir, "/agg_samp_ana_summary.csv"),
                                                  show_col_types = FALSE) %>%
    tidyr::pivot_longer(cols = -statistic, names_to = "patient_id") %>%
    dplyr::mutate(statistic = statistic %>% gsub(paste0(genome, "\\_"), "", .)) %>%
    dplyr::filter(statistic %in% c("median_tscp_per_cell",
                                "median_genes_per_cell",
                                "tso_fraction_in_read1",
                                "fraction_tscp_in_cells")) %>%
    ggplot2::ggplot(ggplot2::aes(x = patient_id, y = value)) +
    ggplot2::geom_col() +
    ggplot2::facet_grid(statistic ~ ., scales = "free") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

  # load parse output to seurat object, save to outdir
  seu <- load_parse_to_seurat(
    dge_dir, min_genes_per_cell = 0, min_cells_per_gene = 0, sample_subset
  )
  saveRDS(seu, file = paste0(out$base, "/seu.rds"))

  # cells per sample
  cat("Plotting cells per sample...\n")
  plots[["cells_per_sample_bar"]] <- dplyr::tibble(
    sample_id = names(table(seu$sample)),
    n_cells = table(seu$sample)) %>%
    ggplot2::ggplot(ggplot2::aes(x = sample_id, y = n_cells)) +
    ggplot2::geom_col() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

  # check distribution of genes per cell and cells per gene
  cat("Plotting distribution of number of cells/genes...\n")
  plots[["genes_per_cell_dist"]] <- seu@meta.data %>%
    dplyr::as_tibble() %>%
    dplyr::transmute(nCount_RNA, nFeature_RNA,
                     `bottom 95%` = 0.95, all = 1, `top 5%` = 0.95) %>%
    tidyr::pivot_longer(c("nCount_RNA", "nFeature_RNA"),
                        names_to = "statistic") %>%
    tidyr::pivot_longer(c("bottom 95%", "all", "top 5%"),
                        names_to = "label", values_to = "quantile") %>%
    dplyr::group_by(label) %>%
    dplyr::filter(grepl("top", label) & value >= quantile(value, quantile) |
                    grepl("bottom", label) & value <= quantile(value, quantile) |
                    label == "all") %>%
    dplyr::mutate(label = paste0(label, " (n = ", dplyr::n(), ")")) %>%
    ggplot2::ggplot(ggplot2::aes(value)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ..density..),
                            bins = 50,
                            fill = "grey") +
    ggplot2::geom_density() +
    ggplot2::facet_wrap(statistic ~ label, scales = "free") +
    ggplot2::theme_classic()

  # get proportions of relevant transcript types
  purrr::pwalk(transcript_types, function(...) {
    tt <- tibble::tibble(...)
    cat("Calculating %", tt$name, "genes...\n")
    cat(tt$message, "\n\n")
    # use `<<` for global assignment
    seu <<- seu %>%
      Seurat::PercentageFeatureSet(
        pattern = tt$pattern,
        col.name = paste0("percent_", tt$name)
      )
  })

  # plot transcript type abundances
  cat("Plotting abundance of different transcript types...\n")
  plots[["genes_per_cell_vln"]] <- Seurat::VlnPlot(
    seu, c("tscp_count", "gene_count"), group.by = "sample", raster = F)
  plots[["transcript_types_per_cell_vln"]] <- Seurat::VlnPlot(
    seu, c("percent_mito", "percent_ribo", "percent_globin"),
    group.by = "sample", pt.size = 0.1, ncol = 3, raster = F)
  plots[["rna_molecules_vs_rna_transcripts_scatter"]] <- Seurat::FeatureScatter(
    seu, "nCount_RNA", "nFeature_RNA", group.by = "sample", raster = F)
  plots[["rna_molecules_vs_pct_mito_scatter"]] <- Seurat::FeatureScatter(
    seu, "nCount_RNA", "percent_mito", group.by = "sample", raster = F)
  plots[["pct_globin_vs_ribo_scatter"]] <- Seurat::FeatureScatter(
    seu, "percent_globin", "percent_ribo", group.by = "sample", raster = F)

  # most highly expressed genes
  cat("Plotting top 20 most highly expressed genes...\n")
  plots[["top_genes_boxplot"]] <- boxplot_top_genes(seu, 20)

  # NORMALISATION AND SCALING (SCTransform) ----
  cat("Running SCTransform() for normalisation and scaling...\n")
  seu <- Seurat::FindVariableFeatures(seu)
  seu <- Seurat::SCTransform(seu)

  # LINEAR DIMENSIONALITY REDUCTION (PCA) ----

  # run PCA
  seu <- Seurat::RunPCA(seu)
  saveRDS(seu, file = paste0(out$base, "/seu_transformed.rds"))

  # assign cell cycle phase scores
  cat("Assigning cell cycle phases...\n")
  seu <- Seurat::CellCycleScoring(
    seu,
    s.features = Seurat::cc.genes$s.genes,
    g2m.features = Seurat::cc.genes$g2m.genes,
    set.indent = T
  )

  # pca plots
  cat("Plotting principle componentss...\n")
  plots[["pca_by_sample"]] <- Seurat::DimPlot(seu, reduction = "pca", group.by = "sample", raster = F)
  plots[["pca_by_phase"]] <- Seurat::DimPlot(seu, reduction = "pca", group.by = "Phase", raster = F)
  plots[["pca_dim_loadings"]] <- Seurat::VizDimLoadings(seu, dims = 1:4, reduction = "pca")
  plots[["pca_dim_heatmap"]] <- Seurat::DimHeatmap(seu, dims = 1:15, cells = 500, balanced = TRUE, raster = F)

  # elbow plot
  cat("Plotting elbow plot...\n")
  plots[["pca_elbow"]] <- Seurat::ElbowPlot(seu)

  # NON-LINEAR DIMENSIONALITY REDUCTION (UMAP) ----

  # run umap, plot
  seu <- Seurat::RunUMAP(seu, dims = 1:n_dims)
  plots[["umap_by_sample"]] <- Seurat::DimPlot(seu, reduction = "umap", group.by = "sample", raster = F)

  # INTEGRATION(?) ----

  if(do_integration == T){
    # integrate samples
    cat("Integrating samples...\n")
    seu_list <- Seurat::SplitObject(seu, split.by = "sample")
    for (i in 1:length(seu_list)) {
      seu_list[[i]] <- Seurat::NormalizeData(seu_list[[i]])
      seu_list[[i]] <- Seurat::FindVariableFeatures(seu_list[[i]], selection.method = "vst", nfeatures = 2000,
                                                    verbose = FALSE)
    }
    seu_anchors <- Seurat::FindIntegrationAnchors(object.list = seu_list, dims = 1:25)
    seu <- Seurat::IntegrateData(anchorset = seu_anchors, dims = 1:25)
    Seurat::DefaultAssay(seu) <- "integrated"
  }

  # CLUSTERING ----

  # cluster umap at different resolutions
  cat("Getting UMAP at different resolutions...\n")
  resolutions <- seq(0.1, 0.8, by = 0.1)
  seu <- Seurat::FindNeighbors(seu, dims = 1:n_dims)
  seu <- Seurat::FindClusters(seu, resolution = resolutions)
  saveRDS(seu, file = paste0(out$base, "/seu_transformed_and_clustered.rds"))

  # cluster tree plot of different resolutions
  cat("Plotting clustering tree for UMAP at different resolutions...\n")
  library(clustree)
  snn_res_prefixes <- paste0(Seurat::DefaultAssay(seu), "_snn_res.")
  plots[["clustering_tree"]] <- clustree::clustree(
    seu@meta.data[, grep(snn_res_prefixes, colnames(seu@meta.data))],
    prefix = snn_res_prefixes
  )

  # umap plot at different resolutions
  cat("Plotting UMAP at different resolutions...\n")
  for(res in resolutions) {
    plots[[paste0("umap_by_cluster_res.", res)]] <- Seurat::DimPlot(
      seu, group.by = paste0(Seurat::DefaultAssay(seu), "_snn_res.", res), raster = F)
  }

  # save plots
  cat("Saving all plots...\n")
  library(gridExtra)
  plots_grob <- marrangeGrob(grobs = plots, nrow = 1, ncol = 1)
  ggsave(paste0(out$base, "/qc_report.pdf"), plots_grob,
         width = 21, height = 20, units = "cm")
}
