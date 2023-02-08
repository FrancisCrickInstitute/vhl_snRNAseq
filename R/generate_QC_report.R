generate_qc_report <- function(
    experiment,
    parse_pipeline_dir = "/camp/project/tracerX/working/VHL_GERMLINE/tidda/parse_pipeline/",
    genome = "hg38",
    sublibrary = "comb",
    parse_analysis_subdir = "/all-well/DGE_unfiltered/",
    n_dims = 50,
    out_dir = NULL,
    sample_subset = NULL,
    do_timestamp = F,
    do_integration = F,
    integration_col = "sample") {

  # testing: # library(devtools);load_all();experiment="221202_A01366_0326_AHHTTWDMXY";genome="hg38";sublibrary="SHE5052A9_S101";parse_analysis_subdir="all-well/DGE_filtered";parse_pipeline_dir = paste0(get_base_dir(), "/parse_pipeline/");n_dims = 50;out_dir = NULL;sample_subset = NULL;do_timestamp = F;do_integration = T;integration_col="sample";
  # testing: # experiment="230127_A01366_0343_AHGNCVDMXY";genome="hg38";sublibrary="comb"

  # cell quality control
  # -> unusually high transcript/gene counts indicate multiplets
  # -> unusually low transcript/gene counts indicate barcoding of cells with
  #    damaged membranes (uninformative)
  # -> high % mito genes indicates death / loss of cytoplasmic RNA / increased
  #    apoptosis (for scRNA-seq, not snRNA-seq)
  # -> high % globin and low % ribo suggests erythrocytes

  # gridExtra must be loaded in the enrivonment
  library(gridExtra)

  # set directories
  data_dir <- paste(parse_pipeline_dir, "analysis", experiment, genome, sublibrary, sep = "/")
  dge_dir <- paste(data_dir, parse_analysis_subdir, sep = "/")
  out <- define_out(experiment, genome, sublibrary, parse_analysis_subdir,
                    out_dir, do_timestamp, do_integration)
  plots <- list()

  # LOAD DATA ----

  # load parse output to seurat object with no cut-offs, remove NAs, save to outdir
  seu <- load_parse_to_seurat(
    dge_dir, min_genes_per_cell = 0, min_cells_per_gene = 0, sample_subset
  )
  seu <- subset(seu, subset = sample %in% seu$sample[!is.na(seu$sample)])
  saveRDS(seu, file = paste0(out$base, "/seu.rds"))

  # sample groupings to check at clustering stage
  groupings <- c("sample", "Phase",
                 "percent_mito", "percent_ribo", "percent_globin")

  # sample stats
  summary_stats <-  readr::read_csv(paste0(data_dir, "/agg_samp_ana_summary.csv"),
                                    show_col_types = FALSE) %>%
    tidyr::pivot_longer(cols = -statistic, names_to = "sample") %>%
    dplyr::mutate(statistic = statistic %>% gsub(paste0(genome, "\\_"), "", .)) %>%
    dplyr::filter(statistic %in% c("median_tscp_per_cell",
                                   "median_genes_per_cell",
                                   "tso_fraction_in_read1",
                                   "fraction_tscp_in_cells"))

  # sample metadata
  sample_metadata_file <- paste(parse_pipeline_dir, "expdata", experiment, "sample_metadata.tsv", sep = "/")
  if(file.exists(sample_metadata_file)) {

    # get sample / batch metadata
    sample_metadata <- sample_metadata_file %>%
      readr::read_tsv() %>%
      janitor::clean_names()

    # add to seu@meta.data
    seu@meta.data <- seu@meta.data%>%
      dplyr::left_join(sample_metadata,
                       by = "sample")
    rownames(seu@meta.data) <- colnames(seu)
    md_groupings <- c("date_prep", "patient_id", "rin", "sample_type", "size")

    # append to summary stats
    summary_stats <- summary_stats %>%
      dplyr::bind_rows(
        sample_metadata %>%
          tidyr::pivot_longer(cols = tidyr::all_of(md_groupings) & where(is.numeric),
                              names_to = "statistic") %>%
          dplyr::select(statistic, sample, value)
      )
    groupings <- c(groupings, md_groupings)

  }

  # CELL AND GENE SUMMARY STATISTICS ----

  # plot summary stats
  cat("Plotting summary stats...\n")
  plots[["run_summary_stats"]] <- summary_stats %>%
    ggplot2::ggplot(ggplot2::aes(x = sample, y = value)) +
    ggplot2::geom_col() +
    ggplot2::facet_grid(statistic ~ ., scales = "free") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

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
    cat("\tCalculating %", tt$name, "genes...\n")
    cat("\t", tt$message, "\n\n")
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

  # highest expressed genes
  cat("Plotting top 20 most highly expressed genes...\n")
  plots[["top_genes_boxplot"]] <- boxplot_top_genes(seu)

  # highest variable genes (HVGs)
  cat("Plotting most highly variable genes (HVGs)...\n")
  seu <- Seurat::FindVariableFeatures(seu)
  plots[["hvg_scatter"]] <- Seurat::LabelPoints(
    plot = Seurat::VariableFeaturePlot(seu, raster = F),
    points = head(Seurat::VariableFeatures(seu), 10),
    repel = T
  ) +
    ggplot2::theme(legend.position = "top")

  if(do_integration == T){

    # INTEGRATION ----
    cat("Integrating samples...\n")

    seu_list <- seu %>%
      # separate samples
      Seurat::SplitObject(split.by = integration_col) %>%
      # normalise and identify variable features for each independent dataset
      purrr::map(function(samp) {
        samp %>%
          Seurat::NormalizeData() %>%
          Seurat::FindVariableFeatures()
      })

    # find anchors
    int_anchors <- seu_list %>%
      Seurat::FindIntegrationAnchors(
        anchor.features = Seurat::SelectIntegrationFeatures(seu_list)
    )

    # integrate, set as default assay
    seu <- Seurat::IntegrateData(anchorset = int_anchors)
    Seurat::DefaultAssay(seu) <- "integrated"

    # SCALING (POST-INTEGRATION) ----
    cat("Scaling integrated dataset...\n")
    seu <- Seurat::ScaleData(seu)

    # save
    saveRDS(seu, paste0(out$base, "seu_integrated.rds"))

  } else {

    # NORMALISATION AND SCALING (SCTransform) ----
    cat("Running SCTransform() for normalisation and scaling...\n")
    seu <- Seurat::SCTransform(seu)

    # save
    saveRDS(seu, file = paste0(out$base, "/seu_transformed.rds"))

  }

  # LINEAR DIMENSIONALITY REDUCTION (PCA) ----

  # run PCA
  cat("Running PCA...\n")
  seu <- Seurat::RunPCA(seu)

  # assign cell cycle phase scores
  cat("Assigning cell cycle phases...\n")
  seu <- Seurat::CellCycleScoring(
    seu,
    s.features = Seurat::cc.genes$s.genes,
    g2m.features = Seurat::cc.genes$g2m.genes,
    set.indent = T
  )

  # plot pca
  cat("Plotting principle components...\n")
  plots[["pca_by_sample"]] <- Seurat::DimPlot(seu, reduction = "pca", group.by = "sample", raster = F)
  plots[["pca_dim_loadings"]] <- Seurat::VizDimLoadings(seu, dims = 1:4, reduction = "pca")
  plots[["pca_dim_heatmap"]] <- Seurat::DimHeatmap(seu, dims = 1:15, cells = 500, balanced = TRUE, raster = F, fast = F)

  # elbow plot
  cat("Plotting elbow plot...\n")
  plots[["pca_elbow"]] <- Seurat::ElbowPlot(seu)

  # NON-LINEAR DIMENSIONALITY REDUCTION (UMAP) ----

  # run and plot umap
  seu <- Seurat::RunUMAP(seu, dims = 1:n_dims)
  plots[["umap_by_sample"]] <- Seurat::DimPlot(seu, reduction = "umap", group.by = "sample", raster = F)

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

  # ANNOTATION ----

  # clustering of different sample/patient/batch-level features in different reductions
  purrr::walk(groupings, function(grouping) {
    cat("\tPlotting PCA and UMAP vs", grouping, "\n")
    purrr::walk(c("pca", "umap"), function(redu) {
      # use `<<` for global assignment
      if (is.numeric(seu@meta.data[,grouping])) {
        p <- Seurat::FeaturePlot(object = seu,
                            reduction = redu,
                            features = grouping,
                            raster = F)
      } else {
        p <- seu %>%
          Seurat::DimPlot(reduction = redu,
                          group.by = grouping,
                          raster = F)
      }
      plots[[paste0(redu, "_vs_", grouping)]] <<- p
    })
  })

  # save plots as pdf and list
  cat("Saving all plots...\n")
  if(length(dev.list()) != 0) { dev.off() }
  plots_grob <- marrangeGrob(grobs = plots, nrow = 1, ncol = 1)
  ggsave(paste0(out$base, "/qc_report_plots.pdf"), plots_grob,
         width = 20, height = 20, units = "cm")
  if(length(dev.list()) != 0) { dev.off() }
  saveRDS(plots, file = paste0(out$base, "/qc_report_plots.rds"))

  cat("DONE!\n\n")
}
