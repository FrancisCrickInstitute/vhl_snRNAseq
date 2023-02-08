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

  # capture arguments
  args <- as.list(environment())

  # testing: # library(devtools);load_all();experiment="221202_A01366_0326_AHHTTWDMXY";genome="hg38";sublibrary="SHE5052A9_S101";parse_analysis_subdir="all-well/DGE_filtered";parse_pipeline_dir = paste0(get_base_dir(), "/parse_pipeline/");n_dims = 50;out_dir = NULL;sample_subset = NULL;do_timestamp = F;do_integration = T;integration_col="sample";
  # testing: # experiment="230127_A01366_0343_AHGNCVDMXY";genome="hg38";sublibrary="SHE5052A11_S164";do_integration=F
  # testing: # setwd(paste0(ifelse(Sys.info()["nodename"]=="Alexs-MacBook-Air-2.local","/Volumes/TracerX/working/VHL_GERMLINE/tidda/","/camp/project/tracerX/working/VHL_GERMLINE/tidda/"),"vhl/"));library(devtools);load_all();args <- dget("out/230127_A01366_0343_AHGNCVDMXY/hg38/comb/all-well/DGE_filtered/args_for_generate_qc_report.R") ; list2env(args,globalenv()); parse_pipeline_dir=paste0(get_base_dir(), "/parse_pipeline/")

  # cell quality control
  # -> unusually high transcript/gene counts indicate multiplets
  # -> unusually low transcript/gene counts indicate barcoding of cells with
  #    damaged membranes (uninformative)
  # -> high % mito genes indicates death / loss of cytoplasmic RNA / increased
  #    apoptosis (for scRNA-seq, not snRNA-seq)
  # -> high % globin and low % ribo suggests erythrocytes

  # gridExtra must be loaded in the environment
  library(gridExtra)

  # set directories
  data_dir <- paste(parse_pipeline_dir, "analysis", experiment, genome, sublibrary, sep = "/")
  dge_dir <- paste(data_dir, parse_analysis_subdir, sep = "/")
  out <- define_out(experiment, genome, sublibrary, parse_analysis_subdir,
                    out_dir,
                    do_timestamp, do_integration)
  plots <- list()

  # save arguments
  dput(args, file = paste0(out$base, "/args_for_generate_qc_report.R"))

  # LOAD DATA ----

  # load parse output to seurat object with no cut-offs, remove NAs, save to outdir
  seu <- load_parse_to_seurat(
    dge_dir, min_genes_per_cell = 0, min_cells_per_gene = 0, sample_subset
  )
  # remove NA samples
  seu <- subset(seu, subset = sample %in% seu$sample[!is.na(seu$sample)])
  saveRDS(seu, file = paste0(out$base, "/seu.rds"))

  # sample groupings to check at clustering stage
  groupings <- c("sample", "percent_mito", "percent_ribo", "percent_globin") %>%
    # if genome is human, do cell cycle scoring
    { if(grepl("hg38", genome)) c(., "Phase") else . }

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
      readr::read_tsv(show_col_types = FALSE) %>%
      janitor::clean_names() %>%
      # add patient x sample type label
      dplyr::mutate(
        sample_type_code = dplyr::case_when(sample_type == "renal_cyst" ~ "c",
                                            sample_type == "solid" ~ "t",
                                            sample_type == "normal_renal" ~ "n",
                                            sample_type == "metastasis" ~ "met",
                                            sample_type == "mixed" ~ "mix",
                                            TRUE ~ sample_type),
        sample_label = paste0(sample, "_", sample_type_code))

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
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90),
      strip.text.y = ggplot2::element_text(angle = 0)) 

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
    ggplot2::facet_wrap(statistic ~ label, scales = "free")

  # get proportions of relevant transcript types
  cat("Calculating abundance of different transcript types...\n")
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
  plots[["rna_molecules_vs_rna_transcripts_scatter"]] <- dittoSeq::dittoScatterPlot(
    seu, "nCount_RNA", "nFeature_RNA", "sample", raster = F)
  plots[["rna_molecules_vs_pct_mito_scatter"]] <- dittoSeq::dittoScatterPlot(
    seu, "nCount_RNA", "percent_mito", "sample", raster = F)
  plots[["pct_globin_vs_ribo_scatter"]] <- dittoSeq::dittoScatterPlot(
    seu, "percent_globin", "percent_ribo", "sample", raster = F)

  # highest expressed genes
  # cat("Plotting top 20 most highly expressed genes...\n")
  # plots[["top_genes_boxplot"]] <- boxplot_top_genes(seu)
  # TODO: fix this - it kills the job, too memory intensive

  # highest variable genes (HVGs)
  cat("Plotting most highly variable genes (HVGs)...\n")
  seu <- Seurat::FindVariableFeatures(seu)
  plots[["hvg_scatter"]] <- Seurat::LabelPoints(
    plot = Seurat::VariableFeaturePlot(seu, raster = F),
    points = head(Seurat::VariableFeatures(seu), 10),
    repel = T
  ) +
    ggplot2::theme(legend.position = "top")
  plots[["hvg_heatmap"]] <- dittoSeq::dittoHeatmap(
    seu, genes = head(Seurat::VariableFeatures(seu), 20),  
    annot.by = "sample", scaled.to.max = T)
  
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

  # assign cell cycle phase scores
  cat("Assigning cell cycle phases...\n")
  seu <- Seurat::CellCycleScoring(
    seu,
    s.features = Seurat::cc.genes$s.genes,
    g2m.features = Seurat::cc.genes$g2m.genes
  )

  # LINEAR DIMENSIONALITY REDUCTION (PCA) ----

  # run PCA
  cat("Running PCA...\n")
  seu <- Seurat::RunPCA(seu)

  # plot pca dims
  cat("Plotting principle components...\n")
  plots[["pca_dim_loadings"]] <- Seurat::VizDimLoadings(seu, dims = 1:4, reduction = "pca")
  plots[["pca_dim_heatmap"]] <- Seurat::DimHeatmap(seu, dims = 1:15, cells = 500, balanced = TRUE, raster = F, fast = F)

  # elbow plot
  cat("Plotting elbow plot...\n")
  plots[["pca_elbow"]] <- Seurat::ElbowPlot(seu, ndims = 50)

  # NON-LINEAR DIMENSIONALITY REDUCTION (UMAP) ----

  # run umap
  seu <- Seurat::RunUMAP(seu, dims = 1:n_dims)

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

  # clustering of different sample/patient/batch-level features in different reductions
  plots[["umap_vs_sample_label"]] <- dittoSeq::dittoDimPlot(
    seu, "sample_label", reduction.use = "umap", raster = F)
  # plot all reduction / grouping projections, side-by-side with sample annots for comparison
  purrr::walk(groupings, function(grouping) {
    cat("\tPlotting PCA and UMAP vs", grouping, "\n")
    purrr::walk(c("pca", "umap"), function(redu) {
      # use `<<` for global assignment
      plots[[paste0(redu, "_vs_", grouping)]] <<- 
        dittoSeq::dittoDimPlot(seu, grouping, reduction.use = redu, raster = F) +
        plots[["umap_vs_sample_label"]] 
    })
  })
  
  # ANNOTATION ----
  
  # gene modules scoring
  gene_modules %>% names %>%
    purrr::walk(function(module) {
      cat(module, "\n")
      seu <<- Seurat::AddModuleScore(
        seu,
        features = list(gene_modules[[module]]),
        name = module,
        nbin = 10
      )
      # umap
      plots[[paste0("umap_vs_", module, "_module_score")]] <<-
        dittoSeq::dittoDimPlot(
          seu, paste0(module, "1"),
          main = paste0(module, " module (n=", length(gene_modules[[module]]), ")")) +
        plots[["umap_vs_sample_label"]] 
      # ridgeplot
      plots[[paste0(module, "_module_score_vs_sample_type_ridge")]] <<- 
        dittoSeq::dittoRidgePlot(seu, paste0(module, "1"), group.by = "sample_type")
    })
  
  # celldex annotations from the human primary cell atlas
  human_primary_ref <- celldex::HumanPrimaryCellAtlasData()
  seu_SingleR <- SingleR::SingleR(
    test = Seurat::GetAssayData(seu, slot = "data"),
    ref = human_primary_ref,
    labels = human_primary_ref$label.main)
  plots[["singleR_annots_heatmap"]] <- SingleR::plotScoreHeatmap(seu_SingleR)
  plots[["singleR_annots_delta_dist"]] <- SingleR::plotDeltaDistribution(seu_SingleR)
  
  # add labels
  seu$SingleR_annot <- seu_SingleR %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(labels) %>%
    dplyr::transmute(
      n = dplyr::n(),
      SingleR_annot = replace(labels, n < 10, "none")
    ) %>%
    dplyr::pull(SingleR_annot)
    
  plots[["umap_vs_singleR_annot"]] <- 
    dittoSeq::dittoDimPlot(seu, "SingleR_annot", size = 0.7) + 
    plots[["umap_vs_sample_type"]] 
  plots[["singleR_annot_bar"]] <- dittoSeq::dittoBarPlot(
    seu, var = "SingleR_annot", group.by = "sample", 
    split.by = c("sample_type"), 
    split.adjust = list(scales = "free"),
    split.nrow = 1)
  
  # save annotated seu
  saveRDS(seu, file = paste0(out$base, "/seu_annotated.rds"))

  # save plots as pdf and list
  cat("Saving all plots to", paste0(out$base, "/qc_report_plots.pdf"), "...\n")
  if(length(dev.list()) != 0) { dev.off() }
  plots_grob <- marrangeGrob(grobs = plots, nrow = 1, ncol = 1)
  ggsave(paste0(out$base, "/qc_report_plots.pdf"), plots_grob,
         width = 20, height = 20, units = "cm")
  if(length(dev.list()) != 0) { dev.off() }
  saveRDS(plots, file = paste0(out$base, "/qc_report_plots.rds"))

  # TODO: GENE ONTOLOGY ----
  # 
  # # choose a cluster resolution
  # clust_res <- 0.3
  # deg_list <- split(rownames(seu), seu[, paste0(snn_res_prefixes, clust_res)])
  
  cat("\nDONE!\n\n")

}
