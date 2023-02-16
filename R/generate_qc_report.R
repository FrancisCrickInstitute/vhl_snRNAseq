# min_nFeature_RNA - passed to SeuratObject::CreateAssayObject() min.features argument
generate_qc_report <- function(experiment,
                               parse_pipeline_dir = "/camp/project/tracerX/working/VHL_GERMLINE/tidda/parse_pipeline/",
                               genome = "hg38",
                               sublibrary = "comb",
                               parse_analysis_subdir = "/all-well/DGE_unfiltered/",
                               do_filtering = T,
                               remove_doublets = T,
                               min_cells_per_gene = 5,
                               max_nCount_RNA = NULL,
                               min_nFeature_RNA = NULL,
                               max_nFeature_RNA = NULL,
                               max_percent_mito = NULL,
                               n_dims = NULL,
                               clustering_resolutions = seq(0.1, 0.8, by = 0.1),
                               final_clustering_resolution = NULL,
                               out_dir = NULL,
                               sample_subset = NULL,
                               do_timestamp = F,
                               do_integration = F,
                               integration_col = "sample") {
  # capture arguments
  args <- as.list(environment())

  # test runs
  # testing: setwd # setwd(paste0(ifelse(Sys.info()["nodename"]=="Alexs-MacBook-Air-2.local","/Volumes/TracerX/working/VHL_GERMLINE/tidda/","/camp/project/tracerX/working/VHL_GERMLINE/tidda/"),"vhl/"));
  # testing: pilot # library(devtools);load_all();experiment="221202_A01366_0326_AHHTTWDMXY";genome="hg38";sublibrary="SHE5052A9_S101";parse_analysis_subdir="all-well/DGE_filtered";parse_pipeline_dir = paste0(get_base_dir(), "/parse_pipeline/");do_filtering=T;remove_doublets = T;min_cells_per_gene=NULL;min_nFeature_RNA=NULL;max_nCount_RNA=NULL;max_nFeature_RNA=NULL;max_percent_mito=NULL;n_dims=NULL;clustering_resolutions = seq(0.1, 0.8, by = 0.1);final_clustering_resolution=NULL;out_dir = NULL;sample_subset = NULL;do_timestamp = F;do_integration = F;integration_col="sample";
  # testing: 2 SLs # experiment="230127_A01366_0343_AHGNCVDMXY";genome="hg38";sublibrary="comb"
  # testing: 8 SLs # experiment="230210_A01366_0351_AHNHCFDSX5";genome="hg38";sublibrary="comb"
  # testing: args  # library(devtools);load_all();args <- dget("out/230127_A01366_0343_AHGNCVDMXY/hg38/comb/all-well/DGE_filtered/args_for_generate_qc_report.R") ; list2env(args,globalenv()); parse_pipeline_dir=paste0(get_base_dir(), "/parse_pipeline/")
  # testing: archived # data_dir="/Volumes/TracerX/working/VHL_GERMLINE/tidda//parse_pipeline/archive/230209////analysis/221202_A01366_0326_AHHTTWDMXY/hg38/SHE5052A9_S101";dge_dir="/Volumes/TracerX/working/VHL_GERMLINE/tidda//parse_pipeline/archive/230209////analysis/221202_A01366_0326_AHHTTWDMXY/hg38/SHE5052A9_S101/all-well/DGE_filtered";sample_metadata$sample<-c("N090_V1027"  ,  "N090_V1024A" ,  "K1026_T2D_CL" ,"Mouse_nuclei")

  # gridExtra and clustree must be loaded in the environment due to bugs
  library(clustree)
  library(gridExtra)

  # set ggplot theme, initiate plots list
  plots <- list()
  ggplot2::theme_set(ggplot2::theme_bw())
  ditto_colours <- list(scale_fill_manual(values = dittoSeq::dittoColors()),
                        scale_colour_manual(values = dittoSeq::dittoColors()))

  # set directories
  data_dir <- paste(parse_pipeline_dir, "analysis", experiment, genome, sublibrary, sep = "/")
  dge_dir <- paste(data_dir, parse_analysis_subdir, sep = "/")
  out <- define_out(experiment, genome, sublibrary, parse_analysis_subdir, out_dir,
                    do_filtering, do_integration, do_timestamp)

  # save arguments
  if("args" %in% ls()) { dput(args, paste0(out$base, "/args_for_generate_qc_report.R")) }

  # LOAD DATA ----

  # load parse pipeline output to seurat object with no cut-offs, remove sample
  # NAs, add sample_metadata, optionally subset samples
  seu <- load_parse_to_seurat(
    dge_dir,
    min_nFeature_RNA = 0,
    min_cells_per_gene = 0,
    sample_subset,
    remove_na_samples = T,
    do_add_sample_metadata = T,
    parse_pipeline_dir,
    experiment
  )
  # save seu object
  saveRDS(seu, file = paste0(out$base, "/seu.rds"))

  # sample groupings to check at clustering stage
  md_groupings <- c("date_prep", "patient_id", "rin", "sample_type", "size")
  groupings <- c("sample", "percent_mito", "percent_ribo", "percent_globin", "doublet",
                 md_groupings) %>%
    # if genome is human, do cell cycle scoring (doesn't work with other genomes)
    { if (grepl("hg38", genome)) c(., "Phase") else . }

  # sample-level summary stats of interest
  summary_stats <- get_summary_stats(data_dir) %>%
    # append numeric sample metadata to summary stats
    dplyr::bind_rows(
      seu@misc$sample_metadata %>%
        tidyr::pivot_longer(
          cols = tidyr::any_of(md_groupings) & where(is.numeric),
          names_to = "statistic"
        ) %>% dplyr::select(statistic, sample, value)
    )

  # plot summary stats
  cat("Plotting summary stats...\n")
  plots[["run_summary_stats"]] <- summary_stats %>%
    ggplot2::ggplot(ggplot2::aes(x = sample, y = value)) +
    ggplot2::geom_col() +
    ggplot2::facet_grid(statistic ~ ., scales = "free") +
    ggplot2::labs(title = "per-sample summary statistics")

  # CELL AND GENE SUMMARY STATISTICS ----
  # gene-level (min_cells_per_gene)
  # -> genes that are present in very few cells of the dataset (<3) are uninformative
  #    and unlikely to play any part in differentiating groups of cells (in general,
  #    most genes removed by this filtering will be those not detected in any cell)
  # cell-level (nCount_RNA, nFeature_RNA, percent_mito, doublet)
  # -> unusually high transcript/gene counts indicate multiplets
  # -> unusually low transcript/gene counts indicate barcoding of cells with
  #    damaged membranes (uninformative)
  # -> high % of mito genes indicates death / loss of cytoplasmic RNA / increased
  #    apoptosis (for scRNA-seq, not snRNA-seq)
  # -> high % of globin and low % of ribo suggests erythrocytes
  # -> high features:counts per cell ratio could be dying cells

  # identify doublets (by creating artificial doublets and looking at their clustering)
  cat("Annotating suspected doublets...\n")
  seu <- annotate_doublets(seu)

  # get proportions of relevant transcript types
  cat("Annotating abundance of different transcript types...\n")
  seu <- annotate_proportions_of_transcript_types(seu)

  # get user-defined filters or designate filters
  # -> nCount/Feature_RNA, percent_mito = median plus 5 x the median absolute deviation
  #    (recommended by https://romanhaa.github.io/projects/scrnaseq_workflow)
  # -> min_nFeature_RNA/min_cells_per_gene = defaults recommended by Parse
  filters <- get_filters(seu,
                         do_filtering,
                         remove_doublets,
                         max_nCount_RNA,
                         min_nFeature_RNA,
                         max_nFeature_RNA,
                         max_percent_mito)

  # create misc table of included/excluded cells, with exclude criteria
  seu@misc$cell_filtering <-
    seu@meta.data %>%
    dplyr::as_tibble(rownames = "cell") %>%
    dplyr::select(cell, sample, tidyr::any_of(names(filters))) %>%
    tidyr::pivot_longer(-c(cell, sample)) %>%
    dplyr::left_join(filters %>%
                       purrr::map(dplyr::as_tibble) %>%
                       dplyr::bind_rows(.id = "name")) %>%
    dplyr::mutate(include = value >= min & value <= max,
                  exclude_criteria = ifelse(include == F, name, NA)) %>%
    dplyr::group_by(cell) %>%
    dplyr::mutate(include_cell = all(include),
                  exclude_criteria_cell = exclude_criteria %>% unique %>% na.omit %>%
                    paste(collapse = ",") %>% dplyr::na_if(., "")
    )

  # create misc table of included/excluded genes
  seu@misc$gene_filtering <- tibble::tibble(
    gene = rownames(seu),
    n_transcripts = Matrix::rowSums(seu),
    n_cells = rowSums(as.matrix(seu@assays$RNA@counts) > 0)
  ) %>%
    dplyr::mutate(min = min_cells_per_gene,
                  include_gene = n_cells >= min)

  # cells per sample
  cat("Plotting cells per sample...\n")
  plots[["cells_per_sample_bar"]] <-
    dplyr::tibble(sample_id = names(table(seu$sample)),
                  n_cells = table(seu$sample)) %>%
    ggplot2::ggplot(ggplot2::aes(x = sample_id, y = n_cells)) +
    ggplot2::geom_col() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

  # plot filters on distributions
  plots[["filters_vln"]] <-
    seu@misc$cell_filtering %>%
    dplyr::group_by(name) %>%
    dplyr::mutate(name = paste0(name,
                                "\nmin=", min, ", max=", max %>% {formatC(signif(., digits=5), digits=5, format="fg", flag="#")},
                                "\nexcluded=", sum(!include), ", included=", sum(include))) %>%
    ggplot2::ggplot(ggplot2::aes(x = sample, y = value)) +
    ggplot2::geom_jitter(ggplot2::aes(colour = include),
                         height = 0, alpha = 0.5, size = 0.6) +
    ggplot2::geom_violin(ggplot2::aes(fill = sample), alpha = 0.8,
                         draw_quantiles = 0.5) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = min), colour = "red") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = max), colour = "red") +
    ggplot2::facet_wrap(~ name, scales = "free") +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = element_text(angle = 90)) +
    ditto_colours +
    ggplot2::scale_colour_manual(values = c("red", "black"))

  # plot doublets
  plots[["doublets_per_sample"]] <-
    seu@misc$cell_filtering %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(n_doublets = sum(name == "doublet" & value == 1)) %>%
    dplyr::group_by(cell) %>%
    dplyr::mutate(doublet = any(name == "doublet" & value == 1)) %>%
    ggplot2::ggplot(ggplot2::aes(y = value, x = sample, colour = doublet)) +
    ggplot2::geom_jitter(data = . %>% dplyr::filter(!doublet), height = 0) +
    ggplot2::geom_jitter(data = . %>% dplyr::filter(doublet), height = 0) +
    ggplot2::labs(title = "detected doublets in each sample (using scDblFinder)") +
    ggplot2::facet_wrap(~ name, scales = "free") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +
    ditto_colours +
    seu@misc$cell_filtering %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(n_singlets = sum(name == "doublet" & value == 0),
                  n_doublets = sum(name == "doublet" & value == 1)) %>%
    dplyr::distinct(sample, n_doublets, n_singlets) %>%
    tidyr::pivot_longer(-sample) %>%
    ggplot2::ggplot(ggplot2::aes(x = sample, y = value, fill = name)) +
    ggplot2::geom_col() +
    ggplot2::theme(axis.text.x = element_text(angle = 90)) +
    ditto_colours

  # plot transcript type abundances
  cat("Plotting abundance of different transcript types...\n")
  plots[["filters_scatter"]] <-
    plot_cell_scatter_with_filters(seu, "nCount_RNA", "percent_mito", filters) +
    plot_cell_scatter_with_filters(seu, "nCount_RNA", "nFeature_RNA", filters)

  # filtering
  if(do_filtering == T) {

    # perform cell filtering
    cat("Filtering cells...\n")
    seu <- subset(seu, cells = unique(dplyr::filter(seu@misc$cell_filtering, include_cell)$cell))

    # perform gene filtering
    cat("Filtering genes...\n")
    seu <- subset(seu, features = unique(dplyr::filter(seu@misc$gene_filtering, include_gene)$gene))

  }

  # highest expressed genes
  # cat("Plotting top 20 most highly expressed genes...\n")
  # plots[["top_genes_boxplot"]] <- boxplot_top_genes(seu)
  # TODO: fix this - it kills the job, too memory intensive

  # highest variable genes (HVGs)
  cat("Finding most highly variable genes (HVGs)...\n")
  seu <- Seurat::FindVariableFeatures(seu)

  # plot hvgs
  cat("Plotting HVGs...\n")
  plots[["hvg_scatter"]] <- Seurat::LabelPoints(
    plot = Seurat::VariableFeaturePlot(seu, raster = F),
    points = head(Seurat::VariableFeatures(seu), 10),
    repel = T
  ) +
    ggplot2::theme(legend.position = "top")
  plots[["hvg_heatmap"]] <- dittoSeq::dittoHeatmap(
    seu,
    genes = head(Seurat::VariableFeatures(seu), 20),
    annot.by = "sample",
    scaled.to.max = T
  )

  if (do_integration == T) {
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
      Seurat::FindIntegrationAnchors(anchor.features = Seurat::SelectIntegrationFeatures(seu_list))

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

  # choose dimensionality (final_n_dims) of the dataset
  if(!is.null(n_dims)) {

    # n_dims passed by the user
    final_n_dims <- n_dims

  } else {

    # test dimensionality of the dataset using intrinsicDimension package
    n_dims_iD <- intrinsicDimension::maxLikGlobalDimEst(
      seu@reductions$pca@cell.embeddings,
      k = 10
    )$dim.est %>% round(0)

    # better to be more generous with the number of dimensions
    # therefore double dimension cutoff suggested by intrinsicDimension
    generous_n_dims <- round(2 * n_dims_iD, 0)
    final_n_dims <- generous_n_dims

  }

  # elbow plot
  cat("Plotting elbow plot...\n")
  plots[["pca_elbow"]] <- Seurat::ElbowPlot(seu, ndims = max(50, final_n_dims)) &
    ggplot2::geom_vline(ggplot2::aes(xintercept = n_dims_iD, colour = "intrisicDimensions"),
                        linetype = "dashed") &
    ggplot2::geom_vline(ggplot2::aes(xintercept = generous_n_dims, colour = "generous (iD x 2)"),
                        linetype = "dashed") &
    ggplot2::scale_colour_manual(
      name = "n_dims cut-off",
      values = c(user_provided = dittoSeq::dittoColors()[[1]],
                 intrisicDimensions = dittoSeq::dittoColors()[[2]],
                 `generous (iD x 2)` = dittoSeq::dittoColors()[[3]])
    ) &
    ggplot2::theme(legend.position = c(0.7, 0.9),
                   legend.background = element_rect(fill = "white"))
  if(!is.null(n_dims)) {
    plots[["pca_elbow"]] <- plots[["pca_elbow"]] &
      ggplot2::geom_vline(ggplot2::aes(xintercept = n_dims, colour = "user_provided"),
                          linetype = "dashed")
  }

  # plot pca dims
  cat("Plotting principle components...\n")
  plots[["pca_dim_loadings"]] <-
    Seurat::VizDimLoadings(seu, dims = 1:4, reduction = "pca")
  plots[["pca_dim_heatmap"]] <-
    Seurat::DimHeatmap(
      seu,
      dims = 1:final_n_dims,
      cells = 500,
      balanced = TRUE,
      raster = F,
      fast = F
    )

  # NON-LINEAR DIMENSIONALITY REDUCTION (UMAP, t-SNE) ----

  # run t-SNE and UMAP
  seu <- Seurat::RunTSNE(seu)
  seu <- Seurat::RunUMAP(seu, dims = 1:final_n_dims)

  # plot projection of different sample features in different reductions
  plots[["umap_vs_sample_label"]] <- dittoSeq::dittoDimPlot(seu,
                                                            "sample_label",
                                                            reduction.use = "umap",
                                                            raster = F)

  # split out sample clusters, side-by-side
  plots[["umap_split_by_patient"]] <-
    dittoSeq::dittoDimPlot(seu,
                           "sample_label",
                           reduction.use = "umap",
                           split.by = "patient_id")

  # plot all reduction / grouping projections, side-by-side with sample annots for comparison
  groupings_vs_reductions <- list()
  purrr::walk(groupings, function(grouping) {
    cat("Plotting reductions vs", grouping, "\n")
    purrr::walk(Seurat::Reductions(seu), function(redu) {
      title <- paste0(redu, "_vs_", grouping)
      p <- dittoSeq::dittoDimPlot(seu,
                                  grouping,
                                  reduction.use = redu,
                                  size = 0.5,
                                  raster = F,
                                  show.axes.numbers = F,
                                  show.grid.lines = F,
                                  main = gsub("_", " ", title))
      groupings_vs_reductions[[paste0(grouping, "_legend")]] <<- lemon::g_legend(p)
      groupings_vs_reductions[[title]] <<- p + ggplot2::theme(legend.position = "none")
    })
  })

  # convert to grob
  groupings_vs_reductions_grob <- marrangeGrob(
    grobs = groupings_vs_reductions,
    nrow = length(groupings),
    ncol = length(Seurat::Reductions(seu)),
    layout_matrix = matrix(1:length(groupings_vs_reductions),
                           length(groupings),
                           length(Seurat::Reductions(seu)) + 1,
                           TRUE)
  )

  # CLUSTERING ----

  # cluster data at different resolutions
  cat("Clustering data at different resolutions...\n")

  # compute k.param nearest neighbours based on euclidean distance in PCA space
  # and check shared nearest neighbours between any 2 cells (Jaccard similarity)
  seu <- Seurat::FindNeighbors(seu, dims = 1:final_n_dims)

  # identify clusters of cells by a shared nearest neighbours modularity
  # optimisation based clustering algorithm
  seu <- Seurat::FindClusters(seu, resolution = clustering_resolutions)
  saveRDS(seu, file = paste0(out$base, "/seu_transformed_and_clustered.rds"))

  # cluster tree plot of increasing resolutions
  cat("Plotting clustering tree at different resolutions...\n")
  snn_res_prefixes <- paste0(Seurat::DefaultAssay(seu), "_snn_res.")
  plots[["clustering_tree"]] <- clustree::clustree(
    seu@meta.data[, grep(snn_res_prefixes, colnames(seu@meta.data))],
    prefix = snn_res_prefixes)

  # plot clustering of reductions at different resolutions
  cat("Plotting reductions of clusters at different resolutions...\n")
  clustered_reductions <- list()
  purrr::walk(clustering_resolutions, function(res) {
    purrr::walk(Seurat::Reductions(seu), function(redu) {
      title <- paste0(redu, "_by_cluster_res_", res)
      clustered_reductions[[title]] <<-
        dittoSeq::dittoDimPlot(
          seu,
          paste0(Seurat::DefaultAssay(seu), "_snn_res.", res),
          reduction.use = redu,
          size = 0.5,
          raster = F,
          legend.show = F,
          show.axes.numbers = F,
          show.grid.lines = F,
          main = gsub("_", " ", title)
        )
    })
  })

  # convert to grob
  clustered_reductions_grob <- marrangeGrob(
    grobs = clustered_reductions,
    nrow = length(clustering_resolutions),
    ncol = length(Seurat::Reductions(seu)),
    layout_matrix = matrix(1:length(clustered_reductions),
                           length(clustering_resolutions),
                           length(Seurat::Reductions(seu)),
                           TRUE)
  )

  # choose a resolution
  if(!is.null(final_clustering_resolution)) { # final_clustering_resolution=0.3

    # define clusters at final resolution
    seu$cluster <- seu@meta.data[, paste0(snn_res_prefixes, final_clustering_resolution)]

    # plot sample/cluster proportions
    samples_by_clusters_bar <-
      dittoSeq::dittoBarPlot(seu,
                             var = "cluster",
                             scale = "count",
                             group.by = "sample",
                             ) +
      ggplot2::theme(legend.position = "left") +
      dittoSeq::dittoBarPlot(seu,
                             var = "sample",
                             scale = "count",
                             group.by = "cluster") +
      patchwork::plot_layout(ncol = 2, widths = c(
        seu$sample %>% unique() %>% length(),
        seu$cluster %>% unique() %>% length()
      ))

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
          dittoSeq::dittoDimPlot(seu,
                                 paste0(module, "1"),
                                 main = paste0(module, " module (n=", length(gene_modules[[module]]), ")")) +
          plots[["umap_vs_sample_label"]]
        # ridgeplot
        plots[[paste0(module, "_module_score_vs_sample_type_ridge")]] <<-
          dittoSeq::dittoRidgePlot(seu, paste0(module, "1"), group.by = "sample_type")
      })

    # celldex annotations from the human primary cell atlas
    human_primary_ref <- celldex::HumanPrimaryCellAtlasData()
    seu_singler <- SingleR::SingleR(
      test = Seurat::GetAssayData(seu, slot = "data"),
      ref = human_primary_ref,
      labels = human_primary_ref$label.main
    )

    # plot heatmap of annotation scores
    plots[["singler_annots_heatmap"]] <-
      SingleR::plotScoreHeatmap(
        seu_singler,
        annotation_col = data.frame(
          patient = seu$patient_id,
          sample_type = seu$sample_type
        )
      )

    # plot delta of cell type annotations
    plots[["singler_annots_delta_dist"]] <-
      SingleR::plotDeltaDistribution(seu_singler)

    # add cell type annotations
    seu$singler_annot <- seu_singler %>%
      dplyr::as_tibble() %>%
      dplyr::group_by(labels) %>%
      # (remove cell types with less than 10 cells)
      dplyr::transmute(n = dplyr::n(),
                       singler_annot = replace(labels, n < 10, "none")) %>%
      dplyr::pull(singler_annot)

    # add cell type annotation labels for plotting
    seu$singler_annot_label <- seu@meta.data %>%
      dplyr::group_by(cluster, singler_annot) %>%
      dplyr::count() %>%
      dplyr::group_by(cluster) %>%
      dplyr::mutate(total = sum(n),
                    prop = round(n / sum(n) * 100, 0),
                    label = paste0(singler_annot, " (", prop, "%)")) %>%
      dplyr::filter(prop == max(prop) & prop > 50 | max(prop) < 50 & rank(-prop) <= 2) %>%
      dplyr::summarise(label = paste0("cluster ", unique(cluster), " (n = ", unique(total), ")\n", paste(label, collapse = "\n"))) %>%
      {dplyr::left_join(seu@meta.data[, "cluster", drop = F], multiple = "all", .)} %>%
      dplyr::pull(label)

    # plot cell types against reduction
    plots[["umap_vs_singler_annot"]] <-
      dittoSeq::dittoDimPlot(seu, "singler_annot", size = 0.7) +
      ggplot2::geom_label(
        data = seu@meta.data %>%
          dplyr::right_join(get_centroids(seu, "umap", cluster)) %>%
          dplyr::distinct(cluster, singler_annot_label, x, y),
        ggplot2::aes(x, y, label = singler_annot_label),
        size = 3, fill = 'white', color = 'black', alpha = 0.5, label.size = 0,
        show.legend = FALSE, fontface = "bold",
        vjust = "inward", hjust = "inward"
      ) +
      ggplot2::labs(title = "majority celldex annotation(s) per cluster") +
      plots[["umap_vs_sample_label"]]

    # plot cell type composition of samples
    plots[["singler_annot_bar"]] <-
      dittoSeq::dittoBarPlot(
        seu,
        var = "singler_annot",
        group.by = "sample",
        scale = "count",
        legend.show = F
      ) +
      dittoSeq::dittoBarPlot(
        seu,
        var = "singler_annot",
        group.by = "cluster",
        scale = "count"
      ) +
      patchwork::plot_layout(
        ncol = 2,
        widths = c(seu$sample %>% dplyr::n_distinct(),
                   seu$cluster %>% dplyr::n_distinct())
      )

    # alluvial plots
    plots[["alluvial_sample_cluster_singler"]] <-
      seu@meta.data %>%
      dplyr::group_by(sample, cluster, singler_annot) %>%
      dplyr::count() %>%
      dplyr::ungroup() %>%
      ggforce::gather_set_data(1:3) %>%
      dplyr::mutate(x = dplyr::case_when(x == 1 ~ "sample",
                                         x == 2 ~ "cluster",
                                         x == 3 ~ "singler_annotation") %>%
                      factor(levels = c("sample", "cluster", "singler_annotation"))) %>%
      ggplot2::ggplot(ggplot2::aes(
                        x = x,
                        id = id,
                        split = y,
                        value = n
                      )) +
      ggforce::geom_parallel_sets(ggplot2::aes(fill = cluster),
                                  axis.width = 0.15,
                                  alpha = 0.75) +
      ggforce::geom_parallel_sets_axes(fill = "lightgrey",
                                       axis.width = 0.15) +
      ggforce::geom_parallel_sets_labels(angle = 0,
                                         nudge_x = c(
                                           rep(-0.3, dplyr::n_distinct(seu$sample)),
                                           rep(0, dplyr::n_distinct(seu$cluster)),
                                           rep(0.3, dplyr::n_distinct(seu$singler_annot))
                                         )) +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "none",
                     axis.text.x = element_text(face = "bold", size = 15)) +
      ditto_colours

    # save annotated seu
    saveRDS(seu, file = paste0(out$base, "/seu_annotated.rds"))

    # TODO: GENE ONTOLOGY ----

    # set identities based on the final clustering resolution
    seu <- Seurat::SetIdent(
      seu, value = seu@meta.data[,paste0(snn_res_prefixes, final_clustering_resolution)]
    )

    # reset default assay to RNA
    Seurat::DefaultAssay(seu) <- "RNA"

    # pathway analysis
    gsva_result <- ReactomeGSA::analyse_sc_clusters(seu, verbose = T)

    # plot heatmap of pathway expression values
    # (ranks pathways based on their expression difference)
    ReactomeGSA::plot_gsva_heatmap(gsva_result, max_pathways = 15)


  }

  # literature markers
  # # checking that all markers have been found
  # unmatched_markers <- unlist(literature_markers)[
  #   !(unlist(literature_markers) %in% unique(rownames(seu@assays$RNA)))
  # ]
  # if(length(unmatched_markers) > 0) {
  #   warning("Marker(s) ", paste(unmatched_markers, collapse = ", "), " not found!")
  # }

  # SAVE PLOTS ----
  # save plots as pdf and list
  cat("Saving all plots to", paste0(out$base, "/qc_report_plots.pdf"), "...\n")
  if (length(dev.list()) != 0) { dev.off() }
  plots_grob <- marrangeGrob(grobs = plots, nrow = 1, ncol = 1)
  ggsave(paste0(out$base, "/qc_report_plots.pdf"), plots_grob,
         width = 20, height = 20, units = "cm")
  ggsave(paste0(out$base, "/groupings_vs_reductions.pdf"), groupings_vs_reductions_grob,
  width = 20, height = 60, units = "cm")
  ggsave(paste0(out$base, "/clustered_reductions.pdf"), clustered_reductions_grob,
         width = 20, height = 60, units = "cm")
  if (length(dev.list()) != 0) { dev.off() }
  saveRDS(plots, file = paste0(out$base, "/qc_report_plots.rds"))
  saveRDS(groupings_vs_reductions_grob, file = paste0(out$base, "/groupings_vs_reductions_grob.rds"))
  saveRDS(clustered_reductions_grob, file = paste0(out$base, "/clustered_reductions_grob.rds"))
  # TODO: add clustered_reductions and groupings_vs_reductions to output

  cat("\nDONE!\n\n")

  return(seu)

}
