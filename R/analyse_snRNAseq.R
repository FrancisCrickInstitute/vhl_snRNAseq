#' analyse snRNAseq output from Parse Biosciences split-pipe
#'
#' @param parse_dir Path to Parse Biosciences split-pipe directory.
#' @param experiment The name of the experiment run via split-pipe.
#' @param genome The reference genome build used by split-pipe for alignment. Default is "hg38".
#' @param sublibrary The sublibrary to analyse. Default is "comb".
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
#' @param cluster_resolutions Optional. A vector of clustering resolutions to test. A higher resolution will result in a larger number of communities. Default is 0.1-0.8.
#' @param final_clustering_resolution Optional. The chosen clustering resolution of the dataset to use for downstream analysis. This is the resolution at which clusters appear to capture true biological groupings of interest in the dataset. This can be decided by consulting the clustering tree. If no value given, downstream analysis will not proceed.
#' @param out_dir Optional. Output directory. If no value given, the output will be saved to `out/{experiment}/{genome}/{sublibrary}/{parse_analysis_subdir}/{integrated,unintegrated}`.
#' @param sample_subset Optional. Vector of sample IDs to subset to.
#' @param do_timestamp If `TRUE`, will save the output to a time-stamped subdirectory (e.g. `20230217_105554/`).
#' @param do_integration If `TRUE`, will integrate the dataset and save the output to `integrated/` subdirectory. If `FALSE` (default), output will save to `unintegrated/` subdirectory.
#' @param integration_col If `do_integration` is `TRUE`, Seurat metadata column upon which to integrate the dataset.
#' @return A Seurat object.
#'
#' @export
analyse_snRNAseq <- function(parse_dir = "/camp/project/tracerX/working/VHL_GERMLINE/tidda/parse_pipeline/",
                             experiment,
                             genome = "hg38",
                             sublibrary = "comb",
                             parse_analysis_subdir = "all-well/DGE_unfiltered/",
                             do_filtering = T,
                             remove_doublets = F,
                             min_nuclei_per_gene = 5,
                             min_nCount_RNA = NULL,
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

  # check arguments
  check_analyse_snRNAseq_args(args)

  # for internal test runs
  # testing: setwd, default params, load_all # base_dir=ifelse(Sys.info()["nodename"]=="Alexs-MacBook-Air-2.local","/Volumes/TracerX/","/camp/project/tracerX/");setwd(paste0(base_dir,"working/VHL_GERMLINE/tidda/vhl/"));library(devtools);load_all();parse_dir=paste0(base_dir,"working/VHL_GERMLINE/tidda/parse_pipeline/");genome="hg38";sublibrary="comb";parse_analysis_subdir="all-well/DGE_filtered/";do_filtering=T;remove_doublets=F;min_nuclei_per_gene=5;min_nFeature_RNA=NULL;min_nCount_RNA=NULL;max_nCount_RNA=NULL;max_nFeature_RNA=NULL;max_percent_mito=NULL;n_dims=NULL;clustering_resolutions = seq(0.1, 0.8, by = 0.1);final_clustering_resolution=0.3;out_dir = NULL;sample_subset = NULL;do_timestamp = F;do_integration = F;integration_col="sample";
  # testing: pilot # experiment="221202_A01366_0326_AHHTTWDMXY";sublibrary="SHE5052A9_S101"
  # testing: 2 SLs # experiment="230127_A01366_0343_AHGNCVDMXY"
  # testing: 8 SLs # experiment="230210_A01366_0351_AHNHCFDSX5"
  # testing: args  # library(devtools);load_all(); args <- dget("out/230127_A01366_0343_AHGNCVDMXY/hg38/comb/all-well/DGE_filtered/args_for_generate_qc_report.R") ; list2env(args,globalenv())

  # gridExtra and clustree must be loaded in the environment (due to package-specific bugs)
  library(clustree)
  library(gridExtra)

  # initiate plots list, set ggplot2 theme
  plots <- list()
  ggplot2::theme_set(ggplot2::theme_bw())
  ditto_colours <- list(scale_fill_manual(values = dittoSeq::dittoColors()),
                        scale_colour_manual(values = dittoSeq::dittoColors()))

  # define output directory
  out <- {
    # if out_dir not given, use same output structure as in the parse analysis/ directory
    if (is.null(out_dir)) "out/" %>%
      paste(experiment, genome, sublibrary, parse_analysis_subdir, sep = "/") %>%
      { if (do_integration) paste0(., "/integrated/") else paste0(., "/unintegrated/")} %>%
      { if (do_timestamp) paste0(., format(Sys.time(), "%Y%m%d_%H%M%S"), "/") else . }
    # if out_dir is given, use out_dir
    else out_dir
    # clean path (// -> /)
  } %>% clean_path() %>% {
    # pre-set file names (TODO: define these once all outputs are finalised)
    purrr::map(list(base = ""), function(x) paste0(., x))
  }

  # create the output directory
  cat("Output will be saved to", out$base, "\n")
  dir.create(out$base, showWarnings = F, recursive = T)

  # save captured arguments
  if ("args" %in% ls()) { dput(args, paste0(out$base, "/args_for_generate_qc_report.R")) }

  # sample groupings to check at clustering stage
  groupings <- c("sample", "percent_mito", "percent_ribo", "percent_globin",
                 "date_prep", "nih_pid", "rin", "lesion_type", "tumour_size", "fuhrman_grade") %>%
    # if genome is human, do cell cycle scoring (doesn't work with other genomes)
    { if (grepl("hg38", genome)) c(., "Phase") else . } %>%
    # if checking for doublets, add to the groupings
    { if (remove_doublets) c(., "doublet") else . }

  # LOAD DATA AND GENE FILTERING ----
  cat("\n\nLOAD DATA AND GENE FILTERING ----\n\n")
  cat("Creating Seurat object and filtering nuclei per gene...\n")
  cat("Minimum number of nuclei per gene =", min_nuclei_per_gene, "\n")
  # load parse pipeline output to seurat object with no cut-offs, remove sample
  # NAs, add sample_metadata, add summary_stats, optionally subset samples
  # gene-level filtering (min_nuclei_per_gene)
  # -> genes that are present in very few nuclei of the dataset (<3) are uninformative
  #    and unlikely to play any part in differentiating groups of nuclei (in general,
  #    most genes removed by this filtering will be those not detected in any nuclei)

  seu <- load_parse_to_seurat(
    parse_dir,
    experiment,
    genome,
    sublibrary,
    parse_analysis_subdir,
    min_nFeature_RNA = 0,
    # perform gene filtering
    min_nuclei_per_gene = min_nuclei_per_gene,
    sample_subset,
    remove_na_samples = T,
    do_add_sample_metadata = T,
    do_add_summary_stats = T,
    groupings
  ) # subset: # seu <- seu[1:1000, 1:1000] ; out$base <- "out/test/" ; dir.create(out$base)

  # save seu object
  saveRDS(seu, file = paste0(out$base, "/seu.rds"))
  # read seu object: # seu <- readRDS(paste0(out$base, "/seu.rds"))

  # NUCLEUS FILTERING ----
  cat("\n\nNUCLEUS FILTERING ----\n\n")
  # nucleus-level (nCount_RNA, nFeature_RNA, percent_mito, doublet)
  # -> unusually high transcript/gene counts indicate multiplets
  # -> unusually low transcript/gene counts indicate barcoding of nuclei with
  #    damaged membranes (uninformative)
  # -> high % of mito genes indicates death / loss of cytoplasmic RNA / increased
  #    apoptosis (for scRNA-seq, not snRNA-seq)
  # -> high % of globin and low % of ribo suggests erythrocytes
  # -> high features:counts per nucleus ratio could be dying cells
  plots[["nucleus_and_gene_filtering"]] <- list()

  # identify doublets (by creating artificial doublets and looking at their clustering)
  # remove_doublets = F
  if (remove_doublets == T) {
    cat("Annotating suspected doublets...\n")
    seu <- annotate_doublets(seu)
  } else {
    seu$doublet <- 0
  }

  # get proportions of relevant transcript types
  cat("Annotating abundance of different transcript types...\n")
  seu <- annotate_proportions_of_transcript_types(seu)

  # get user-defined filters or designate filters
  # -> nCount/Feature_RNA, percent_mito = median plus 5 x the median absolute deviation
  #    (recommended by https://romanhaa.github.io/projects/scrnaseq_workflow)
  # -> min_nFeature_RNA = defaults recommended by Parse
  if (do_filtering) { cat("\nGetting nucleus and gene filters...\n") }
  seu <- get_filters(
    seu,
    do_filtering,
    remove_doublets,
    min_nCount_RNA,
    max_nCount_RNA,
    min_nFeature_RNA,
    max_nFeature_RNA,
    max_percent_mito
  )

  # save filters
  seu@misc$nucleus_filtering %>%
    readr::write_tsv(paste0(out$base, "seu_nucleus_filtering.tsv"))

  # nuclei per sample
  cat("Plotting nuclei per sample...\n")
  plots[["nucleus_and_gene_filtering"]][["nuclei_per_sample_bar"]] <-
    dplyr::tibble(sample_id = names(table(seu$sample)),
                  n_nuclei = table(seu$sample)) %>%
    ggplot2::ggplot(ggplot2::aes(x = sample_id, y = n_nuclei)) +
    ggplot2::geom_col() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

  # plot summary stats
  cat("Plotting summary stats...\n")
  plots[["nucleus_and_gene_filtering"]][["run_summary_stats"]] <-
    seu@misc$summary_stats %>%
    ggplot2::ggplot(ggplot2::aes(x = sample, y = value)) +
    ggplot2::geom_col() +
    ggplot2::facet_grid(statistic ~ ., scales = "free") +
    ggplot2::labs(title = "per-sample summary statistics") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

  # plot filters on distributions
  cat("Plotting filters onto distributions...\n")
  plots[["nucleus_and_gene_filtering"]][["filters_vln"]] <-
    seu@misc$nucleus_filtering %>%
    tidyr::pivot_longer(names(seu@misc$filters)) %>%
    dplyr::left_join(seu@misc$filters %>%
                       purrr::map(~dplyr::as_tibble(.)) %>%
                       dplyr::bind_rows(.id = "name")) %>%
    dplyr::group_by(fail_criteria) %>%
    ggplot2::ggplot(ggplot2::aes(x = sample, y = value)) +
    ggplot2::geom_jitter(ggplot2::aes(colour = pass),
                         height = 0, alpha = 0.5, size = 0.6) +
    ggplot2::geom_violin(ggplot2::aes(fill = sample), alpha = 0.8,
                         draw_quantiles = 0.5) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = min), colour = "red") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = max), colour = "red") +
    ggplot2::facet_wrap(~ name, scales = "free") +
    ggplot2::theme(axis.text.x = element_text(angle = 90)) +
    ditto_colours +
    ggplot2::scale_colour_manual(values = c("darkgrey", "black"))

  plots[["nucleus_and_gene_filtering"]][["filters_bar"]] <-
    seu@misc$nucleus_filtering %>%
    dplyr::group_by(sample) %>%
    dplyr::count(fail_criteria) %>%
    ggplot2::ggplot(ggplot2::aes(x = sample, y = n, fill = fail_criteria)) +
    ggplot2::geom_col() +
    ditto_colours

  # plot doublets
  cat("Plotting detected doublets...\n")
  plots[["nucleus_and_gene_filtering"]][["doublets_per_sample"]] <-
    seu@misc$nucleus_filtering %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(n_doublets = sum(doublet),
                  is_doublet = doublet == 1) %>%
    tidyr::pivot_longer(names(seu@misc$filters)) %>%
    ggplot2::ggplot(ggplot2::aes(y = value, x = sample, colour = is_doublet)) +
    ggplot2::geom_jitter(data = . %>% dplyr::filter(!is_doublet), height = 0) +
    ggplot2::geom_jitter(data = . %>% dplyr::filter(is_doublet), height = 0) +
    ggplot2::labs(title = "detected doublets in each sample (using scDblFinder)") +
    ggplot2::facet_wrap(~ name, scales = "free") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +
    ditto_colours

  # plot cell scatters
  cat("Plotting abundance of different transcript types...\n")
  plots[["nucleus_and_gene_filtering"]][["filters_scatter"]] <-
    plot_nucleus_scatter_with_filters(seu, "nCount_RNA", "percent_mito", log_y = T) +
    ggplot2::theme(legend.position = "none") +
    plot_nucleus_scatter_with_filters(seu, "nCount_RNA", "nFeature_RNA")

  # filtering
  if (do_filtering == T) {

    # perform nucleus filtering
    cat("Filtering nuclei...\n")
    n_retained <- nrow(dplyr::filter(seu@misc$nucleus_filtering, pass))
    cat(round((n_retained / ncol(seu)) * 100, 1), "% (", n_retained, "/", ncol(seu), ") of nuclei retained\n")
    seu <- subset(seu, cells = unique(dplyr::filter(seu@misc$nucleus_filtering, pass)$nucleus))

    # save
    saveRDS(seu, paste0(out$base, "seu_filtered.rds"))

  }

  # highest expressed genes
  # cat("Plotting top 20 most highly expressed genes...\n")
  # plots[["nucleus_and_gene_filtering"]][["top_genes_boxplot"]] <- boxplot_top_genes(seu)
  # TODO: fix this - it kills the job, too memory intensive

  # highest variable genes (HVGs)
  cat("Finding most highly variable genes (HVGs)...\n")
  seu <- Seurat::FindVariableFeatures(seu)

  # plot hvgs
  cat("Plotting HVGs...\n")
  plots[["nucleus_and_gene_filtering"]][["hvg_scatter"]] <- Seurat::LabelPoints(
    plot = Seurat::VariableFeaturePlot(seu, raster = F),
    points = head(Seurat::VariableFeatures(seu), 10),
    repel = T
  ) +
    ggplot2::theme(legend.position = "top")
  plots[["nucleus_and_gene_filtering"]][["hvg_heatmap"]] <- dittoSeq::dittoHeatmap(
    seu,
    genes = head(Seurat::VariableFeatures(seu), 20),
    annot.by = "sample",
    scaled.to.max = T
  )

  if (do_integration == T) {

    # INTEGRATION ----
    cat("\n\nINTEGRATION ----\n\n")
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
    cat("\n\nNORMALISATION AND SCALING (SCTransform) ----\n\n")
    cat("Running SCTransform() for normalisation and scaling...\n")
    seu <- Seurat::SCTransform(seu, vars.to.regress = c("percent_mito", "nCount_RNA"))

    # save
    saveRDS(seu, file = paste0(out$base, "/seu_transformed.rds"))

  }

  # assign nucleus cycle phase scores
  cat("Assigning cell cycle phases...\n")
  seu <- Seurat::CellCycleScoring(
    seu,
    s.features = Seurat::cc.genes$s.genes,
    g2m.features = Seurat::cc.genes$g2m.genes
  )

  # LINEAR DIMENSIONALITY REDUCTION (PCA) ----
  cat("\n\nLINEAR DIMENSIONALITY REDUCTION (PCA) ----\n\n")
  plots[["linear_dimensionality_reduction"]] <- list()

  # run PCA
  cat("Running PCA...\n")
  seu <- Seurat::RunPCA(seu)

  # choose dimensionality (final_n_dims) of the dataset
  if (!is.null(n_dims)) {

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
  plots[["linear_dimensionality_reduction"]][["pca_elbow"]] <-
    Seurat::ElbowPlot(seu, ndims = max(50, final_n_dims)) &
    ggplot2::geom_vline(ggplot2::aes(xintercept = final_n_dims, colour = "chosen n dims"),
                        alpha = 0.5) &
    ggplot2::geom_vline(ggplot2::aes(xintercept = n_dims_iD, colour = "intrisicDimensions"),
                        linetype = "dashed") &
    ggplot2::geom_vline(ggplot2::aes(xintercept = generous_n_dims, colour = "generous (iD x 2)"),
                        linetype = "dashed") &
    ggplot2::scale_colour_manual(
      name = "n_dims cut-off",
      values = c(`chosen n dims` = "red",
                 user_provided = dittoSeq::dittoColors()[[1]],
                 intrisicDimensions = dittoSeq::dittoColors()[[2]],
                 `generous (iD x 2)` = dittoSeq::dittoColors()[[3]])
    ) &
    ggplot2::theme(legend.position = c(0.7, 0.9),
                   legend.background = element_rect(fill = "white"))
  if (!is.null(n_dims)) {
    plots[["linear_dimensionality_reduction"]][["pca_elbow"]] <-
      plots[["linear_dimensionality_reduction"]][["pca_elbow"]] &
      ggplot2::geom_vline(ggplot2::aes(xintercept = n_dims, colour = "user_provided"),
                          linetype = "dashed")
  }

  # plot pca dims
  cat("Plotting principle components...\n")
  plots[["linear_dimensionality_reduction"]][["pca_dim_loadings"]] <-
    Seurat::VizDimLoadings(seu, dims = 1:4, reduction = "pca")
  plots[["linear_dimensionality_reduction"]][["pca_dim_heatmap"]] <-
    Seurat::DimHeatmap(
      seu,
      dims = 1:final_n_dims,
      cells = 500,
      balanced = TRUE,
      raster = F,
      fast = F
    )

  # NON-LINEAR DIMENSIONALITY REDUCTION (UMAP, t-SNE) ----
  cat("\n\nNON-LINEAR DIMENSIONALITY REDUCTION (UMAP, t-SNE) ----\n\n")
  plots[["nonlinear_dimensionality_reduction"]] <- list()

  # run t-SNE and UMAP
  seu <- Seurat::RunTSNE(seu)
  seu <- Seurat::RunUMAP(seu, dims = 1:final_n_dims)

  # plot projection of different sample features in different reductions
  plots[["nonlinear_dimensionality_reduction"]][["umap_vs_sample"]] <-
    dittoSeq::dittoDimPlot(seu, "sample", reduction.use = "umap", raster = F)

  # split out sample clusters, side-by-side
  plots[["nonlinear_dimensionality_reduction"]][["umap_split_by_patient"]] <-
    dittoSeq::dittoDimPlot(seu, "sample", reduction.use = "umap", split.by = "nih_pid")

  # amend groupings
  groupings <- groupings[groupings %in% colnames(seu@meta.data)]

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

  # CLUSTERING ----
  cat("\n\nCLUSTERING ----\n\n")
  plots[["clustering"]] <- list()

  # cluster data at different resolutions
  cat("Clustering data at different resolutions...\n")

  # compute k.param nearest neighbours based on euclidean distance in PCA space
  # and check shared nearest neighbours between any 2 nuclei (Jaccard similarity)
  seu <- Seurat::FindNeighbors(seu, dims = 1:final_n_dims)

  # identify clusters of nuclei by a shared nearest neighbours modularity
  # optimisation based clustering algorithm
  seu <- Seurat::FindClusters(seu, resolution = clustering_resolutions, verbose = F)
  saveRDS(seu, file = paste0(out$base, "/seu_transformed_and_clustered.rds"))

  # cluster tree plot of increasing resolutions
  cat("Plotting clustering tree at different resolutions...\n")
  snn_res_prefixes <- paste0(Seurat::DefaultAssay(seu), "_snn_res.")
  plots[["clustering"]][["clustering_tree"]] <-
    clustree::clustree(seu@meta.data[, grep(snn_res_prefixes, colnames(seu@meta.data))],
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

  # choose a resolution
  if (!is.null(final_clustering_resolution)) {

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

    # CELLTYPE ANNOTATION ----
    cat("\n\nCELLTYPE ANNOTATION ----\n\n")
    plots[["celltype_annotation"]] <- list()

    # gene modules scoring
    gene_modules %>% names %>%
      purrr::walk(function(module) {
        cat("Scoring gene module", module, "\n")
        seu <<- Seurat::AddModuleScore(
          seu,
          features = list(gene_modules[[module]]),
          name = module,
          nbin = 10
        )
        # umap
        plots[["celltype_annotation"]][[paste0("umap_vs_", module, "_module_score")]] <<-
          dittoSeq::dittoDimPlot(seu,
                                 paste0(module, "1"),
                                 main = paste0(module, " module (n=", length(gene_modules[[module]]), ")")) +
          plots[["celltype_annotation"]][["umap_vs_sample"]]
        # ridgeplot
        plots[["celltype_annotation"]][[paste0(module, "_module_score_vs_lesion_type_ridge")]] <<-
          dittoSeq::dittoRidgePlot(seu, paste0(module, "1"), group.by = "lesion_type")
      })

    # celldex annotations from the human primary cell atlas
    human_primary_ref <- celldex::HumanPrimaryCellAtlasData()
    seu_singler <- SingleR::SingleR(
      test = Seurat::GetAssayData(seu, slot = "data"),
      ref = human_primary_ref,
      labels = human_primary_ref$label.main
    )

    # plot heatmap of annotation scores
    plots[["celltype_annotation"]][["singler_annots_heatmap"]] <-
      SingleR::plotScoreHeatmap(
        seu_singler,
        annotation_col = data.frame(
          patient = seu$nih_pid,
          lesion_type = seu$lesion_type
        )
      )

    # plot delta of celltype annotations
    plots[["celltype_annotation"]][["singler_annots_delta_dist"]] <-
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
      {dplyr::left_join(seu@meta.data[, "cluster", drop = F], multiple = "all", ., by = "cluster")} %>%
      dplyr::pull(label)

    # plot cell types against reduction
    plots[["celltype_annotation"]][["umap_vs_singler_annot"]] <-
      dittoSeq::dittoDimPlot(seu, "singler_annot", size = 0.7) +
      ggplot2::geom_label(
        data = seu@meta.data %>%
          dplyr::right_join(get_centroids(seu, "umap", cluster), by = "cluster") %>%
          dplyr::distinct(cluster, singler_annot_label, x, y),
        ggplot2::aes(x, y, label = singler_annot_label),
        size = 3, fill = 'white', color = 'black', alpha = 0.5, label.size = 0,
        show.legend = FALSE, fontface = "bold",
        vjust = "inward", hjust = "inward"
      ) +
      ggplot2::labs(title = "majority celldex annotation(s) per cluster") +
      plots[["celltype_annotation"]][["umap_vs_sample"]]

    # plot cell type composition of samples
    plots[["celltype_annotation"]][["singler_annot_bar"]] <-
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
    plots[["celltype_annotation"]][["alluvial_sample_cluster_singler"]] <-
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

    # literature markers
    # # checking that all markers have been found
    # unmatched_markers <- unlist(literature_markers)[
    #   !(unlist(literature_markers) %in% unique(rownames(seu@assays$RNA)))
    # ]
    # if (length(unmatched_markers) > 0) {
    #   warning("Marker(s) ", paste(unmatched_markers, collapse = ", "), " not found!")
    # }

    # save annotated seu
    saveRDS(seu, file = paste0(out$base, "/seu_annotated.rds"))

    # TODO: PATHWAY ANALYSIS ----
    cat("\n\nPATHWAY ANALYSIS ----\n\n")
    plots[["pathway_analysis"]] <- list()

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
    plots[["pathway_analysis"]][["reactome_gsa_heatmap"]] <-
      ReactomeGSA::plot_gsva_heatmap(gsva_result, max_pathways = 15)

  }

  # SAVE PLOTS ----
  cat("\n\nSAVE PLOTS ----\n\n")
  dev_off_if()

  # plot lists of plots
  names(plots) %>%
    purrr::walk(function(section) {
      names(plots[[section]]) %>%
        purrr::walk(function(p) {
          section_plots_file <- paste0(out$base, "/", section, "_", p, "_plots.pdf")
          cat("-> Plotting", section, "plot", p, "\n")
          print(plots[[section]][[p]])
        })
    })
  dev_off_if()
  #cat("\nPlotting section", section, "plots...\nSaving to", section_plots_file, "\n")
  # ggsave(paste0(out$base, "/qc_report_plots.pdf"), plots_grob,
  #        width = 20, height = 20, units = "cm")
  # saveRDS(plots, file = paste0(out$base, "/qc_report_plots.rds"))

  # plot groupings vs reductions grob
  cat(paste0(out$base, "/groupings_vs_reductions.pdf\n"))
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
  ggsave(paste0(out$base, "/groupings_vs_reductions.pdf"), groupings_vs_reductions_grob,
         width = 20, height = 60, units = "cm")
  dev_off_if()

  # plot clustered reductions grob
  rmd_script <- system.file("rmd", "generate_qc_report.rmd", package = "vhl")
  # render(input = rmd_script, params = list())

  cat(paste0(out$base, "/clustered_reductions.pdf\n"))
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
  ggsave(paste0(out$base, "/clustered_reductions.pdf"), clustered_reductions_grob,
         width = 20, height = 60, units = "cm")
  dev_off_if()

  cat("\n\nDONE!\n\n")

  return(seu)

}
