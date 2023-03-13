# Not %in% function
`%ni%` <- Negate('%in%')

# Convert all double slashes in path to single
clean_path <- function(path) {path %>% gsub("//", "/", .)}

# dev.off() if there is a device
dev_off_if <- function() {if (length(dev.list()) != 0) { dev.off() }}

# Find matches to any of a vector of character strings (for filtering)
greplany <- function(patterns, v) {
  match <- rep(FALSE, length(v))
  for (pattern in patterns) {
    match <- match | grepl(pattern, v)
  }
  return(match)
}

# Check args
check_args <- function(args) {
  list2env(args,globalenv())
  dge_dir <- paste(parse_dir, "analysis", experiment, genome, sublibrary, parse_analysis_subdir, sep = "/")
  if(!all(c("DGE.mtx", "all_genes.csv", "cell_metadata.csv") %in% list.files(dge_dir))) {
    stop("Provided directory:\n\n", dge_dir, "\n\ndoes not exist or is not a Parse Biosciences split-pipe analysis directory.",
         " Must contain DGE.mtx, all_genes.csv, and cell_metadata.csv files.")
  }
  if(length(experiment) != length(sublibrary)) {
    stop("Length of the experiment argument does not equal the length of the sublibrary argument!")
  }
}

# Get base_dir
get_base_dir <- function() {
  if (Sys.info()["nodename"]=="Alexs-MacBook-Air-2.local") {
    "/Volumes/TracerX/working/VHL_GERMLINE/tidda/"
  } else {
    "/camp/project/tracerX/working/VHL_GERMLINE/tidda/"
  }
}

# Get cluster centroids from a Seurat object
get_centroids <- function(seu, reduction, ...) {
  dplyr::tibble(
    x = seu@reductions[[reduction]]@cell.embeddings[,1],
    y = seu@reductions[[reduction]]@cell.embeddings[,2],
    seu@meta.data
  ) %>%
    dplyr::group_by(...) %>%
    dplyr::summarise(x = median(x), y = median(y))
}

# Get filters (anything below the median plus 5 x the median absolute deviation)
get_max <- function(x, n_mads = 5) { median(x) + n_mads * mad(x) }

# Annotate proportions of relevant transcript types
annotate_proportions_of_transcript_types <- function(seu) {
  purrr::pwalk(transcript_types, function(...) {
    tt <- tibble::tibble(...)
    # use `<<` for global assignment
    seu <<- seu %>%
      Seurat::PercentageFeatureSet(pattern = tt$pattern,
                                   col.name = paste0("percent_", tt$name))
  })

  return(seu)
}

# Annotate doublets using the scDblFinder package
annotate_doublets <- function(seu) {
  seu_dbl <- scDblFinder::scDblFinder(
    Seurat::as.SingleCellExperiment(seu), samples = "sample"
  )
  seu$doublet <- as.numeric(seu_dbl$scDblFinder.class == "doublet")
  return(seu)
}

# Check distribution of genes per cell
plot_genes_per_nucleus_dist <- function(seu) {
  seu@meta.data %>%
    dplyr::as_tibble() %>%
    dplyr::transmute(
      nCount_RNA,
      nFeature_RNA,
      `bottom 95%` = 0.95,
      all = 1,
      `top 5%` = 0.95
    ) %>%
    tidyr::pivot_longer(c("nCount_RNA", "nFeature_RNA"),
                        names_to = "statistic") %>%
    tidyr::pivot_longer(c("bottom 95%", "all", "top 5%"),
                        names_to = "label",
                        values_to = "quantile") %>%
    dplyr::group_by(label) %>%
    dplyr::filter(
      grepl("top", label) & value >= quantile(value, quantile) |
        grepl("bottom", label) &
        value <= quantile(value, quantile) |
        label == "all"
    ) %>%
    dplyr::mutate(label = paste0(label, " (n = ", dplyr::n(), ")")) %>%
    ggplot2::ggplot(ggplot2::aes(value)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ..density..),
                            bins = 50,
                            fill = "grey") +
    ggplot2::geom_density() +
    ggplot2::facet_wrap(statistic ~ label, scales = "free")
}

# Plot nucleus scatter with filters
plot_nucleus_scatter_with_filters <- function(seu, x, y, ditto_colours, log_x = F, log_y = F) {
  dat <- seu@misc$nucleus_filtering %>%
    dplyr::mutate(x_var = get(x) %>% { if(log_x) log(.) else . },
                  y_var = get(y) %>% { if(log_y) log(.) else . },
                  title = paste0(x, " n fail = ", sum(grepl(x, fail_criteria)), " / ", dplyr::n(), "\n",
                                 y, " n fail = ", sum(grepl(y, fail_criteria)), " / ", dplyr::n()))
  dat %>%
    ggplot2::ggplot(ggplot2::aes(x = x_var, y = y_var, colour = fail_criteria)) +
    ggplot2::geom_vline(xintercept = unlist(seu@misc$filters[[x]]) %>% { if(log_x) log(.) else . },
                        colour = "red") +
    ggplot2::geom_hline(yintercept = unlist(seu@misc$filters[[y]]) %>% { if(log_y) log(.) else . },
                        colour = "red") +
    ggplot2::geom_point(size = 0.5, alpha = 0.7) +
    ggplot2::geom_point(data = . %>% dplyr::filter(pass == F), size = 0.5, alpha = 0.7) +
    ggplot2::labs(title = unique(dat$title),
                  x = paste0(ifelse(log_x, "log ", ""), x),
                  y = paste0(ifelse(log_y, "log ", ""), y)) +
    ditto_colours +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10))

}

# Get filters list object
get_filters <- function(seu,
                        do_filtering,
                        remove_doublets,
                        min_nCount_RNA,
                        max_nCount_RNA,
                        min_nFeature_RNA,
                        max_nFeature_RNA,
                        max_percent_mito) {

  # get filters (user-provided or median + n_mads * MAD)
  if(do_filtering == T) {
    filters <- list(
      "doublet" = list(
        min = 0,
        max = ifelse(remove_doublets == T, 0, 1)
      ),
      "nCount_RNA" = list(
        min = ifelse(!is.null(min_nCount_RNA), min_nCount_RNA, 0),
        max = ifelse(!is.null(max_nCount_RNA), max_nCount_RNA, get_max(seu$nCount_RNA, n_mads = 3))
      ),
      "nFeature_RNA" = list(
        min = ifelse(!is.null(min_nFeature_RNA), min_nFeature_RNA, 100),
        max = ifelse(!is.null(max_nFeature_RNA), max_nFeature_RNA, get_max(seu$nFeature_RNA, n_mads = 3))
      ),
      "percent_mito" = list(
        min = 0,
        max = ifelse(!is.null(max_percent_mito), max_percent_mito, get_max(seu$percent_mito, n_mads = 3))
      )
    )
  } else {
    filters <- list(
      "doublet" = list(min = -Inf, max = Inf),
      "nCount_RNA" = list(min = -Inf, max = Inf),
      "nFeature_RNA" = list(min = -Inf, max = Inf),
      "percent_mito" = list(min = -Inf, max = Inf)
    )
  }

  # add to seu
  seu@misc$filters <- filters

  # create misc table of pass/fail nuclei, with fail criteria
  seu@misc$nucleus_filtering <-
    seu@meta.data[, colnames(seu@meta.data) %in% c("sample", names(filters))] %>%
    dplyr::as_tibble(rownames = "nucleus") %>%
    tidyr::pivot_longer(-c(nucleus, sample)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      fail_criteria = ifelse(
        value < filters[[name]]$min, paste0("min_", name),
        ifelse(
          value > filters[[name]]$max, paste0("max_", name),
          NA
        )
      ),
      pass = is.na(fail_criteria)) %>%
    dplyr::group_by(nucleus, sample) %>%
    dplyr::mutate(pass = all(pass),
                  fail_criteria = fail_criteria %>% unique %>% na.omit %>% paste(collapse = ",") %>% dplyr::na_if(., "")) %>%
    tidyr::pivot_wider(id_cols = c("nucleus", "sample", "pass", "fail_criteria")) %>%
    dplyr::ungroup()

  return(seu)
}

# Generate a boxplot of the top n genes
boxplot_top_genes <- function(seu, n_genes = 20) {
  # matrix of raw counts
  cts <- seu %>%
    Seurat::GetAssayData(assay = "RNA", slot = "counts") %>%
    as.matrix()

  # get percentage/nucleus
  cts <- t(cts) / colSums(cts) * 100
  medians <- apply(cts, 2, median)

  # get top n genes
  most_expressed <- order(medians, decreasing = T)[n_genes:1]
  most_exp_matrix <- as.matrix((cts[, most_expressed]))

  # prepare for plotting
  most_exp_df <- stack(as.data.frame(most_exp_matrix))
  colnames(most_exp_df) <- c("perc_total", "gene")

  # boxplot with ggplot2
  most_exp_df %>%
    ggplot2::ggplot(ggplot2::aes(x = gene, y = perc_total)) +
    ggplot2::geom_boxplot() +
    ggplot2::coord_flip()
}

# Get available markers (marker_list is a named list of marker modules)
get_available_markers <- function(seu, marker_list) {
  ml <- marker_list %>% purrr::map(intersect, rownames(seu))
  ml <- Filter(function(x) length(x) > 0, ml)
  return(ml)
}

# Plot all markers in the marker list (marker_list is a named list of marker modules)
plot_markers_on_umap <- function(seu, ml, seu_final_clusters_umap, umap_void_theme) {
  # list of plots
  p <- list()
  p[["clusters"]] <- seu_final_clusters_umap
  purrr::walk(ml, function(g) {
    p[[g]] <<-
      dittoSeq::dittoDimPlot(seu, g,
                             size = 0.3, xlab = NULL, ylab = NULL) +
      ggplot2::scale_colour_gradientn(colours = c("lightgrey", "#ECE147", "#2374B0"),
                                      values = scales::rescale(
                                        c(0,
                                          min(seu@assays$RNA@counts[g,seu@assays$RNA@counts[g,] > 0],
                                              na.rm = T),
                                          max(seu@assays$RNA@counts[g,], na.rm = T)))) +
      umap_void_theme
    })
  # create grob layout
  p <-
    gridExtra::marrangeGrob(grobs = p,
                            ncol = 4,
                            nrow = ceiling(length(p) / 4),
                            top = NULL)
  return(p)
}



# Calculate plot height for grids
get_fig_dims <- function(n_plots, n_cols = NULL, grid_width = 10, height_to_width_ratio = 1) {
  # get dims
  dims <- ggplot2::wrap_dims(n_plots, n_cols)
  # get scale
  scale <- grid_width / dims[1]
  # get height
  grid_height <- (dims[2] * scale) * height_to_width_ratio
  # return dims
  c(grid_width, grid_height)
}

# Function to plot markers for seu
plot_seu_markers <- function(chunk_name,
                             ct,
                             n_ct_markers,
                             seu) {
  knit1 <- c()
  if(n_ct_markers > 1) {
    knit1 <- c(knit1,
      '',
      paste0('```{r ', chunk_name, '_module, echo = F, warning = F, message = F, fig.dim = c(10,5)}'),
      'seu <- Seurat::AddModuleScore(seu, features = list(ct_markers), name = chunk_name)',
      paste0('p <- seu_final_clusters_umap + dittoSeq::dittoDimPlot(seu, paste0(chunk_name, "1"), size = 0.5, main = "', ct, ' module score") + umap_void_theme'),
      'p',
      '```',
      ''
    )
  }
  if(n_ct_markers < 50) {
    knit2 <- c(
      '',
      paste0('```{r ', chunk_name, ', echo = F, warning = F, message = F, fig.dim = get_fig_dims((n_ct_markers + 1), n_cols = 4, height_to_width_ratio = 1.2)}'),
      'p <- plot_markers_on_umap(seu, ml = ct_markers, seu_final_clusters_umap, umap_void_theme)',
      'p',
      '```',
      ''
    )
  } else {
    knit2 <- paste(n_ct_markers, "markers in this set - too many to print!\n")
  }
  knit <- c(knit1, knit2)
  knit
}

# Function to plot marker modules for cds
plot_cds_marker_modules <- function(chunk_name,
                                    cds) {
  knit <- c()
  for(group in c("clusters", "partitions")) {
    if (dplyr::n_distinct(cds@clusters$UMAP[[group]]) > 1) {
      knit <- c(knit,
                 '',
                 paste0('```{r ', chunk_name, '_', group, 'module_heatmap, echo = F, warning = F, message = F, fig.dim = c(10, 5)}'),
                 'gene_group_df <- avail_markers %>% purrr::map(dplyr::as_tibble) %>% dplyr::bind_rows(.id = "group") %>% dplyr::select(value, group)',
                 paste0('cell_group_df <- tibble::tibble(cell = row.names(SummarizedExperiment::colData(cds)), cell_group = monocle3::', group, '(cds))'),
                 'agg_mat <- monocle3::aggregate_gene_expression(cds, gene_group_df, cell_group_df)',
                 'p <- pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")',
                 'p',
                 '```',
                 ''
      )
    }
  }
  knit
}

# Function to plot markers for cds
plot_cds_markers <- function(chunk_name, n_ct_markers) {
  knit <- c(
    '',
    paste0('```{r ', chunk_name, ', echo = F, warning = F, message = F, fig.dim = get_fig_dims(n_ct_markers)}'),
    'p <- monocle3::plot_cells(cds, genes = c(avail_markers[[ct]]), show_trajectory_graph = FALSE, label_cell_groups = FALSE, label_leaves = FALSE)',
    'p',
    '```',
    ''
  )
  knit
}
