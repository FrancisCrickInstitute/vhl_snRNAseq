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
  if (!is.null(dge_mtx_dir)) {
    dge_dir <- dge_mtx_dir
  } else {
    dge_dir <- paste(data_dir, experiment, genome, sublibrary, parse_analysis_subdir, sep = "/")
  }
  if(!all(c("DGE.mtx", "all_genes.csv", "cell_metadata.csv") %in% list.files(dge_dir))) {
    stop("Provided directory:\n\n", dge_dir, "\n\ndoes not exist or is not a Parse Biosciences split-pipe analysis directory.",
         " Must contain DGE.mtx, all_genes.csv, and cell_metadata.csv files.")
  }
  if(!is.null(experiment) & length(experiment) != length(sublibrary)) {
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

# Get cluster centroids from a Seurat or cell_data_set object
get_centroids <- function(object, reduction, lvl) {

  if (class(object)[1] == "Seurat") {
    embeddings <- object@reductions[[reduction]]@cell.embeddings
    metadata <- object@meta.data
  } else if (class(object)[1] == "cell_data_set") {
    embeddings <- SingleCellExperiment::reducedDims(object)[[reduction]]
    metadata <- SummarizedExperiment::colData(object) %>% dplyr::as_tibble()
  }

  dplyr::tibble(
    x = embeddings[,1],
    y = embeddings[,2],
    metadata
  ) %>%
    dplyr::group_by(get(lvl)) %>%
    dplyr::summarise(x = median(x), y = median(y)) %>%
    dplyr::rename(!!lvl := `get(lvl)`)

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

# Get available markers for a given object
get_available_markers <- function(object, markers) {
  markers %>%
    purrr::map(purrr::map, purrr::map, intersect, rownames(object)) %>%
    purrr::map(function(x) {
      purrr::map(x, function(y) {
        Filter(function(z) length(z) > 0, y)
      })})
}

# Plot all markers in the marker list (marker_list is a named list of marker modules)
plot_markers_on_umap <- function(object, ml, final_umap) {

  p <- final_umap

  purrr::walk(ml, function(g) {

    # get ranges for colour gradient
    if (class(object)[1] == "Seurat") {
      dat <- object@assays$RNA@counts[g, object@assays$RNA@counts[g,] > 0]
    } else if (class(object)[1] == "cell_data_set") {
      dat <- object@assays@data$counts[g, object@assays@data$counts[g,] > 0]
    }
    minmax <- range(dat, na.rm = T)

    # plot
    p[[g]] <<-
      dittoSeq::dittoDimPlot(object, g, size = 0.3, xlab = NULL, ylab = NULL) +
      ggplot2::geom_point(data =  dplyr::tibble(x = SingleCellExperiment::reducedDims(object)$UMAP[names(dat), 1],
                                                y = SingleCellExperiment::reducedDims(object)$UMAP[names(dat), 2]) %>%
                            dplyr::mutate(count = dat) %>%
                            dplyr::as_tibble(),
                          ggplot2::aes(x, y, colour = count), size = 0.4) +
      ggplot2::scale_colour_gradientn(colours = c("grey", "#2374B0", "#ECE147"),
                                      values = scales::rescale(c(0, minmax[1], minmax[2]))) +
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

# Get singler annot labels for plotting
get_singler_annot_label <- function(cds, annot, lvl) {
  level = lvl
  dplyr::as_tibble(SummarizedExperiment::colData(cds)) %>%
    dplyr::mutate(lvl = get(lvl), annot = get(annot)) %>%
    dplyr::group_by(lvl, annot) %>%
    dplyr::count() %>%
    dplyr::group_by(lvl) %>%
    dplyr::mutate(total = sum(n),
                  prop = n / sum(n) * 100,
                  label = paste0(annot, " (", round(prop, 0), "%)")) %>%
    dplyr::filter(prop == max(prop) & prop > 50 | max(prop) < 50 & rank(-prop, ties.method = "min") <= 2) %>%
    dplyr::arrange(-prop) %>%
    dplyr::summarise(label = paste0(level, " ", unique(lvl), " (n = ", unique(total), ")\n", paste(label, collapse = "\n"))) %>%
    {dplyr::left_join(dplyr::as_tibble(SummarizedExperiment::colData(cds)) %>% dplyr::transmute(lvl = get(lvl)), multiple = "all", ., by = "lvl")} %>%
    dplyr::pull(label)
}

# get patient ids from sample ids
get_patients_from_samples <- function(samples) {
  gsub("\\_.*", "", gsub("K891", "N23", samples))
}
