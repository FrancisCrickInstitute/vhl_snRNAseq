# Not %in% function
`%ni%` <- Negate('%in%')

# Find matches to any of a vector of character strings (for filtering)
greplany <- function(patterns, v) {
  match <- rep(FALSE, length(v))
  for (pattern in patterns) {
    match <- match | grepl(pattern, v)
  }
  return(match)
}

# Check analyse_snRNAseq args
check_analyse_snRNAseq_args <- function(args) {
  list2env(args,globalenv())
  dge_dir <- paste(parse_dir, "analysis", experiment, genome, sublibrary, parse_analysis_subdir, sep = "/")
  if(!all(c("DGE.mtx", "all_genes.csv", "cell_metadata.csv") %in% list.files(dge_dir))) {
    stop("Provided directory:\n\n", dge_dir, "\n\ndoes not exist or is not a Parse Biosciences split-pipe analysis directory.",
         " Must contain DGE.mtx, all_genes.csv, and cell_metadata.csv files.")
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
    cat("Calculating %", tt$name, "genes...\n")
    cat(tt$message, "\n")
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
    Seurat::as.SingleCellExperiment(seu), samples = "sample",
    BPPARAM = BiocParallel::MulticoreParam(3)
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
plot_nucleus_scatter_with_filters <- function(seu, x, y, ...) {
  #filters <- seu@misc$filters
  dat <- seu@meta.data %>%
    tibble::as_tibble(rownames = "nucleus") %>%
    dplyr::mutate(x_var = get(x), y_var = get(y),
                  include = x_var >= seu@misc$filters[[x]]$min & x_var <= seu@misc$filters[[x]]$max &
                    y_var >= seu@misc$filters[[y]]$min & y_var <= seu@misc$filters[[y]]$max,
                  n_include = sum(include),
                  n_exclude = dplyr::n() - include,
                  title = paste0("include=", n_include, ", exclude=", n_exclude))
  dat %>%
    ggplot2::ggplot(ggplot2::aes(x = x_var, y = y_var, colour = include)) +
    ggplot2::geom_vline(xintercept = unlist(seu@misc$filters[[x]]), colour = "red") +
    ggplot2::geom_hline(yintercept = unlist(seu@misc$filters[[y]]), colour = "red") +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::labs(title = unique(dat$title),
                  x = x,
                  y = y) +
    ggplot2::scale_colour_manual(values = c("grey", "black")) +
    ggplot2::theme(legend.position = "none")
}

# Get filters list object
get_filters <- function(seu,
                        do_filtering,
                        remove_doublets,
                        max_nCount_RNA,
                        min_nFeature_RNA,
                        max_nFeature_RNA,
                        max_percent_mito) {

  # get filters
  if(do_filtering == T) {
    filters <- list(
      "doublet" = list(
        min = 0,
        max = ifelse(remove_doublets == T, 0, 1)
      ),
      "nCount_RNA" = list(
        min = 0,
        max = ifelse(!is.null(max_nCount_RNA), max_nCount_RNA, get_max(seu$nCount_RNA))
      ),
      "nFeature_RNA" = list(
        min = ifelse(!is.null(min_nFeature_RNA), min_nFeature_RNA, 100),
        max = ifelse(!is.null(max_nFeature_RNA), max_nFeature_RNA, get_max(seu$nFeature_RNA))
      ),
      "percent_mito" = list(
        min = 0,
        max = ifelse(!is.null(max_percent_mito), max_percent_mito, get_max(seu$percent_mito))
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

  # create misc table of included/excluded nuclei, with exclude criteria
  seu@misc$nucleus_filtering <-
    seu@meta.data %>%
    dplyr::as_tibble(rownames = "nucleus") %>%
    dplyr::select(nucleus, sample, tidyr::any_of(names(filters))) %>%
    tidyr::pivot_longer(-c(nucleus, sample)) %>%
    dplyr::left_join(filters %>%
                       purrr::map(dplyr::as_tibble) %>%
                       dplyr::bind_rows(.id = "name"),
                     by = "name") %>%
    dplyr::mutate(include = value >= min & value <= max,
                  exclude_criteria = ifelse(include == F, name, NA)) %>%
    dplyr::group_by(nucleus) %>%
    dplyr::mutate(include_nucleus = all(include),
                  exclude_criteria_nucleus = exclude_criteria %>% unique %>% na.omit %>%
                    paste(collapse = ",") %>% dplyr::na_if(., "")
    )

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
