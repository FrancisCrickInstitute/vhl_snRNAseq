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

# Get summary stats from the parse-pipeline output directory
get_summary_stats <- function(
    data_dir,
    statistics = c(
      "median_tscp_per_cell",
      "median_genes_per_cell",
      "tso_fraction_in_read1",
      "fraction_tscp_in_cells"
    )) {
  readr::read_csv(paste0(data_dir, "/agg_samp_ana_summary.csv"),
                  show_col_types = FALSE) %>%
    tidyr::pivot_longer(cols = -statistic, names_to = "sample") %>%
    dplyr::mutate(statistic = statistic %>% gsub(paste0(genome, "\\_"), "", .)) %>%
    dplyr::filter(
      statistic %in% statistics
    )
}

# Get sample metadata for the experiment
add_sample_metadata <- function(seu, parse_pipeline_dir, experiment) {

  # load sample metadata
  sample_metadata <- paste(parse_pipeline_dir,
                           "expdata",
                           experiment,
                           "sample_metadata.tsv",
                           sep = "/") %>%
    readr::read_tsv(show_col_types = FALSE) %>%
    janitor::clean_names() %>%
    # add patient x sample type label
    dplyr::mutate(
      sample_type_code = dplyr::case_when(
        sample_type == "renal_cyst" ~ "c",
        sample_type == "solid" ~ "t",
        sample_type == "normal_renal" ~ "n",
        sample_type == "metastasis" ~ "met",
        sample_type == "mixed" ~ "mix",
        TRUE ~ sample_type
      ),
      sample_label = paste0(sample, "_", sample_type_code)
    )

  # add to Seurat object meta.data
  seu@meta.data <- seu@meta.data %>% dplyr::left_join(sample_metadata, by = "sample")
  rownames(seu@meta.data) <- colnames(seu)

  # add sample-level table to misc slot
  seu@misc$sample_metadata <- sample_metadata

  return(seu)

}

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
    Seurat::as.SingleCellExperiment(seu), samples = "sample"
  )
  seu$doublet <- as.numeric(seu$multiplet_class == "doublet")
  return(seu)
}

# Check distribution of genes per cell
plot_genes_per_cell_dist <- function(seu) {
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

# Plot cell scatter with filters
plot_cell_scatter_with_filters <- function(seu, x, y, filters, ...) {
  dat <- seu@meta.data %>%
    tibble::as_tibble(rownames = "cell") %>%
    dplyr::mutate(x_var = get(x), y_var = get(y),
                  include = x_var >= filters[[x]]$min & x_var <= filters[[x]]$max &
                    y_var >= filters[[y]]$min & y_var <= filters[[y]]$max,
                  n_include = sum(include),
                  n_exclude = dplyr::n() - include,
                  title = paste0("include=", n_include, ", exclude=", n_exclude))
  dat %>%
    ggplot2::ggplot(ggplot2::aes(x = x_var, y = y_var, colour = include)) +
    ggplot2::geom_vline(xintercept = unlist(filters[[x]]), colour = "red") +
    ggplot2::geom_hline(yintercept = unlist(filters[[y]]), colour = "red") +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::labs(title = unique(dat$title),
                  x = x,
                  y = y) +
    ggplot2::scale_colour_manual(values = c("grey", "black")) +
    ggplot2::theme(legend.position = "none")
}

# Get filters list object
get_filters <- function(seu,
                        remove_doublets,
                        max_nCount_RNA,
                        min_nFeature_RNA,
                        max_nFeature_RNA,
                        max_percent_mito) {
  if(do_filtering == T) {
    list(
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
    list(
      "doublet" = list(min = -Inf, max = Inf),
      "nCount_RNA" = list(min = -Inf, max = Inf),
      "nFeature_RNA" = list(min = -Inf, max = Inf),
      "percent_mito" = list(min = -Inf, max = Inf)
    )
  }
}

