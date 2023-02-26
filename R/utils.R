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
    dplyr::mutate(pass = value >= filters[[name]]$min & value <= filters[[name]]$max,
                  fail_criteria = ifelse(pass == F, name, NA)) %>%
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
