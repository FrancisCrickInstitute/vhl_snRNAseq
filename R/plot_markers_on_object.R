plot_markers_on_object <- function(object, object_type) {
  # testing: # object=seu;object_type="Seurat"
  # testing: # object=cds;object_type="Monocle3"

  # get available markers
  avail_markers <- get_available_markers(object, markers)

  to_knit <- c()

  for(study in names(avail_markers)) {

    to_knit <- c(to_knit, get_subheading(study))

    for(source in names(avail_markers[[study]])) {

      to_knit <- c(to_knit, get_subsubheading(study, source, avail_markers))
      author <- get_author(study)
      chunk_name <- paste(object_type, author, substr(source, 1, 3)) %>% janitor::make_clean_names()

      if (object_type == "Monocle3") {

        knit <- plot_cds_marker_modules(chunk_name, cds)
        knit <- knitr::knit_child(text = knit, envir = environment(), quiet = TRUE)
        to_knit <- c(to_knit, knit)

      }

      for(ct in names(avail_markers[[study]][[source]])) {

        to_knit <- c(to_knit, paste0("##### Population: **", ct, "**\n"))
        chunk_name <- paste(object_type, author, substr(source, 1, 3), ct) %>% janitor::make_clean_names()
        ct_markers <- avail_markers[[study]][[source]][[ct]]
        n_ct_markers <- length(ct_markers)

        if (object_type == "Seurat") {

          knit <- plot_seu_markers(chunk_name, ct, n_ct_markers, seu)

        } else if (object_type == "Monocle3") {

          knit <- plot_cds_markers(chunk_name, n_ct_markers)

        }

        knit <- knitr::knit_child(text = knit, envir = environment(), quiet = TRUE)
        to_knit <- c(to_knit, knit)

      }
    }
  }

  cat(unlist(to_knit), sep = '\n')

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
                                    cds,
                                    group = "clusters") {
  knit <- c()
  if (dplyr::n_distinct(cds@clusters$UMAP[[group]]) > 1) {
    knit <- c(knit,
              '',
              paste0('```{r ', chunk_name, '_', group, '_module_heatmap, echo = F, warning = F, message = F, fig.dim = c(10, 5)}'),
              'annotations <- avail_markers[[study]][[source]] %>% purrr::map(dplyr::as_tibble) %>% dplyr::bind_rows(.id = "annotation") %>% dplyr::group_by(value) %>% dplyr::summarise(annotation = paste(annotation, collapse = ", ")) %>% dplyr::arrange(annotation) %>% tibble::column_to_rownames("value")',
              'annotation_cols <- list(annotation = c(dittoSeq::dittoColors()[1:dplyr::n_distinct(annotations$annotation)]) %>% setNames(unique(annotations$annotation)))',
              paste0('cell_group_df <- tibble::tibble(cell = row.names(SummarizedExperiment::colData(cds)), cell_group = monocle3::', group, '(cds))'),
              'agg_mat <- monocle3::aggregate_gene_expression(cds[rownames(annotations)], cell_group_df = cell_group_df)',
              'p <- pheatmap::pheatmap(agg_mat, scale = "row", cluster_rows = F, clustering_method = "ward.D2", annotation_row = annotations, annotation_colors = annotation_cols, annotation_names_row = F, gaps_row = which(annotations$annotation != dplyr::lag(annotations$annotation)) - 1, fontsize = 3)',
              'p',
              '```',
              ''
    )
  }
  knit
}

# Function to plot markers for cds
plot_cds_markers <- function(chunk_name, n_ct_markers) {
  knit <- c()
  if(n_ct_markers < 50) {
    knit <- c(
      '',
      paste0('```{r ', chunk_name, ', echo = F, warning = F, message = F, fig.dim = get_fig_dims(n_ct_markers)}'),
      'p <- monocle3::plot_cells(cds, genes = ct_markers, show_trajectory_graph = FALSE, label_cell_groups = FALSE, label_leaves = FALSE)',
      'p',
      '```',
      ''
    )
  }
  knit
}

get_author <- function(study) {
  gsub(".*\\_", "", study)
}

get_subheading <- function(study) {
  author <- get_author(study)
  info <- stringr::str_split(study, "_")[[1]]
  paste(
    "####",
    ifelse(
      all(info == author),
      author, paste0(author, " et al. (", info[2], ")")),
    "markers\n"
  )
}

get_subsubheading <- function(study, source, avail_markers) {
  author <- get_author(study)
  info <- stringr::str_split(study, "_")[[1]]
  paste0(
    "**", length(unlist(markers[[study]][[source]])),
    "** markers sourced from **", source, "** by ",
    ifelse(
      all(info == author),
      author,
      paste0("[", author, " et al. (", info[2], ")", "](https://pubmed.ncbi.nlm.nih.gov/", info[1], ")")), ".\n",
    "**", length(unlist(avail_markers[[study]][[source]])), "** markers are present in the dataset.\n"
  )
}

