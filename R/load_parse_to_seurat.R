load_parse_to_seurat <-
  function(data_dir = NULL,
           experiment = NULL,
           genome = NULL,
           sublibrary = NULL,
           parse_analysis_subdir = NULL,
           dge_mtx_dir = NULL,
           sample_metadata_file = NULL,
           min_n_feature_RNA = 0,
           min_nuclei_per_gene,
           sample_subset = NULL,
           cell_subset = NULL,
           remove_na_samples = F,
           do_add_sample_metadata = T,
           do_add_summary_stats = T,
           groupings,
           statistics) {

    # testing: # remove_na_samples = F;do_add_sample_metadata = T;do_add_summary_stats = T

    if(is.null(experiment) & is.null(sublibrary)) {
      experiment = ""
      sublibrary = ""
    }

    seu_ls <- purrr::map2(c(experiment), c(sublibrary), function(exp, sublib) {

      cat(exp, sublib, "\n")
      # testing: # exp=experiment[2];sublib=sublibrary[2]

      # read in DGE matrix
      if (!is.null(dge_mtx_dir)) {
        dge_dir <- dge_mtx_dir
      } else {
        sublib_dir <- paste(data_dir, exp, genome, sublib, sep = "/")
        dge_dir <- paste(sublib_dir, parse_analysis_subdir, sep = "/")
      }
      dge_mat <- Seurat::ReadParseBio(dge_dir)

      # read in cell metadata
      cell_metadata <-
        paste0(dge_dir, "/cell_metadata.csv") %>%
        readr::read_csv(show_col_types = F) %>%
        tibble::column_to_rownames("bc_wells") %>%
        as.data.frame()

      # create Seurat DGE object
      seu_i <-
        dge_mat %>%
        Seurat::CreateSeuratObject(
          names.field = 0,
          meta.data = cell_metadata,
          min_genes = min_n_feature_RNA,
          min_cells = min_nuclei_per_gene
        )

      # add experiment and sublibrary col
      seu_i$experiment <- exp
      seu_i$sublibrary <- sublib

      # remove NA samples
      if(remove_na_samples == T) {
        seu_i <- subset(seu_i, subset = sample %in% seu_i$sample[!is.na(seu_i$sample)])
      }

      # add sample metadata
      if(do_add_sample_metadata == T) {

        # if sample metadata path given, read it
        if (is.null(sample_metadata_file)) {
          sm_file <- paste0(data_dir, "/../expdata/", exp, "/sample_metadata.tsv")
        } else {
          sm_file <- sample_metadata_file
        }

        # read in sample metadata and add to misc slot
        seu_i@misc$sample_metadata <-
          readr::read_tsv(sm_file, show_col_types = F)

        # add to Seurat object meta.data
        seu_i@meta.data <-
          dplyr::left_join(
            seu_i@meta.data,
            seu_i@misc$sample_metadata %>%
              dplyr::select(sample, dplyr::any_of(groupings)),
            by = "sample")
        rownames(seu_i@meta.data) <- colnames(seu_i)

      }

      # add summary statistics
      if(do_add_summary_stats == T) {

        # get summary stats from parse analysis
        seu_i@misc$summary_stats <-
          paste0(dge_dir, "/../../agg_samp_ana_summary.csv") %>%
          readr::read_csv(show_col_types = FALSE) %>%
          dplyr::rename(statistic = 1) %>%
          tidyr::pivot_longer(cols = -c("statistic"), names_to = "sample") %>%
          dplyr::mutate(statistic = statistic %>% gsub(paste0(genome, "\\_"), "", .)) %>%
          dplyr::filter(
            statistic %in% statistics
          )  %>%
          # append numeric sample metadata to summary stats
          dplyr::bind_rows(
            seu_i@misc$sample_metadata %>%
              tidyr::pivot_longer(
                cols = tidyr::any_of(groupings) & where(is.numeric),
                names_to = "statistic"
              ) %>% dplyr::select(statistic, sample, value)
          )

      }

      # subset to cell subset
      if (!is.null(cell_subset)) {
        seu_i <- seu_i[, intersect(cell_subset, colnames(seu_i))]
      }

      # subset to sample subset
      if (!is.null(sample_subset)) {

        if (length(intersect(sample_subset, seu_i$sample)) == 0) {

          seu_i <- NULL

        } else {

          seu_i <- subset(x = seu_i, subset = sample %in% sample_subset)

          seu_i@misc$sample_metadata <-
            seu_i@misc$sample_metadata %>%
            dplyr::filter(sample %in% sample_subset)

          seu_i@misc$summary_stats <-
            seu_i@misc$summary_stats %>%
            dplyr::filter(sample %in% sample_subset)

        }

      }

      seu_i

    })

    # drop empty elements
    seu_ls <- Filter(Negate(is.null), seu_ls)

    # merge if multiple runs, combine misc slot
    if (length(seu_ls) > 1) {
      seu <- Reduce(merge, seu_ls)
      seu@misc <- do.call(Map, c(rbind, purrr::map(seu_ls, ~ .x@misc))) %>% purrr::map(dplyr::distinct)
    } else {
      seu <- seu_ls[[1]]
    }

    # remove genes with no counts
    counts <- as.matrix(seu@assays$RNA@counts)
    seu <- subset(seu, features = rownames(counts[rowSums(counts) > 0, ]))

    return(seu)
  }



