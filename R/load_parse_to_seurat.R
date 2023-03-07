load_parse_to_seurat <-
  function(parse_dir,
           experiment,
           genome,
           sublibrary,
           parse_analysis_subdir,
           min_nFeature_RNA,
           min_nuclei_per_gene,
           sample_subset,
           remove_na_samples = F,
           do_add_sample_metadata = T,
           do_add_summary_stats = T,
           groupings) {

    # testings: # remove_na_samples = F;do_add_sample_metadata = T;do_add_summary_stats = T

    seu_ls <- purrr::map2(c(experiment), c(sublibrary), function(exp, sublib) {

      # read in DGE matrix
      sublib_dir <- paste(parse_dir, "/analysis/", exp, genome, sublib, sep = "/")
      dge_dir <- paste(sublib_dir, parse_analysis_subdir, sep = "/")
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
          min_genes = min_genes_per_nucleus,
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

        # read in sample metadata and add to misc slot
        seu_i@misc$sample_metadata <-
          paste0(parse_dir, "/expdata/", exp, "/sample_metadata.tsv") %>%
          readr::read_tsv(show_col_types = F)

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
          paste(parse_dir, "/analysis/", exp, genome, sublib, "/agg_samp_ana_summary.csv", sep = "/") %>%
          readr::read_csv(show_col_types = FALSE) %>%
          tidyr::pivot_longer(cols = -statistic, names_to = "sample") %>%
          dplyr::mutate(statistic = statistic %>% gsub(paste0(genome, "\\_"), "", .)) %>%
          dplyr::filter(
            statistic %in% c(
              "median_tscp_per_cell",
              "median_genes_per_cell",
              "fraction_tscp_in_cells"
            )
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

      # subset to sample subset
      if(!is.null(sample_subset)) {

        seu_i <- subset(x = seu_i, subset = sample %in% sample_subset)

        seu_i@misc$sample_metadata <-
          seu_i@misc$sample_metadata %>%
          dplyr::filter(sample %in% sample_subset)

        seu_i@misc$summary_stats <-
          seu_i@misc$summary_stats %>%
          dplyr::filter(sample %in% sample_subset)

      }


      seu_i

    })

    # merge if multiple runs, combine misc slot
    if (length(seu_ls) > 1) {
      seu <- Reduce(merge, seu_ls)
      seu@misc <- do.call(Map, c(rbind, purrr::map(seu_ls, ~ .x@misc))) %>% purrr::map(dplyr::distinct)
    } else {
      seu <- seu_ls[[1]]
    }

    return(seu)
  }



