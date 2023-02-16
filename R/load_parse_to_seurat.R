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
           do_add_sample_metadata = F,
           do_add_summary_stats = T,
           groupings) {

    # read in DGE matrix
    sublib_dir <- paste(parse_dir, "/analysis/", experiment, genome, sublibrary, sep = "/")
    dge_dir <- paste(sublib_dir, parse_analysis_subdir, sep = "/")
    dge_mat <- Seurat::ReadParseBio(dge_dir)

    # read in cell metadata
    cell_metadata <- paste0(dge_dir, "/cell_metadata.csv") %>%
      readr::read_csv(show_col_types = F) %>%
      tibble::column_to_rownames("bc_wells") %>%
      as.data.frame()

    # create Seurat DGE object
    seu <- dge_mat %>%
      Seurat::CreateSeuratObject(
        names.field = 0,
        meta.data = cell_metadata,
        min_genes = min_genes_per_nucleus,
        min_cells = min_nuclei_per_gene
      )

    # subset to sample subset
    if(!is.null(sample_subset)) {
      seu <- subset(x = seu, subset = sample %in% sample_subset)
    }

    # remove NA samples
    if(remove_na_samples == T) {
      seu <- subset(seu, subset = sample %in% seu$sample[!is.na(seu$sample)])
    }

    # add sample metadata
    if(do_add_sample_metadata == T) {

      # read in sample metadata
      sample_metadata <- paste0(parse_dir, "/expdata/", experiment, "/sample_metadata.tsv") %>%
        readr::read_tsv(show_col_types = F)

      # add to Seurat object meta.data
      seu@meta.data <-
        dplyr::left_join(seu@meta.data,
                         sample_metadata %>% dplyr::select(dplyr::any_of(groupings)),
                         by = "sample")
      rownames(seu@meta.data) <- colnames(seu)

      # add sample-level table to misc slot
      seu@misc$sample_metadata <- sample_metadata

    }

    # add summary statistics
    if(do_add_summary_stats == T) {

      # get summary stats from parse analysis
      seu@misc$summary_stats <-
        paste(parse_dir, "/analysis/", experiment, genome, sublibrary, "/agg_samp_ana_summary.csv", sep = "/") %>%
        readr::read_csv(show_col_types = FALSE) %>%
        tidyr::pivot_longer(cols = -statistic, names_to = "sample") %>%
        dplyr::mutate(statistic = statistic %>% gsub(paste0(genome, "\\_"), "", .)) %>%
        dplyr::filter(
          statistic %in% c(
            "median_tscp_per_cell",
            "median_genes_per_cell",
            "tso_fraction_in_read1",
            "fraction_tscp_in_cells"
          )
        )  %>%
        # append numeric sample metadata to summary stats
        dplyr::bind_rows(
          seu@misc$sample_metadata %>%
            tidyr::pivot_longer(
              cols = tidyr::any_of(groupings) & where(is.numeric),
              names_to = "statistic"
            ) %>% dplyr::select(statistic, sample, value)
        )

    }

    return(seu)
  }
