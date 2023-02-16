load_parse_to_seurat <-
  function(parse_dir,
           experiment,
           genome,
           sublibrary,
           parse_analysis_subdir,
           min_nFeature_RNA,
           min_cells_per_gene,
           sample_subset,
           remove_na_samples = F,
           do_add_sample_metadata = F) {

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
        min_genes = min_genes_per_cell,
        min_cells = min_cells_per_gene
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
      seu@meta.data <- seu@meta.data %>% dplyr::left_join(sample_metadata, by = "sample")
      rownames(seu@meta.data) <- colnames(seu)

      # add sample-level table to misc slot
      seu@misc$sample_metadata <- sample_metadata

    }

    return(seu)
  }
