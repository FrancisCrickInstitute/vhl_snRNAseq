load_parse_to_seurat <-
  function(experiment,
           genome,
           parse_dir,
           parse_analysis_subdir,
           min_genes,
           min_cells,
           sample_subset) {

    # get data dir
    data_dir <- paste(parse_dir, experiment, genome, parse_analysis_subdir,
                       sep = "/")

    # read in DGE matrix
    dge_mat <- Seurat::ReadParseBio(data_dir)

    # read in cell metadata
    cell_metadata <- paste0(data_dir, "cell_metadata.csv") %>%
      readr::read_csv(show_col_types = F) %>%
      tibble::column_to_rownames("bc_wells") %>%
      as.data.frame()

    # create Seurat DGE object
    seu <- dge_mat %>%
      Seurat::CreateSeuratObject(
        names.field = 0,
        meta.data = cell_metadata,
        # keep cells that have at least n genes
        min_genes = min_genes,
        # keep genes expressed in at least n cells
        min_cells = min_cells
      )

    # subset to sample subset
    if(!is.null(sample_subset)) {
      seu <- subset(x = seu, subset = sample %in% sample_subset)
    }

    return(seu)
  }
