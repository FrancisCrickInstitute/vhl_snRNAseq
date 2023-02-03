load_parse_to_seurat <-
  function(dge_dir,
           min_genes_per_cell,
           min_cells_per_gene,
           sample_subset) {

    # read in DGE matrix
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
        # keep cells that have at least n genes
        min_genes = min_genes_per_cell,
        # keep genes expressed in at least n cells
        min_cells = min_cells_per_gene
      )

    # subset to sample subset
    if(!is.null(sample_subset)) {
      seu <- subset(x = seu, subset = sample %in% sample_subset)
    }

    return(seu)
  }
