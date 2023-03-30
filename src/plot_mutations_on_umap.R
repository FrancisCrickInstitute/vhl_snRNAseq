base_dir <- ifelse(Sys.info()["nodename"] == "Alexs-MacBook-Air-2.local",
                   "/Volumes/TracerX/",
                   "/camp/project/tracerX/")
setwd(paste0(base_dir,"working/VHL_GERMLINE/tidda/vhl/"))
devtools::load_all()

out_dir <- "out/230210_A01366_0351_AHNHCFDSX5_x_221202_A01366_0326_AHHTTWDMXY/hg38/comb_x_SHE5052A9_S101/all-well/DGE_filtered/unintegrated/"
out <- get_out(out_dir)

# load cds
cds <- readRDS(paste0(out$cache, "cds_annotated.rds"))

# get infercnv metadata
curr_analysis_mode <- "samples"
infercnv_metadata <- list()
SummarizedExperiment::colData(cds)$nih_pid %>%
  unique() %>%
  purrr::walk(function(nih_pid) {
    print(nih_pid)
    infercnv_dir <- paste(out$base, nih_pid, "infercnv", curr_analysis_mode, sep = "/")

    if(file.exists(paste0(infercnv_dir, "/run.final.infercnv_obj"))) {
      # write metadata file
      infercnv::add_to_seurat(infercnv_output_path = infercnv_dir)
      infercnv_metadata[[nih_pid]] <<-
        read.table(paste0(infercnv_dir, "/map_metadata_from_infercnv.txt"),
                   row.names = 1) %>%
        tibble::as_tibble(rownames = "cell")
    } else {
      cat("No run.final.infercnv_obj file for", nih_pid, "!\n")
    }
  })

# add to cds object
SummarizedExperiment::colData(cds) <-
  cbind(
    SummarizedExperiment::colData(cds),
    dplyr::left_join(
      tibble::tibble(cell = colnames(cds)),
      infercnv_metadata %>% dplyr::bind_rows(),
      by = "cell"
    ) %>% dplyr::select(-cell)
  )

# save cds
saveRDS(cds, paste0(out$cache, "cds_infercnv.rds"))

# plot driver losses on umap
dittoSeq::dittoDimPlot(cds, "has_loss_chr3", size = 0.1)

# plot 3p loss frequency
prop_cn <-
  SummarizedExperiment::colData(cds) %>%
  tibble::as_tibble(rownames = "cell") %>%
  tidyr::pivot_longer(cols = dplyr::starts_with("proportion_scaled_")) %>%
  dplyr::group_by(partition_annot, name) %>%
  dplyr::summarise(prop_cn = mean(value, na.rm = T))

prop_cn %>%
  dplyr::filter(!grepl("_cnv_", name)) %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(event = gsub("proportion\\_scaled\\_", "", name),
                direction = gsub("\\_.*", "", event),
                chr = as.numeric(gsub(".*chr", "", event)),
                n = dplyr::n()) %>%
  dplyr::filter(partition_annot %ni% c("myeloid", "lymphoid")) %>%
  ggplot2::ggplot(ggplot2::aes(x = reorder(event, -prop_cn), y = prop_cn, fill = direction)) +
  ggplot2::geom_col() +
  ggplot2::facet_grid(rows = "partition_annot") +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))



