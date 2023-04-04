# config
base_dir <- ifelse(Sys.info()["nodename"] == "Alexs-MacBook-Air-2.local",
                   "/Volumes/TracerX/",
                   "/camp/project/tracerX/")
setwd(paste0(base_dir,"working/VHL_GERMLINE/tidda/vhl/"))
devtools::load_all()
out_dir <- "out/230210_A01366_0351_AHNHCFDSX5_x_221202_A01366_0326_AHHTTWDMXY/hg38/comb_x_SHE5052A9_S101/all-well/DGE_filtered/unintegrated/"
out <- get_out(out_dir)

# run tumour-specific markers analysis from Wu et al.
source("src/get_tumour_specific_markers.R")

# read in degs
malignant_specific_DEGs <-
  read.table(paste0(out_dir, "tumour_specific_markers/malignant_vs_others_pairwise_DE.txt")) %>%
  tibble::as_tibble(rownames = "sample.gene") %>% dplyr::arrange(-avg_log2FC)
DE_genes_filtered <-
  read.table(paste0(out_dir, "tumour_specific_markers/malignant_specific_DEG.txt"))$V1

# plot degs
p_dat <-
  tumor_vs_rest_DE_total %>%
  dplyr::filter(gene_symbol %in% DE_genes_filtered,
                cell_type %in% c("endothelial", "epithelial")) %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::summarise(neg_log_p_val =  mean(-log10(p_val)),
                   avg_log2FC = mean(avg_log2FC))
p_dat %>%
  ggplot2::ggplot(ggplot2::aes(x = avg_log2FC, y = neg_log_p_val)) +
  ggplot2::geom_point() +
  ggrepel::geom_text_repel(data = p_dat %>% dplyr::slice_max(neg_log_p_val, n = 20),
                           ggplot2::aes(label = gene_symbol))



# load cds
cds <- readRDS(paste0(out$cache, "cds_annotated.rds"))

# get cell annotations
annots <-
  readr::read_tsv(paste0(out$base, "cell_annotations.tsv")) %>%
  dplyr::select(-cluster, -partition)
patients <- unique(SummarizedExperiment::colData(cds)$nih_pid)

top_lin_markers <- list()
for(curr_nih_pid in patients) {
  curr_cds <-
    readRDS(list.files(path = paste0(out$base, "/", curr_nih_pid, "/cache/"),
                       pattern = "cds\\_clustered", full.names = T))
  # get partition lineages
  SummarizedExperiment::colData(curr_cds)$partition_lineage <-
    dplyr::left_join(
      SummarizedExperiment::colData(curr_cds) %>%
        tibble::as_tibble(rownames = "cell"),
      SummarizedExperiment::colData(curr_cds) %>%
        tibble::as_tibble(rownames = "cell") %>%
        dplyr::left_join(annots) %>%
        dplyr::filter(!is.na(cluster_annot)) %>%
        dplyr::group_by(partition, partition_lineage) %>%
        dplyr::count(partition_lineage) %>%
        dplyr::group_by(partition) %>%
        dplyr::filter(n == max(n)) %>%
        dplyr::distinct(partition, partition_lineage)
    ) %>% dplyr::pull(partition_lineage)
  top_lin_markers[[curr_nih_pid]] <-
    top_markers_fix(curr_cds, group_cells_by = "partition_lineage")
}

# get common markers
top_lin_markers_df <-
  top_lin_markers %>%
  dplyr::bind_rows(.id = "nih_pid") %>%
  tibble::as_tibble() %>%
  dplyr::filter(
    # sig
    marker_test_q_value < 0.05,
    # malig
    cell_group == "malignant"
    ) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(n_patients = dplyr::n()) %>%
  dplyr::filter(n_patients > 1) %>%
  dplyr::summarise(dplyr::across(where(is.numeric), ~ mean(.x)),
                   patients = paste(nih_pid, collapse = ",")) %>%
  dplyr::arrange(desc(pseudo_R2))
top_lin_markers_df %>%
  readr::write_tsv(paste0(out$dea, "top_malignant_lineage_markers_across_patients.tsv"))

# # get differentially expressed genes across groupings in each patient
#
# top_markers <-
#   patients %>%
#   purrr::map(function(curr_nih_pid) {
#     print(curr_nih_pid)
#     # get partition lineages
#     curr_partition_lineages <-
#       SummarizedExperiment::colData(
#         readRDS(list.files(path = paste0(out$base, "/", curr_nih_pid, "/cache/"),
#                            pattern = "cds\\_clustered", full.names = T))) %>%
#       tibble::as_tibble(rownames = "cell") %>%
#       dplyr::left_join(annots) %>%
#       dplyr::filter(!is.na(cluster_annot)) %>%
#       dplyr::group_by(partition, partition_lineage) %>%
#       dplyr::count(partition_lineage) %>%
#       dplyr::group_by(partition) %>%
#       dplyr::filter(n == max(n))
#     # get lineage annotations
#     curr_markers <-
#       readRDS(list.files(path = paste0(out$base, "/", curr_nih_pid, "/cache/"),
#                          pattern = "cds\\_markers", full.names = T))$partition %>%
#       tibble::as_tibble() %>%
#       dplyr::rename(partition = cell_group) %>%
#       dplyr::left_join(curr_partition_lineages)
#   }) %>%
#   setNames(patients) %>%
#   dplyr::bind_rows(.id = "nih_pid")




