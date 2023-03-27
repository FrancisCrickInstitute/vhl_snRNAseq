base_dir <- ifelse(Sys.info()["nodename"] == "Alexs-MacBook-Air-2.local",
                   "/Volumes/TracerX/",
                   "/camp/project/tracerX/")
setwd(paste0(base_dir,"working/VHL_GERMLINE/tidda/vhl/"))
devtools::load_all()

patient <- "N045"
sample_pattern <- ifelse(patient == "N23", "K891", patient)

out_dir <- "out/230210_A01366_0351_AHNHCFDSX5_x_221202_A01366_0326_AHHTTWDMXY/hg38/comb_x_SHE5052A9_S101/all-well/DGE_filtered/unintegrated/"
patient_dir <- paste0(out_dir, "/", patient, "/")

out <- get_out(patient_dir)
dir.create(out$dea, showWarnings = F)

# epithelial and tumour

# load cds
cds <- readRDS(list.files(out$cache, pattern = "cds\\_singler\\_annotated", full.names = T))

# isolate only epithelial and tumour lineages (remove immune / endothelial)
renal_annotations <-
  readr::read_tsv(paste0(out_dir, "cell_annotations.tsv")) %>%
  dplyr::filter(partition_annot %in% c("epithelial", "malignant"),
                grepl(sample_pattern, sample),
                cell %in% colnames(cds))

# subset cds
cds <- cds[, renal_annotations$cell]

# add annotations to colData
SummarizedExperiment::colData(cds) <-
  cbind(
    SummarizedExperiment::colData(cds),
    renal_annotations %>% dplyr::select(-intersect(colnames(renal_annotations), colnames(SummarizedExperiment::colData(cds))))
  )

# re-cluster cells
cds <- monocle3::cluster_cells(cds)

# fit a principal graph
cds <- monocle3::learn_graph(cds)

# plot graph
monocle3::plot_cells(cds,
           color_cells_by = "cluster_annot",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_principal_points = T,
           graph_label_size=1.5)

# order cells in pseudotime
root_pr_nodes <- c("Y_97", "Y_51", "Y_1")
cds <- monocle3::order_cells(
  cds, root_pr_nodes = root_pr_nodes)

# plot graph
monocle3::plot_cells(cds,
                     color_cells_by = "pseudotime",
                     label_cell_groups=FALSE,
                     label_leaves=FALSE,
                     label_branch_points=FALSE, label_principal_points = T,
                     graph_label_size=1.5)

# dea
cds_deg <- monocle3::graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
degs <- cds_deg %>%
  tibble::as_tibble(rownames = "gene") %>%
  dplyr::arrange(plyr::desc(morans_test_statistic), plyr::desc(-q_value))

# plot top 9 degs
p_degs <- monocle3::plot_cells(cds, genes=head(degs$gene, 9),
                     show_trajectory_graph=FALSE,
                     label_cell_groups=FALSE,
                     label_leaves=FALSE)
png(paste0(out$dea, "cds_top_pseudotimed_degs.png"))
p_degs
dev.off()

# save cds, degs, and plots
saveRDS(cds, paste0(out$base, "cds_pseudotimed.rds"))
readr::write_tsv(cds_deg, paste0(out$base, "cds_pseudotimed_deg.tsv"))

# epithelial ----

# load cds
cds <- readRDS(list.files(out$cache, pattern = "cds\\_singler\\_annotated", full.names = T))

# isolate only epithelial lineages (tumour + normal)
epi_annotations <-
  readr::read_tsv(paste0(out_dir, "cell_annotations.tsv")) %>%
  dplyr::filter(partition_annot == "epithelial",
                grepl(sample_pattern, sample),
                cell %in% colnames(cds))

# subset cds
cds <- cds[, epi_annotations$cell]

# add annotations to colData
SummarizedExperiment::colData(cds) <-
  cbind(
    SummarizedExperiment::colData(cds),
    epi_annotations %>% dplyr::select(-intersect(colnames(epi_annotations), colnames(SummarizedExperiment::colData(cds))))
  )

# re-cluster cells
cds <- monocle3::cluster_cells(cds)

# fit a principal graph
cds <- monocle3::learn_graph(cds)

# plot graph
monocle3::plot_cells(cds,
                     color_cells_by = "cluster_annot",
                     label_cell_groups=FALSE,
                     label_leaves=TRUE,
                     label_principal_points = T,
                     graph_label_size=1.5)


cds_subset <- monocle3::choose_cells(cds)
cds_subset <- monocle3::order_cells(cds_subset)
subset_pr_test_res <- monocle3::graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
degs <- subset_pr_test_res %>%
  tibble::as_tibble(rownames = "gene") %>%
  dplyr::arrange(plyr::desc(morans_test_statistic), plyr::desc(-q_value))
monocle3::plot_cells(cds_subset, genes=head(degs$gene, 9),
                     show_trajectory_graph=FALSE,
                     label_cell_groups=FALSE,
                     label_leaves=FALSE)


# order cells in pseudotime
root_pr_nodes <- c("Y_97", "Y_51", "Y_1")
cds <- monocle3::order_cells(
  cds, root_pr_nodes = root_pr_nodes)

# plot graph
monocle3::plot_cells(cds,
                     color_cells_by = "pseudotime",
                     label_cell_groups=FALSE,
                     label_leaves=FALSE,
                     label_branch_points=FALSE, label_principal_points = T,
                     graph_label_size=1.5)

# deg
cds_deg <- monocle3::graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
degs <- cds_deg %>%
  tibble::as_tibble(rownames = "gene") %>%
  dplyr::arrange(plyr::desc(morans_test_statistic), plyr::desc(-q_value))

monocle3::plot_cells(cds, genes=head(degs$gene, 6),
                     show_trajectory_graph=FALSE,
                     label_cell_groups=FALSE,
                     label_leaves=FALSE)



