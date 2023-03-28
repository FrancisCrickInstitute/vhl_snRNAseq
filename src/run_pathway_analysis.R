base_dir <- ifelse(Sys.info()["nodename"] == "Alexs-MacBook-Air-2.local",
                   "/Volumes/TracerX/",
                   "/camp/project/tracerX/")
setwd(paste0(base_dir,"working/VHL_GERMLINE/tidda/vhl/"))
devtools::load_all()

out_dir <- "out/230210_A01366_0351_AHNHCFDSX5_x_221202_A01366_0326_AHHTTWDMXY/hg38/comb_x_SHE5052A9_S101/all-well/DGE_filtered/unintegrated/"
# out_dir <- "out/221202_A01366_0326_AHHTTWDMXY/hg38/SHE5052A9_S101/all-well/DGE_filtered/unintegrated/"

out <- get_out(out_dir)

# load cds
cds <- readRDS(paste0(out$cache, "cds_annotated.rds"))

# patient subset
patient <- "N045"

# subset cds
cds <- cds[,
           SummarizedExperiment::colData(cds)$nih_pid == patient &
             SummarizedExperiment::colData(cds)$partition_lineage != "immune"]

# differential expression between tumour vs non-tumour
marker_test_res <-
  top_markers_develop(cds,
                      group_cells_by = "partition",
                      reference_cells = 1000)
top_specific_markers <- marker_test_res %>%
  dplyr::filter(fraction_expressing >= 0.10) %>%
  dplyr::group_by(cell_group) %>%
  dplyr::top_n(1, pseudo_R2)
top_markers %>% dplyr::as_tibble()

monocle3::plot_genes_by_group(
  cds[, SummarizedExperiment::colData(cds)$partition_lineage != "immune"],
  top_specific_markers$gene_id,
  group_cells_by = "partition_lineage",
  ordering_type = "maximal_on_diag",
  max.size = 3)

# regression analysis of normal vs tumour
nvt_cds <- cds[, SummarizedExperiment::colData(cds)$partition_lineage != "immune"]
gene_fits <-
  monocle3::fit_models(nvt_cds, model_formula_str = "~ partition_lineage + nih_pid")

ube <-
  tibble::tibble(expr = as.matrix(cds["UBE2D2",]@assays@data$counts)[1,],
                 lin = SummarizedExperiment::colData(cds)$partition_lineage)
ggplot2::ggplot(ube, ggplot2::aes(x = lin, y = expr)) + ggplot2::geom_boxplot()

# regression analysis, correcting for patient-level effects
gene_fits <-
  monocle3::fit_models(cds, model_formula_str = "~ partition_lineage + nih_pid")
