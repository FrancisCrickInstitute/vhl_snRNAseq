base_dir <- ifelse(Sys.info()["nodename"] == "Alexs-MacBook-Air-2.local",
                   "/Volumes/TracerX/",
                   "/camp/project/tracerX/")
setwd(paste0(base_dir,"working/VHL_GERMLINE/tidda/vhl/"))
devtools::load_all()

out_dir <- "out/230210_A01366_0351_AHNHCFDSX5_x_221202_A01366_0326_AHHTTWDMXY/hg38/comb_x_SHE5052A9_S101/all-well/DGE_filtered/unintegrated/"
out <- get_out(out_dir)

# load cds
cds <- readRDS(paste0(out$cache, "cds_annotated.rds"))
col_data <- SummarizedExperiment::colData(cds) %>% tibble::as_tibble(rownames = "cell")

# get mutation data
muts <-
  readr::read_tsv(paste0(out$base, "mutect2_mutations.tsv"))

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

# add to cnv and muts and hif activation to cds object
SummarizedExperiment::colData(cds) <-
  cbind(
    SummarizedExperiment::colData(cds),
    tibble::tibble(cell = colnames(cds),
                   sample = SummarizedExperiment::colData(cds)$sample) %>%
      dplyr::left_join(
        infercnv_metadata %>% dplyr::bind_rows(),
        by = "cell") %>%
      dplyr::left_join(
        muts %>%
          dplyr::select(-nih_pid) %>%
          dplyr::rename_with(.cols = -c("sample"), ~ paste0(.x, "_mut")),
        by = "sample") %>%
      dplyr::select(-sample, -cell)
  )

# add hif expression in each cell to coldata
hif_activ <-
  monocle3::aggregate_gene_expression(
    cds,
    gene_group_df =
      tibble::enframe(markers[c("HIF", "proximal tubule")]) %>%
      tidyr::unnest() %>%
      dplyr::select(2,1))
SummarizedExperiment::colData(cds)$hif <-
  hif_activ["HIF",]

# add vhl status to coldata
SummarizedExperiment::colData(cds)$vhl_status <-
  SummarizedExperiment::colData(cds) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(vhl_status = dplyr::case_when(
    partition_lineage == "immune" ~ NA,
    has_loss_chr3 == 1 & VHL_mut == 1 ~ "VHL mut + 3p loss",
    has_loss_chr3 == 1 ~ "3p loss",
    VHL_mut == 1 ~ "VHL mut",
    TRUE ~ "WT")) %>%
  dplyr::pull(vhl_status)

# save cds
saveRDS(cds, paste0(out$cache, "cds_infercnv.rds"))
# cds <- readRDS(paste0(out$cache, "cds_infercnv.rds"))
col_data <- SummarizedExperiment::colData(cds) %>% tibble::as_tibble(rownames = "cell")

# get hif activation
hif_genes <-
  monocle3::aggregate_gene_expression(
    cds[rownames(cds) %in% c(markers$HIF),
        col_data$partition_lineage != "immune" & !is.na(col_data$has_loss_chr3)],
    cell_group_df = col_data %>%
      dplyr::transmute(cell,
                       lin_chr = paste(partition_annot, has_loss_chr3, sep = "_")))
pheatmap::pheatmap(hif_genes, scale = "row", clustering_method = "ward.D2")

agg_hif <-
  monocle3::aggregate_gene_expression(
    cds[rownames(cds) %in% c(markers$HIF, markers$`proximal tubule`),
        col_data$partition_lineage != "immune" & !is.na(col_data$has_loss_chr3)],
    gene_group_df =
      tibble::enframe(markers[c("HIF", "proximal tubule")]) %>% tidyr::unnest() %>% dplyr::select(2,1))


# get col_data
kruskal.test(proportion_scaled_loss_chr3 ~ cluster_annot,
             data = SummarizedExperiment::colData(cds) %>%
               tibble::as_tibble())

col_data %>%
  dplyr::filter(partition_lineage != "immune") %>%
  {
    pairwise.wilcox.test(.$proportion_scaled_loss_chr3,
                         .$partition_lineage,
                         p.adjust.method = "BH")
  }
col_data %>%
  dplyr::filter(partition_lineage != "immune") %>%
  dplyr::group_by(partition_lineage) %>%
  dplyr::summarise(mean(proportion_scaled_loss_chr3, na.rm = T))

# ridge plot
pdf(paste0(out$base, "3p_loss_ridge.pdf"))
ridge_3p <- list()
purrr::walk(c("malignant", "normal"), function(lin) {
  p_dat <-
    cds[, SummarizedExperiment::colData(cds)$partition_lineage == lin]
  ridge_3p[[lin]] <<-
    dittoSeq::dittoRidgePlot(
      p_dat,
      "proportion_scaled_loss_chr3",
      group.by = "cluster_annot",
      ylab = "scaled proportion of chr3 loss",
      ridgeplot.lineweight = 0.5,
      legend.show = F, main = NULL) +
    ggplot2::lims(x = c(-0.15, 0.75)) +
    ggplot2::theme(axis.title = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size = 15))
})
patchwork::wrap_plots(
  list(
    ridge_3p$malignant +
      ggplot2::scale_fill_manual(values = "darkred"),
    ridge_3p$normal),
  ncol = 1, heights = c(1, 3))
dev.off()

# plot mutations
muts_to_plot <-
  colnames(SummarizedExperiment::colData(cds))[
    grepl("_mut", colnames(SummarizedExperiment::colData(cds)))
  ]
monocle3::plot_cells(cds, color_cells_by = "VHL_mut")
monocle3::plot_cells(cds, color_cells_by = "PBRM1_mut")

pdf(paste0(out$base, "chr3_status_by_lesion_type.pdf"), width = 16, height = 13)
unique(col_data$lesion_type) %>%
  purrr::map(function(lt) {
    cds[,col_data$lesion_type == lt] %>%
      monocle3::plot_cells(color_cells_by = "has_loss_chr3") +
      ggplot2::labs(title = lt)
  }) %>%
  patchwork::wrap_plots(ncol = 2)
dev.off()

vhl_status_colours <- c("3p loss" = "#F0E442",
                        "VHL mut" = "#ff82da",
                        "VHL mut + 3p loss" = "#98DFD6",
                        "WT" = "#00235B"
)

pdf(paste0(out$base, "vhl_status.pdf"))
monocle3::plot_cells(cds, color_cells_by = "vhl_status",
                     label_cell_groups = F) +
  ggplot2::scale_colour_manual(values = vhl_status_colours) +
  ggplot2::coord_equal()
dev.off()

cds[,SummarizedExperiment::colData(cds)$partition_lineage == "normal"] %>%
  {
    list(
      monocle3::plot_cells(., genes = "VCAM1"),
      monocle3::plot_cells(., color_cells_by = "vhl_status", label_cell_groups = F) +
        ggplot2::scale_colour_manual(values = vhl_status_colours),
      monocle3::plot_cells(., color_cells_by = "lesion_type", label_cell_groups = F),
      monocle3::plot_cells(., color_cells_by = "hif")
    )
  } %>%
  patchwork::wrap_plots()


c(0, 1) %>%
  purrr::map(function(chr3) {
    chr3_loss_subset <-
      SummarizedExperiment::colData(cds)$has_loss_chr3 == chr3 &
      SummarizedExperiment::colData(cds)$partition_annot == "epithelial"
    chr3_loss_subset[is.na(chr3_loss_subset)] <- FALSE
    SummarizedExperiment::colData(cds)$has_loss_chr3 <-
      as.character(SummarizedExperiment::colData(cds)$has_loss_chr3)
    cds[,chr3_loss_subset] %>%
      monocle3::plot_cells(color_cells_by = "has_loss_chr3") +
      ggplot2::scale_colour_manual(values = c("0" = "#00235B", "1" = "#F0E442")) +
      ggplot2::coord_equal() +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "none")
  }) %>%
  patchwork::wrap_plots(ncol = 2)


cds[,SummarizedExperiment::colData(cds)$partition_lineage %in% c("normal", "malignant")] %>%
  monocle3::plot_cells(color_cells_by = "partition_lineage", label_cell_groups = F) +
  ggplot2::scale_colour_manual(values = c("malignant" = "#56B4E9",
                                          "normal" = "#009E73")) +
  ggplot2::coord_equal() +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none")


monocle3::plot_cells(cds, color_cells_by = "partition") +
  ggplot2::coord_equal() +
  ggplot2::theme_void() +
  ditto_colours

pdf(paste0(out$base, "chr3_status_by_lesion_type.pdf"), width = 16, height = 13)
unique(col_data$lesion_type) %>%
  purrr::map(function(lt) {
    cds[,col_data$lesion_type == lt] %>%
      monocle3::plot_cells(color_cells_by = "vhl_status") +
      ggplot2::labs(title = lt) +
      ggplot2::scale_colour_manual(values = vhl_status_colours)
  }) %>%
  patchwork::wrap_plots(ncol = 2)
dev.off()

pdf(paste0(out$base, "vhl_status_by_partition_annot.pdf"))
col_data %>%
  dplyr::filter(partition_lineage != "immune") %>%
  dplyr::inner_join(
    col_data %>%
      dplyr::count(partition_annot, vhl_status) %>%
      dplyr::group_by(partition_annot) %>%
      dplyr::mutate(prop = (n / sum(n)) * 100) %>%
      dplyr::filter(vhl_status == "3p loss") %>%
      dplyr::select(partition_annot, prop, n) %>%
      dplyr::filter(n > 10)
    ) %>%
  dplyr::mutate(lesion_type = factor(lesion_type,
                                         levels = c("solid",
                                                    "renal_cyst",
                                                    "metastasis",
                                                    "normal_renal"))) %>%
  ggplot2::ggplot(ggplot2::aes(x = reorder(partition_annot, -n),
                               fill = vhl_status)) +
  ggplot2::geom_bar() +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0)) +
  ggplot2::scale_fill_manual(values = vhl_status_colours) +
  ggplot2::facet_grid(~ lesion_type, scales = "free", space = "free")
dev.off()

col_data %>%
  dplyr::filter(partition_lineage != "immune") %>%
  dplyr::left_join(
    col_data %>%
      dplyr::count(partition, vhl_status) %>%
      dplyr::group_by(partition) %>%
      dplyr::mutate(prop = (n / sum(n)) * 100) %>%
      dplyr::filter(vhl_status == "3p loss") %>%
      dplyr::select(partition, prop, n)
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = reorder(partition, -prop),
                               fill = vhl_status)) +
  ggplot2::geom_bar(position = "fill") +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 270, hjust = 0)) +
  ggplot2::scale_fill_manual(values = vhl_status_colours) +
  ggplot2::facet_grid(~partition_lineage, scales = "free", space = "free")

# plot driver losses on umap
pdf(paste0(out$base, "chr3_status.pdf"), width = 16, height = 13)
c("has_loss_chr3", "VHL_mut", "partition_annot", "sample") %>%
  purrr::map(function(var) {
    p <- monocle3::plot_cells(cds, color_cells_by = var, label_cell_groups = F)
    if(var %in% c("partition_annot", "sample")) { p + ditto_colours } else { p }
  }) %>%
  patchwork::wrap_plots(ncol = 2)
dev.off()

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



