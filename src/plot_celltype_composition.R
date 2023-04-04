# config
base_dir <- ifelse(Sys.info()["nodename"] == "Alexs-MacBook-Air-2.local",
                   "/Volumes/TracerX/",
                   "/camp/project/tracerX/")
setwd(paste0(base_dir,"working/VHL_GERMLINE/tidda/vhl/"))
devtools::load_all()
out_dir <- "out/230210_A01366_0351_AHNHCFDSX5_x_221202_A01366_0326_AHHTTWDMXY/hg38/comb_x_SHE5052A9_S101/all-well/DGE_filtered/unintegrated/"
out <- get_out(out_dir)


# plot annotations on umap
cds <- readRDS(paste0(out$cache, "cds_annotated.rds"))
pdf(paste0(out$base, "annotations_on_umap.pdf"), width = 15, height = 6.5)
dittoSeq::dittoDimPlot(cds, "partition_lineage") +
  umap_void_theme + ggplot2::theme(legend.position = "bottom") +
  dittoSeq::dittoDimPlot(cds, "partition_annot") +
    umap_void_theme + ggplot2::theme(legend.position = "bottom") +
  dittoSeq::dittoDimPlot(cds, "cluster_annot") +
    umap_void_theme + ggplot2::theme(legend.position = "bottom")
dev.off()

# plot heatmap of tumour interesting markers
tum_markers <-
  markers[c("HIF",
            "proximal tubule", "proximal tubule (progenitor)", "proximal tubule (injured)",
            "tumour suppressors")]
agg_mat <-
  monocle3::aggregate_gene_expression(
    cds[rownames(cds) %in% unlist(tum_markers),],
    gene_group_df = tibble::enframe(tum_markers) %>% tidyr::unnest(value) %>% dplyr::select(value, name),
    cell_group_df = SummarizedExperiment::colData(cds)[, "cluster", drop = F] %>%
      tibble::as_tibble(rownames = "cell"))
agg_mat_CA9 <-
  monocle3::aggregate_gene_expression(
    cds[rownames(cds) == "CA9",],
    cell_group_df = SummarizedExperiment::colData(cds)[, "cluster", drop = F] %>%
      tibble::as_tibble(rownames = "cell"))
agg_mat <- rbind(agg_mat, agg_mat_CA9)
agg_mat %>%
  tibble::as_tibble(rownames = "module") %>%
  dplyr::select(dplyr::all_of(order))

order <-
  c("3", "5", "6", "7", "9", "10", "11", "12", "15", "17", "18", "19", "21", "23", "1", "2", "14", "16", "20", "24", "25", "4", "8", "13", "22")
pheatmap::pheatmap(agg_mat[,order], cluster_cols = F, border_color = NA, fontsize = 12)

agg_mat <-
  monocle3::aggregate_gene_expression(
    cds[rownames(cds) == "CA9",,drop=F],
    cell_group_df = SummarizedExperiment::colData(cds)[, "cluster", drop = F] %>%
      tibble::as_tibble(rownames = "cell"))


# hif de
p_dat <-
  cbind(
    as.matrix(
      cds[rownames(cds) %in% markers$HIF,
          SummarizedExperiment::colData(cds)$partition_lineage == "malignant",
          drop = F]@assays@data$counts) %>%
      rowMeans(),
    as.matrix(
      cds[rownames(cds) %in% markers$HIF,
          SummarizedExperiment::colData(cds)$partition_lineage == "normal",
          drop = F]@assays@data$counts) %>%
      rowMeans()
  ) %>%
  tibble::as_tibble(rownames = "gene") %>%
  dplyr::rename(malignant = V1, normal = V2) %>%
  tidyr::pivot_longer(-gene)

# Shapiro-Wilk normality test for the differences
d <- with(p_dat,
          value[name == "normal"] - value[name == "malignant"])
shapiro.test(d)

# paired samples t test
t.test(
  dplyr::filter(p_dat, name == "normal")$value, dplyr::filter(p_dat, name == "malignant")$value,
  paired = T, alternative = "two.sided"
) -> hif_paired_t

p_dat %>%
  dplyr::mutate(order = dplyr::case_when(name == "normal" ~ 1, TRUE ~ 2)) %>%
  ggplot2::ggplot(ggplot2::aes(x = reorder(name, order), y = log(value))) +
  ggplot2::geom_jitter(height = 0, width = 0.2) +
  ggplot2::geom_boxplot(alpha = 0.5) +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 12))

# get celltype annotations
cts <- readr::read_tsv(paste0(out$base, "cell_annotations.tsv")) %>%
  dplyr::mutate(patient = gsub("\\_.*", "", ifelse(grepl("K891", sample), "N23", sample))) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(n = dplyr::n())

cluster_annot_counts <-
  SummarizedExperiment::colData(cds) %>% tibble::as_tibble() %>%
  dplyr::count(cluster_annot) %>%
  dplyr::mutate(cluster_annot_counts = paste0(cluster_annot, " (", n, ")")) %>%
  dplyr::arrange(-n)
SummarizedExperiment::colData(cds)$cluster_annot_counts <-
  dplyr::left_join(
    SummarizedExperiment::colData(cds) %>% tibble::as_tibble(),
    cluster_annot_counts
  ) %>% dplyr::pull(cluster_annot_counts)

partition_lineage_counts <-
  SummarizedExperiment::colData(cds) %>% tibble::as_tibble() %>%
  dplyr::count(partition_lineage) %>%
  dplyr::mutate(partition_lineage_counts = paste0(partition_lineage, " (", n, ")")) %>%
  dplyr::arrange(-n)
SummarizedExperiment::colData(cds)$partition_lineage_counts <-
  dplyr::left_join(
    SummarizedExperiment::colData(cds) %>% tibble::as_tibble(),
    partition_lineage_counts
  ) %>% dplyr::pull(partition_lineage_counts)

pdf(paste0(out$base, "cluster_annots_umap.pdf"), height = 6, width = 20)
dittoSeq::dittoDimPlot(cds, "cluster_annot_counts", size = 0.1) +
  ggplot2::coord_equal() +
  ggplot2::scale_colour_manual(ct_colours)
dittoSeq::dittoDimPlot(cds, "partition_lineage_counts", size = 0.1) +
  ggplot2::coord_equal()
dev.off()

# get counts (patient, sample, lineage, annot levels)
counts <-
  cts %>%
  dplyr::left_join(
    cts %>% dplyr::count(patient, sample, partition_lineage, cluster_annot, name = "annot_n")
  ) %>%
  dplyr::left_join(
    cts %>% dplyr::count(patient, sample, partition_lineage, name = "lineage_n")
  ) %>%
  dplyr::left_join(
    cts %>% dplyr::count(patient, name = "patient_n")
  ) %>%
  dplyr::left_join(
    cts %>% dplyr::count(patient, sample, name = "sample_n")
  ) %>%
  dplyr::left_join(
    cts %>%
      dplyr::group_by(patient, sample) %>%
      dplyr::mutate(n = dplyr::n(),
                    partition_lineage = paste0(partition_lineage, "_lineage")) %>%
      dplyr::group_by(partition_lineage, .add = T) %>%
      dplyr::summarise(n = dplyr::n()) %>%
      tidyr::pivot_wider(names_from = "partition_lineage",
                         values_from = "n",
                         names_prefix = "n_",
                         values_fill = 0) %>%
      dplyr::mutate(n_total = sum(dplyr::across(dplyr::starts_with("n_"))),
                    dplyr::across(dplyr::starts_with("n_"), ~ .x / n_total,
                                  .names = "prop_{.col}"))
  ) %>%
  dplyr::left_join(
    cts %>%
      dplyr::group_by(patient, sample) %>%
      dplyr::mutate(n = dplyr::n(),
                    cluster_annot = paste0(cluster_annot, "_cluster")) %>%
      dplyr::group_by(cluster_annot, .add = T) %>%
      dplyr::summarise(n = dplyr::n()) %>%
      tidyr::pivot_wider(names_from = "cluster_annot",
                         values_from = "n",
                         names_prefix = "n_",
                         values_fill = 0) %>%
      dplyr::mutate(n_total = sum(dplyr::across(dplyr::starts_with("n_"))),
                    dplyr::across(dplyr::starts_with("n_"), ~ .x / n_total,
                                  .names = "prop_{.col}"))
  ) %>%
  dplyr::mutate(
    sample_id = gsub(".*\\_", "", sample),
    partition_lineage = factor(partition_lineage, levels = c("immune", "normal", "malignant")),
    annot_prop = annot_n / sample_n)
counts$cluster_annot <-
  factor(counts$cluster_annot,
         levels = unique(
           counts$cluster_annot[order(counts$partition_lineage, counts$cluster_annot)]),
         ordered = TRUE)

# prepare data with ordering by unique sample
p_dat <-
  counts %>%
  dplyr::select(sample, sample_id, patient,
                dplyr::starts_with("n_"),
                dplyr::starts_with("prop_"),
                cluster_annot,
                partition_lineage) %>%
  dplyr::arrange(sample, n_malignant_lineage)

# plot all function
plot_all <- function(counts, y_title, no_legend = F, ...) {
  counts %>%
    ggplot2::ggplot(ggplot2::aes(x = reorder(sample, -as.numeric(prop_n_malignant_lineage)), fill = cluster_annot)) +
    ggplot2::geom_bar(...) +
    ditto_colours +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.95, vjust = 0.2),
                   panel.grid = ggplot2::element_blank(),
                   legend.position = ifelse(no_legend, "none", "right")) +
    ggplot2::labs(x = "sample", y = y_title, fill = "population")
}

# plot lineages function
plot_by_lineage <- function(p_dat, y_title, no_legend = F, ...) {
  p_dat %>%
    ggplot2::ggplot(ggplot2::aes(x = tidytext::reorder_within(sample_id, -n_malignant_lineage, patient), fill = cluster_annot)) +
    ggplot2::geom_bar(...) +
    ggplot2::facet_grid(partition_lineage ~ patient,
                        scales = "free", space = "free_x", switch = "x") +
    ditto_colours +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.95, vjust = 0.2),
                   panel.grid = ggplot2::element_blank(),
                   legend.position = ifelse(no_legend, "none", "right")) +
    ggplot2::labs(x = "sample", fill = "population", y = y_title)
}

# plots
p <- list()
p[["all"]] <-
  (plot_all(counts, "proportion of cells", no_legend = T, position = "fill") +
     plot_all(counts, "n cells"))
p[["lineages"]] <-
  (plot_by_lineage(p_dat, "proportion of cells", no_legend = T, position = "fill") +
  plot_by_lineage(p_dat, "n cells"))

# save proportions to output
counts %>%
  dplyr::distinct(
    patient, sample,
    dplyr::across(dplyr::starts_with("prop_")),
    dplyr::across(dplyr::starts_with("n_"))
  ) %>%
  readr::write_tsv(paste0(out$base, "sample_composition.tsv"))


count_annot_lvl <- function(counts, annot_i) {
  counts %>%
    dplyr::rename(annot = dplyr::matches(annot_i)) %>%
    dplyr::group_by(annot) %>%
    dplyr::summarise(n = dplyr::n(),
                     clusters = paste(unique(sort(cluster)), collapse = ","),
                     partitions = paste(unique(sort(partition)), collapse = ","),
                     annot_level = annot_i) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(prop = paste0(round(n / sum(n) * 100, 1), "%")) %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::select(annot_level, annot, n, prop, clusters, partitions)
}

# save annotation counts at different levels
c("partition_annot", "partition_lineage", "cluster_annot") %>%
  purrr::map(function(annot_i) {count_annot_lvl(counts, annot_i)}) %>%
  dplyr::bind_rows() %>%
  readr::write_tsv(paste0(out$base, "annot_counts.tsv"))

# save plots
pdf(paste0(out$base, "sample_composition.pdf"), height = 8, width = 12, onefile = T)
purrr::walk(p, print)
dev.off()

# colours <-
#   dittoSeq::dittoColors()[1:length(unique(SummarizedExperiment::colData(cds)$sample))]
# samples <- sort(unique(SummarizedExperiment::colData(cds)$sample))
# names(colours) <- c(samples[2:24], samples[1])
# dittoSeq::dittoDimPlot(cds, "sample") + ggplot2::scale_colour_manual(values = colours)
