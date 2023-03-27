base_dir <- ifelse(Sys.info()["nodename"] == "Alexs-MacBook-Air-2.local",
                   "/Volumes/TracerX/",
                   "/camp/project/tracerX/")
setwd(paste0(base_dir,"working/VHL_GERMLINE/tidda/vhl/"))
devtools::load_all()

out_dir <- "out/230210_A01366_0351_AHNHCFDSX5_x_221202_A01366_0326_AHHTTWDMXY/hg38/comb_x_SHE5052A9_S101/all-well/DGE_filtered/unintegrated/"

out <- get_out(out_dir)

# get celltype annotations
cts <- readr::read_tsv(paste0(out$base, "cell_annotations.tsv")) %>%
  dplyr::mutate(patient = gsub("\\_.*", "", ifelse(grepl("K891", sample), "N23", sample))) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(n = dplyr::n())

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

# plot annotations, by proportion, ordered by prop_n_malignant
counts %>%
  ggplot2::ggplot(ggplot2::aes(x = reorder(sample, -as.numeric(prop_n_malignant)), fill = cluster_annot)) +
  ggplot2::geom_bar(position = "fill") +
  ditto_colours +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.95, vjust = 0.2),
                 panel.grid = ggplot2::element_blank()) +
  ggplot2::labs(x = "sample", y = "", fill = "population")

# plot annotations, by count, ordered by n_malignant
counts %>%
  ggplot2::ggplot(ggplot2::aes(x = reorder(sample, -as.numeric(n_malignant)), fill = cluster_annot)) +
  ggplot2::geom_bar() +
  ditto_colours +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.95, vjust = 0.2),
                 panel.grid = ggplot2::element_blank()) +
  ggplot2::labs(x = "sample", y = "", fill = "population")

# prepare data with ordering by unique sample
p_dat <-
  counts %>%
  dplyr::select(sample, sample_id,
                dplyr::starts_with("n_"),
                dplyr::starts_with("prop_"),
                cluster_annot,
                partition_lineage, patient) %>%
  dplyr::arrange(sample, n_malignant_lineage)
p_labels <- unique(p_dat$sample_id)

# plot lineages and annotations, by proportion, ordered by n_malignant within patients
p_dat %>%
  dplyr::filter(partition_lineage != "malignant") %>%
  ggplot2::ggplot(ggplot2::aes(x = reorder(sample, - n_malignant_lineage), fill = cluster_annot)) +
  ggplot2::geom_bar(position = "fill") +
  ggplot2::facet_grid(partition_lineage ~ patient,
                      scales = "free", space = "free") +
  ditto_colours +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.95, vjust = 0.2),
                 panel.grid = ggplot2::element_blank()) +
  ggplot2::labs(x = "sample", y = "", fill = "population") +
  ggplot2::scale_x_discrete(labels = p_labels)

# plot lineages and annotations, by count, ordered by prop_n_malignant within patients
p_dat %>%
  ggplot2::ggplot(ggplot2::aes(x = reorder(sample, -p), fill = cluster_annot)) +
  ggplot2::geom_bar() +
  ggplot2::facet_grid(partition_lineage ~ patient,
                      scales = "free", space = "free") +
  ditto_colours +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.95, vjust = 0.2),
                 panel.grid = ggplot2::element_blank()) +
  ggplot2::labs(x = "sample", y = "", fill = "population") +
  ggplot2::scale_x_discrete(labels = p_labels)

# save proportions to output
counts %>%
  dplyr::distinct(
    patient, sample,
    dplyr::across(dplyr::starts_with("prop_")),
    dplyr::across(dplyr::starts_with("n_"))
  ) %>%
  readr::write_tsv(paste0(out$base, "sample_composition.tsv"))

