base_dir <- ifelse(Sys.info()["nodename"] == "Alexs-MacBook-Air-2.local",
                   "/Volumes/TracerX/",
                   "/camp/project/tracerX/")
setwd(paste0(base_dir,"working/VHL_GERMLINE/tidda/vhl/"))
devtools::load_all()

out <- get_out("out/230210_A01366_0351_AHNHCFDSX5_x_221202_A01366_0326_AHHTTWDMXY/hg38/comb_x_SHE5052A9_S101/all-well/DGE_filtered/unintegrated/")

# get sample subset
sample_subset <-
  readr::read_tsv(paste0(base_dir, "/working/VHL_GERMLINE/tidda/parse_pipeline/expdata/sample_subset.tsv")) %>%
  dplyr::mutate(lesion_id = gsub(".*\\_", "", sample))

# sample metadata
sample_metadata <-
  readr::read_csv(paste0(base_dir, '/working/CMELA/alex/work/ucl/data/vhl/sampledb/20230218_sample_data.csv')) %>%
  janitor::clean_names() %>%
  dplyr::right_join(sample_subset) %>%
  # get n samples per patient
  dplyr::group_by(nih_pid) %>%
  dplyr::mutate(n_samples = dplyr::n()) %>%
  dplyr::ungroup()

# clinical metadata
clinical_metadata <-
  readr::read_csv(paste0(base_dir, '/working/VHL_GERMLINE/tidda/parse_pipeline/expdata/clinical_metadata.csv')) %>%
  # fix N023 -> N23
  dplyr::mutate(nih_pid = dplyr::case_when(nih_id == "N023" ~ "N23",
                                           TRUE ~ nih_id),
                vhl_germline_mutation = vhl_mutation_type) %>%
  dplyr::filter(nih_pid %in% sample_subset$nih_pid) %>%
  dplyr::select(-nih_id)

# pp (purity, ploidy, wgII)
pp <-
  readr::read_tsv(paste0(base_dir, '/working/CMELA/alex/work/ucl/data/vhl/absolute/master/2023.all.summary.unjoinedv4.w.wgiiv3.tsv')) %>%
  janitor::clean_names() %>%
  dplyr::mutate(nih_pid = gsub("NIH\\_", "", gsub("K891", "N23", patient)),
                lesion_id = gsub("\\.|d.*", "", sample),
                tumour_sample_barcode = gsub("\\.", "", sample)) %>%
  dplyr::select(nih_pid, lesion_id, tumour_sample_barcode, purity, ploidy, wgii) %>%
  dplyr::inner_join(sample_subset)

# mutations (with CCF)
muts <-
  readr::read_tsv(paste0(base_dir, '/working/CMELA/alex/work/ucl/data/vhl/mutect2/master/2023.mutations.w.cn.and.ccfv4.w.rescue.statusv4.tsv')) %>%
  janitor::clean_names() %>%
  dplyr::mutate(nih_pid = gsub("K891", "N23", gsub("NIH\\_", "", patient)),
                tumour_sample_barcode = gsub("\\-", "", tumor_sample_barcode),
                lesion_id = gsub("\\.|d.*", "", tumour_sample_barcode)) %>%
  dplyr::filter(filter == "PASS") %>%
  dplyr::select(mut = funco_gene, tumour_sample_barcode,
                nih_pid, lesion_id,
                mut_ccf = report_ccf) %>%
  dplyr::distinct() %>%
  dplyr::inner_join(sample_subset)
# get repeat muts / muts from Turajlic 2018 Fig1
turajlic_2018_muts <-
  c("VHL", "PBRM1", "SETD2", "PIK3CA", "MTOR", "PTEN", "KDM5C",
    "CSMD3", "BAP1", "TP53", "TSC1", "TSC2", "ARID1A", "TCEB1")
muts_to_plot <-
  muts %>%
  dplyr::count(mut) %>%
  dplyr::filter(n > 1 | n == 1 & mut %in% turajlic_2018_muts) %>%
  dplyr::pull(mut)
muts <-
  muts %>%
  dplyr::filter(mut %in% muts_to_plot) %>%
  # >1 ARID1A mut entry for N090_V127, take most clonal
  dplyr::arrange(desc(mut_ccf)) %>%
  dplyr::group_by(lesion_id, mut) %>%
  dplyr::slice_max(mut_ccf) %>%
  tidyr::pivot_wider(names_from = mut, values_from = mut_ccf)

# arm-level CNV
scnas <-
  readr::read_csv(paste0(base_dir, '/working/VHL_GERMLINE/SCOTT/Out/20230213_CN/20230301_arm_SCNA.csv')) %>%
  dplyr::rename(scna = `...1`) %>%
  tidyr::pivot_longer(-scna, names_to = "sample") %>%
  dplyr::mutate(sample = gsub("NIH\\_", "", sample)) %>%
  dplyr::inner_join(sample_subset) %>%
  dplyr::filter(value == 1) %>%
  dplyr::select(-value)

# CNV drivers
scna_drivers <-
  readr::read_csv(paste0(base_dir, '/working/VHL_GERMLINE/SCOTT/Out/20230213_CN/20230301_seg_tab_with_drivers.csv')) %>%
  dplyr::select(-`...1`) %>%
  janitor::clean_names() %>%
  dplyr::transmute(nih_pid = gsub("NIH\\_", "", patient),
                   tumour_sample_barcode = gsub("\\.", "", sample),
                   lesion_id = gsub("\\.|d.*", "", tumour_sample_barcode),
                   ccf = cancer_cell_frac,
                   scna_driver = driver_scna) %>%
  dplyr::filter(!is.na(scna_driver))

# CNV drivers by cytoband
scna_cyto_drivers <-
  readr::read_csv(paste0(base_dir, '/working/VHL_GERMLINE/SCOTT/Out/20230213_CN/20230301_driver_cytobands.csv')) %>%
  dplyr::rename(scna_cyto_driver = `...1`) %>%
  tidyr::pivot_longer(-scna_cyto_driver, names_to = "sample") %>%
  dplyr::mutate(sample = gsub("NIH\\_", "", sample),
                scna_cyto_driver_direction = dplyr::case_when(grepl("Gain", scna_cyto_driver) ~ "gain",
                                                              grepl("Loss", scna_cyto_driver) ~ "loss")) %>%
  dplyr::inner_join(sample_subset) %>%
  dplyr::filter(value == 1) %>%
  dplyr::select(-value)

# merge
genotypes <-
  clinical_metadata %>%
  # add sample metadata
  dplyr::full_join(sample_metadata, multiple = "all") %>%
  # add pp
  dplyr::full_join(pp, multiple = "all") %>%
  # add muts
  dplyr::full_join(muts, multiple = "all") %>%
  # add scna_cyto_drivers
  dplyr::full_join(scna_cyto_drivers, multiple = "all")

# heatmap
p_dat <-
  genotypes %>%
  dplyr::select(
    # patient-level variables
    nih_pid, sex, vhl_germline_mutation,
    # sample-level variables
    sample, age_at_surgery, fuhrman_grade, tumour_size, lesion_type,
    # whole genome-level variables
    purity, ploidy, wgii,
    # mutation-level variables
    dplyr::all_of(muts_to_plot),
    # scna cyto drivers-level variables
    scna_cyto_driver, scna_cyto_driver_direction
  ) %>%
  dplyr::mutate(sample = factor(sample, ordered = T))

plot_variable <- function(p_dat, variable, show_col_names = F) {
  p <-
    p_dat %>%
    dplyr::select(sample, dplyr::all_of(variable)) %>%
    tidyr::pivot_longer(-sample) %>%
    ggplot2::ggplot(
      ggplot2::aes(x = sample, y = name, fill = value)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::theme_void() +
    ggplot2::labs(fill = variable) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(hjust = 1),
                   legend.position = "none")
  if (show_col_names == T) {
    p <- p +
      ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 0, angle = 90)) +
      ggplot2::scale_x_discrete(position = "top")
  }
  if (is.numeric(unlist(p_dat[, variable , drop=T]))) {
    p <- p + ggplot2::scale_fill_gradient(
      low = "#fae3f5",
      high = "#611d52",
      na.value = 'white')
  } else {
    p <- p + ditto_colours
  }
  p
}


# variable levels
variable_lvls <- list(
  patient = list("sex", "vhl_germline_mutation"),
  sample = list("age_at_surgery", "fuhrman_grade", "tumour_size", "lesion_type"),
  wg = list("purity", "ploidy", "wgii"),
  mut = list(muts = muts_to_plot)
  )

# patient-level variables
p_patients <-
  plot_variable(p_dat, "nih_pid", show_col_names = T)
p_variables <-
  names(variable_lvls) %>%
  purrr::map(function(lvl) {
    variable_lvls[[lvl]] %>%
      purrr::map(function(variable) {
        plot_variable(p_dat, variable)
      }) %>% patchwork::wrap_plots(ncol = 1)
  }) %>% patchwork::wrap_plots(ncol = 1, heights = c(2,4,3,8))
p_cn <-
  p_dat %>%
  dplyr::select(sample, dplyr::all_of(c("scna_cyto_driver", "scna_cyto_driver_direction"))) %>%
  dplyr::mutate(scna_cyto_driver = tidyr::replace_na(scna_cyto_driver, "")) %>%
  ggplot2::ggplot(
    ggplot2::aes(x = sample, y = scna_cyto_driver, fill = scna_cyto_driver_direction)) +
  ggplot2::geom_tile(colour = "white") +
  ggplot2::theme_void() +
  ggplot2::labs(fill = variable) +
  ggplot2::theme(legend.position = "none",
                 axis.text.y = ggplot2::element_text(hjust = 1)) +
  ggplot2::scale_fill_manual(values = c("loss" = "#1c429c", "gain" = "#a31b0f"),
                             na.value = "white")

pdf(paste0(out$base, "sample_heatmap.pdf"), width = 5, height = 7)
list(p_patients, p_variables, p_cn) %>%
  patchwork::wrap_plots(ncol = 1, heights = c(1, 16, 12))
dev.off()

