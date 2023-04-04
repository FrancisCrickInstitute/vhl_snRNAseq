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
  paste0(base_dir, "working/CMELA/alex/work/ucl/data/vhl/sampledb/20230218_sample_data.csv") %>%
  readr::read_csv() %>%
  janitor::clean_names() %>%
  # fix lesion id
  dplyr::mutate(lesion_id = region_id) %>%
  dplyr::inner_join(sample_subset) %>%
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

# growth rates
growth_rates <-
  readr::read_csv(paste0(base_dir, '/working/VHL_GERMLINE/SCOTT/CLINICAL_ANNOTATION/20230329_growth_rate.csv')) %>%
  janitor::clean_names() %>%
  # get nih_sam_id
  dplyr::mutate(nih_sam_id = substr(ut_number, 1, 13)) %>%
  # get nih_pid and sample
  dplyr::inner_join(
   sample_metadata %>%
     dplyr::select(nih_sam_id, nih_pid, sample), .
  ) #%>%
  # # convert to long format
  # dplyr::mutate(lesion_types = gsub("\\[|\\]", "", lesion_types)) %>%
  # tidyr::separate_wider_delim(lesion_types, ", ", names = paste0("lesion_types_", 1:100),
  #                             too_few = "align_start") %>%
#   dplyr::rename_with(.cols = dplyr::starts_with("size"),
#                      .fn = ~ gsub("size", "size\\_", .x)) %>%
#   dplyr::rename_with(.cols = dplyr::starts_with("rate"),
#                      .fn = ~ gsub("rate", "rate\\_", .x)) %>%
#   dplyr::mutate(across(tidyselect::starts_with(c("size_",
#                                                  "days_since_first_",
#                                                  "rate_",
#                                                  "lesion_type_")), as.character)) %>%
#   tidyr::pivot_longer(cols = tidyselect::starts_with(c("size_",
#                                                        "days_since_first_",
#                                                        "rate_",
#                                                        "lesion_type_")),
#                       names_to = c("name", "timepoint"),
#                       names_pattern = "(.*)\\_(.*)") %>%
#   dplyr::filter(!is.na(value)) %>%
#   tidyr::pivot_wider()  %>%
#   readr::type_convert() %>%
#   dplyr::group_by(ut_number) %>%
#   dplyr::select(-total_growth_rate, -crick_id)
# growth_rates <- growth_rates %>%
#   dplyr::left_join(
#     growth_rates %>%
#       dplyr::filter(timepoint == 1 | timepoint == max(timepoint)) %>%
#       dplyr::transmute(ut_number,
#                        name = dplyr::case_when(timepoint == 1 ~ "size_start",
#                                             TRUE ~ "size_end"),
#                        value = size) %>%
#       tidyr::pivot_wider()
#   ) %>%
#   dplyr::mutate(annual_growth_rate = (size_end - size_start) / max(days_since_first) * 365)
growth_rates <-
  growth_rates %>%
  dplyr::distinct(nih_pid, sample, total_growth_rate)

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
  dplyr::filter(mut %in% turajlic_2018_muts) %>%
  dplyr::pull(mut)
muts <-
  muts %>%
  dplyr::filter(mut %in% muts_to_plot) %>%
  # >1 ARID1A mut entry for N090_V127, take most clonal
  dplyr::arrange(desc(mut_ccf)) %>%
  dplyr::group_by(lesion_id, mut) %>%
  dplyr::slice_max(mut_ccf) %>%
  tidyr::pivot_wider(names_from = mut, values_from = mut_ccf)
muts %>%
  dplyr::select(-tumour_sample_barcode, -lesion_id) %>%
  readr::write_tsv(paste0(out$base, "mutect2_mutations.tsv"))

# CNV drivers
scna_drivers <-
  readr::read_csv(paste0(base_dir, '/working/VHL_GERMLINE/SCOTT/Out/20230213_CN/20230301_seg_tab_with_drivers.csv')) %>%
  dplyr::select(-`...1`) %>%
  janitor::clean_names() %>%
  dplyr::transmute(nih_pid = gsub("NIH\\_", "", patient),
                   tumour_sample_barcode = gsub("\\.", "", sample),
                   lesion_id = gsub("\\.|d.*", "", tumour_sample_barcode),
                   driver_scna_ccf = cancer_cell_frac,
                   scna_driver = driver_scna) %>%
  dplyr::filter(!is.na(scna_driver))
all_scna_drivers <-
  readr::read_csv(paste0(base_dir, '/working/VHL_GERMLINE/SCOTT/Out/20230213_CN/20230301_seg_tab_with_drivers.csv')) %>%
  dplyr::select(-`...1`) %>%
  janitor::clean_names() %>%
  dplyr::filter(!is.na(driver_scna)) %>%
  dplyr::mutate(
    scna_driver = driver_scna %>% tolower(),
    x = gsub(".*\\_", "", scna_driver)
  ) %>%
  dplyr::transmute(chr = as.numeric(gsub("p.*|q.*", "", x)),
                   arm = gsub("[0-9]|\\.", "", x),
                   direction = gsub("\\_.*", "", scna_driver),
                   driver = paste0(direction, "_", chr, arm)) %>%
  dplyr::distinct()

# arm-level CNV
scnas_all <-
  readr::read_csv(paste0(base_dir, '/working/VHL_GERMLINE/SCOTT/Out/20230213_CN/20230301_arm_SCNA.csv')) %>%
  dplyr::rename(scna_arm = `...1`) %>%
  tidyr::pivot_longer(-scna_arm, names_to = "sample") %>%
  dplyr::mutate(sample = gsub("NIH\\_", "", sample),
                scna_arm_direction = gsub("\\_.*", "", scna_arm)) %>%
  dplyr::inner_join(sample_subset)
scnas <-
  scnas_all %>%
  dplyr::filter(value == 1) %>%
  dplyr::select(-value) %>%
  dplyr::group_by(scna_arm) %>%
  dplyr::mutate(scna_arm_n = dplyr::n()) %>%
  dplyr::ungroup()

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
  dplyr::select(-value) %>%
  dplyr::group_by(scna_cyto_driver) %>%
  dplyr::mutate(scna_cyto_driver_n = dplyr::n()) %>%
  dplyr::arrange(desc(scna_cyto_driver_n), scna_cyto_driver_direction)

# merge
genotypes <-
  clinical_metadata %>%
  # add sample metadata
  dplyr::full_join(sample_metadata, multiple = "all") %>%
  # add growth rates
  dplyr::full_join(growth_rates, multiple = "all") %>%
  # add pp
  dplyr::full_join(pp) %>%
  # add muts
  dplyr::full_join(muts, multiple = "all") %>%
  # add scna_cyto_drivers
  dplyr::full_join(scna_cyto_drivers, multiple = "all") %>%
  # add scnas
  dplyr::full_join(scnas, multiple = "all")

# write full clinical annotations file
genotypes %>%
  readr::write_tsv(paste0(out$base, "clinical_and_genotype_annotations.tsv"))
growth_rates %>%
  readr::write_tsv(paste0(out$base, "growth_rates.tsv"))

# heatmap
p_dat <-
  genotypes %>%
  dplyr::select(
    # patient-level variables
    nih_pid, sex, vhl_germline_mutation,
    # sample-level variables
    sample, age_at_surgery, fuhrman_grade, tumour_size, lesion_type, annual_growth_rate,
    # whole genome-level variables
    purity, ploidy, wgii,
    # mutation-level variables
    dplyr::all_of(muts_to_plot),
    # scna cyto drivers-level variables
    dplyr::starts_with("scna_cyto_driver"),
    # scna arm-level variabls
    dplyr::starts_with("scna_arm")
  ) %>%
  dplyr::mutate(sample = factor(sample, ordered = T))

plot_variable <- function(p_dat, variable, lvl, variable_colours,
                          show_col_names = F, return_legend = F) {
  p <-
    p_dat %>%
    dplyr::select(sample, dplyr::all_of(variable)) %>%
    tidyr::pivot_longer(-sample) %>%
    ggplot2::ggplot(
      ggplot2::aes(x = sample, y = name, fill = value)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::theme_void() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(hjust = 1))
  if (show_col_names == T) {
    p <- p +
      ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 0, angle = 90)) +
      ggplot2::scale_x_discrete(position = "top")
  }
  if (all(variable %in% names(variable_colours))) {
    p <- p +
      ggplot2::scale_fill_manual(values = variable_colours[[variable]],
                                 na.value = "white")
  } else if (is.numeric(unlist(p_dat[, variable , drop=T]))) {
    p <- p + ggplot2::scale_fill_gradient(
      low = "#fae3f5",
      high = "#611d52",
      na.value = 'white')
  } else {
    p <- p + ditto_colours
  }
  if (return_legend == T) {
    lemon::g_legend(
      p +
        ggplot2::labs(fill = ifelse(lvl == "mut", "mut ccf", variable)) +
        ggplot2::theme(legend.position = "bottom")
    )
  } else {
    p + ggplot2::theme(legend.position = "none")
  }
}

# variable levels
variable_lvls <- list(
  patient = list("sex", "vhl_germline_mutation"),
  sample = list("fuhrman_grade", "lesion_type", "age_at_surgery", "tumour_size",
                "annual_growth_rate"
                ),
  wg = list("purity", "ploidy", "wgii"),
  mut = list(muts = muts_to_plot)
  )

# get colours
avail_colours <- dittoSeq::dittoColors()
variable_colours <- list()
p_leg <- list()
for(variable in c("nih_pid", unlist(variable_lvls))) {
  if (!is.numeric(unlist(p_dat[,variable]))) {
    cat(variable, "categorical\n")
    n_colours <- length(unique(p_dat[,variable,drop=T]))
    variable_colours[[variable]] <-
      avail_colours[1:n_colours]
    avail_colours <- avail_colours[-c(1:n_colours)]
  }
}
variable_colours[["scna_cyto_driver"]] <-
  variable_colours[["scna_arm"]] <-
  c("loss" = "#1c429c", "gain" = "#a31b0f")
variable_colours[["sex"]] <-
  c("f" = "#ffb3fc", "m" = "#7eaaed")

# patient level variables
p_patients <-
  plot_variable(p_dat, "nih_pid", variable_colours, "patient", show_col_names = T)
p_leg[["patients"]] <-
  plot_variable(p_dat, "nih_pid", variable_colours, "patient", return_legend = T)

# patient / sample / wg / mut level variables
p_variables <-
  names(variable_lvls) %>%
  purrr::map(function(lvl) {
    variable_lvls[[lvl]] %>%
      purrr::map(function(variable) {
        cat(lvl, variable, "\n")
        p_leg[[ifelse(lvl == "mut", lvl, variable)]] <<-
          plot_variable(p_dat, variable, lvl, variable_colours, return_legend = T)
        plot_variable(p_dat, variable, lvl, variable_colours)
      }) %>% patchwork::wrap_plots(ncol = 1)
  }) %>% patchwork::wrap_plots(ncol = 1, heights = c(2,4,3,4.5))

plot_cn <-function(p_dat, prefix) {
  p_dat %>%
    dplyr::select(sample,
                  y = dplyr::all_of(prefix),
                  y_direction = dplyr::all_of(paste0(prefix, "_direction")),
                  y_n = dplyr::all_of(paste0(prefix, "_n"))) %>%
    dplyr::mutate(y_no_na = tidyr::replace_na(as.character(y), "")) %>%
    dplyr::arrange(desc(y_direction), desc(y_n)) %>%
    dplyr::mutate(order = dplyr::row_number()) %>%
    ggplot2::ggplot(
      ggplot2::aes(x = sample, y = reorder(y_no_na, -order), fill = y_direction)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::theme_void() +
    ggplot2::labs(fill = "direction") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(hjust = 1),
                   legend.position = "bottom") +
    ggplot2::scale_fill_manual(values = variable_colours[[prefix]],
                               na.value = "white")
}

# scna arm level variables
p_scna_arm <-
  plot_cn(p_dat, "scna_arm") +
  ggplot2::theme(legend.position = "none")

# chromosome lengths
chr_len <-
  readr::read_tsv("data/gencode/chromosomes.tsv") %>%
  dplyr::transmute(y_n = as.numeric(gsub("chr", "", chr)),
                   ratio = length / max(length),
                   y_cumsum = cumsum(ratio),
                   y_pos = y_cumsum - 0.5 * ratio) %>%
  dplyr::filter(!is.na(y_n)) %>%
  dplyr::left_join(tibble::tibble(
    y = apply(expand.grid(1:22, c("p", "q")), 1, paste, collapse="") %>% gsub(" ", "", .),
    y_n = rep(1:22, 2)),
    multiple = "all"
  ) %>%
  dplyr::mutate(y_pos = dplyr::case_when(grepl("p", y) ~ y_cumsum - 0.25*ratio,
                                         grepl("q", y) ~ y_cumsum + 0.25*ratio))

# infercnv sample order
infercnv_sample_order <-
  tibble::tibble(
    sample = c(
      "N059_V103", "N059_V102A", "N059_V003", "N059_M001",
      "N088_V108", "N088_V106", "N088_V008", "N088_V006", "N088_V004",
      "N045_V008C", "N045_V004", "N045_V003",
      "N090_V127", "N090_V126", "N090_V124D", "N090_V124A", "N090_V116",
      "K891_V014")
  ) %>% dplyr::mutate(sample_order = dplyr::row_number())


# full chrom arm level events
p_scna_all_dat <-
  tibble::tibble(
    y = apply(expand.grid(1:22, c("p", "q")), 1, paste, collapse="") %>% gsub(" ", "", .),
    y_n = rep(1:22, 2),
    A = 1
  ) %>%
  dplyr::cross_join(
    p_dat %>% dplyr::distinct(nih_pid, sample)
  ) %>%
  dplyr::inner_join(
    infercnv_sample_order
  ) %>%
  dplyr::left_join(
    genotypes %>%
      dplyr::transmute(sample, nih_pid,
                       y = gsub(".*\\_", "", scna_arm),
                       y_direction = scna_arm_direction,
                       y_n = as.numeric(gsub("q|p", "", y)),
                       arm = gsub("[0-9]", "", y),
                       value = 1) %>%
      dplyr::filter(!is.na(y)) %>%
      dplyr::distinct()
  ) %>%
  dplyr::left_join(
    all_scna_drivers %>%
      dplyr::transmute(y_n = chr, arm, y_direction = direction, driver = "*")
  ) %>%
  dplyr::distinct() %>%
  dplyr::arrange(y_n) %>%
  dplyr::mutate(order = dplyr::row_number(),
                y_no_na = tidyr::replace_na(y, "")) %>%
  dplyr::left_join(chr_len)
p_scna_all <-
  p_scna_all_dat %>%
  dplyr::filter(!is.na(sample)) %>%
  ggplot2::ggplot(ggplot2::aes(x = reorder(sample, -sample_order),
                               y = reorder(y_no_na, -order),
                               height = 1)) +
  ggplot2::geom_tile(ggplot2::aes(fill = y_direction), colour = "white") +
  ggplot2::theme_void() +
  ggplot2::labs(fill = "direction") +
  ggplot2::theme(axis.text.y = ggplot2::element_text(hjust = 1, size = 7),
                 legend.position = "none") +
  ggplot2::scale_fill_manual(values = variable_colours[["scna_arm"]],
                             na.value = "white")  +
  ggplot2::geom_line(data = data.frame(x = c(0, length(unique(p_scna_all_dat$sample))) + 0.5,
                                       y = rep(seq(1, length(unique(p_scna_all_dat$y)), by = 2),
                                               each = 2) - 0.5),
                     ggplot2::aes(x, y, group = y), colour = "grey") +
  ggplot2::geom_text(data = p_scna_all_dat,
                     ggplot2::aes(x = sample, y = reorder(y_no_na, -order), label = driver),
                     colour = "white", vjust = 0.8)

# scna cyto driver level variables
p_scna_cyto_driver <-
  plot_cn(p_dat, "scna_cyto_driver")
p_leg[["scna_cyto_driver"]] <- lemon::g_legend(p_scna_cyto_driver)
p_scna_cyto_driver <- p_scna_cyto_driver + ggplot2::theme(legend.position = "none")

# save plots
pdf(paste0(out$base, "sample_heatmap.pdf"), width = 6.6, height = 8)
list(p_patients, p_variables, p_scna_cyto_driver) %>%
  patchwork::wrap_plots(ncol = 1, heights = c(1, 16, 12))
dev.off()

pdf(paste0(out$base, "sample_heatmap_legends.pdf"), width = 6, height = 10)
p_leg %>% patchwork::wrap_plots(ncol = 1)
dev.off()

pdf(paste0(out$base, "sample_heatmap_scna_arm.pdf"), width = 6, height = 8)
list(p_patients, p_scna_arm) %>%
  patchwork::wrap_plots(ncol = 1, heights = c(1, 25))
dev.off()

pdf(paste0(out$base, "sample_heatmap_scna_arm_full.pdf"), width = 6, height = 8)
list(p_patients, p_scna_all) %>%
  patchwork::wrap_plots(ncol = 1, heights = c(1, 30))
dev.off()

pdf(paste0(out$base, "sample_heatmap_scna_arm_full_labelled.pdf"), width = 20, height = 8)
p_scna_all + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, size = 7))
dev.off()
