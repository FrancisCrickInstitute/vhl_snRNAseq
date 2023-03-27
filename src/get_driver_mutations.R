base_dir <- ifelse(Sys.info()["nodename"] == "Alexs-MacBook-Air-2.local",
                   "/Volumes/TracerX/",
                   "/camp/project/tracerX/")
setwd(paste0(base_dir,"working/VHL_GERMLINE/tidda/vhl/"))
devtools::load_all()

# get sample subset
sample_subset <- readr::read_tsv(paste0(base_dir, "/working/VHL_GERMLINE/tidda/parse_pipeline/expdata/sample_subset.tsv")) %>%
  dplyr::mutate(lesion_id = gsub(".*\\_", "", sample))

# sample metadata
sample_metadata <-
  readr::read_csv(paste0(base_dir, '/working/CMELA/alex/work/ucl/data/vhl/sampledb/20230218_sample_data.csv')) %>%
  janitor::clean_names() %>%
  dplyr::right_join(sample_subset)

# clinical metadata
clinical_metadata <-
  readr::read_csv(paste0(base_dir, '/working/VHL_GERMLINE/tidda/parse_pipeline/expdata/clinical_metadata.csv')) %>%
  dplyr::rename() %>%
  # fix N023 -> N23
  dplyr::mutate(nih_pid = dplyr::case_when(nih_id == "N023" ~ "N23",
                                           TRUE ~ nih_id)) %>%
  dplyr::filter(nih_pid %in% sample_subset$nih_pid)

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
mut <-
  readr::read_tsv(paste0(base_dir, '/working/CMELA/alex/work/ucl/data/vhl/mutect2/master/2023.mutations.w.cn.and.ccfv4.w.rescue.statusv4.tsv')) %>%
  janitor::clean_names() %>%
  dplyr::mutate(nih_pid = gsub("K891", "N23", gsub("NIH\\_", "", patient)),
                tumour_sample_barcode = gsub("\\-", "", tumor_sample_barcode),
                lesion_id = gsub("\\.|d.*", "", tumour_sample_barcode)) %>%
  dplyr::filter(filter == "PASS") %>%
  dplyr::select(gene = funco_gene, tumour_sample_barcode,
                nih_pid, lesion_id, abs_ccf, report_ccf, purity, ploidy, wgii) %>%
  dplyr::inner_join(sample_subset, multiple = "all")

# arm-level CNV
scna_events <-
  readr::read_csv(paste0(base_dir, '/working/VHL_GERMLINE/SCOTT/Out/20230213_CN/20230301_arm_SCNA.csv')) %>%
  dplyr::rename(event = `...1`) %>%
  tidyr::pivot_longer(-event, names_to = "sample") %>%
  dplyr::mutate(sample = gsub("NIH\\_", "", sample)) %>%
  dplyr::inner_join(sample_subset)

# CNV drivers
scna_drivers <-
  readr::read_csv(paste0(base_dir, '/working/VHL_GERMLINE/SCOTT/Out/20230213_CN/20230301_seg_tab_with_drivers.csv')) %>%
  dplyr::select(-`...1`) %>%
  janitor::clean_names() %>%
  dplyr::transmute(nih_pid = gsub("NIH\\_", "", patient),
                   tumour_sample_barcode = gsub("\\.", "", sample),
                   lesion_id = gsub("\\.|d.*", "", tumour_sample_barcode),
                   ccf = cancer_cell_frac,
                   driver_scna) %>%
  dplyr::filter(!is.na(driver_scna))

# CNV drivers by cytoband
scna_cyto_drivers <-
  readr::read_csv(paste0(base_dir, '/working/VHL_GERMLINE/SCOTT/Out/20230213_CN/20230301_driver_cytobands.csv')) %>%
  dplyr::rename(event = `...1`) %>%
  tidyr::pivot_longer(-event, names_to = "sample") %>%
  dplyr::mutate(sample = gsub("NIH\\_", "", sample)) %>%
  dplyr::inner_join(sample_subset)

