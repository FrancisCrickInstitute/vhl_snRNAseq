base_dir <- ifelse(Sys.info()["nodename"] == "Alexs-MacBook-Air-2.local",
                   "/Volumes/TracerX/",
                   "/camp/project/tracerX/")
setwd(paste0(base_dir,"working/VHL_GERMLINE/tidda/vhl/"))
devtools::load_all()

out_dir <- "out/230210_A01366_0351_AHNHCFDSX5_x_221202_A01366_0326_AHHTTWDMXY/hg38/comb_x_SHE5052A9_S101/all-well/DGE_filtered/unintegrated/"
out <- get_out(out_dir)

# get tumour cells
cell_annotations <- readr::read_tsv(paste0(out$base, "cell_annotations.tsv"))
tumour_cells <-
  cell_annotations %>%
  dplyr::filter(cluster_annot == "tumour") %>%
  dplyr::pull(cell)

# get sample subset
sample_subset <-
  readr::read_tsv(paste0(base_dir, "/working/VHL_GERMLINE/tidda/parse_pipeline/expdata/sample_subset.tsv"))$sample

# split samples by patient
patient_subsets <- sample_subset %>% split(get_patients_from_samples(sample_subset))

# run on each patient
purrr::map(names(patient_subsets), function(patient) {
  print(patient)
  generate_qc_report(
    experiment = c("221202_A01366_0326_AHHTTWDMXY",
                   "230210_A01366_0351_AHNHCFDSX5"),
    sublibrary = c("SHE5052A9_S101", "comb"),
    sample_subset = patient_subsets[[patient]],
    remove_doublets = F,
    do_cell_cycle_scoring = F,
    out_subdir = paste0("tumour_cells/", patient, "/"),
    rerun = F
  )
})

# run on tumour cells of each patient
purrr::map(names(patient_subsets), function(patient) {
  print(patient)
  generate_qc_report(
    experiment = c("221202_A01366_0326_AHHTTWDMXY",
                   "230210_A01366_0351_AHNHCFDSX5"),
    sublibrary = c("SHE5052A9_S101", "comb"),
    cell_subset = tumour_cells,
    sample_subset = patient_subsets[[patient]],
    remove_doublets = F,
    do_cell_cycle_scoring = F,
    out_subdir = paste0("tumour_cells/", patient, "/"),
    rerun = F
  )
})

# run on all tumour cells
generate_qc_report(
  experiment = c("221202_A01366_0326_AHHTTWDMXY",
                 "230210_A01366_0351_AHNHCFDSX5"),
  sublibrary = c("SHE5052A9_S101", "comb"),
  cell_subset = tumour_cells,
  sample_subset = sample_subset,
  remove_doublets = F,
  out_subdir = "tumour_cells"
)

