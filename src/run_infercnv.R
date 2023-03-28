base_dir <- ifelse(Sys.info()["nodename"] == "Alexs-MacBook-Air-2.local",
                   "/Volumes/TracerX/",
                   "/camp/project/tracerX/")
setwd(paste0(base_dir,"working/VHL_GERMLINE/tidda/vhl/"))
devtools::load_all()

out_dir <- "out/230210_A01366_0351_AHNHCFDSX5_x_221202_A01366_0326_AHHTTWDMXY/hg38/comb_x_SHE5052A9_S101/all-well/DGE_filtered/unintegrated/"
out <- get_out(out_dir)

# define included lineages
excl_lin <- "immune"
ref_lin <- "normal"
query_lin <- "malignant"

# do_each_tumour?
do_each_tumour = F
do_all_tumours = T

# run infercnv on matched normal and tumour cells from each patient

# get coldata
col_data <-
  SummarizedExperiment::colData(readRDS(paste0(out$cache, "cds_annotated.rds"))) %>%
  tibble::as_tibble(rownames = "cell")
patients <-
  col_data %>%
  dplyr::filter(lesion_type != "normal_renal") %>%
  dplyr::filter(nih_pid %ni% c("N045", "N059")) %>%
  dplyr::distinct(nih_pid, sample) %>%
  {split(.$sample, .$nih_pid)}

purrr::map(names(patients), function(patient) {

  if (do_all_tumours == T) {

  # all tumours per patient
  print(patient)

  # get outdir
  patient_out <- get_out(out_dir = paste0(out_dir, "/", patient))

  # create infercnv outdir
  dir.create(patient_out$infercnv, showWarnings = F, recursive = T)

  # write annotations_file for infercnv
  infercnv_annotations_file <- paste0(patient_out$infercnv, "infercnv_annotations.tsv")
  infercnv_annotations <-
   col_data %>%
   dplyr::filter(partition_lineage != excl_lin,
                 nih_pid == patient) %>%
   dplyr::transmute(
     cell,
     infercnv_lineage = dplyr::case_when(partition_lineage == "malignant" ~ paste0("malignant_", sample),
                                         TRUE ~ cluster_annot))
  ref_group_names <-
   infercnv_annotations %>%
   dplyr::filter(!grepl("malignant", infercnv_lineage)) %>%
   dplyr::pull(infercnv_lineage) %>%
   unique()
  write.table(infercnv_annotations, infercnv_annotations_file,
             row.names = F, col.names = F, quote = F, sep = "\t")

  # create infercnv object
  infercnv_obj <-
   infercnv::CreateInfercnvObject(
     raw_counts_matrix = as.matrix(readRDS(paste0(out$cache, "cds_annotated.rds"))@assays@data$counts),
     annotations_file = infercnv_annotations_file,
     gene_order_file = "data/gencode/gencode.v43.basic.annotation_clean.bed",
     ref_group_names = ref_group_names
   )

  options(scipen = 100)

  # perform infercnv operations to reveal cnv signal
  infercnv_obj <-
   infercnv::run(
     infercnv_obj,
     cutoff = 0.1,
     out_dir = patient_out$infercnv,
     cluster_by_groups = T,
     denoise = T,
     HMM = T, resume_mode = T
   )

  }

  if (do_each_tumour == T) {
    purrr::map(patients[[patient]], function(tumour) {

      # each tumour per patient
      print(patient) ; print(tumour)

      # get outdir
      tumour_out <- get_out(out_dir = paste0(out_dir, "/", patient, "/", tumour))

      # create infercnv outdir
      dir.create(tumour_out$infercnv, showWarnings = F, recursive = T)

      # write annotations_file for infercnv
      infercnv_annotations_file <- paste0(tumour_out$infercnv, "infercnv_annotations.tsv")
      infercnv_annotations <-
        col_data %>%
        dplyr::filter(partition_lineage != excl_lin,
                      nih_pid == patient,
                      sample == tumour | partition_lineage == "normal") %>%
        dplyr::transmute(
          cell,
          infercnv_lineage = dplyr::case_when(partition_lineage == "malignant" ~ paste0("malignant_", sample),
                                              TRUE ~ cluster_annot))
      ref_group_names <-
        infercnv_annotations %>%
        dplyr::filter(!grepl("malignant", infercnv_lineage)) %>%
        dplyr::pull(infercnv_lineage) %>%
        unique()
      write.table(infercnv_annotations, infercnv_annotations_file,
                  row.names = F, col.names = F, quote = F, sep = "\t")

      # create infercnv object
      infercnv_obj <-
        infercnv::CreateInfercnvObject(
          raw_counts_matrix = as.matrix(readRDS(paste0(out$cache, "cds_annotated.rds"))@assays@data$counts),
          annotations_file = infercnv_annotations_file,
          gene_order_file = "data/gencode/gencode.v43.basic.annotation_clean.bed",
          ref_group_names = ref_group_names
        )

      options(scipen = 100)

      # perform infercnv operations to reveal cnv signal
      infercnv_obj <-
        infercnv::run(
          infercnv_obj,
          cutoff = 0.1,
          out_dir = tumour_out$infercnv,
          cluster_by_groups = T,
          denoise = T,
          HMM = T, resume_mode = T
        )

    })
  }

})

