base_dir <- ifelse(Sys.info()["nodename"] == "Alexs-MacBook-Air-2.local",
                   "/Volumes/TracerX/",
                   "/camp/project/tracerX/")
setwd(paste0(base_dir,"working/VHL_GERMLINE/tidda/vhl/"))
devtools::load_all()

out_dir <- "out/230210_A01366_0351_AHNHCFDSX5_x_221202_A01366_0326_AHHTTWDMXY/hg38/comb_x_SHE5052A9_S101/all-well/DGE_filtered/unintegrated/"
out <- get_out(out_dir)

# load cds
cds <- readRDS(list.files(out$cache, pattern = "cds\\_singler\\_annotated", full.names = T))

# define partitions and lineages (for infercnv) (based on consensus_modules_vs_partitions_heatmap-1.png)
partition_annotations <-
  readr::read_tsv(paste0(out$base, "partition_annotations.tsv")) %>%
  dplyr::mutate(partition = factor(partition))

# encode in colData
SummarizedExperiment::colData(cds) <-
  cbind(
    SummarizedExperiment::colData(cds),
    monocle3::partitions(cds) %>%
      tibble::as_tibble(rownames = "cell") %>%
      dplyr::left_join(partition_annotations, by = c("value" = "partition")) %>%
      dplyr::select(partition_annot = annotation,
                    partition_lineage = lineage)
  )

# define clusters (from marker_notes.xlsx)
cluster_annotations <-
  readr::read_tsv(paste0(out$base, "cluster_annotations.tsv")) %>%
  dplyr::mutate(cluster = factor(cluster))

# encode in colData
SummarizedExperiment::colData(cds) <-
  cbind(
    SummarizedExperiment::colData(cds),
    monocle3::clusters(cds) %>%
      tibble::as_tibble(rownames = "cell") %>%
      dplyr::left_join(cluster_annotations, by = c("value" = "cluster")) %>%
      dplyr::select(cluster_annot = annotation)
  )

# save cell_annotations df
SummarizedExperiment::colData(cds) %>%
  tibble::as_tibble(rownames = "cell") %>%
  dplyr::select(cell, sample,
                partition, partition_annot, partition_lineage,
                cluster, cluster_annot) %>%
  readr::write_tsv(paste0(out$base, "cell_annotations.tsv"))

# save seu_annotated
saveRDS(cds, paste0(out$cache, "cds_annotated.rds"))
