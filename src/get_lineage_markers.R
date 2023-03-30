
gene_fits <-
  monocle3::fit_models(cds, "~partition_lineage + nih_pid")
